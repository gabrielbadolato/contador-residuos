#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
contar_residuos.py - Conta os residuos de aminoacidos de uma estrutura, por tipo.

Uso:
    python contar_residuos.py entrada.pdb
    python contar_residuos.py entrada.cif -o contagem.txt

Formatos aceitos (reconhecidos pela extensao e pelo conteudo):
    PDB   (.pdb .ent .brk)  Protein Data Bank e saidas de GROMACS, AMBER,
                            CHARMM/NAMD, CHARMM-GUI, VMD, PyMOL
    PQR   (.pqr)            PDB2PQR / APBS (colunas separadas por espaco)
    mmCIF (.cif .mmcif)     PDBx/mmCIF do RCSB/PDBe
    GRO   (.gro)            coordenadas do GROMACS
    Qualquer um deles compactado em .gz (ex.: 1abc.cif.gz).

Saida: quantidade de cada um dos 20 aminoacidos, total, nomes convertidos e a
lista de residuos nao identificados (fora dos 20: modificados, ligantes, ions,
lipidios, acucares, nucleotideos...). A mesma saida e gravada em um arquivo
de texto (por padrao, <entrada>_contagem.txt).

Regras de contagem:
- Um residuo novo comeca quando cadeia, segmento, numero ou codigo de
  insercao mudam de um atomo para o seguinte (ou apos TER). Assim funcionam
  arquivos sem cadeia e numeracao que reinicia apos 9999 ou 99999.
- Conformacoes alternativas (altloc) do mesmo residuo contam uma vez.
- Variantes de campo de forca (HID, HSD, HISE, CYX, CYS2, ASH, GLH, LYN...)
  e terminais do AMBER/GROMACS (NALA, CLYS...) contam como o aminoacido padrao.
- Aguas sao ignoradas. Em arquivos com varios modelos ou quadros, so o
  primeiro e lido.
"""

import argparse
import gzip
import os
import re
import sys
from collections import OrderedDict

AMINOACIDOS = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS",
               "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP",
               "TYR", "VAL"]

SINONIMOS = {
    # histidina (AMBER, CHARMM, GROMACS, GROMOS)
    "HID": "HIS", "HIE": "HIS", "HIP": "HIS", "HIN": "HIS", "HSD": "HIS",
    "HSE": "HIS", "HSP": "HIS", "HISD": "HIS", "HISE": "HIS", "HISH": "HIS",
    "HISP": "HIS", "HISA": "HIS", "HISB": "HIS", "HIS1": "HIS", "HIS2": "HIS",
    # cisteina (ponte dissulfeto, desprotonada, protonada)
    "CYX": "CYS", "CYM": "CYS", "CYN": "CYS", "CYS2": "CYS", "CYSH": "CYS",
    # acidos protonados
    "ASH": "ASP", "ASPH": "ASP", "ASPP": "ASP",
    "GLH": "GLU", "GLUH": "GLU", "GLUP": "GLU",
    # lisina e arginina neutras/protonadas
    "LYN": "LYS", "LYSH": "LYS", "LYSN": "LYS", "LSN": "LYS", "ARGN": "ARG",
}

AGUAS = {"HOH", "WAT", "H2O", "DOD", "D2O", "SOL", "TIP", "TIP2", "TIP3",
         "TIP4", "TIP5", "TP3", "TP4", "T3P", "T4P", "T5P", "SPC", "SPCE",
         "SWM4", "OPC", "OPC3"}


def padronizar(nome):
    """Aminoacido padrao correspondente ao nome, ou None."""
    if nome in AMINOACIDOS:
        return nome
    if nome in SINONIMOS:
        return SINONIMOS[nome]
    if len(nome) == 4 and nome[0] in "NC":      # NALA, CLYS, NHIE...
        resto = nome[1:]
        if resto in AMINOACIDOS:
            return resto
        if resto in SINONIMOS:
            return SINONIMOS[resto]
    return None


# ---------------------------------------------------------------- leitura

def abrir_linhas(caminho):
    with open(caminho, "rb") as f:
        inicio = f.read(2)
    if inicio == b"\x1f\x8b":                   # gzip
        dados = gzip.open(caminho, "rb").read()
    else:
        dados = open(caminho, "rb").read()
    return dados.decode("latin-1").splitlines()


def detectar_formato(caminho, linhas):
    nome = caminho.lower()
    if nome.endswith(".gz"):
        nome = nome[:-3]
    ext = os.path.splitext(nome)[1]
    if ext in (".cif", ".mmcif", ".pdbx"):
        return "mmCIF"
    if ext == ".gro":
        return "GRO"
    if ext == ".pqr":
        return "PQR"
    if ext in (".pdb", ".ent", ".brk", ".pdbqt"):
        return "PDB"
    # sem extensao conhecida: olha o conteudo
    for l in linhas[:200]:
        if l.startswith("data_") or l.startswith("_atom_site."):
            return "mmCIF"
        if l.startswith(("ATOM", "HETATM")):
            return "PDB"
    if len(linhas) > 2 and linhas[1].strip().isdigit():
        return "GRO"
    return "PDB"


def ler_pdb(linhas):
    """PDB de colunas fixas. Devolve lista de (cadeia, numero, nome)."""
    residuos, anterior = [], None
    for l in linhas:
        reg = l[:6]
        if reg == "ENDMDL":
            break                                # so o primeiro modelo
        if reg.startswith("TER"):
            anterior = None
            continue
        if reg not in ("ATOM  ", "HETATM") or len(l) < 22:
            continue
        # cadeia (22) + numero (23-26) + insercao (27) + segmento (73-76)
        chave = l[21:27] + "|" + l[72:76]
        if chave == anterior:
            continue
        anterior = chave
        cadeia = l[21].strip() or l[72:76].strip()
        numero = l[22:27].strip()
        nome = l[17:21].strip().upper()          # aceita nomes de 4 letras
        residuos.append((cadeia, numero, nome))
    return residuos


def ler_pqr(linhas):
    """PQR (campos separados por espaco, cadeia opcional)."""
    residuos, anterior = [], None
    for l in linhas:
        if l.startswith("ENDMDL"):
            break
        if l.startswith("TER"):
            anterior = None
            continue
        if not l.startswith(("ATOM", "HETATM")):
            continue
        c = l.split()
        if len(c) < 10:
            continue
        # ATOM serial atomo residuo [cadeia] numero x y z carga raio
        if len(c) >= 11 and not re.match(r"^-?\d", c[4]):
            cadeia, numero = c[4], c[5]
        else:
            cadeia, numero = "", c[4]
        chave = (cadeia, numero)
        if chave == anterior:
            continue
        anterior = chave
        residuos.append((cadeia, numero, c[3].upper()))
    return residuos


def tokens_cif(linha):
    """Separa uma linha de mmCIF em campos, respeitando aspas."""
    return [a or b or c for a, b, c in
            re.findall(r"'(.*?)'(?=\s|$)|\"(.*?)\"(?=\s|$)|(\S+)", linha)]


def ler_mmcif(linhas):
    """mmCIF: le o bloco loop_ _atom_site."""
    campos, valores = [], []
    i, n = 0, len(linhas)
    while i < n:                                 # acha o loop de _atom_site
        if linhas[i].strip() == "loop_" and i + 1 < n \
                and linhas[i + 1].strip().startswith("_atom_site."):
            i += 1
            while i < n and linhas[i].strip().startswith("_atom_site."):
                campos.append(linhas[i].strip().split()[0][len("_atom_site."):])
                i += 1
            break
        i += 1
    if not campos:
        return []
    while i < n:
        l = linhas[i].strip()
        if l.startswith(("_", "loop_", "data_", "#")) and valores \
                and len(valores) % len(campos) == 0:
            break
        if l and not l.startswith("#"):
            valores.extend(tokens_cif(l))
        i += 1

    def col(*nomes):
        for nm in nomes:
            if nm in campos:
                return campos.index(nm)
        return None

    k_nome = col("label_comp_id", "auth_comp_id")
    k_cad = col("auth_asym_id", "label_asym_id")
    k_num = col("auth_seq_id", "label_seq_id")
    k_ins = col("pdbx_PDB_ins_code")
    k_mod = col("pdbx_PDB_model_num")
    k_ent = col("label_entity_id")
    if k_nome is None or k_num is None:
        return []
    residuos, anterior, primeiro_modelo = [], None, None
    nc = len(campos)
    for j in range(0, len(valores) - nc + 1, nc):
        r = valores[j:j + nc]
        if k_mod is not None:
            if primeiro_modelo is None:
                primeiro_modelo = r[k_mod]
            elif r[k_mod] != primeiro_modelo:
                break
        cadeia = r[k_cad] if k_cad is not None else ""
        numero = r[k_num]
        ins = r[k_ins] if k_ins is not None else ""
        ins = "" if ins in ("?", ".") else ins
        ent = r[k_ent] if k_ent is not None else ""
        if numero in ("?", "."):                 # sem numero: usa a entidade
            numero = ""
        chave = (cadeia, numero, ins, ent if numero == "" else "")
        if chave == anterior:
            continue
        anterior = chave
        residuos.append((cadeia, numero + ins, r[k_nome].upper()))
    return residuos


def ler_gro(linhas):
    """GRO: primeira linha titulo, segunda numero de atomos."""
    try:
        natomos = int(linhas[1].split()[0])
    except (IndexError, ValueError):
        return []
    residuos, anterior = [], None
    for l in linhas[2:2 + natomos]:             # so o primeiro quadro
        chave = l[0:10]                          # numero (1-5) + nome (6-10)
        if chave == anterior:
            continue
        anterior = chave
        residuos.append(("", l[0:5].strip(), l[5:10].strip().upper()))
    return residuos


LEITORES = {"PDB": ler_pdb, "PQR": ler_pqr, "mmCIF": ler_mmcif, "GRO": ler_gro}


# ---------------------------------------------------------------- contagem

def contar(residuos):
    contagem = OrderedDict((aa, 0) for aa in AMINOACIDOS)
    convertidos = OrderedDict()
    nao_identificados = []
    for cadeia, numero, nome in residuos:
        aa = padronizar(nome)
        if aa:
            contagem[aa] += 1
            if aa != nome:
                k = f"{nome}->{aa}"
                convertidos[k] = convertidos.get(k, 0) + 1
        elif nome not in AGUAS:
            nao_identificados.append((cadeia, numero, nome))
    return contagem, convertidos, nao_identificados


def relatorio(arquivo, formato, contagem, convertidos, nao_identificados):
    out = [f"Arquivo: {arquivo}", f"Formato: {formato}", "",
           "Residuo  Quantidade", "-------------------"]
    for aa, q in contagem.items():
        out.append(f"{aa:<8}{q:>11d}")
    out += ["-------------------",
            f"{'TOTAL':<8}{sum(contagem.values()):>11d}", ""]
    if convertidos:
        out.append("Nomes convertidos para o aminoacido padrao (ja incluidos acima):")
        out.append("  " + ", ".join(f"{k} x{v}" for k, v in convertidos.items()))
        out.append("")
    if nao_identificados:
        resumo = OrderedDict()
        for _, _, nome in nao_identificados:
            resumo[nome] = resumo.get(nome, 0) + 1
        out.append(f"Residuos nao identificados: {len(nao_identificados)} "
                   "(fora da contagem acima)")
        out.append("  Resumo: " + ", ".join(f"{k} x{v}" for k, v in resumo.items()))
        out.append("  Nome  Cadeia  Numero")
        for cadeia, numero, nome in nao_identificados:
            out.append(f"  {nome:<5} {cadeia or '-':<7} {numero}")
    else:
        out.append("Residuos nao identificados: nenhum")
    return "\n".join(out) + "\n"


def main():
    p = argparse.ArgumentParser(description="Conta os residuos de uma estrutura "
                                "(PDB, PQR, mmCIF, GRO, tambem .gz).")
    p.add_argument("arquivo", help="arquivo de entrada")
    p.add_argument("-o", "--saida", help="arquivo de saida (texto)")
    a = p.parse_args()

    if not os.path.isfile(a.arquivo):
        sys.exit(f"Arquivo nao encontrado: {a.arquivo}")
    linhas = abrir_linhas(a.arquivo)
    formato = detectar_formato(a.arquivo, linhas)
    residuos = LEITORES[formato](linhas)
    if not residuos:
        sys.exit(f"Nenhum atomo encontrado (formato lido: {formato}).")

    texto = relatorio(os.path.basename(a.arquivo), formato, *contar(residuos))
    print(texto, end="")
    base = a.arquivo[:-3] if a.arquivo.lower().endswith(".gz") else a.arquivo
    saida = a.saida or os.path.splitext(base)[0] + "_contagem.txt"
    try:
        with open(saida, "w", encoding="utf-8") as f:
            f.write(texto)
        print(f"\nSaida gravada em: {saida}")
    except OSError as erro:
        print(f"\nNao foi possivel gravar {saida}: {erro}")


if __name__ == "__main__":
    main()
