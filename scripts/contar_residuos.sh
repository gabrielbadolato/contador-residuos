#!/bin/bash
# contar_residuos.sh - Conta os residuos de aminoacidos de uma estrutura, por tipo.
#
# Uso:
#     bash contar_residuos.sh entrada.pdb
#     bash contar_residuos.sh entrada.cif contagem.txt
#
# Formatos aceitos (reconhecidos pela extensao e pelo conteudo):
#     PDB   (.pdb .ent .brk)  Protein Data Bank e saidas de GROMACS, AMBER,
#                             CHARMM/NAMD, CHARMM-GUI, VMD, PyMOL
#     PQR   (.pqr)            PDB2PQR / APBS
#     mmCIF (.cif .mmcif)     PDBx/mmCIF do RCSB/PDBe
#     GRO   (.gro)            coordenadas do GROMACS
#     Qualquer um deles compactado em .gz.
#
# Saida: quantidade de cada um dos 20 aminoacidos, total, nomes convertidos e
# lista dos residuos nao identificados. Tambem gravada em arquivo (padrao:
# <entrada>_contagem.txt).
#
# Regras: um residuo novo comeca quando cadeia, segmento, numero ou codigo de
# insercao mudam de um atomo para o seguinte; altlocs contam uma vez;
# variantes de campo de forca (HID, HSD, HISE, CYX, CYS2...) e terminais
# (NALA, CLYS...) contam como o aminoacido padrao; aguas sao ignoradas; so o
# primeiro modelo/quadro e lido.
# Usa apenas bash, awk e gzip, presentes em qualquer Linux ou macOS.

if [ $# -lt 1 ]; then
    echo "Uso: bash $0 entrada [saida.txt]"
    exit 1
fi
arq="$1"
if [ ! -f "$arq" ]; then
    echo "Arquivo nao encontrado: $arq"
    exit 1
fi

# Descompacta se for gzip
ler() {
    if [ "$(head -c 2 "$arq" | od -An -tx1 | tr -d ' \n')" = "1f8b" ]; then
        gzip -dc "$arq"
    else
        cat "$arq"
    fi
}

# Formato pela extensao (ignorando .gz) e, se preciso, pelo conteudo
base="${arq%.gz}"; base="${base%.GZ}"
ext="$(echo "${base##*.}" | tr '[:upper:]' '[:lower:]')"
[ "$base" = "${base%.*}" ] && ext=""
case "$ext" in
    cif|mmcif|pdbx)      fmt="mmCIF" ;;
    gro)                 fmt="GRO" ;;
    pqr)                 fmt="PQR" ;;
    pdb|ent|brk|pdbqt)   fmt="PDB" ;;
    *)
        if ler | head -200 | grep -qE '^(data_|_atom_site\.)'; then fmt="mmCIF"
        elif ler | head -200 | grep -qE '^(ATOM|HETATM)'; then fmt="PDB"
        elif ler | sed -n 2p | grep -qE '^ *[0-9]+ *$'; then fmt="GRO"
        else fmt="PDB"; fi ;;
esac

saida="${2:-${base%.*}_contagem.txt}"

ler | awk -v arquivo="$(basename "$arq")" -v fmt="$fmt" '
function padronizar(nome,   resto, p) {
    if (nome in padrao) return nome
    if (nome in sinonimo) return sinonimo[nome]
    p = substr(nome, 1, 1)
    if (length(nome) == 4 && (p == "N" || p == "C")) {
        resto = substr(nome, 2)
        if (resto in padrao) return resto
        if (resto in sinonimo) return sinonimo[resto]
    }
    return ""
}
function trim(s) { gsub(/^[ \t\r]+|[ \t\r]+$/, "", s); return s }
function registrar(cadeia, numero, nome,   p, c) {
    nome = toupper(nome)
    p = padronizar(nome)
    if (p != "") {
        cont[p]++
        if (p != nome) {
            c = nome "->" p
            if (!(c in conv)) { nc++; ordc[nc] = c }
            conv[c]++
        }
        return
    }
    if (nome in agua) return
    if (cadeia == "") cadeia = "-"
    nd++; dnome[nd] = nome; dcad[nd] = cadeia; dnum[nd] = numero
    if (!(nome in resumo)) { nr++; ordem[nr] = nome }
    resumo[nome]++
}
# separa uma linha de mmCIF em campos, respeitando aspas
function tokens_cif(s,   i, n, ch, q, tok, dentro) {
    n = length(s); i = 1
    while (i <= n) {
        ch = substr(s, i, 1)
        if (ch == " " || ch == "\t" || ch == "\r") { i++; continue }
        if (ch == "\047" || ch == "\"") {
            q = ch; tok = ""; i++
            while (i <= n) {
                ch = substr(s, i, 1)
                if (ch == q && (i == n || substr(s, i + 1, 1) ~ /[ \t\r]/)) { i++; break }
                tok = tok ch; i++
            }
        } else {
            tok = ""
            while (i <= n) {
                ch = substr(s, i, 1)
                if (ch == " " || ch == "\t" || ch == "\r") break
                tok = tok ch; i++
            }
        }
        nv++; val[nv] = tok
    }
}
function linha_cif(   j, base, cad, num, ins, ent, chave) {
    while (nv - usado >= ncampos) {
        base = usado
        usado += ncampos
        if (k_mod) {
            if (primeiro_modelo == "") primeiro_modelo = val[base + k_mod]
            else if (val[base + k_mod] != primeiro_modelo) { fim = 1; return }
        }
        cad = k_cad ? val[base + k_cad] : ""
        num = val[base + k_num]
        ins = k_ins ? val[base + k_ins] : ""
        if (ins == "?" || ins == ".") ins = ""
        ent = k_ent ? val[base + k_ent] : ""
        if (num == "?" || num == ".") num = ""; else ent = ""
        chave = cad SUBSEP num SUBSEP ins SUBSEP ent
        if (chave != anterior) { anterior = chave; registrar(cad, num ins, val[base + k_nome]) }
        for (j = base + 1; j <= usado; j++) delete val[j]
    }
}
BEGIN {
    n = split("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL", aa, " ")
    for (i = 1; i <= n; i++) { padrao[aa[i]] = 1; cont[aa[i]] = 0 }
    m = split("HID HIE HIP HIN HSD HSE HSP HISD HISE HISH HISP HISA HISB HIS1 HIS2", t, " ")
    for (i = 1; i <= m; i++) sinonimo[t[i]] = "HIS"
    m = split("CYX CYM CYN CYS2 CYSH", t, " ");   for (i = 1; i <= m; i++) sinonimo[t[i]] = "CYS"
    m = split("ASH ASPH ASPP", t, " ");           for (i = 1; i <= m; i++) sinonimo[t[i]] = "ASP"
    m = split("GLH GLUH GLUP", t, " ");           for (i = 1; i <= m; i++) sinonimo[t[i]] = "GLU"
    m = split("LYN LYSH LYSN LSN", t, " ");       for (i = 1; i <= m; i++) sinonimo[t[i]] = "LYS"
    sinonimo["ARGN"] = "ARG"
    m = split("HOH WAT H2O DOD D2O SOL TIP TIP2 TIP3 TIP4 TIP5 TP3 TP4 T3P T4P T5P SPC SPCE SWM4 OPC OPC3", t, " ")
    for (i = 1; i <= m; i++) agua[t[i]] = 1
    nd = 0; nr = 0; nc = 0; anterior = ""; fim = 0
    estado = 0; ncampos = 0; nv = 0; usado = 0; primeiro_modelo = ""
}
fim { next }

# ---------- PDB (colunas fixas)
fmt == "PDB" {
    if ($0 ~ /^ENDMDL/) { fim = 1; next }
    if ($0 ~ /^TER/) { anterior = ""; next }
    if ($0 !~ /^(ATOM  |HETATM)/ || length($0) < 22) next
    chave = substr($0, 22, 6) "|" substr($0, 73, 4)
    if (chave == anterior) next
    anterior = chave
    cadeia = trim(substr($0, 22, 1)); if (cadeia == "") cadeia = trim(substr($0, 73, 4))
    registrar(cadeia, trim(substr($0, 23, 5)), trim(substr($0, 18, 4)))
    next
}

# ---------- PQR (campos separados por espaco)
fmt == "PQR" {
    if ($0 ~ /^ENDMDL/) { fim = 1; next }
    if ($0 ~ /^TER/) { anterior = ""; next }
    if ($0 !~ /^(ATOM|HETATM)/ || NF < 10) next
    if (NF >= 11 && $5 !~ /^-?[0-9]/) { cadeia = $5; numero = $6 } else { cadeia = ""; numero = $5 }
    chave = cadeia SUBSEP numero
    if (chave == anterior) next
    anterior = chave
    registrar(cadeia, numero, $4)
    next
}

# ---------- GRO (so o primeiro quadro)
fmt == "GRO" {
    if (NR == 2) { natomos = $1 + 0; next }
    if (NR < 3) next
    if (NR > natomos + 2) { fim = 1; next }
    chave = substr($0, 1, 10)
    if (chave == anterior) next
    anterior = chave
    registrar("", trim(substr($0, 1, 5)), trim(substr($0, 6, 5)))
    next
}

# ---------- mmCIF (loop _atom_site)
fmt == "mmCIF" {
    s = trim($0)
    if (estado == 0) { if (s == "loop_") estado = 1; next }
    if (estado == 1) {
        if (s ~ /^_atom_site\./) {
            split(s, partes, " ")
            ncampos++; nomecampo = substr(partes[1], 12)
            if (nomecampo == "label_comp_id") k_nome = ncampos
            if (nomecampo == "auth_comp_id")  k_auth_nome = ncampos
            if (nomecampo == "auth_asym_id")  k_cad = ncampos
            if (nomecampo == "label_asym_id") k_lcad = ncampos
            if (nomecampo == "auth_seq_id")   k_num = ncampos
            if (nomecampo == "label_seq_id")  k_lnum = ncampos
            if (nomecampo == "pdbx_PDB_ins_code") k_ins = ncampos
            if (nomecampo == "pdbx_PDB_model_num") k_mod = ncampos
            if (nomecampo == "label_entity_id") k_ent = ncampos
            next
        }
        if (ncampos == 0) { estado = (s == "loop_") ? 1 : 0; next }
        if (!k_nome) k_nome = k_auth_nome
        if (!k_cad) k_cad = k_lcad
        if (!k_num) k_num = k_lnum
        if (!k_nome || !k_num) { fim = 1; next }
        estado = 2
    }
    if (estado == 2) {
        if (s ~ /^(_|loop_|data_|#)/ && (nv - usado) == 0 && usado > 0) { fim = 1; next }
        if (s != "" && s !~ /^#/) { tokens_cif(s); linha_cif() }
    }
    next
}

END {
    print "Arquivo: " arquivo
    print "Formato: " fmt
    print ""
    print "Residuo  Quantidade"
    print "-------------------"
    total = 0
    for (i = 1; i <= n; i++) { printf "%-8s%11d\n", aa[i], cont[aa[i]]; total += cont[aa[i]] }
    print "-------------------"
    printf "%-8s%11d\n", "TOTAL", total
    print ""
    if (nc > 0) {
        print "Nomes convertidos para o aminoacido padrao (ja incluidos acima):"
        linha = "  "
        for (i = 1; i <= nc; i++) linha = linha (i > 1 ? ", " : "") ordc[i] " x" conv[ordc[i]]
        print linha
        print ""
    }
    if (nd > 0) {
        print "Residuos nao identificados: " nd " (fora da contagem acima)"
        linha = "  Resumo: "
        for (i = 1; i <= nr; i++) linha = linha (i > 1 ? ", " : "") ordem[i] " x" resumo[ordem[i]]
        print linha
        print "  Nome  Cadeia  Numero"
        for (i = 1; i <= nd; i++) printf "  %-5s %-7s %s\n", dnome[i], dcad[i], dnum[i]
    } else {
        print "Residuos nao identificados: nenhum"
    }
}' | tee "$saida"

echo
echo "Saida gravada em: $saida"
