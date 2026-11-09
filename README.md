# Contador de Resíduos

Conta os aminoácidos de uma estrutura de proteína por tipo e lista o que não é aminoácido padrão (resíduos modificados, ligantes, íons, lipídios, açúcares).

**Página:** https://gabrielbadolato.github.io/contador-residuos/

Na página você pode enviar um arquivo ou buscar uma estrutura pelo código PDB (ex.: `7D9Z`). Na busca por código, a página também mostra as informações da estrutura vindas do RCSB PDB: título, método, resolução, datas, moléculas e organismos, quantos resíduos da sequência têm coordenadas e a publicação. Tudo roda no navegador; o arquivo não é enviado a nenhum servidor. Também dá para abrir direto uma estrutura pelo endereço: `https://gabrielbadolato.github.io/contador-residuos/?pdb=7D9Z`.

## Formatos aceitos

| Formato | Extensões | Observações |
|---|---|---|
| PDB | `.pdb` `.ent` | Inclusive saídas de GROMACS, AMBER, CHARMM/NAMD e CHARMM-GUI |
| PQR | `.pqr` | PDB2PQR / APBS |
| mmCIF | `.cif` `.mmcif` | PDBx/mmCIF do RCSB e PDBe |
| GRO | `.gro` | Coordenadas do GROMACS |
| gzip | `.gz` | Qualquer um dos anteriores compactado |

## Scripts

A pasta [`scripts`](scripts) tem a mesma contagem em quatro versões, que dão o mesmo resultado:

```bash
python contar_residuos.py entrada.pdb             # Python 3.6+, sem pacotes extras
bash contar_residuos.sh entrada.pdb               # Linux/macOS, só bash, awk e gzip
.\contar_residuos.ps1 entrada.pdb                 # Windows PowerShell 5.1 ou PowerShell 7
```

No VMD, com a molécula já aberta, no Tk Console (Extensions > Tk Console):

```tcl
source /caminho/completo/contar_residuos.tcl
contar_residuos top "chain A" contagem.txt
```

Se o Windows bloquear a execução de scripts do PowerShell:
`powershell -ExecutionPolicy Bypass -File .\contar_residuos.ps1 entrada.pdb`

## Regras de contagem (resumo)

- Um resíduo novo começa quando cadeia, segmento, número ou código de inserção mudam de um átomo para o seguinte. Isso funciona em arquivos sem cadeia e com numeração que reinicia após 9999.
- Conformações alternativas (altloc) contam uma vez.
- Nomes de campos de força (HID, HSD, HISE, CYX, CYS2, LYSH...) e terminais (NALA, CLYS...) contam como o aminoácido padrão.
- Águas são ignoradas. Só o primeiro modelo ou quadro é lido.
- Resíduos sem coordenadas no arquivo não são contados.

A explicação completa, com a validação, está na própria página.
