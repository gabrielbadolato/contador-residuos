<#
contar_residuos.ps1 - Conta os residuos de aminoacidos de uma estrutura, por tipo.

Uso (PowerShell no Windows, Linux ou macOS):
    .\contar_residuos.ps1 entrada.pdb
    .\contar_residuos.ps1 entrada.cif contagem.txt

Se o Windows bloquear a execucao de scripts, rode uma vez:
    powershell -ExecutionPolicy Bypass -File .\contar_residuos.ps1 entrada.pdb

Formatos aceitos (reconhecidos pela extensao e pelo conteudo):
    PDB   (.pdb .ent .brk)  Protein Data Bank e saidas de GROMACS, AMBER,
                            CHARMM/NAMD, CHARMM-GUI, VMD, PyMOL
    PQR   (.pqr)            PDB2PQR / APBS
    mmCIF (.cif .mmcif)     PDBx/mmCIF do RCSB/PDBe
    GRO   (.gro)            coordenadas do GROMACS
    Qualquer um deles compactado em .gz.

Saida: quantidade de cada um dos 20 aminoacidos, total, nomes convertidos e
lista dos residuos nao identificados. Tambem gravada em arquivo (padrao:
<entrada>_contagem.txt).

Regras: um residuo novo comeca quando cadeia, segmento, numero ou codigo de
insercao mudam de um atomo para o seguinte; altlocs contam uma vez;
variantes de campo de forca (HID, HSD, HISE, CYX, CYS2...) e terminais
(NALA, CLYS...) contam como o aminoacido padrao; aguas sao ignoradas; so o
primeiro modelo/quadro e lido. Funciona no Windows PowerShell 5.1 e no
PowerShell 7.
#>
param(
    [Parameter(Mandatory = $true, Position = 0)][string]$Arquivo,
    [Parameter(Position = 1)][string]$Saida = ""
)

$ErrorActionPreference = "Stop"

$AMINOACIDOS = @("ALA","ARG","ASN","ASP","CYS","GLN","GLU","GLY","HIS","ILE",
                 "LEU","LYS","MET","PHE","PRO","SER","THR","TRP","TYR","VAL")
$SINONIMOS = @{
    "HID"="HIS"; "HIE"="HIS"; "HIP"="HIS"; "HIN"="HIS"; "HSD"="HIS"; "HSE"="HIS";
    "HSP"="HIS"; "HISD"="HIS"; "HISE"="HIS"; "HISH"="HIS"; "HISP"="HIS";
    "HISA"="HIS"; "HISB"="HIS"; "HIS1"="HIS"; "HIS2"="HIS";
    "CYX"="CYS"; "CYM"="CYS"; "CYN"="CYS"; "CYS2"="CYS"; "CYSH"="CYS";
    "ASH"="ASP"; "ASPH"="ASP"; "ASPP"="ASP"; "GLH"="GLU"; "GLUH"="GLU"; "GLUP"="GLU";
    "LYN"="LYS"; "LYSH"="LYS"; "LYSN"="LYS"; "LSN"="LYS"; "ARGN"="ARG"
}
$AGUAS = @("HOH","WAT","H2O","DOD","D2O","SOL","TIP","TIP2","TIP3","TIP4","TIP5",
           "TP3","TP4","T3P","T4P","T5P","SPC","SPCE","SWM4","OPC","OPC3")

function Padronizar([string]$nome) {
    if ($AMINOACIDOS -ccontains $nome) { return $nome }
    if ($SINONIMOS.ContainsKey($nome)) { return $SINONIMOS[$nome] }
    if ($nome.Length -eq 4 -and ($nome[0] -ceq 'N' -or $nome[0] -ceq 'C')) {
        $resto = $nome.Substring(1)
        if ($AMINOACIDOS -ccontains $resto) { return $resto }
        if ($SINONIMOS.ContainsKey($resto)) { return $SINONIMOS[$resto] }
    }
    return ""
}

function Corte([string]$s, [int]$ini, [int]$tam) {
    # substring segura (colunas a partir de 0)
    if ($ini -ge $s.Length) { return "" }
    if ($ini + $tam -gt $s.Length) { $tam = $s.Length - $ini }
    return $s.Substring($ini, $tam)
}

# ---------------------------------------------------------------- leitura
if (-not (Test-Path -LiteralPath $Arquivo -PathType Leaf)) {
    Write-Host "Arquivo nao encontrado: $Arquivo"; exit 1
}
$caminho = (Resolve-Path -LiteralPath $Arquivo).Path
$bytes = [System.IO.File]::ReadAllBytes($caminho)
if ($bytes.Length -ge 2 -and $bytes[0] -eq 0x1f -and $bytes[1] -eq 0x8b) {
    $entrada = New-Object System.IO.MemoryStream(,$bytes)
    $gz = New-Object System.IO.Compression.GZipStream($entrada, [System.IO.Compression.CompressionMode]::Decompress)
    $saidaMem = New-Object System.IO.MemoryStream
    $gz.CopyTo($saidaMem); $gz.Close()
    $bytes = $saidaMem.ToArray()
}
$texto = [System.Text.Encoding]::GetEncoding(28591).GetString($bytes)
$linhas = $texto -split "\r?\n"

$nomeArq = [System.IO.Path]::GetFileName($caminho)
$base = $nomeArq -replace '\.gz$', ''
$ext = [System.IO.Path]::GetExtension($base).ToLower()
switch ($ext) {
    { $_ -in ".cif", ".mmcif", ".pdbx" } { $fmt = "mmCIF"; break }
    ".gro" { $fmt = "GRO"; break }
    ".pqr" { $fmt = "PQR"; break }
    { $_ -in ".pdb", ".ent", ".brk", ".pdbqt" } { $fmt = "PDB"; break }
    default {
        $fmt = ""
        foreach ($l in ($linhas | Select-Object -First 200)) {
            if ($l.StartsWith("data_") -or $l.StartsWith("_atom_site.")) { $fmt = "mmCIF"; break }
            if ($l.StartsWith("ATOM") -or $l.StartsWith("HETATM")) { $fmt = "PDB"; break }
        }
        if ($fmt -eq "") {
            if ($linhas.Count -gt 2 -and $linhas[1].Trim() -match '^\d+$') { $fmt = "GRO" } else { $fmt = "PDB" }
        }
    }
}

# Cada item: @(cadeia, numero, nome)
$residuos = New-Object System.Collections.Generic.List[object]

if ($fmt -eq "PDB") {
    $anterior = $null
    foreach ($l in $linhas) {
        if ($l.StartsWith("ENDMDL")) { break }
        if ($l.StartsWith("TER")) { $anterior = $null; continue }
        if (-not ($l.StartsWith("ATOM  ") -or $l.StartsWith("HETATM")) -or $l.Length -lt 22) { continue }
        $chave = (Corte $l 21 6) + "|" + (Corte $l 72 4)
        if ($chave -ceq $anterior) { continue }
        $anterior = $chave
        $cadeia = (Corte $l 21 1).Trim()
        if ($cadeia -eq "") { $cadeia = (Corte $l 72 4).Trim() }
        $residuos.Add(@($cadeia, (Corte $l 22 5).Trim(), (Corte $l 17 4).Trim().ToUpper()))
    }
}
elseif ($fmt -eq "PQR") {
    $anterior = $null
    foreach ($l in $linhas) {
        if ($l.StartsWith("ENDMDL")) { break }
        if ($l.StartsWith("TER")) { $anterior = $null; continue }
        if (-not ($l.StartsWith("ATOM") -or $l.StartsWith("HETATM"))) { continue }
        $c = $l.Trim() -split '\s+'
        if ($c.Count -lt 10) { continue }
        if ($c.Count -ge 11 -and $c[4] -notmatch '^-?\d') { $cadeia = $c[4]; $numero = $c[5] }
        else { $cadeia = ""; $numero = $c[4] }
        $chave = "$cadeia|$numero"
        if ($chave -ceq $anterior) { continue }
        $anterior = $chave
        $residuos.Add(@($cadeia, $numero, $c[3].ToUpper()))
    }
}
elseif ($fmt -eq "GRO") {
    $natomos = 0
    if ($linhas.Count -gt 1) { [int]::TryParse(($linhas[1].Trim() -split '\s+')[0], [ref]$natomos) | Out-Null }
    $anterior = $null
    $fimGro = [Math]::Min($linhas.Count - 1, $natomos + 1)
    for ($i = 2; $i -le $fimGro; $i++) {
        $l = $linhas[$i]
        $chave = Corte $l 0 10
        if ($chave -ceq $anterior) { continue }
        $anterior = $chave
        $residuos.Add(@("", (Corte $l 0 5).Trim(), (Corte $l 5 5).Trim().ToUpper()))
    }
}
else {
    # mmCIF: bloco loop_ _atom_site
    $campos = New-Object System.Collections.Generic.List[string]
    $i = 0; $n = $linhas.Count
    while ($i -lt $n) {
        if ($linhas[$i].Trim() -eq "loop_" -and ($i + 1) -lt $n -and $linhas[$i + 1].Trim().StartsWith("_atom_site.")) {
            $i++
            while ($i -lt $n -and $linhas[$i].Trim().StartsWith("_atom_site.")) {
                $campos.Add((($linhas[$i].Trim() -split '\s+')[0]).Substring(11)); $i++
            }
            break
        }
        $i++
    }
    function Col([string[]]$nomes) {
        foreach ($nm in $nomes) { $k = $campos.IndexOf($nm); if ($k -ge 0) { return $k } }
        return -1
    }
    $kNome = Col @("label_comp_id", "auth_comp_id")
    $kCad = Col @("auth_asym_id", "label_asym_id")
    $kNum = Col @("auth_seq_id", "label_seq_id")
    $kIns = Col @("pdbx_PDB_ins_code")
    $kMod = Col @("pdbx_PDB_model_num")
    $kEnt = Col @("label_entity_id")
    $nc = $campos.Count
    if ($nc -gt 0 -and $kNome -ge 0 -and $kNum -ge 0) {
        $tok = New-Object System.Collections.Generic.List[string]
        $re = [regex]"'(.*?)'(?=\s|$)|""(.*?)""(?=\s|$)|(\S+)"
        $anterior = $null; $primeiroModelo = $null; $parar = $false; $usado = 0
        while ($i -lt $n -and -not $parar) {
            $l = $linhas[$i].Trim(); $i++
            if (($l.StartsWith("_") -or $l.StartsWith("loop_") -or $l.StartsWith("data_") -or $l.StartsWith("#")) -and $usado -gt 0 -and $tok.Count -eq 0) { break }
            if ($l -eq "" -or $l.StartsWith("#")) { continue }
            foreach ($m in $re.Matches($l)) {
                if ($m.Groups[1].Success) { $tok.Add($m.Groups[1].Value) }
                elseif ($m.Groups[2].Success) { $tok.Add($m.Groups[2].Value) }
                else { $tok.Add($m.Groups[3].Value) }
            }
            while ($tok.Count -ge $nc) {
                $r = $tok.GetRange(0, $nc); $tok.RemoveRange(0, $nc); $usado++
                if ($kMod -ge 0) {
                    if ($null -eq $primeiroModelo) { $primeiroModelo = $r[$kMod] }
                    elseif ($r[$kMod] -ne $primeiroModelo) { $parar = $true; break }
                }
                $cadeia = ""; if ($kCad -ge 0) { $cadeia = $r[$kCad] }
                $numero = $r[$kNum]
                $ins = ""; if ($kIns -ge 0) { $ins = $r[$kIns] }
                if ($ins -eq "?" -or $ins -eq ".") { $ins = "" }
                $ent = ""; if ($kEnt -ge 0) { $ent = $r[$kEnt] }
                if ($numero -eq "?" -or $numero -eq ".") { $numero = "" } else { $ent = "" }
                $chave = "$cadeia|$numero|$ins|$ent"
                if ($chave -ceq $anterior) { continue }
                $anterior = $chave
                $residuos.Add(@($cadeia, "$numero$ins", $r[$kNome].ToUpper()))
            }
        }
    }
}

if ($residuos.Count -eq 0) {
    Write-Host "Nenhum atomo encontrado (formato lido: $fmt)."; exit 1
}

# ---------------------------------------------------------------- contagem
$cont = [ordered]@{}
foreach ($aa in $AMINOACIDOS) { $cont[$aa] = 0 }
$conv = [ordered]@{}
$resumo = [ordered]@{}
$naoId = New-Object System.Collections.Generic.List[object]
foreach ($r in $residuos) {
    $nome = $r[2]
    $p = Padronizar $nome
    if ($p -ne "") {
        $cont[$p]++
        if ($p -cne $nome) {
            $k = "$nome->$p"
            if ($conv.Contains($k)) { $conv[$k]++ } else { $conv[$k] = 1 }
        }
    }
    elseif ($AGUAS -cnotcontains $nome) {
        $naoId.Add($r)
        if ($resumo.Contains($nome)) { $resumo[$nome]++ } else { $resumo[$nome] = 1 }
    }
}

$out = New-Object System.Collections.Generic.List[string]
$out.Add("Arquivo: $nomeArq"); $out.Add("Formato: $fmt"); $out.Add("")
$out.Add("Residuo  Quantidade"); $out.Add("-------------------")
$total = 0
foreach ($aa in $AMINOACIDOS) { $out.Add(("{0,-8}{1,11}" -f $aa, $cont[$aa])); $total += $cont[$aa] }
$out.Add("-------------------")
$out.Add(("{0,-8}{1,11}" -f "TOTAL", $total))
$out.Add("")
if ($conv.Count -gt 0) {
    $out.Add("Nomes convertidos para o aminoacido padrao (ja incluidos acima):")
    $out.Add("  " + (($conv.Keys | ForEach-Object { "$_ x$($conv[$_])" }) -join ", "))
    $out.Add("")
}
if ($naoId.Count -gt 0) {
    $out.Add("Residuos nao identificados: $($naoId.Count) (fora da contagem acima)")
    $out.Add("  Resumo: " + (($resumo.Keys | ForEach-Object { "$_ x$($resumo[$_])" }) -join ", "))
    $out.Add("  Nome  Cadeia  Numero")
    foreach ($r in $naoId) {
        $cad = $r[0]; if ($cad -eq "") { $cad = "-" }
        $out.Add(("  {0,-5} {1,-7} {2}" -f $r[2], $cad, $r[1]))
    }
}
else {
    $out.Add("Residuos nao identificados: nenhum")
}

$relatorio = ($out -join "`n") + "`n"
Write-Host -NoNewline $relatorio
if ($Saida -eq "") {
    $pasta = [System.IO.Path]::GetDirectoryName($caminho)
    $Saida = Join-Path $pasta (([System.IO.Path]::GetFileNameWithoutExtension($base)) + "_contagem.txt")
}
try {
    [System.IO.File]::WriteAllText($Saida, $relatorio, (New-Object System.Text.UTF8Encoding($false)))
    Write-Host ""
    Write-Host "Saida gravada em: $Saida"
}
catch {
    Write-Host ""
    Write-Host "Nao foi possivel gravar ${Saida}: $($_.Exception.Message)"
}
