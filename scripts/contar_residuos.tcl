# contar_residuos.tcl - Conta os residuos da molecula ja aberta no VMD.
#
# Uso: com a molecula carregada, digite no Tk Console (Extensions > Tk Console):
#
#     source /caminho/completo/contar_residuos.tcl
#
# A contagem da molecula "top" aparece na hora. Depois disso o comando
# contar_residuos fica disponivel para repetir:
#
#     contar_residuos                          ;# molecula top
#     contar_residuos 0                        ;# molecula de ID 0
#     contar_residuos top "chain A"            ;# so uma parte (selecao do VMD)
#     contar_residuos top all contagem.txt     ;# grava tambem em arquivo
#
# Saida: quantidade de cada um dos 20 aminoacidos, total, nomes convertidos e
# lista dos residuos nao identificados (qualquer nome fora dos 20). Aguas sao
# ignoradas.
#
# Formatos: qualquer um que o VMD abra (PDB, PQR, mmCIF/PDBx, GRO, PSF+DCD,
# MOL2, XYZ com nomes de residuo...). Quem separa os residuos e o proprio VMD
# (atributo "residue"), entao trajetorias e varios quadros nao mudam a contagem.
#
# Funciona tambem com estruturas de simulacao (GROMACS, AMBER, CHARMM/NAMD):
# nomes de 4 letras (HISE, CYS2...), variantes de protonacao (HID, HSD,
# CYX...) e terminais (NALA, CLYS...) contam como o aminoacido padrao, e as
# conversoes feitas aparecem na saida.

# Devolve o aminoacido padrao correspondente ao nome, ou "" se nao houver
proc padronizar_residuo {nome} {
    set padrao [list ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL]
    set sinonimos [list \
        HID HIS HIE HIS HIP HIS HIN HIS HSD HIS HSE HIS HSP HIS \
        HISD HIS HISE HIS HISH HIS HISP HIS HISA HIS HISB HIS HIS1 HIS HIS2 HIS \
        CYX CYS CYM CYS CYN CYS CYS2 CYS CYSH CYS \
        ASH ASP ASPH ASP ASPP ASP GLH GLU GLUH GLU GLUP GLU \
        LYN LYS LYSH LYS LYSN LYS LSN LYS ARGN ARG]
    array set sin $sinonimos
    if {[lsearch -exact $padrao $nome] >= 0} {
        return $nome
    }
    if {[info exists sin($nome)]} {
        return $sin($nome)
    }
    # terminais do AMBER/GROMACS: NALA, CLYS, NHIE...
    set primeira [string index $nome 0]
    if {[string length $nome] == 4 && ($primeira eq "N" || $primeira eq "C")} {
        set resto [string range $nome 1 end]
        if {[lsearch -exact $padrao $resto] >= 0} {
            return $resto
        }
        if {[info exists sin($resto)]} {
            return $sin($resto)
        }
    }
    return ""
}

proc contar_residuos {{molid top} {selecao all} {saida ""}} {
    set aminoacidos [list ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL]

    if {[molinfo num] == 0} {
        puts "Nenhuma molecula carregada no VMD."
        return
    }
    if {$molid eq "top"} {
        set molid [molinfo top]
    }

    foreach aa $aminoacidos {
        set cont($aa) 0
    }

    # aguas: macro "water" do VMD e os nomes usados por campos de forca
    set aguas "HOH WAT H2O DOD D2O SOL TIP TIP2 TIP3 TIP4 TIP5 TP3 TP4 T3P T4P T5P SPC SPCE SWM4 OPC OPC3"
    set sel [atomselect $molid "($selecao) and not water and not resname $aguas"]
    set natomos [$sel num]
    if {$natomos == 0} {
        $sel delete
        puts "Nenhum atomo encontrado na selecao."
        return
    }
    set lres   [$sel get residue]
    set lnome  [$sel get resname]
    set lcad   [$sel get chain]
    set lnum   [$sel get resid]
    set lins   [$sel get insertion]
    $sel delete

    # Cada residuo conta uma unica vez (indice interno "residue" do VMD)
    set desconhecidos [list]
    set ordem [list]
    set ordemconv [list]
    for {set i 0} {$i < $natomos} {incr i} {
        set res [lindex $lres $i]
        if {[info exists visto($res)]} {
            continue
        }
        set visto($res) 1
        set nome [string toupper [string trim [lindex $lnome $i]]]
        set aa [padronizar_residuo $nome]
        if {$aa ne ""} {
            incr cont($aa)
            if {$aa ne $nome} {
                set c "$nome->$aa"
                if {![info exists conv($c)]} {
                    set conv($c) 0
                    lappend ordemconv $c
                }
                incr conv($c)
            }
        } else {
            set cad [string trim [lindex $lcad $i]]
            if {$cad eq ""} {
                set cad "-"
            }
            set num [lindex $lnum $i]
            append num [string trim [lindex $lins $i]]
            lappend desconhecidos [list $nome $cad $num]
            if {![info exists resumo($nome)]} {
                set resumo($nome) 0
                lappend ordem $nome
            }
            incr resumo($nome)
        }
    }

    set linhas [list]
    lappend linhas "Molecula: [molinfo $molid get name]"
    lappend linhas ""
    lappend linhas "Residuo  Quantidade"
    lappend linhas "-------------------"
    set total 0
    foreach aa $aminoacidos {
        lappend linhas [format "%-8s%11d" $aa $cont($aa)]
        incr total $cont($aa)
    }
    lappend linhas "-------------------"
    lappend linhas [format "%-8s%11d" TOTAL $total]
    lappend linhas ""
    if {[llength $ordemconv] > 0} {
        lappend linhas "Nomes convertidos para o aminoacido padrao (ja incluidos acima):"
        set partes [list]
        foreach c $ordemconv {
            lappend partes "$c x$conv($c)"
        }
        lappend linhas "  [join $partes ", "]"
        lappend linhas ""
    }
    set nd [llength $desconhecidos]
    if {$nd > 0} {
        lappend linhas "Residuos nao identificados: $nd (fora da contagem acima)"
        set partes [list]
        foreach n $ordem {
            lappend partes "$n x$resumo($n)"
        }
        lappend linhas "  Resumo: [join $partes ", "]"
        lappend linhas "  Nome  Cadeia  Numero"
        foreach d $desconhecidos {
            lappend linhas [format "  %-5s %-7s %s" [lindex $d 0] [lindex $d 1] [lindex $d 2]]
        }
    } else {
        lappend linhas "Residuos nao identificados: nenhum"
    }

    set texto [join $linhas "\n"]
    puts $texto
    if {$saida ne ""} {
        set f [open $saida w]
        puts $f $texto
        close $f
        puts ""
        puts "Saida gravada em: [file normalize $saida]"
    }
    return
}

# Ao dar "source", ja conta a molecula top
contar_residuos
