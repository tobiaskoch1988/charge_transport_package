#!/bin/bash
# V3 mit Reiterationen der Geometrieoptimierung, wenn eine versagt, wird eine "benachbare" Geometrie verwendet also z.B. statt der neutrale Geometrie, wird aus der optimierten geladenen Geometrie (el,lo) begonnen um wieder eine neutrale Geometrie zu erzeugen.
# run_lambda_in.sh  monomerA monomerB foldername lambda_in_results.dat oniom_ok=[oniom/false]
# Ausfuehrung im Basisordner
 if [ $# -le 3 ]; then
    echo "Usage: ${0} monomerA monomerB foldername lambda_in_results.dat oniom_ok=[oniom/false] use_reiterations=[true/false] sub_g09=[sub_g09/false] "
    echo "Execute ${0} in foldername/${0}"
    exit 1
 fi

monomerA=${1}
monomerB=${2}
foldername=${3}
zieldatei=${4}
oniom_ok='false'
use_reiterations='false'
sub_g09='false' ## calculate g09 here!
CURDIR=$( pwd)

if [ $# -ge 5 ]; then # submit g09
    oniom_ok=${5}
    if [[ "${oniom_ok}" == "oniom" ]] ; then    	
        oniom_ok="true"
        echo "Oniom-Rechnung ausgewählt in ${0} "
    else
        oniom_ok="false"
    fi
fi	

if [ $# -ge 6 ]; then # use reiterations
    use_reiterations=${6}
    if [[ "${use_reiterations}" == "true" ]] ; then    	
        use_reiterations="true"
        iteration_max=8     ## Anzahl der Iterationen fuer Resubmitierung
        echo " Reiterationen ausgewaehlt ${iteration_max} "
    else
        use_reiterations="false"
        iteration_max=0
    fi
fi	# use reiterations


if [ $# -ge 7 ]; then # sub_g09
    sub_g09=${7}
    if [[ "${sub_g09}" == "sub_g09" ]] ; then    	
        sub_g09="true"
    else
        sub_g09='false'
    fi
fi	# sub_g09


if [[ ${foldername} == *'UFF'* ]] ; then #UFF_Abfrage
        extra_SP='true'
fi
    
 extra_SP='false'
 Au_to_eV=$(echo "scale=10; 27.211396132" | bc -l)
 
echo start ${0}  ${monomerA} ${monomerB} ${foldername} ${zieldatei} "oniom_ok=${oniom_ok}" " sub_g09=${sub_g09} " "extra_SP=${extra_SP}"


# function zur Multiplikation zweier Zahlen in Bash ueber ein python skript
function mult {
        if [ $# -lt 2 ] ; then
                echo "> Argument missing. Expected: 2"
                exit 1
        fi

        if [ -f script_mult.py ] ; then
                rm script_mult.py
        fi

        echo "#! /usr/bin/python"                  >> script_mult.py
        echo "import sys"                          >> script_mult.py
        echo "in1 = float(sys.argv[1])"            >> script_mult.py
        echo "in2 = float(sys.argv[2])"            >> script_mult.py
        echo "print \"{0:15.8f}\".format(in1*in2)" >> script_mult.py

        chmod a+x script_mult.py

        RESULT=$(./script_mult.py ${1} ${2})
        echo $RESULT

        rm script_mult.py
}



## Funktion um Optimierte Geometrien in SCF-Rechnungen zu setzen
## Opt_check_to_SCF  Opt  SCF   mol  ${foldername} ${zieldatei}
## Opt_check_to_SCF ${R1} ${R9} molA ${foldername} ${zieldatei}
##  Number to read    1     2     3       4              5 
##  Opt_check_to_SCF ${R1} ${R7} molA ${foldername} ${zieldatei}
function Opt_check_to_SCF(){
 Optfilename=${1}  # Optimierungsrechnungsname
 SCFfilename=${2}
 molA=${3}
 charge=${4}
 foldername=${5}
 zieldatei=${6}
 oniom_ok=${7}
 sub_g09=${8}
 
 XYZ=""
 Freeze_list=""
if [ -e ${molA}/${Optfilename}/${Optfilename}'.log' ]; then
    cd  ${molA}/${Optfilename}/
    normalR1=$(tail -1 ${Optfilename}.log | cut -c 1-19)        
    if [[ ${normalR1} == *"Normal termination"* ]]; then
        OptR1=$(grep "Optimized Parameters" ${Optfilename}.log | cut -c 32-40)
        OptMM=$(grep "Stationary point found" ${Optfilename}.log | awk '{print $2}')
        Max_Force=$(grep "Maximum Force" ${Optfilename}.log | cut -c 53-56 | tail -1)
        RMS_Force=$(grep "RMS     Force" ${Optfilename}.log | cut -c 53-56 | tail -1)
        Max_Displace=$(grep "Maximum Displacement" ${Optfilename}.log | cut -c 53-56 | tail -1)
        RMS_Displace=$(grep "RMS     Displacement" ${Optfilename}.log | cut -c 53-56 | tail -1)
        if [ "${OptR1}" == "Optimized" -o "${OptMM}" == "Stationary" ] && [ "${Max_Force}" == "YES" ] && [ "${Max_Displace}" == "YES" ]; then
            echo " Optimierung in ${Optfilename}.log OK"  
            if [ -e ${Optfilename}'.chk' ]; then
                        # get optimized coordinates from ${Optfilename}.xyz file
                        # input .ichk output -ixyz ###Ausgabe als R1. xyz  File
                        
                        
                        if [ ! -s ${Optfilename}'_conv.xyz' ]; then
                            newzmat -ichk -oxyz -step 99999 ${Optfilename}.chk ${Optfilename}.xyz
                        	Lines=$(wc -l ${Optfilename}.xyz) ; Lines=${Lines%% *};
                            echo ${Lines} >> ${Optfilename}_conv.xyz
                            echo " " >> ${Optfilename}_conv.xyz
                        	
                            XYZ=$(head -${Lines} ${Optfilename}.xyz)
                        	echo "${XYZ}" >> ${Optfilename}_conv.xyz
                        	rm ${Optfilename}.xyz
                        	echo " Datei ${Optfilename}_conv.xyz erstellt"
                        else
                            Lines=$(wc -l ${Optfilename}_conv.xyz )
                            Lines=$((${Lines%% *}-2))
                            XYZ=$(tail -${Lines} ${Optfilename}_conv.xyz)
                        fi # Erstellung _conv.xyz
                        # Einsetzen der optimierten Koordinaten in SCF-file und SCF-Berechnung 
                 
                        if [ -e ../${SCFfilename} ] ; then
                            cd ../${SCFfilename} 
                        else ## Erstellen neuen SCFfilename Ordner- und *.inp file !
                            mkdir ../${SCFfilename}
                            head -9 ../${SCFfilename%_V*}/${SCFfilename%_V*}.inp >> ../${SCFfilename}/${SCFfilename}.inp
                            cd ../${SCFfilename} 
                            echo " ${SCFfilename} erstellt" 
                        fi # Existiert Ordner ${SCFfilename} 
                        cp ../${Optfilename}/${Optfilename}'_conv.xyz' .
                     

                        ###echo ${SCFfilename} 
                        ###echo ${Optfilename}
                        ###tail -${Lines} ${Optfilename}_conv.xyz >> ${SCFfilename}.inp
                        ###echo " " >> ${SCFfilename}.inp
                        ##echo "oniom_ok=${oniom_ok}"
                        ##if [[ "${oniom_ok}" == "true" ]] ; then 
                        ##    if [ -e ../../Freeze_list_${molA}_${charge}.dat ]; then
                        ##        Freeze_list_A=$(tail -${Lines} ../../Freeze_list_${molA}_${charge}.dat)
                        ##        Freeze_list_A=$(echo ${Freeze_list_A})
                        ##        OIFS="$IFS"
                        ##        IFS=' '
                        ##        read -a Freeze_list <<< "${Freeze_list_A}"
                        ##        IFS="$OIFS"
                        ##        
                        ##    else
                        ##        echo "Fehler: oniom_ok=${oniom_ok}_but_file_../Freeze_list_${molA}.dat_is_not_available"  >> ../../${zieldatei}
                        ##    fi
                        ##fi
                        ##echo ${Freeze_list[@]}
                        ##echo ${Freeze_list[4]} ${Freeze_list[5]}
                        ##
                        ##XYZ=$(echo ${XYZ})
                        ##OIFS="$IFS"
                        ##IFS=' '
                        ##read -a XYZ_a <<< "${XYZ}"
                        ##IFS="$OIFS"
                        ##
                        ##i=0;j=1;k=2;l=3;m=1;f=0
                        ##echo ${SCFfilename}.inp
                        ##while [ ${m} -le ${Lines} ]
                        ##do  
                        ##    if [[ "${oniom_ok}" == "true" ]] ; then
                        ##        echo "${XYZ_a[${i}]} ${Freeze_list[${f}]} ${XYZ_a[${j}]} ${XYZ_a[${k}]} ${XYZ_a[${l}]} ${Freeze_list[$((${f}+1))]}" >> "${SCFfilename}.inp"
                        ##        echo "${XYZ_a[${i}]} ${Freeze_list[${f}]} ${XYZ_a[${j}]} ${XYZ_a[${k}]} ${XYZ_a[${l}]} ${Freeze_list[$((${f}+1))]}"
                        ##    else
                        ##        echo "${XYZ_a[${i}]} ${XYZ_a[${j}]} ${XYZ_a[${k}]} ${XYZ_a[${l}]}" >> "${SCFfilename}.inp"
                        ##    fi
                        ##    
                        ##    m=$((${m}+1))
                        ##    i=$((${i}+4))
                        ##    j=$((${j}+4))
                        ##    k=$((${k}+4))
                        ##    l=$((${l}+4))
                        ##    f=$((${f}+2))
                        ##done
                        if [[ "${oniom_ok}" == "true" ]] ; then 
                            Freeze_list_A="Freeze_list_${molA}_${charge}.dat" 
                            if [ -e ../../${Freeze_list_A} ]; then
                                tail -n +3 ${Optfilename}'_conv.xyz' > tmpA.xyz
                                paste tmpA.xyz ../../${Freeze_list_A} | awk '{ print " "$5" "$6"  "$2" "$3" "$4"  "$7" "$8" "$9}' >> ${SCFfilename}.inp ; rm tmpA.xyz 
                            else
                                echo "Fehler: oniom_ok=${oniom_ok}_but_file_../Freeze_list_${molA}.dat_is_not_available"  >> ../../${zieldatei}
                            fi
                        else
                            XYZ=$(echo ${XYZ})
                            OIFS="$IFS"
                            IFS=' '
                            read -a XYZ_a <<< "${XYZ}"
                            IFS="$OIFS"
                            
                            i=0;j=1;k=2;l=3;m=1;f=0
                            while [ ${m} -le ${Lines} ]
                            do  
                                echo "${XYZ_a[${i}]} ${XYZ_a[${j}]} ${XYZ_a[${k}]} ${XYZ_a[${l}]}" >> "${SCFfilename}.inp"
                                
                                m=$((${m}+1))
                                i=$((${i}+4))
                                j=$((${j}+4))
                                k=$((${k}+4))
                                l=$((${l}+4))
                                f=$((${f}+2))
                            done
                        fi
                        

                        echo " "  >> ${SCFfilename}.inp
                        if [ "${sub_g09}" == "true" ]; then
                            g09 ${SCFfilename}.inp	   
                        fi # sub_g09                        
            else
                echo "Error: Die Datei ${Optfilename}.chk existiert nicht " 
                echo "Fehler: Die_Datei_${Optfilename}.chk_existiert_nicht ${foldername} "  >> ../../${zieldatei}
                echo " EXIT" ; exit 1
            fi
        else
            echo " No Optimized Geometry in ${foldername}/${molA}/${Optfilename}/${Optfilename}.log "
            echo " Optimized Geometry  : ${OptR1}     "
            echo " Maximum Force       : ${Max_Force} "
            echo " RMS     Force       : ${RMS_Force} " 
            echo " Maximum Displacement: ${Max_Displace}"
            echo " RMS     Displacement: ${RMS_Displace}"
            echo " Fehler: No_Optimized_Geometry_in_${foldername}/${molA}/${Optfilename}/${Optfilename}.log ${foldername} "  >> ../../${zieldatei}
            echo " EXIT" ; exit 1
        fi
    else
        echo "Error: No Normal termination in ${foldername}/${molA}/${Optfilename}/${Optfilename}.log " 
        echo "Fehler: No_Normal_termination_in_${foldername}/${molA}/${Optfilename}/${Optfilename}.log ${foldername}"  >> ../../${zieldatei}
        echo " EXIT" ; exit 1
    fi #Normal termination
    cd ../../ # Zurück in Ausgangsordner
else
    echo "Error: file does not exist: ${foldername}/${molA}/${Optfilename}/${Optfilename}.log " 
    echo "Fehler: file does not exist: ${foldername}/${molA}/${Optfilename}/${Optfilename}.log ${foldername}" >> ${zieldatei}
    echo " EXIT" ; exit 1
fi
   return   
  } # Endefunction Opt_check_to_SCF  

  
    #Create file structure for Optimization
    R1=OptS0_n_${monomerA}
    R2=OptS0_el_${monomerA}
    R3=OptS0_lo_${monomerA}
    R4=OptS0_n_${monomerB}
    R5=OptS0_el_${monomerB}
    R6=OptS0_lo_${monomerB}

    #Create file structure for SCF-Calculations molA
    R7=SCF_el_OptS0_n_${monomerA}
    R8=SCF_lo_OptS0_n_${monomerA}
    R9=SCF_n_OptS0_el_${monomerA}
    R10=SCF_n_OptS0_lo_${monomerA}
    # SCF file structure for molB
    R11=SCF_el_OptS0_n_${monomerB}
    R12=SCF_lo_OptS0_n_${monomerB}
    R13=SCF_n_OptS0_el_${monomerB}
    R14=SCF_n_OptS0_lo_${monomerB}
 
    if [ "${extra_SP}" == 'true' ]; then
        R15=SCF_n_OptS0_n_${monomerA}
        R16=SCF_el_OptS0_el_${monomerA}
        R17=SCF_lo_OptS0_lo_${monomerA}
    
        R18=SCF_n_OptS0_n_${monomerB}
        R19=SCF_el_OptS0_el_${monomerB}
        R20=SCF_lo_OptS0_lo_${monomerB}
    fi
 #################
 # Berechne molA #
 #################
 cd molA/${R1}
 if [ -e "${R1}.inp" ]; then
        if [ ! -e "${R1}.log" ]; then
            if [ "${sub_g09}" == "true" ]; then
                g09 ${R1}.inp
            fi # sub_g09
        else
            echo " ${R1}.log existiert bereits"
        fi # *.log existiert nicht
 else
    echo "Fehler: Die Datei ${R1}.inp existiert nicht. ENDE"
    exit 1
 fi # existier input
 cd ../../
 
 Opt_check_to_SCF ${R1} ${R7} molA el ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
 Opt_check_to_SCF ${R1} ${R8} molA lo ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
 
 cd molA/${R2}
 if [ -e "${R2}.inp" ]; then
        if [ ! -e "${R2}.log" ]; then
            if [ "${sub_g09}" == "true" ]; then
                g09 ${R2}.inp
            fi # sub_g09
        else
            echo " ${R2}.log existiert bereits"
        fi # *.log existiert nicht
 else
    echo "Fehler: Die Datei ${R2}.inp existiert nicht. ENDE"
    exit 1
 fi # existier input
 cd ../../
 
 Opt_check_to_SCF ${R2} ${R9} molA n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
 
 cd molA/${R3}
 if [ -e "${R3}.inp" ]; then
        if [ ! -e "${R3}.log" ]; then
            if [ "${sub_g09}" == "true" ]; then
                g09 ${R3}.inp
            fi # sub_g09
        else
            echo " ${R3}.log existiert bereits"
        fi # *.log existiert nicht
 else
    echo "Fehler: Die Datei ${R3}.inp existiert nicht. ENDE"
    exit 1
 fi # existier input
 cd ../../
 
 
 Opt_check_to_SCF ${R3} ${R10} molA n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
 
  if [ "${extra_SP}" == 'true' ]; then
        Opt_check_to_SCF ${R1} ${R15} molA n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
        Opt_check_to_SCF ${R2} ${R16} molA el ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
        Opt_check_to_SCF ${R3} ${R17} molA lo ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
  fi
################### 
#  Berechne molB  #
###################
cd molB/${R4}
 if [ -e "${R4}.inp" ]; then
        if [ ! -e "${R4}.log" ]; then
            if [ "${sub_g09}" == "true" ]; then
                g09 ${R4}.inp
            fi # sub_g09
        else
            echo " ${R4}.log existiert bereits"
        fi # *.log existiert nicht
 else
    echo "Fehler: Die Datei ${R4}.inp existiert nicht. ENDE"
    exit 1
 fi # existier input
cd ../../

Opt_check_to_SCF ${R4} ${R11} molB el ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
Opt_check_to_SCF ${R4} ${R12} molB lo ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}

cd molB/${R5}
 if [ -e "${R5}.inp" ]; then
        if [ ! -e "${R5}.log" ]; then
            if [ "${sub_g09}" == "true" ]; then
                g09 ${R5}.inp
            fi # sub_g09
        else
            echo " ${R5}.log existiert bereits"
        fi # *.log existiert nicht
 else
    echo "Fehler: Die Datei ${R5}.inp existiert nicht. ENDE"
    exit 1
 fi # existier input
cd ../../

Opt_check_to_SCF ${R5} ${R13} molB n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}

cd molB/${R6}
 if [ -e "${R6}.inp" ]; then
        if [ ! -e "${R6}.log" ]; then
            if [ "${sub_g09}" == "true" ]; then
                g09 ${R6}.inp
            fi # sub_g09
        else
            echo " ${R6}.log existiert bereits"
        fi # *.log existiert nicht
 else
    echo "Fehler: Die Datei ${R6}.inp existiert nicht. ENDE"
    exit 1
 fi # existier input
cd ../../


Opt_check_to_SCF ${R6} ${R14} molB n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
  if [ "${extra_SP}" == 'true' ]; then
        Opt_check_to_SCF ${R4} ${R18} molB n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
        Opt_check_to_SCF ${R5} ${R19} molB el ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
        Opt_check_to_SCF ${R6} ${R20} molB lo ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
  fi

cd ../ 

if [ "${extra_SP}" == 'true' ]; then
  # lambda Auswertung mit erneuter SCF-Berechnung nach der Optimierungsrechnung
        echo "*** Starte Auswertung ***" 
        counter=0
        molA='molA'
        j=0
        for SCF in ${R7} ${R8} ${R9} ${R10} ${R15} ${R16} ${R17} ${R11} ${R12} ${R13} ${R14} ${R18} ${R19} ${R20}
        do
            normal="A"
            let j=${j}+1
            if [ ${j} -ge 8 ] ; then
                molA='molB'
            fi
            if [ -e ${foldername}/${molA}/${SCF}/${SCF}.log ]; then
                            normal=$( grep "Normal termination" ${foldername}/${molA}/${SCF}/${SCF}.log | cut -c 1-19)
                            if [[ ${normal} == *"Normal termination"* ]]; then
                                let counter=${counter}+1
                            else
                                echo "Fehler: Rechnung ${SCF}.log nicht mit Normal termination beendet"
                                echo "Fehler: Rechnung_${SCF}.log_nicht_mit_Normal_termination_beendet ${foldername}" >> ${foldername}/${zieldatei}
                                echo "Die Auswertung wird beendet"
                                exit  1
                            fi
            else
                echo "Die Datei ${SCF}.log existiert nicht"
                echo "Fehler: Die_Datei_${SCF}.log_existiert_nicht ${foldername}" >> ${foldername}/${zieldatei}
                echo "Die Auswertung wird beendet"
                #exit 1   
            fi
        done
        # Beginne Auswertung
        j=0
        molA='molA'
                     
        if [[ ${counter} == 14 ]]; then
                echo "Alle SCF-Rechnungen normal beendet."
                for SCF in ${R15} ${R16} ${R17} ${R18} ${R19} ${R20} ${R7} ${R8} ${R9} ${R10} ${R11} ${R12} ${R13} ${R14} ${R1} ${R2} ${R3} ${R4} ${R5} ${R6}
                do
                    let j=${j}+1
                    SCF_Done=$( grep "SCF Done" ${foldername}/${molA}/${SCF}/${SCF}.log | tail -1 )
                    SCF_Done=${SCF_Done##*=} ; SCF_Done=${SCF_Done%%A*} # Vorne und hinten einkuerzen
                    if [[ "${oniom_ok}" == "true" ]] ; then
                        SCF_Done=$(grep "ONIOM: extrapolated energy =" ${foldername}/${molA}/${SCF}/${SCF}.log | tail -1) 
                        SCF_Done=${SCF_Done##*=} # Vorne kürzen,ONIOM: extrapolated energy
                    elif [[ ${foldername} == *'UFF'* ]] ; then
                        SCF_Done=$( grep "Energy=" ${foldername}/${molA}/${SCF}/${SCF}.log | head -1)
                        SCF_Done=${SCF_Done%NIter*} ; SCF_Done=${SCF_Done##*=} #Vorne und hinten einkuerzen
                    fi
                    Au_to_eV=$(echo "scale=10; 27.211396132" | bc -l)
                    SCF_Energie[${j}]=$(mult ${SCF_Done} ${Au_to_eV})
                    echo  " ${SCF} ${SCF_Energie[${j}]} eV ${SCF_Done} AU " ## ${Au_to_eV}
                    
                    if [ ${j} == 3 ]; then
                            molA='molB'
                    elif [ ${j} == 6 ]; then
                            molA='molA'
                    elif [ ${j} == 10 ]; then
                            molA='molB'
                    elif [ ${j} == 14 ]; then
                            molA='molA'
                    elif [ ${j} == 17 ]; then
                            molA='molB'
                    fi
                done
                lambda_int_el_A_B=$( echo "scale=10; ${SCF_Energie[9]}  - ${SCF_Energie[1]} + ${SCF_Energie[11]} - ${SCF_Energie[5]} " | bc -l)
                lambda_int_el_B_A=$( echo "scale=10; ${SCF_Energie[13]} - ${SCF_Energie[4]} + ${SCF_Energie[7]}  - ${SCF_Energie[2]} " | bc -l)

                lambda_int_lo_A_B=$( echo "scale=10; ${SCF_Energie[10]} - ${SCF_Energie[1]} + ${SCF_Energie[12]} - ${SCF_Energie[6]} " | bc -l)
                lambda_int_lo_B_A=$( echo "scale=10; ${SCF_Energie[14]} - ${SCF_Energie[4]} + ${SCF_Energie[8]}  - ${SCF_Energie[3]} " | bc -l)

                dE_el_AB=$( echo "scale=10; ${SCF_Energie[2]} - ${SCF_Energie[1]} - ${SCF_Energie[5]}  + ${SCF_Energie[4]} " | bc -l)
                dE_lo_AB=$( echo "scale=10; ${SCF_Energie[3]} - ${SCF_Energie[1]} - ${SCF_Energie[6]}  + ${SCF_Energie[4]} " | bc -l)
                # Schreibe daten in file
                # Test ob Daten bereits in der Datei stehen
                echo "SCF Energien in eV."

                echo "1"  ${R15}  ${SCF_Energie[1]}
                echo "2"  ${R16}  ${SCF_Energie[2]}
                echo "3"  ${R17}  ${SCF_Energie[3]}
                echo "4"  ${R18}  ${SCF_Energie[4]}
                echo "5"  ${R19}  ${SCF_Energie[5]}
                echo "6"  ${R20}  ${SCF_Energie[6]}
                echo "7"  ${R7}  ${SCF_Energie[7]}
                echo "8"  ${R8}  ${SCF_Energie[8]}
                echo "9"  ${R9}  ${SCF_Energie[9]}
                echo "10" ${R10} ${SCF_Energie[10]}
                echo "11" ${R11} ${SCF_Energie[11]}
                echo "12" ${R12} ${SCF_Energie[12]}
                echo "13" ${R13} ${SCF_Energie[13]}
                echo "14" ${R14} ${SCF_Energie[14]}
                echo " --------------------------  "
                echo "15" ${R1} ${SCF_Energie[15]}
                echo "16" ${R2} ${SCF_Energie[16]}
                echo "17" ${R3} ${SCF_Energie[17]}
                echo "18" ${R4} ${SCF_Energie[18]}
                echo "19" ${R5} ${SCF_Energie[19]}
                echo "20" ${R6} ${SCF_Energie[20]}
                
                echo "---------------------------------"
                echo "lambda_int_el_A_B = $lambda_int_el_A_B "
                echo "lambda_int_el_B_A = $lambda_int_el_B_A "          
                echo "dE_el_AB = $dE_el_AB"
                echo "lambda_int_lo_A_B = $lambda_int_lo_A_B "
                echo "lambda_int_lo_B_A = $lambda_int_lo_B_A "
                echo "dE_lo_AB = $dE_lo_AB"
		if [ "${use_reiterations}" != "true" ]; then 
                	if [ -e ${foldername}/${zieldatei} ]; then
                        	str=$( grep "${foldername}" ${zieldatei} )
                        	if [ -z "$str" ]; then
                                	echo  ${lambda_int_el_A_B} ${lambda_int_el_B_A} ${dE_el_AB} ${lambda_int_lo_A_B} ${lambda_int_lo_B_A} ${dE_lo_AB} ${SCF_Energie[1]}	${SCF_Energie[2]}	${SCF_Energie[3]}	${SCF_Energie[7]}	${SCF_Energie[8]}	${SCF_Energie[9]}	${SCF_Energie[10]}	${SCF_Energie[4]}	${SCF_Energie[5]}	${SCF_Energie[6]}	${SCF_Energie[11]}	${SCF_Energie[12]}	${SCF_Energie[13]}	${SCF_Energie[14]}  ${SCF_Energie[15]} ${SCF_Energie[16]} ${SCF_Energie[17]} ${SCF_Energie[18]} ${SCF_Energie[19]} ${SCF_Energie[20]} ${foldername} >> ${foldername}/${zieldatei}
                                	echo "Daten in file ${zieldatei} geschrieben."
                                	#exit 1
                        	else
                                	echo "Daten in file ${zieldatei} sind bereits vorhanden." 
                                	#exit 1
                        	fi # data already in file ?
                	else
                    			echo  ${lambda_int_el_A_B} ${lambda_int_el_B_A} ${dE_el_AB} ${lambda_int_lo_A_B} ${lambda_int_lo_B_A} ${dE_lo_AB} ${SCF_Energie[1]}	${SCF_Energie[2]}	${SCF_Energie[3]}	${SCF_Energie[7]}	${SCF_Energie[8]}	${SCF_Energie[9]}	${SCF_Energie[10]}	${SCF_Energie[4]}	${SCF_Energie[5]}	${SCF_Energie[6]}	${SCF_Energie[11]}	${SCF_Energie[12]}	${SCF_Energie[13]}	${SCF_Energie[14]} ${SCF_Energie[15]} ${SCF_Energie[16]} ${SCF_Energie[17]} ${SCF_Energie[18]} ${SCF_Energie[19]} ${SCF_Energie[20]} ${foldername} >> ${foldername}/${zieldatei}
                                	echo "Datei ${zieldatei} neu erstellt. "
                                	if [ "${use_reiterations}" != "true" ]; then
                                	    exit 1
                        	        fi #
                	fi # zieldatei existiert - Print daten to file
		fi # use_reiterations = true ?
        fi # Ende Auswertung




else # if "${extra_SP}" == 'false'
      # lambda Auswertung ohne erneute SCF-Berechnung
        echo "*** Starte Auswertung ***" 
        counter=0
        molA='molA'
        j=0
        for SCF in ${R7} ${R8} ${R9} ${R10} ${R11} ${R12} ${R13} ${R14}
        do
           normal="A"
           let j=${j}+1
           if [ ${j} -ge 5 ] ; then
                molA='molB'
           fi
           if [ -e ${foldername}/${molA}/${SCF}/${SCF}.log ]; then
                        normal=$( grep "Normal termination" ${foldername}/${molA}/${SCF}/${SCF}.log | cut -c 1-19)
                        if [[ ${normal} == *"Normal termination"* ]]; then
                            let counter=${counter}+1
                        else
                            echo "Fehler: Rechnung ${SCF}.log nicht mit Normal termination beendet"
                            echo "Fehler: Rechnung_${SCF}.log_nicht_mit_Normal_termination_beendet ${foldername}" >> ${foldername}/${zieldatei}
                            echo "Die Auswertung wird beendet"
                            exit  1
                        fi
           else
              echo "Die Datei ${SCF}.log existiert nicht"
              echo "Fehler: Die_Datei_${SCF}.log_existiert_nicht ${foldername}" >> ${foldername}/${zieldatei}
              echo "Die Auswertung wird beendet"
              #exit 1   
           fi
        done
        # Beginne Auswertung
        j=0
        molA='molA'
                     
        if [[ ${counter} == 8 ]]; then
                 echo "Alle SCF-Rechnungen normal beendet."
                 for SCF in ${R1} ${R2} ${R3} ${R4} ${R5} ${R6} ${R7} ${R8} ${R9} ${R10} ${R11} ${R12} ${R13} ${R14}
                 do
                        let j=${j}+1
                        SCF_Done=$( grep "SCF Done" ${foldername}/${molA}/${SCF}/${SCF}.log | tail -1 )
                        SCF_Done=${SCF_Done##*=} ; SCF_Done=${SCF_Done%%A*} # Vorne und hinten einkuerzen
                        if [[ "${oniom_ok}" == "true" ]] ; then
                            SCF_Done=$(grep "ONIOM: extrapolated energy =" ${foldername}/${molA}/${SCF}/${SCF}.log | tail -1) 
                            SCF_Done=${SCF_Done##*=} # Vorne kürzen,ONIOM: extrapolated energy
                        elif [[ ${foldername} == *'UFF'* ]]; then
                            SCF_Done=$( grep "Energy=" ${foldername}/${molA}/${SCF}/${SCF}.log | head -1)
                            SCF_Done=${SCF_Done%NIter*} ; SCF_Done=${SCF_Done##*=} #Vorne und hinten einkuerzen
                        fi
                        Au_to_eV=$(echo "scale=10; 27.211396132" | bc -l)
                        SCF_Energie[${j}]=$(mult ${SCF_Done} ${Au_to_eV})
                        echo  "${SCF} ${SCF_Energie[${j}]} eV ${SCF_Done} Au"  ## ${Au_to_eV}
                        
                        if [ ${j} == 3 ]; then
                                molA='molB'
                        elif [ ${j} == 6 ]; then
                                molA='molA'
                        elif [ ${j} == 10 ]; then
                                molA='molB'
                        fi
                done
                lambda_int_el_A_B=$( echo "scale=10; ${SCF_Energie[9]}  - ${SCF_Energie[1]} + ${SCF_Energie[11]} - ${SCF_Energie[5]} " | bc -l)
                lambda_int_el_B_A=$( echo "scale=10; ${SCF_Energie[13]} - ${SCF_Energie[4]} + ${SCF_Energie[7]}  - ${SCF_Energie[2]} " | bc -l)

                lambda_int_lo_A_B=$( echo "scale=10; ${SCF_Energie[10]} - ${SCF_Energie[1]} + ${SCF_Energie[12]} - ${SCF_Energie[6]} " | bc -l)
                lambda_int_lo_B_A=$( echo "scale=10; ${SCF_Energie[14]} - ${SCF_Energie[4]} + ${SCF_Energie[8]}  - ${SCF_Energie[3]} " | bc -l)

                dE_el_AB=$( echo "scale=10; ${SCF_Energie[2]} - ${SCF_Energie[1]} - ${SCF_Energie[5]}  + ${SCF_Energie[4]} " | bc -l)
                dE_lo_AB=$( echo "scale=10; ${SCF_Energie[3]} - ${SCF_Energie[1]} - ${SCF_Energie[6]}  + ${SCF_Energie[4]} " | bc -l)
                # Schreibe daten in file
                # Test ob Daten bereits in der Datei stehen
                echo "SCF Energien in eV."

                echo "1"  ${R1}  ${SCF_Energie[1]}
                echo "2"  ${R2}  ${SCF_Energie[2]}
                echo "3"  ${R3}  ${SCF_Energie[3]}
                echo "4"  ${R4}  ${SCF_Energie[4]}
                echo "5"  ${R5}  ${SCF_Energie[5]}
                echo "6"  ${R6}  ${SCF_Energie[6]}
                echo "7"  ${R7}  ${SCF_Energie[7]}
                echo "8"  ${R8}  ${SCF_Energie[8]}
                echo "9"  ${R9}  ${SCF_Energie[9]}
                echo "10" ${R10} ${SCF_Energie[10]}
                echo "11" ${R11} ${SCF_Energie[11]}
                echo "12" ${R12} ${SCF_Energie[12]}
                echo "13" ${R13} ${SCF_Energie[13]}
                echo "14" ${R14} ${SCF_Energie[14]}
                echo "---------------------------------"
                echo "lambda_int_el_A_B = $lambda_int_el_A_B "
                echo "lambda_int_el_B_A = $lambda_int_el_B_A "          
                echo "dE_el_AB = $dE_el_AB"
                echo "lambda_int_lo_A_B = $lambda_int_lo_A_B "
                echo "lambda_int_lo_B_A = $lambda_int_lo_B_A "
                echo "dE_lo_AB = $dE_lo_AB"

                if [ -e ${foldername}/${zieldatei} ]; then
                        str=$( grep "${foldername}" ${zieldatei} )
                        if [ -z "$str" ]; then
                                echo  ${lambda_int_el_A_B} ${lambda_int_el_B_A} ${dE_el_AB} ${lambda_int_lo_A_B} ${lambda_int_lo_B_A} ${dE_lo_AB} ${SCF_Energie[1]}	${SCF_Energie[2]}	${SCF_Energie[3]}	${SCF_Energie[7]}	${SCF_Energie[8]}	${SCF_Energie[9]}	${SCF_Energie[10]}	${SCF_Energie[4]}	${SCF_Energie[5]}	${SCF_Energie[6]}	${SCF_Energie[11]}	${SCF_Energie[12]}	${SCF_Energie[13]}	${SCF_Energie[14]} ${foldername} >> ${foldername}/${zieldatei}
                                echo "Daten in file ${zieldatei} geschrieben."
                                #exit 1
                        else
                                echo  ${lambda_int_el_A_B} ${lambda_int_el_B_A} ${dE_el_AB} ${lambda_int_lo_A_B} ${lambda_int_lo_B_A} ${dE_lo_AB} ${SCF_Energie[1]}     ${SCF_Energie[2]}       ${SCF_Energie[3]}       ${SCF_Energie[7]}       ${SCF_Energie[8]}       ${SCF_Energie[9]}       ${SCF_Energie[10]}      ${SCF_Energie[4]}       ${SCF_Energie[5]}       ${SCF_Energie[6]}       ${SCF_Energie[11]}      ${SCF_Energie[12]}      ${SCF_Energie[13]}      ${SCF_Energie[14]} ${foldername} >> ${foldername}/${zieldatei}

				echo "Daten in file ${zieldatei} sind bereits vorhanden." 
                                #exit 1
                        fi
                else
                    echo  ${lambda_int_el_A_B} ${lambda_int_el_B_A} ${dE_el_AB} ${lambda_int_lo_A_B} ${lambda_int_lo_B_A} ${dE_lo_AB} ${SCF_Energie[1]}	${SCF_Energie[2]}	${SCF_Energie[3]}	${SCF_Energie[7]}	${SCF_Energie[8]}	${SCF_Energie[9]}	${SCF_Energie[10]}	${SCF_Energie[4]}	${SCF_Energie[5]}	${SCF_Energie[6]}	${SCF_Energie[11]}	${SCF_Energie[12]}	${SCF_Energie[13]}	${SCF_Energie[14]} ${foldername} >> ${foldername}/${zieldatei}
                                echo "Datei ${zieldatei} neu erstellt. "
                                #exit 1
                fi
        fi # Ende Auswertung

fi # Ende if "${extra_SP}" == 'true'



### Resubmit if negativ lambda occures.
#
if [ "${use_reiterations}" != "true" ]; then
    echo " Use of Reiterations not selected: ${use_reiterations} ; ENDE"
    exit 1
fi # exit if not use_reiterations
echo "Start resubmit program"
######### RRSUBBMIT PROGRAM!


#lambda_int_el_A_B=12.0
#lambda_int_el_B_A=-2.0
#lambda_int_lo_A_B=3.0
#lambda_int_lo_B_A=4.0
null=0.0

counter=0 # counts the iteration in the loop to select SCF-Calculation, which need to be corrected

CURDIR=$( pwd )

function get_geometry_and_calc_opt {
SCF_Geo=${1}  # SCF-file nach der Optimierung die energetisch unterhalb der Energie von Old_Geo-Energie liegt.
Old_Geo=${2}  # Geometrie der alten Optimierung möglich R1-R6
iteration_V=${3}
molX=${4}
Optfilename=${5} # Optimierungsrechnung vor der SCF-Rechnung SCF_${charge}_..._${Optfilename}
foldername=${6}
CURDIR=${7}
charge=${8} # Ladung der Optimierung fuer das neue file / entspricht Ladung von Old_Geo
oniom_ok=${9}
extra_SP=${10}
zieldatei=${11}
monomer=${12}
sub_g09=${13}

New_Opt_Geo=${Old_Geo}_${iteration_V}

echo ${SCF_Geo}
echo ${New_Opt_Geo}
echo ${molX}
echo ${CURDIR}

## Auswahl um die Geometrie aus dem anderen File zu holen, falls es zweimal nicht konvergiert ist.
## Es sollen so, im Wechsel jeweils die anderen beiden Ausgangsgeometrien verwendet werden. Bsp OptS0_n -> el,lo
Number_calculated_optimizations=$( ls -d ${CURDIR}/${foldername}/${molX}/${New_Opt_Geo%%_V*}* | wc -l )
SCF_Geo_charge=${SCF_Geo#*OptS0_} ; SCF_Geo_charge=${SCF_Geo_charge%%_*} ## Schneide Ladung der alten Geometrie aus

if ((  $(echo "  ${Number_calculated_optimizations}  > 2 " |bc -l) )) ; then
    for oldcharge in n el lo ; do
        if [ "${oldcharge}" != "${charge}" ] && [ "${oldcharge}" != "${SCF_Geo_charge}" ] ; then
            vorne=${Old_Geo%_${charge}_*} ; hinten=${Old_Geo#*_${charge}_} 
            echo "$oldcharge"
            if [ $(($Number_calculated_optimizations % 2)) -eq 0 ] ; then  ## hier gerade   ## Wechsel nur bei ungerader Anzahl an Durchlaeufen
                Old_Geo=$( ls -d ${CURDIR}/${foldername}/${molX}/${vorne}_${SCF_Geo_charge}_${hinten%%_V*}* | tail -1 ) 
                Old_Geo=${Old_Geo##*/}
            else # ungerade -> wechsel 
                Old_Geo=$( ls -d ${CURDIR}/${foldername}/${molX}/${vorne}_${oldcharge}_${hinten%%_V*}* | tail -1 ) 
                Old_Geo=${Old_Geo##*/}
            fi # check even, uneven calculation
        fi # check charges
    done # Schleife ueber Ladungstypen
fi # adapt selected Old_Geo

echo " use: ${Old_Geo} "
echo " Number_calculated_optimizations: ${Number_calculated_optimizations} "
echo " charge: ${charge}"
echo " SCF_Geo_charge: ${SCF_Geo_charge} "


    if [ -e ${CURDIR}/${foldername}/${molX}/${New_Opt_Geo} ]; then
        echo " Der Ordner ${New_Opt_Geo} existiert bereits"
        #exit 1
    else	
        mkdir ${CURDIR}/${foldername}/${molX}/${New_Opt_Geo}
    fi # existiert New_Opt_Geo bereits ?
    
   if [ -e  ${CURDIR}/${foldername}/${molX}/${SCF_Geo} ]; then
    if [ -e ${CURDIR}/${foldername}/${molX}/${SCF_Geo}/${Optfilename}'_conv.xyz' ]; then
        cp ${CURDIR}/${foldername}/${molX}/${SCF_Geo}/${Optfilename}'_conv.xyz' ${CURDIR}/${foldername}/${molX}/${New_Opt_Geo}
        cd ${CURDIR}/${foldername}/${molX}/${New_Opt_Geo}
        head -9 ${CURDIR}/${foldername}/${molX}/${Old_Geo}/${Old_Geo}.inp > ${New_Opt_Geo}.inp # Erstelle neuen Head
        sed -i s/${Old_Geo}/${New_Opt_Geo}/g ${New_Opt_Geo}.inp # Namen im file durch neuen Namen ersetzen
        ### Add new coords to file
        
        Lines=$(wc -l ${Optfilename}_conv.xyz )
        Lines=$((${Lines%% *}-2))
        XYZ=$(tail -${Lines} ${Optfilename}_conv.xyz)
        echo 'start geo'
        
        if [[ "${oniom_ok}" == "true" ]] ; then 
            Freeze_list_A="Freeze_list_${molX}_${charge}.dat" 
            echo "use: ${Freeze_list_A}"
            if [ -e ../../${Freeze_list_A} ]; then
                tail -n +3 ${Optfilename}'_conv.xyz' > tmpA.xyz
                paste tmpA.xyz ../../${Freeze_list_A} | awk '{ print " "$5" "$6"  "$2" "$3" "$4"  "$7" "$8" "$9}' >> ${New_Opt_Geo}.inp ; rm tmpA.xyz 
            else
                echo "Fehler: oniom_ok=${oniom_ok}_but_file_../Freeze_list_${molA}.dat_is_not_available"  >> ../../${zieldatei}
            fi
        else # ohne Oniom
            XYZ=$(echo ${XYZ})
            OIFS="$IFS"
            IFS=' '
            read -a XYZ_a <<< "${XYZ}"
            IFS="$OIFS"
            
            i=0;j=1;k=2;l=3;m=1;f=0
            while [ ${m} -le ${Lines} ]
            do  
                echo "${XYZ_a[${i}]} ${XYZ_a[${j}]} ${XYZ_a[${k}]} ${XYZ_a[${l}]}" >> "${New_Opt_Geo}.inp"
                
                m=$((${m}+1))
                i=$((${i}+4))
                j=$((${j}+4))
                k=$((${k}+4))
                l=$((${l}+4))
                f=$((${f}+2))
            done
        fi
        
        
        echo " "  >> ${New_Opt_Geo}.inp 
        
        if [ "${sub_g09}" == "true" ]; then
            echo " g09 ${New_Opt_Geo}.inp calculate once more"
            g09 ${New_Opt_Geo}.inp	# submit g09
        fi # sub_g09
        
       
                
        cd ${CURDIR} # Wechsel zurück 
        
    else  # *.conv.xyz exists 
        echo " Error ${CURDIR}/${foldername}/${molX}/${SCF_Geo}/${Optfilename}_conv.xyz does not exist "
        echo ' Ende' ; exit 1
    fi # *.conv.xyz exists 
    else 
        echo " Error  folder does ${CURDIR}/${foldername}/${molX}/${SCF_Geo}/ does not exist "
        echo ' Ende' ; exit 1
   fi # folder for ${SCF_Geo} exists

    
 } # end function  get_geometry_and_calc_opt


#### Function holt Energie Daten aus Datei bei _V{iteration_A} -index ##
function get_SCF_Energy_from_calculation {
## SCF_Done_A = get_SCF_Energy_from_calculation R7 iteration molA ${foldername} ${oniom_ok}

SCF_A=${1}  # SCF-file nach der Optimierung die energetisch unterhalb der Energie von Old_Geo-Energie liegt.
iteration_A=${2}
molX=${3}
foldername=${4}
oniom_ok=${5}
if [[ ${iteration_A} == 1 ]] ; then
    SCF_A=${SCF_A%_V*}
else
    SCF_A=${SCF_A%_V*}_V${iteration_A}
fi
#echo "${SCF_A}"
    if [[ "${oniom_ok}" == "true" ]] ; then
        SCF_Done_A=$(grep "ONIOM: extrapolated energy =" ${foldername}/${molX}/${SCF_A}/${SCF_A}.log | tail -1) 
        SCF_Done_A=${SCF_Done_A##*=} # Vorne kürzen,ONIOM: extrapolated energy
    elif [[ ${foldername} == *'UFF'* ]]; then
        SCF_Done_A=$( grep "Energy=" ${foldername}/${molX}/${SCF_A}/${SCF_A}.log | head -1)
        SCF_Done_A=${SCF_Done_A%NIter*} ; SCF_Done_A=${SCF_Done_A##*=} #Vorne und hinten einkuerzen
    else # default SCF-Fall
        SCF_Done_A=$( grep "SCF Done" ${foldername}/${molX}/${SCF_A}/${SCF_A}.log | tail -1 )
        SCF_Done_A=${SCF_Done_A##*=} ; SCF_Done_A=${SCF_Done_A%%A*} # Vorne und hinten einkuerzen
    fi
    #echo " SCF_Done_A: ${SCF_Done_A}"
    echo ${SCF_Done_A}
    return
} #### Function holt Energie Daten aus Datei ##




# Iterationen zur neuberechnung
for iteration in `seq 2 ${iteration_max}`
do # Iterationen der Neuberechnung 
    counter_positiv_lambda=0 # Zaehlt die Anzahl der positiven lambdas, bei 4 positiven lambdas ist alles ok!
    counter=0
    molA_n_done='false' # Reset, pro Iteration soll nur einmal neutral optimiert werden!
    molB_n_done='false' # Reset
    for lambda in ${lambda_int_el_A_B} ${lambda_int_el_B_A} ${lambda_int_lo_A_B} ${lambda_int_lo_B_A} # check negativ lambda 
    do
        let counter=${counter}+1
        echo "Counter lambda: ${counter}"
        echo "lambda ${lambda}"
        if (( $(echo " ${lambda} < ${null} " |bc -l) )); then # lambda negativ ???
            echo " lambda ${lambda} < 0.0 eV  start resubmit."
            ########################################### Start  Block um die Berechnung der negativen lambda zu pruefen und neue Optimierungen mit anschliessenden SCF-Rechnungen durchzufuerhen ###
            if [[ ${counter} == 1 ]]; then      
                if (( $(echo "${SCF_Energie[9]} < ${SCF_Energie[1]}" |bc -l) )) && [ "${molA_n_done}" != "True" ]; then
                echo "  get_geometry_and_calc_opt ${R9} ${R1%_V*} V${iteration} molA ${R2} ${foldername} ${CURDIR} n ${oniom_ok}"
                    get_geometry_and_calc_opt ${R9} ${R1%_V*} V${iteration}  molA ${R2} ${foldername} ${CURDIR} n ${oniom_ok} ${extra_SP}  ${zieldatei}  ${monomerA} ${sub_g09}
                    if [ -e  ${CURDIR}/${foldername}/'molB'/${R1%_V*}_V${iteration}/${R1%_V*}_V${iteration}.log ]; then
                        R1=${R1%_V*}_V${iteration}
                        cd ${CURDIR}/${foldername}
                        Opt_check_to_SCF ${R1} ${R7%_V*}_V${iteration} molA el ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                        Opt_check_to_SCF ${R1} ${R8%_V*}_V${iteration} molA lo ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}  
                        R7=${R7%_V*}_V${iteration}
                        R8=${R8%_V*}_V${iteration}
                        if [ "${extra_SP}" == 'true' ]; then
                            Opt_check_to_SCF ${R1%_V*}_V${iteration} ${R15%_V*}_V${iteration} molA n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                            R15=${R15%_V*}_V${iteration}
                        fi #extra_SP
                        cd ${CURDIR}/
                        molA_n_done='True'
                    fi # *.log nach Iterationsrechnung existiert    
                        
                    echo " recalculated Data in ${R1} "
                fi
                if (( $(echo "${SCF_Energie[11]} < ${SCF_Energie[5]}" |bc -l) )); then
                    echo " get_geometry_and_calc_opt ${R11} ${R5%_V*} V${iteration}  molB ${R4} ${foldername} ${CURDIR} el ${oniom_ok} "
                    get_geometry_and_calc_opt ${R11} ${R5%_V*} V${iteration}  molB ${R4} ${foldername} ${CURDIR} el ${oniom_ok} ${extra_SP}  ${zieldatei}  ${monomerB} ${sub_g09}
                    if [ -e ${CURDIR}/${foldername}/'molB'/${R5%_V*}_V${iteration}/${R5%_V*}_V${iteration}.log ]; then
                        R5=${R5%_V*}_V${iteration}
                        cd ${CURDIR}/${foldername}
                        Opt_check_to_SCF ${R5} ${R13%_V*}_V${iteration} molB n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                        R13=${R13%_V*}_V${iteration}
                        if [ "${extra_SP}" == 'true' ]; then
                            Opt_check_to_SCF ${R5%_V*}_V${iteration}  ${R19%_V*}_V${iteration} molB el ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                            R19=${R19%_V*}_V${iteration}
                        fi #exttra_SP
                        cd ${CURDIR}/
                    fi # *.log nach Iterationsrechnung existiert
                    echo " recalculated Data in ${R5} "
                fi # energy difference <
            ########################################
            elif [[ ${counter} == 2 ]]; then  
                if (( $(echo "${SCF_Energie[13]} < ${SCF_Energie[4]}" |bc -l) )) && [ "${molB_n_done}" != "True" ] ; then
                    echo "  get_geometry_and_calc_opt ${R13} ${R4%_V*} V${iteration} molB ${R5} ${foldername} ${CURDIR} n ${oniom_ok}"
                    get_geometry_and_calc_opt ${R13} ${R4%_V*} V${iteration}  molB ${R5} ${foldername} ${CURDIR} n ${oniom_ok} ${extra_SP}  ${zieldatei}  ${monomerB} ${sub_g09}
                    
                    if [ -e ${CURDIR}/${foldername}/'molB'/${R4%_V*}_V${iteration}/${R4%_V*}_V${iteration}.log ]; then
                        R4=${R4%_V*}_V${iteration}
                        cd ${CURDIR}/${foldername}
                        Opt_check_to_SCF ${R4} ${R11%_V*}_V${iteration} molB el ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                        Opt_check_to_SCF ${R4} ${R12%_V*}_V${iteration} molB lo ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}  
                        R11=${R11%_V*}_V${iteration}
                        R12=${R12%_V*}_V${iteration}
                        if [ "${extra_SP}" == 'true' ]; then
                            Opt_check_to_SCF ${R4%_V*}_V${iteration} ${R18%_V*}_V${iteration} molB n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                            R18=${R18%_V*}_V${iteration}
                        fi #extra_SP
                        cd ${CURDIR}/
                        molB_n_done='True'
                    fi # *.log nach Iterationsrechnung existiert    
                        
                    echo " recalculated Data in ${R4} "
                fi            
            
            
                if (( $(echo "${SCF_Energie[7]} < ${SCF_Energie[2]}" |bc -l) )); then
                    echo " get_geometry_and_calc_opt ${R7} ${R2%_V*} V${iteration}  molA ${R1} ${foldername} ${CURDIR} el ${oniom_ok}"
                    get_geometry_and_calc_opt ${R7} ${R2%_V*} V${iteration}  molA ${R1} ${foldername} ${CURDIR} el ${oniom_ok} ${extra_SP}  ${zieldatei}  ${monomerA} ${sub_g09}
                    if [ -e ${CURDIR}/${foldername}/'molA'/${R2%_V*}_V${iteration}/${R2%_V*}_V${iteration}.log ]; then
                        R2=${R2%_V*}_V${iteration}
                        cd ${CURDIR}/${foldername}
                        Opt_check_to_SCF ${R2} ${R9%_V*}_V${iteration} molA n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                        R9=${R9%_V*}_V${iteration}
                        if [ "${extra_SP}" == 'true' ]; then
                            Opt_check_to_SCF ${R2%_V*}_V${iteration}  ${R16%_V*}_V${iteration} molA el ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                            R16=${R16%_V*}_V${iteration}
                        fi # extra_SP
                        cd ${CURDIR}/
                    fi # *.log nach Iterationsrechnung existiert
                    echo " recalculated Data in ${R2} "
                fi # energy difference <  
                
            elif [[ ${counter} == 3 ]]; then
               
                if (( $(echo "${SCF_Energie[10]} < ${SCF_Energie[1]}" |bc -l) )) && [ "${molA_n_done}" != "True" ]; then
                    echo "  get_geometry_and_calc_opt ${R10} ${R1%_V*} V${iteration} molA ${R3} ${foldername} ${CURDIR} n ${oniom_ok}"
                    get_geometry_and_calc_opt ${R10} ${R1%_V*} V${iteration}  molA ${R3} ${foldername} ${CURDIR} n ${oniom_ok} ${extra_SP}  ${zieldatei}  ${monomerA} ${sub_g09}
                    echo "hier:" $(pwd)
                    echo ${CURDIR}
                    
                    if [ -e ${CURDIR}/${foldername}/'molA'/${R1%_V*}_V${iteration}/${R1%_V*}_V${iteration}.log ]; then
                        R1=${R1%_V*}_V${iteration}
                        cd ${CURDIR}/${foldername}
                        Opt_check_to_SCF ${R1} ${R7%_V*}_V${iteration} molA el ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                        Opt_check_to_SCF ${R1} ${R8%_V*}_V${iteration} molA lo ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}  
                        R7=${R7%_V*}_V${iteration}
                        R8=${R8%_V*}_V${iteration}
                        if [ "${extra_SP}" == 'true' ]; then
                            Opt_check_to_SCF ${R1%_V*}_V${iteration} ${R15%_V*}_V${iteration} molA n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                            R15=${R15%_V*}_V${iteration}                           
                            echo  " Opt_check_to_SCF ${R1%_V*}_V${iteration} ${R15%_V*}_V${iteration} molA n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}"
                        fi
                        cd ${CURDIR}/
                        molA_n_done='True'
                    fi # *.log nach Iterationsrechnung existiert         
                    echo " recalculated Data in ${R1} "
                fi             

                if (( $(echo "${SCF_Energie[12]} < ${SCF_Energie[6]}" |bc -l) )); then
                    echo "get_geometry_and_calc_opt ${R12} ${R6%_V*} V${iteration}  molB ${R4} ${foldername} ${CURDIR} lo ${oniom_ok} "
                    get_geometry_and_calc_opt ${R12} ${R6%_V*} V${iteration}  molB ${R4} ${foldername} ${CURDIR} lo ${oniom_ok} ${extra_SP}  ${zieldatei}  ${monomerB} ${sub_g09}
                    if [ -e ${CURDIR}/${foldername}/'molB'/${R6%_V*}_V${iteration}/${R6%_V*}_V${iteration}.log ]; then
                        R6=${R6%_V*}_V${iteration}
                        cd ${CURDIR}/${foldername}
                        Opt_check_to_SCF ${R6} ${R14%_V*}_V${iteration} molB n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                        R14=${R14%_V*}_V${iteration}
                        if [ "${extra_SP}" == 'true' ]; then
                            Opt_check_to_SCF ${R6%_V*}_V${iteration}  ${R20%_V*}_V${iteration} molB lo ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                            R20=${R20%_V*}_V${iteration}
                        fi #extra_SP
                        cd ${CURDIR}/
                    fi # *.log nach Iterationsrechnung existiert
                    echo " recalculated Data in ${R6} "
                fi # energy difference < 
            

            elif [[ ${counter} == 4 ]]; then 
                if (( $(echo "${SCF_Energie[14]} < ${SCF_Energie[4]}" |bc -l) )) && [ "${molB_n_done}" != "True" ]; then
                    echo "  get_geometry_and_calc_opt ${R14} ${R4%_V*} V${iteration} molB ${R6} ${foldername} ${CURDIR} n ${oniom_ok}"
                    get_geometry_and_calc_opt ${R14} ${R4%_V*} V${iteration}  molB ${R6} ${foldername} ${CURDIR} n ${oniom_ok} ${extra_SP}  ${zieldatei}  ${monomerB} ${sub_g09}
                    if [ -e ${CURDIR}/${foldername}/'molB'/${R4%_V*}_V${iteration}/${R4%_V*}_V${iteration}.log ]; then
                        R4=${R4%_V*}_V${iteration}
                        cd ${CURDIR}/${foldername}
                        Opt_check_to_SCF ${R4} ${R11%_V*}_V${iteration} molB el ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                        Opt_check_to_SCF ${R4} ${R12%_V*}_V${iteration} molB lo ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}  
                        R11=${R11%_V*}_V${iteration}
                        R12=${R12%_V*}_V${iteration}
                        if [ "${extra_SP}" == 'true' ]; then
                            Opt_check_to_SCF ${R4%_V*}_V${iteration} ${R18%_V*}_V${iteration} molB n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                            R18=${R18%_V*}_V${iteration}
                        fi #extra_SP
                        cd ${CURDIR}/
                        molB_n_done='True'
                    fi # *.log nach Iterationsrechnung existiert    
                        
                    echo " recalculated Data in ${R4} "
                fi              
                
                if (( $(echo "${SCF_Energie[7]} < ${SCF_Energie[3]}" |bc -l) )); then
                    echo "get_geometry_and_calc_opt ${R7} ${R3%_V*} V${iteration}  molA ${R1} ${foldername} ${CURDIR} lo ${oniom_ok} "
                    get_geometry_and_calc_opt ${R7} ${R3%_V*} V${iteration}  molA ${R1} ${foldername} ${CURDIR} lo ${oniom_ok} ${extra_SP}  ${zieldatei}  ${monomerA} ${sub_g09}
                    if [ -e ${CURDIR}/${foldername}/'molA'/${R3%_V*}_V${iteration}/${R3%_V*}_V${iteration}.log ]; then
                        R3=${R3%_V*}_V${iteration}
                        cd ${CURDIR}/${foldername}
                        Opt_check_to_SCF ${R3} ${R10%_V*}_V${iteration} molA n ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                        R10=${R10%_V*}_V${iteration}
                        if [ "${extra_SP}" == 'true' ]; then
                            Opt_check_to_SCF ${R3%_V*}_V${iteration}  ${R17%_V*}_V${iteration} molA lo ${foldername} ${zieldatei} ${oniom_ok} ${sub_g09}
                            R17=${R17%_V*}_V${iteration}
                        fi #extra_SP
                        cd ${CURDIR}/
                    fi # *.log nach Iterationsrechnung existiert
                    echo " recalculated Data in ${R3} "
                fi # energy difference <             


        ##if [[ ${counter} == 2 ]]; then 
        ##if (( $(echo "${SCF_Energie[13]} < ${SCF_Energie[4]}" |bc -l) )); then
        ##    get_geometry_and_calc_opt ${R13} ${R4}_V2 molB ${R5} ${foldername} ${CURDIR} n ${oniom_ok} ${sub_g09}
        ##fi
        ##    if (( $(echo "${SCF_Energie[7]} < ${SCF_Energie[2]}" |bc -l) )); then
        ##        get_geometry_and_calc_opt ${R7} ${R2}_V2 molA ${R1}  ${foldername} ${CURDIR} el ${oniom_ok}
        ##    fi   
        ##
        ##if [[ ${counter} == 3 ]]; then 
        ##    if (( $(echo "${SCF_Energie[10]} < ${SCF_Energie[1]}" |bc -l) )); then
        ##        get_geometry_and_calc_opt ${R10} ${R1}_V2 molA ${R3} ${foldername} ${CURDIR} n ${oniom_ok}
        ##    fi
        ##    if (( $(echo "${SCF_Energie[12]} < ${SCF_Energie[6]}" |bc -l) )); then
        ##        get_geometry_and_calc_opt ${R12} ${R6}_V2 molB ${R4} ${foldername} ${CURDIR} lo ${oniom_ok}
        ##    fi
        ##fi         
        ##
        ##if [[ ${counter} == 4 ]]; then     
        ##    if (( $(echo "${SCF_Energie[14]} < ${SCF_Energie[4]}" |bc -l) )); then
        ##        get_geometry_and_calc_opt ${R14} ${R4}_V2 molB ${R6} ${foldername} ${CURDIR} n ${oniom_ok}
        ##    fi
        ##    if (( $(echo "${SCF_Energie[8]} < ${SCF_Energie[3]}" |bc -l) )); then
        ##        get_geometry_and_calc_opt ${R8} ${R3}_V2 molA ${R1} ${foldername} ${CURDIR} lo ${oniom_ok}
        ##    fi
        ##fi      
            
            fi ## Anzahl der counter == 1 bis 4

 
            
        else # wenn lambda positiv ist
            echo " lambda ${lambda} > 0.0 eV ; ok " 
            let counter_positiv_lambda=${counter_positiv_lambda}+1 ## hochzaehlen
        fi # Energy < 0
    done #lambda
    
    
    if [[ ${counter_positiv_lambda} == 4 ]]; then # Wenn alle 4 lambdas positiv sind, schreibe alle Daten in die Zieldatei und Ende!
        if [ -e ${zieldatei} ]; then
                str=$( grep "${foldername}" ${zieldatei} )
                if [ -z "$str" ]; then
                        echo  ${lambda_int_el_A_B} ${lambda_int_el_B_A} ${dE_el_AB} ${lambda_int_lo_A_B} ${lambda_int_lo_B_A} ${dE_lo_AB} ${SCF_Energie[1]}	${SCF_Energie[2]}	${SCF_Energie[3]}	${SCF_Energie[7]}	${SCF_Energie[8]}	${SCF_Energie[9]}	${SCF_Energie[10]}	${SCF_Energie[4]}	${SCF_Energie[5]}	${SCF_Energie[6]}	${SCF_Energie[11]}	${SCF_Energie[12]}	${SCF_Energie[13]}	${SCF_Energie[14]} ${foldername} >> ${foldername}/${zieldatei}
                        echo "Daten in file ${zieldatei} geschrieben."
                        exit 1
                else
                        echo "Daten in file ${zieldatei} sind bereits vorhanden." 
                        exit 1
                fi
        else
            echo  ${lambda_int_el_A_B} ${lambda_int_el_B_A} ${dE_el_AB} ${lambda_int_lo_A_B} ${lambda_int_lo_B_A} ${dE_lo_AB} ${SCF_Energie[1]}	${SCF_Energie[2]}	${SCF_Energie[3]}	${SCF_Energie[7]}	${SCF_Energie[8]}	${SCF_Energie[9]}	${SCF_Energie[10]}	${SCF_Energie[4]}	${SCF_Energie[5]}	${SCF_Energie[6]}	${SCF_Energie[11]}	${SCF_Energie[12]}	${SCF_Energie[13]}	${SCF_Energie[14]} ${foldername} >> ${foldername}/${zieldatei}
                        echo "Datei ${zieldatei} neu erstellt. "
                        exit 1
        fi # existiert Zieldatei ?
    else # Alle lambda_in sind nicht positiv !
     
            echo "*** Starte Auswertung with Recalculation ***" 
 
            counter=0
            molA='molA'
            j=0
            
            if [ "${extra_SP}" == 'true' ]; then # Umordnung der Elemente bei SCF_extra
                scf_swap=${R1}  ;   R1=${R15} ; R15=${scf_swap}
                scf_swap=${R2}  ;   R2=${R16} ; R16=${scf_swap}
                scf_swap=${R3}  ;   R3=${R17} ; R17=${scf_swap}
                scf_swap=${R4}  ;   R4=${R18} ; R18=${scf_swap}
                scf_swap=${R5}  ;   R5=${R19} ; R19=${scf_swap}
                scf_swap=${R6}  ;   R6=${R20} ; R20=${scf_swap}
            fi # Umordnung der Elemente bei SCF_extra
            
            for SCF in ${R7} ${R8} ${R9} ${R10} ${R11} ${R12} ${R13} ${R14}
            do
            normal="A"
            let j=${j}+1
            if [ ${j} -ge 5 ] ; then
                    molA='molB'
            fi
            if [ -e ${foldername}/${molA}/${SCF}/${SCF}.log ]; then
                            normal=$( grep "Normal termination" ${foldername}/${molA}/${SCF}/${SCF}.log | cut -c 1-19)
                            if [[ ${normal} == *"Normal termination"* ]]; then
                                let counter=${counter}+1
                            else
                                echo "Fehler: Rechnung ${SCF}.log nicht mit Normal termination beendet"
                                echo "Fehler: Rechnung_${SCF}.log_nicht_mit_Normal_termination_beendet ${foldername}" >> ${foldername}/${zieldatei}
                                echo "Die Auswertung wird beendet"
                                exit  1
                            fi
            else
                echo "Die Datei ${SCF}.log existiert nicht"
                echo "Fehler: Die_Datei_${SCF}.log_existiert_nicht ${foldername}" >> ${foldername}/${zieldatei}
                echo "Die Auswertung wird beendet"
                exit 1   
            fi
            done
            # Beginne Auswertung
            j=0
            molA='molA'
                        
            if [[ ${counter} == 8 ]]; then
                    echo "Alle SCF-Rechnungen normal beendet."
                    for SCF in ${R1} ${R2} ${R3} ${R4} ${R5} ${R6} ### ${R7} ${R8} ${R9} ${R10} ${R11} ${R12} ${R13} ${R14}
                    do
                            let j=${j}+1
                            ###########SCF-Block ohne Minimalsuche sondern nur letztes Element ####
                            ###SCF_Done=$( grep "SCF Done" ${foldername}/${molA}/${SCF}/${SCF}.log | tail -1 )
                            ###SCF_Done=${SCF_Done##*=} ; SCF_Done=${SCF_Done%%A*} # Vorne und hinten einkuerzen
                            ###if [[ "${oniom_ok}" == "true" ]] ; then
                            ###    SCF_Done=$(grep "ONIOM: extrapolated energy =" ${foldername}/${molA}/${SCF}/${SCF}.log | tail -1) 
                            ###    SCF_Done=${SCF_Done##*=} # Vorne kürzen,ONIOM: extrapolated energy
                            ###elif echo "${foldername}" | grep 'UFF'; then
                            ###    SCF_Done=$( grep "Energy=" ${foldername}/${molA}/${SCF}/${SCF}.log | head -1)
                            ###    SCF_Done=${SCF_Done%NIter*} ; SCF_Done=${SCF_Done##*=} #Vorne und hinten einkuerzen
                            ###fi
                            ###########SCF-Block ohne Minimalsuche sondern nur letztes Element ####
                            
                            ########### Start Block zum herausziehen der minimalen Energie bei mehreren Iterationen ####
                            for i in ${!Energy_Arr[@]} ; do # leere Array ; iterate over array indexes
                                Energy_Arr[${i}]=''       # set the string as empty at this specific index
                            done
                            
                            if [[ "${oniom_ok}" == "true" ]] ; then ## Oniom
                                for Vx in ` seq 1 ${iteration}` ; do ## grep Array aller Energien
                                    if [[ ${Vx} == 1 ]] ; then
                                        Energy_Arr+=$(grep "ONIOM: extrapolated energy =" ${foldername}/${molA}/${SCF%_V*}/${SCF%_V*}.log | awk '{print $5}' |tail -1)               
                                    else
                                        Energy_Arr+=$(grep "ONIOM: extrapolated energy =" ${foldername}/${molA}/${SCF%_V*}_V${Vx}/${SCF%_V*}_V${Vx}.log | awk '{print $5}' |tail -1) 
                                    fi
                                echo "${Energy_Arr[${Vx}]}"
                            done # get all posible energies
                            elif [[ ${foldername} == *'UFF'* ]] ; then ## UFF 
                                Energy_Arr=$( grep "Energy= " ${foldername}/${molA}/${SCF%_V*}*/${SCF%_V*}*.log | awk '{print $3}' )
                                Energy_Arr=( ${Energy_Arr} ) ## Umwandlung String zu Array
                            else ## Standard SCF
                                for Vx in ` seq 1 ${iteration}` ; do
                                    if [[ ${Vx} == 1 ]] ; then
                                        Energy_Arr+=$(grep "SCF Done" ${foldername}/${molA}/${SCF%_V*}/${SCF%_V*}.log | awk '{print $5}' |tail -1)               
                                    else
                                        Energy_Arr+=$(grep "SCF Done" ${foldername}/${molA}/${SCF%_V*}_V${Vx}/${SCF%_V*}_V${Vx}.log | awk '{print $5}' |tail -1) 
                                    fi
                                done # get all posible energies
                            fi
                            
                            ##echo ${Energy_Arr[@]} 
                            
                            SCF_min=${Energy_Arr[0]}  ## Uebergabe erstes Element
                            N_count_elements=0
                            index_SCF_min=0
                            for Energy in ${Energy_Arr[@]}; do  ## Suche kleinestes Element SCF_min in Array Energy_Arr
                                    if (( $(echo "${Energy} <= ${SCF_min}" |bc -l) )); then
                                            SCF_min=${Energy}
                                            index_SCF_min=${N_count_elements}
                                    fi      #  kleinestes Element SCF_min in Array Energy_Arr
                                    ((N_count_elements++))  # Hochzaehlen
                            done # Durchlaufe Energien
                            echo " index_SCF_min: ${index_SCF_min}"
                            
                            ## Erhalte R_min(Dateinamen) zu der SCF_min den Dateinamen
                            if [[ "${oniom_ok}" == "true" ]] ; then ## Oniom 
                                if [[ ${index_SCF_min} == 0 ]] ; then # grep Dateinamen
                                        R_min=${SCF%_V*}
                                else # index_SCF_min nicht 0
                                        R_min=${SCF%_V*}_V${index_SCF_min}
                                fi # index_SCF_min nicht 0
                                R_min=${R_min##*/} # Schneide Pfadrest ab
                            elif [[ ${foldername} == *'UFF'* ]] ; then ## UFF 
                                Energy_Arr=$( grep "Energy= " ${foldername}/${molA}/${SCF%_V*}*/${SCF%_V*}*.log | awk '{print $1}' )
                                Energy_Arr=( ${Energy_Arr} ) ## Umwandlung String zu Array
                                R_min=${Energy_Arr[${index_SCF_min}]##*/}  ## Uebergabe des Namens mit dem kleinesten Energiewert!
                                echo " {Energy_Arr[{index_SCF_min}]##*/}:  ${Energy_Arr[${index_SCF_min}]##*/}  "
                                echo " {Energy_Arr[{index_SCF_min}]##*/}:  ${Energy_Arr[@]##*/}  "
                                R_min=${R_min%%.*} # Endung weg!
                            else ## Default SCF
                                if [[ ${index_SCF_min} == 0 ]] ; then # grep Dateinamen
                                    R_min=${SCF%_V*}
                                else # index_SCF_min nicht 0
                                    R_min=${SCF%_V*}_V${index_SCF_min}
                                fi # index_SCF_min nicht 0
                                R_min=${R_min##*/} # Schneide Pfadrest ab
                            fi # Erhalte SCF
                            
                            
                                echo "SCF_min: ${SCF_min} index: ${index_SCF_min} R_min: ${R_min}  N_counter: ${N_count_elements}" 
                            SCF_Energie[${j}]=$(mult ${SCF_min} ${Au_to_eV})
                            echo  " ${SCF} ${SCF_Energie[${j}]} eV ${SCF_Done} Au" ##${Au_to_eV}

                                SCF_Done=${SCF_min}
                                SCF=${R_min}
                                if [[ ${j} == 1 ]]; then ### Schreibweise ist so sehr schlecht; Array einführen!
                                    R1=${R_min}
                                    get_SCF_Energy_from_calculation ${R7} ${index_SCF_min} molA ${foldername} ${oniom_ok}
                                    SCF_Done_A=$( get_SCF_Energy_from_calculation ${R7} ${index_SCF_min} molA ${foldername} ${oniom_ok} )
                                    SCF_Energie[7]=$(mult ${SCF_Done_A} ${Au_to_eV})
                                    SCF_Done_A=$( get_SCF_Energy_from_calculation ${R8} ${index_SCF_min} molA ${foldername} ${oniom_ok} )
                                    SCF_Energie[8]=$(mult ${SCF_Done_A} ${Au_to_eV})
                                    echo ${SCF_Energie[7]} ${SCF_Energie[8]} 
                                elif [[ ${j} == 2 ]]; then 
                                    R2=${R_min}
                                    SCF_Done_A=$( get_SCF_Energy_from_calculation ${R9} ${index_SCF_min} molA ${foldername} ${oniom_ok} )
                                    SCF_Energie[9]=$(mult ${SCF_Done_A} ${Au_to_eV})
                                elif [[ ${j} == 3 ]]; then   
                                    R3=${R_min}
                                    SCF_Done_A=$( get_SCF_Energy_from_calculation ${R10} ${index_SCF_min} molA ${foldername} ${oniom_ok} )
                                    SCF_Energie[10]=$(mult ${SCF_Done_A} ${Au_to_eV})
                                elif [[ ${j} == 4 ]]; then   
                                    R4=${R_min} 
                                    SCF_Done_A=$( get_SCF_Energy_from_calculation ${R11} ${index_SCF_min} molB ${foldername} ${oniom_ok} )
                                    SCF_Energie[11]=$(mult ${SCF_Done_A} ${Au_to_eV})     
                                    SCF_Done_A=$( get_SCF_Energy_from_calculation ${R12} ${index_SCF_min} molB ${foldername} ${oniom_ok} )
                                    SCF_Energie[12]=$(mult ${SCF_Done_A} ${Au_to_eV})
                                elif [[ ${j} == 5 ]]; then   
                                    R5=${R_min}      
                                    SCF_Done_A=$( get_SCF_Energy_from_calculation ${R13} ${index_SCF_min} molB ${foldername} ${oniom_ok} )
                                    SCF_Energie[13]=$(mult ${SCF_Done_A} ${Au_to_eV})
                                elif [[ ${j} == 6 ]]; then   
                                    R6=${R_min}  
                                    SCF_Done_A=$( get_SCF_Energy_from_calculation ${R14} ${index_SCF_min} molB ${foldername} ${oniom_ok} )
                                    SCF_Energie[14]=$(mult ${SCF_Done_A} ${Au_to_eV})                                                   
                                fi # Array zuweisung der Namen zu den Rechnungen 
                                
                                echo " R_min: ${R_min} "
                                echo " SCF: ${SCF} "
                                echo " ${R1} ${R2} ${R3} ${R4} ${R5} ${R6} ${R7} ${R8} ${R9} ${R10} ${R11} ${R12} ${R13} ${R14} " 
                        ########### Ende Block zum herausziehen der minimalen Energie bei mehreren Iterationen ####                     
                            
                        
                            
                            if [[ ${j} == 3 ]]; then
                                    molA='molB'
                            elif [[ ${j} == 6 ]]; then
                                    molA='molA'
                            elif [[ ${j} == 10 ]]; then
                                    molA='molB'
                            fi
                    done
                    lambda_int_el_A_B=$( echo "scale=10; ${SCF_Energie[9]}  - ${SCF_Energie[1]} + ${SCF_Energie[11]} - ${SCF_Energie[5]} " | bc -l)
                    lambda_int_el_B_A=$( echo "scale=10; ${SCF_Energie[13]} - ${SCF_Energie[4]} + ${SCF_Energie[7]}  - ${SCF_Energie[2]} " | bc -l)
    
                    lambda_int_lo_A_B=$( echo "scale=10; ${SCF_Energie[10]} - ${SCF_Energie[1]} + ${SCF_Energie[12]} - ${SCF_Energie[6]} " | bc -l)
                    lambda_int_lo_B_A=$( echo "scale=10; ${SCF_Energie[14]} - ${SCF_Energie[4]} + ${SCF_Energie[8]}  - ${SCF_Energie[3]} " | bc -l)
    
                    dE_el_AB=$( echo "scale=10; ${SCF_Energie[2]} - ${SCF_Energie[1]} - ${SCF_Energie[5]}  + ${SCF_Energie[4]} " | bc -l)
                    dE_lo_AB=$( echo "scale=10; ${SCF_Energie[3]} - ${SCF_Energie[1]} - ${SCF_Energie[6]}  + ${SCF_Energie[4]} " | bc -l)
                    # Schreibe daten in file
                    # Test ob Daten bereits in der Datei stehen
                    echo "SCF Energien in eV."
    
                    echo "1"  ${R1}  ${SCF_Energie[1]}
                    echo "2"  ${R2}  ${SCF_Energie[2]}
                    echo "3"  ${R3}  ${SCF_Energie[3]}
                    echo "4"  ${R4}  ${SCF_Energie[4]}
                    echo "5"  ${R5}  ${SCF_Energie[5]}
                    echo "6"  ${R6}  ${SCF_Energie[6]}
                    echo "7"  ${R7}  ${SCF_Energie[7]}
                    echo "8"  ${R8}  ${SCF_Energie[8]}
                    echo "9"  ${R9}  ${SCF_Energie[9]}
                    echo "10" ${R10} ${SCF_Energie[10]}
                    echo "11" ${R11} ${SCF_Energie[11]}
                    echo "12" ${R12} ${SCF_Energie[12]}
                    echo "13" ${R13} ${SCF_Energie[13]}
                    echo "14" ${R14} ${SCF_Energie[14]}
                    echo "---------------------------------"
                    echo "lambda_int_el_A_B = $lambda_int_el_A_B "
                    echo "lambda_int_el_B_A = $lambda_int_el_B_A "          
                    echo "dE_el_AB = $dE_el_AB"
                    echo "lambda_int_lo_A_B = $lambda_int_lo_A_B "
                    echo "lambda_int_lo_B_A = $lambda_int_lo_B_A "
                    echo "dE_lo_AB = $dE_lo_AB"
            else 
                echo "Error: Nicht Alle SCF-Rechnungen normal beendet. Iteration: ${iteration} "
                echo "ENDE" ; exit 1
            fi # Alle SCF-Rechnungen normal beendet. ?  

            if [ "${extra_SP}" == 'true' ]; then # Umordnung der Elemente bei SCF_extra
                scf_swap=${R1}  ;   R1=${R15} ; R15=${scf_swap}
                scf_swap=${R2}  ;   R2=${R16} ; R16=${scf_swap}
                scf_swap=${R3}  ;   R3=${R17} ; R17=${scf_swap}
                scf_swap=${R4}  ;   R4=${R18} ; R18=${scf_swap}
                scf_swap=${R5}  ;   R5=${R19} ; R19=${scf_swap}
                scf_swap=${R6}  ;   R6=${R20} ; R20=${scf_swap}
            fi # Umordnung der Elemente bei SCF_extra


    fi # Ende Alle lambda_in sind nicht positiv !
    
 done # interationen 2 bis iteration_max - Anzahl der Iterationen der Neuberechnungen    
 
