#!/bin/bash
# bash-script creates g09-input-files and folder structure to calculate lambda for monomerA and monomerB for
# geometry optimization for a neutral_n ; negative_el ; positive_lo charged molecule.
# skript reades *.xyz-files, the calculation methode and the used foldername
# 
# $1 provides molA.xyz
# $2 provides molB.xyz
# $3 provides functional/basisset
# $4 provides the foldername

#dimer run; prepare orbitals from monomers
#Prepare Optimization and SCF-Input-files to calculate the reorganisation energy lambda.
#### Auslesen von Optimierten Geometrien aus einer *.chk file 
# input .ichk output -ixyz ###Ausgabe als xyz File
# newzmat -ichk -oxyz -step 999 OptS0_n_DIPBI_49.chk Dipbi_opt.xyz
#
if [ $# -le 7 ]; then
	echo "Usage: $0 < monomerA.xyz monomerB.xyz inputline foldername zieldatei run_lambda_in subg09_logical oniom  >"
	echo "Bsp. ${0} MolA.xyz MolB.xyz \"UFF  Opt(Tight,MaxMicroiterations=30000) \" MolA_MolB_B3LYP_6-31G+ results_lambda_in.dat run_lambda_in.sh true oniom"
	echo " logical=true fuer sofortiges Submittieren" 
	exit 1
fi

monomerA=`basename $1 .xyz`
monomerB=`basename $2 .xyz`
sub_ok='false'
oniom_ok='false'
extra_SP='false'
Opt_input="Opt(Tight)"  ##,MaxMicroiterations=50000)"


    if [ -n "$3" ]; then
        input=${3}
    else
        input="UFF"
        Opt_input="Opt(Tight,MaxMicroiterations=50000)"
    fi
	
  
    if [ -n "${4}" ]; then	
        foldername=${4}
    else   	
        foldername="ordner"
    fi		
	
    if [ $# -ge 5 ]; then	
        zieldatei=${5}
    else   	
        zieldatei='results_lambda_in.dat' #Ergebnisse der lambda Berechnung
    fi	


    if [ $# -ge 6 ]; then	
        run_lambda_in=${6}
    else   	
        run_lambda_in='run_lambda_in.sh' #Ergebnisse der lambda Berechnung
    fi
    if [ -e ${run_lambda_in} ] ; then
        echo "use ${run_lambda_in} to calculate ${foldername}"
    else
        echo "Fehler: ${run_lambda_in} ist nicht im Startordner verfuegbar." 
        echo "Fehler: ${run_lambda_in}_ist_nicht_im_Startordner_verfuegbar ${foldername}" >> ${zieldatei}
        echo "EXIT"
        exit 1
    fi


    if [ $# -ge 7 ]; then # submit g09
        sub_ok=${7}
        if [[ "${sub_ok}" == "true" ]] ; then    	
            sub_ok="true"
        else
            sub_ok="false"
        fi
    fi	
    
    function check_file_exists(){     # Schaut ob Datei existiert 
        filename=${1}
        foldername=${2}
        zieldatei=${3}
        if [ -e ${foldername}/${filename} ] ; then
            return 1
        else
            echo "Fehler: ${filename} ist nicht im Startordner verfuegbar." 
            echo "Fehler: ${filename} _ist_nicht_im_Startordner_verfuegbar ${foldername}" >> ${zieldatei}
            echo ${filename}
            echo ${foldername}
            echo "EXIT"
            exit 1
        fi
    } # file_exists
    
    if [ $# -ge 8 ]; then # submit g09
        oniom_ok=${8}
        if [[ "${oniom_ok}" == "oniom" ]] ; then    	
            oniom_ok="true"
            echo "Oniom-Rechnung ausgewählt in ${0} "
                Geo_R1=${1}  ##OptS0_n_${monomerA}
                Geo_R2=${9}  ##OptS0_el_${monomerA}
                Geo_R3=${10} ##OptS0_lo_${monomerA}
                Geo_R4=${2}  ##OptS0_n_${monomerB}
                Geo_R5=${11} ##OptS0_el_${monomerB}
                Geo_R6=${12} ##OptS0_lo_${monomerB}

                check_file_exists ${Geo_R1} ${foldername} ${zieldatei}
                check_file_exists ${Geo_R2} ${foldername} ${zieldatei}
                check_file_exists ${Geo_R3} ${foldername} ${zieldatei}
                check_file_exists ${Geo_R4} ${foldername} ${zieldatei}
                check_file_exists ${Geo_R5} ${foldername} ${zieldatei}
                check_file_exists ${Geo_R6} ${foldername} ${zieldatei}
        else
            oniom_ok="false"
        fi
    fi	
    
    
    if [[ ${input} == *'DFTB'* ]] ; then #DFTB Abfrage
        DFTB_ok='true'
        echo "DFTB_ok = ${DFTB_ok} "
        if [[ ${foldername} == *"DIPBI"* ]] ; then
            echo "Error: DFTB MIO-1-1 does not support Cl in DIPBI"
            echo "EXIT"; exit 1	
        elif [[ ${foldername} == *"DiHPBI"* ]] ; then
            sed -i 's/Cl/H/g' ${foldername}/${monomerA}.xyz
            sed -i 's/Cl/H/g' ${foldername}/${monomerB}.xyz
            echo " DFTB: Cl atoms replaced by H atoms. DIPBI to DiHPBI"
        fi
    else
            DFTB_ok='false'
    fi
   
    function make_dir(){     # Schaut ob Ordner existiert und wird, sonst erstellt
        # make_dir ${foldername}
        foldername2=${1}
        if [ ! -d ${foldername2} ] ; then
            mkdir -p ${foldername2}
        else
            echo " ${foldername} allready exists."
            return 1
        fi
    } # make_dir
   
   
   
   
    echo start: ${0} ${monomerA} ${monomerB} ${input} ${foldername} ${zieldatei}  "sub_ok=${sub_ok}" "oniom_ok=${oniom_ok}" "extra_SP=${extra_SP}"
    nproc=72	
    mem=72	
    NODE=1
    QUEUE='default'
    TIME=160
    CURDIR=`pwd`

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

    if  [[ ${foldername} == *'UFF'* ]] || [[ ${input} == *'UFF'* ]] ; then #UFF_Abfrage
        extra_SP='true'
        R15=SCF_n_OptS0_n_${monomerA}
        R16=SCF_el_OptS0_el_${monomerA}
        R17=SCF_lo_OptS0_lo_${monomerA}
    
        R18=SCF_n_OptS0_n_${monomerB}
        R19=SCF_el_OptS0_el_${monomerB}
        R20=SCF_lo_OptS0_lo_${monomerB}
    fi

    
    
# create functions

       # Erstellt SCF Input files mit Parametern
       # makeSCF Optimierung        SCF        molA   Methode Ladung  Multiplizitaet Ordnername $nproc $mem ${CURDIR} DFTB_ok
       #             1              2           3       4         5             6            7           8        9      10       11     
       # makeSCF ${Optfilame} ${SCFfilename} ${molA} ${input} ${charge} ${multiplizity} ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok}
 function makeSCF(){
	Optfilname=${1} 
	SCFfilename=${2}
 	molA=${3}
	inputline_g09=${4}
	charge=${5}
	multiplizity=${6}
	foldername=${7} 
	nproc=${8}
	mem=${9}
	CURDIR=${10}	
	DFTB_ok=${11}
    zieldatei=${12}

    if [ -e ${CURDIR}/${foldername}/${molA}/${SCFfilename}/${SCFfilename}.log ]; then
       	 echo " Die Datei ${SCFfilename}.log existiert bereits"
    else			
        if [ -e ${CURDIR}/${foldername}/${molA}/${SCFfilename}/${SCFfilename}.inp ]; then
                    echo " Die Datei ${SCFfilename}.inp existiert bereits"
        else # SCF_${charge}_{7} ${SCFfilename}.inp file for ${molA}
            if [ -e ${CURDIR}/${foldername}/${molA}/${SCFfilename}/ ]; then
                    echo "Der Ordner ${foldername}/${molA}/${SCFfilename}/ existiert bereits"
            else
                    make_dir ${CURDIR}/${foldername}/${molA}/${SCFfilename}/
            fi # Ordner existiert
            ### Erstelle SCF-Inputfile
            cd ${CURDIR}/${foldername}/${molA}/${SCFfilename}/
            echo "%NProcShared="${nproc} > ${SCFfilename}.inp
            echo "%Mem="${mem}"GB" >> ${SCFfilename}.inp
            echo "%Chk="${Optfilname}".chk" >> ${SCFfilename}.inp
            echo "#P "${inputline_g09}" NoSymm " >> ${SCFfilename}.inp
            echo "# GFINPUT IOP(6/7=3) 6D 10F " >> ${SCFfilename}.inp
            echo "" >> ${SCFfilename}.inp
            echo "${SCFfilename} ${inputline_g09}  monomer ${molA} nach Optimierung in ${Optfilename} generated with lambda_in_setup_on_the_fly.sh" >> ${SCFfilename}.inp
            echo "" >> ${SCFfilename}.inp
            echo "${charge} ${multiplizity}" >> ${SCFfilename}.inp                        
            ###echo "XYZoptimierteGeometrie " >> ${SCFfilename}.inp
            ###echo "" >> ${SCFfilename}.inp
            if [ "${DFTB_ok}" == 'true' ]; then
                DFTB_skf_Anhang ${SCFfilename}'.inp'
            fi #DFTB
            
            if [ -s ${SCFfilename}.inp ]; then
                echo " Datei ${SCFfilename}.inp erstellt. "
            else
                echo "Fehler: Datei ${SCFfilename}.inp wurde nicht erstellt!"
                echo "Fehler: Datei_${SCFfilename}.inp_wurde_nicht_erstellt  ${foldername}"  >> ../../${zieldatei}
                exit 1
            fi
            cd ${CURDIR}  
        fi # Existiert *.inp
    fi #Existiert *.log     
    return
  } # Ende function makeSCF


function makeSCF_oniom(){
	Optfilname=${1} 
	SCFfilename=${2}
 	molA=${3}
	inputline_g09=${4}
	chargeline=${5}
	foldername=${6} 
	nproc=${7}
	mem=${8}
	CURDIR=${9}	
	DFTB_ok=${10}
    zieldatei=${11}

    if [ -e ${CURDIR}/${foldername}/${molA}/${SCFfilename}/${SCFfilename}.log ]; then
       	 echo " Die Datei ${SCFfilename}.log existiert bereits"
    else			
        if [ -e ${CURDIR}/${foldername}/${molA}/${SCFfilename}/${SCFfilename}.inp ]; then
                    echo " Die Datei ${SCFfilename}.inp existiert bereits"
        else # SCF_${charge}_{7} ${SCFfilename}.inp file for ${molA}
            if [ -e ${CURDIR}/${foldername}/${molA}/${SCFfilename}/ ]; then
                    echo "Der Ordner ${foldername}/${molA}/${SCFfilename}/ existiert bereits"
            else
                    make_dir ${CURDIR}/${foldername}/${molA}/${SCFfilename}/
            fi # Ordner existiert
            ### Erstelle SCF-Inputfile
            cd ${CURDIR}/${foldername}/${molA}/${SCFfilename}/
            echo "%NProcShared="${nproc} > ${SCFfilename}.inp
            echo "%Mem="${mem}"GB" >> ${SCFfilename}.inp
            echo "%Chk="${Optfilname}".chk" >> ${SCFfilename}.inp
            echo "#P "${inputline_g09}" NoSymm " >> ${SCFfilename}.inp
            echo "# GFINPUT IOP(6/7=3) 6D 10F " >> ${SCFfilename}.inp
            echo "" >> ${SCFfilename}.inp
            echo "${SCFfilename} ${inputline_g09}  monomer ${molA} nach Optimierung in ${Optfilename} generated with lambda_in_setup_on_the_fly.sh" >> ${SCFfilename}.inp
            echo "" >> ${SCFfilename}.inp
            echo "${chargeline} " >> ${SCFfilename}.inp                        
            ###echo "XYZoptimierteGeometrie " >> ${SCFfilename}.inp
            ###echo "" >> ${SCFfilename}.inp
            if [ "${DFTB_ok}" == 'true' ]; then
                DFTB_skf_Anhang ${SCFfilename}'.inp'
            fi #DFTB
            
            if [ -s ${SCFfilename}.inp ]; then
                echo " Datei ${SCFfilename}.inp erstellt. "
            else
                echo "Fehler: Datei ${SCFfilename}.inp wurde nicht erstellt!"
                echo "Fehler: Datei_${SCFfilename}.inp_wurde_nicht_erstellt  ${foldername}"  >> ../../${zieldatei}
                exit 1
            fi
            cd ${CURDIR}  
        fi # Existiert *.inp
    fi #Existiert *.log     
    return
  } # Ende function makeSCF_oniom






# function fuer  DFTB Anhang mit skf-Files
# Angabe des Pfades zu den skf files in: skf_file_pfad
# DFTB_skf_Anhang ${filename}'.inp'
function DFTB_skf_Anhang(){
    filename_DFTB=${1}
    skf_file_pfad='/opt/files/DFTB_skf_files/mio-1-1'
    echo "@${skf_file_pfad}/H-C.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/C-H.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/C-C.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/H-H.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/S-S.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/C-S.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/S-C.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/H-S.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/S-H.skf/N" >> ${filename_DFTB}
    
    if [[ ${filename_DFTB} == *"DiHPBI"* ]] ; then
        echo "@${skf_file_pfad}/O-O.skf/N" >> ${filename_DFTB}
        echo "@${skf_file_pfad}/H-O.skf/N" >> ${filename_DFTB}
        echo "@${skf_file_pfad}/O-H.skf/N" >> ${filename_DFTB}
        echo "@${skf_file_pfad}/O-C.skf/N" >> ${filename_DFTB}
        echo "@${skf_file_pfad}/C-O.skf/N" >> ${filename_DFTB}
        echo "@${skf_file_pfad}/N-O.skf/N" >> ${filename_DFTB}
        echo "@${skf_file_pfad}/O-N.skf/N" >> ${filename_DFTB}
        echo "@${skf_file_pfad}/H-N.skf/N" >> ${filename_DFTB}
        echo "@${skf_file_pfad}/N-H.skf/N" >> ${filename_DFTB}
        echo "@${skf_file_pfad}/N-N.skf/N" >> ${filename_DFTB}
        echo "@${skf_file_pfad}/N-C.skf/N" >> ${filename_DFTB}
        echo "@${skf_file_pfad}/C-N.skf/N" >> ${filename_DFTB}
    fi
    
    
    
    echo " " >> ${filename_DFTB}
} # Ende function DFTB_skf_Anhang



#
# Start Actions
#
if [ -e ${foldername} ]; then # Check if direktory exists  
    cp ${run_lambda_in} ${CURDIR}/${foldername}/  # copy run_lambda_in
    cd ${CURDIR}/${foldername}/    
    if [ -e ${monomerA}.xyz ] && [ -e ${monomerB}.xyz ]; then    
        echo " Erstelle Ordnersystem : ${foldername}"			
    	make_dir molA/
    	make_dir molB/
    	make_dir molA/${R1}/
    	make_dir molA/${R2}/
    	make_dir molA/${R3}/    
    	make_dir molB/${R4}/
    	make_dir molB/${R5}/
    	make_dir molB/${R6}/
    else
        echo " Fehler: Die Datei ${monomerA}.xyz oder ${monomerB}.xyz existiert nicht"
        echo " Fehler: Die_Datei_${monomerA}.xyz_oder_${monomerB}.xyz_existiert_nicht ${foldername}" >> ${foldername}/${zieldatei}
        echo " ENDE " ; exit 1  
    fi
     
   # Create Optimization *.inp files
   #Prepare {R1}.inp file for molA neutral
    if [ -e ${monomerA}.xyz ]; then 
        if [ "${oniom_ok}" == "true" ]; then
            mv ${Geo_R1} molA/${R1}/ 
        else
            cp ${monomerA}.xyz molA/${R1}/ 
        fi	
           
        cd molA/${R1}/
    else
        echo "Fehler: Die Datei ${monomerA}.xyz existiert nicht; ENDE! " 
        echo "Fehler: Die_Datei_${monomerA}.xyz_existiert_nicht ${foldername}" >> ${foldername}/${zieldatei}
        exit 1
    fi 
    XYZ=`ls *.xyz` 

    echo "%NProcShared="$nproc > ${R1}.inp
    echo "%Mem="$mem"GB" >> ${R1}.inp
    echo "%Chk="${R1}".chk" >> ${R1}.inp
    echo "#P "$input" ${Opt_input} NoSymm " >> ${R1}.inp
    echo "# GFINPUT IOP(6/7=3) 6D 10F " >> ${R1}.inp
    echo "" >> ${R1}.inp
    echo "$R1 $input  monomer $monomerA neutral generated with lambda_all_setup.sh" >> ${R1}.inp
    echo "" >> ${R1}.inp
    if [ "${oniom_ok}" == "true" ]; then
        echo "0 1  0 1  0 1" >> ${R1}.inp
    else
        echo "0 1" >> ${R1}.inp
    fi
    sed -e "1,2d" $XYZ >> ${R1}.inp
    echo "" >> ${R1}.inp
    if [ "${DFTB_ok}" == 'true' ]; then
       DFTB_skf_Anhang ${R1}'.inp'
    fi

    cd ../../	

   #Prepare {R2}.inp file for molA negative
    if [ "${oniom_ok}" == "true" ]; then
        mv ${Geo_R2} molA/${R2}/ 
    else
        cp ${monomerA}.xyz molA/${R2}/   
    fi
    cd molA/${R2}/   
    XYZ=`ls *.xyz` 
    
    echo "%NProcShared="$nproc > ${R2}.inp
    echo "%Mem="$mem"GB" >> ${R2}.inp
    echo "%Chk="${R2}".chk" >> ${R2}.inp
    echo "#P "$input" ${Opt_input} NoSymm " >> ${R2}.inp
    echo "# GFINPUT IOP(6/7=3) 6D 10F " >> ${R2}.inp  
    echo "" >> ${R2}.inp
    echo "$R2 $input  monomer $monomerA negative generated with lambda_all_setup.sh" >> ${R2}.inp
    echo "" >> ${R2}.inp
    if [ "${oniom_ok}" == "true" ]; then
        echo "-1 2  -1 2  0 1" >> ${R2}.inp
    else
        echo "-1 2" >> ${R2}.inp
    fi
    sed -e "1,2d" $XYZ >> ${R2}.inp
    echo "" >> ${R2}.inp
    if [ "${DFTB_ok}" == 'true' ]; then
       DFTB_skf_Anhang ${R2}'.inp'
    fi
    cd ../../	

   #Prepare {R3}.inp file for molA positive
    if [ "${oniom_ok}" == "true" ]; then
        mv ${Geo_R3} molA/${R3}/ 
    else
        mv ${monomerA}.xyz molA/${R3}/  
    fi
    
    cd molA/${R3}/    
    XYZ=`ls *.xyz` 
    
    echo "%NProcShared="$nproc > ${R3}.inp
    echo "%Mem="$mem"GB" >> ${R3}.inp
    echo "%Chk="${R3}".chk" >> ${R3}.inp
    echo "#P "$input" ${Opt_input} NoSymm " >> ${R3}.inp
    echo "# GFINPUT IOP(6/7=3) 6D 10F " >> ${R3}.inp
    echo "" >> ${R3}.inp
    echo "$R3 $input monomer $monomerA positive generated with lambda_all_setup.sh" >> ${R3}.inp
    echo "" >> ${R3}.inp
    if [ "${oniom_ok}" == "true" ]; then
        echo "+1 2  +1 2  0 1" >> ${R3}.inp
    else
        echo "+1 2" >> ${R3}.inp
    fi
    sed -e "1,2d" $XYZ >> ${R3}.inp
    echo "" >> ${R3}.inp
    if [ "${DFTB_ok}" == 'true' ]; then
       DFTB_skf_Anhang ${R3}'.inp'
    fi
    cd ../../

   #Prepare {R4}.inp file for molB neutral
    if [ -e ${monomerB}.xyz ]; then
        if [ "${oniom_ok}" == "true" ]; then
            mv ${Geo_R4} molB/${R4}/ 
        else
            cp ${monomerB}.xyz molB/${R4}/  
        fi
    	cd molB/${R4}/    
    else
        echo "Die Datei ${monomerB}.xyz existiert nicht; ENDE! " 
         echo "Fehler: Die_Datei_${monomerB}.xyz_existiert_nicht ${foldername}" >> ${foldername}/${zieldatei}
        exit 1
    fi	    
    XYZ=`ls *.xyz` 

    echo "%NProcShared="$nproc > ${R4}.inp
    echo "%Mem="$mem"GB" >> ${R4}.inp
    echo "%Chk="${R4}".chk" >> ${R4}.inp
    echo "#P "$input" ${Opt_input} NoSymm " >> ${R4}.inp
    echo "# GFINPUT IOP(6/7=3) 6D 10F " >> ${R4}.inp
    echo "" >> ${R4}.inp
    echo "$R4 $input  monomer $monomerB neutral generated with lambda_all_setup.sh" >> ${R4}.inp
    echo "" >> ${R4}.inp
    if [ "${oniom_ok}" == "true" ]; then
        echo "0 1  0 1  0 1" >> ${R4}.inp
    else
        echo "0 1" >> ${R4}.inp
    fi
    sed -e "1,2d" $XYZ >> ${R4}.inp
    echo "" >> ${R4}.inp
    if [ "${DFTB_ok}" == 'true' ]; then
       DFTB_skf_Anhang ${R4}'.inp'
    fi

    cd ../../

   #Prepare {R5}.inp file for molB negative
    if [ "${oniom_ok}" == "true" ]; then
        mv ${Geo_R5} molB/${R5}/ 
    else
        cp ${monomerB}.xyz molB/${R5}/  
    fi
    cd molB/${R5}/ 
    XYZ=`ls *.xyz` 
    
    echo "%NProcShared="${nproc} > ${R5}.inp
    echo "%Mem="${mem}"GB" >> ${R5}.inp
    echo "%Chk="${R5}".chk" >> ${R5}.inp
    echo "#P "${input}" ${Opt_input}  NoSymm " >> ${R5}.inp
    echo "# GFINPUT IOP(6/7=3) 6D 10F " >> ${R5}.inp
    echo "" >> ${R5}.inp
    echo "$R5 $input monomer $monomerB negative generated with lambda_all_setup.sh" >> ${R5}.inp
    echo "" >> ${R5}.inp
    if [ "${oniom_ok}" == "true" ]; then
        echo "-1 2  -1 2  0 1" >> ${R5}.inp
    else
        echo "-1 2" >> ${R5}.inp
    fi
    sed -e "1,2d" ${XYZ} >> ${R5}.inp
    echo "" >> ${R5}.inp
    if [ "${DFTB_ok}" == 'true' ]; then
       DFTB_skf_Anhang ${R5}'.inp'
    fi

    cd ../../


   #Prepare {R6}.inp file for molB positive
    if [ "${oniom_ok}" == "true" ]; then
        mv ${Geo_R6} molB/${R6}/ 
    else
        mv ${monomerB}.xyz molB/${R6}/
    fi
    cd molB/${R6}/  
    XYZ=`ls *.xyz`    
    
    echo "%NProcShared="$nproc > ${R6}.inp
    echo "%Mem="$mem"GB" >> ${R6}.inp
    echo "%Chk="${R6}".chk" >> ${R6}.inp
    echo "#P "$input" ${Opt_input} NoSymm " >> ${R6}.inp
    echo "# GFINPUT IOP(6/7=3) 6D 10F " >> ${R6}.inp
    echo "" >> ${R6}.inp
    echo "$R6 $input monomer $monomerB positive generated with lambda_all_setup.sh" >> ${R6}.inp
    echo "" >> ${R6}.inp
    if [ "${oniom_ok}" == "true" ]; then
        echo "+1 2  +1 2  0 1" >> ${R6}.inp
    else
        echo "+1 2" >> ${R6}.inp
    fi
    sed -e "1,2d" $XYZ >> ${R6}.inp
    echo "" >> ${R6}.inp
    if [ "${DFTB_ok}" == 'true' ]; then
       DFTB_skf_Anhang ${R6}'.inp'
    fi

    cd ../../ 

    echo 'Optimierungsfiles erstellt'
   ## Prepare SCF-Files
       # Erstellt SCF Input files mit Parametern
       # makeSCF Optimierung        SCF        molA   Methode Ladung  Multiplizitaet Ordnername $nproc $mem ${CURDIR} 
       #             1              2            3       4        5             6             7       8     9      10           
       # makeSCF ${Optfilame} ${SCFfilename} ${molA} ${input} ${charge} ${multiplizity} ${foldername} ${Optfilename} ${nproc}
       
        cd ${CURDIR}
    if [ "${oniom_ok}" == "false" ]; then
        makeSCF ${R1} ${R7}  molA "${input}" -1 2 ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        makeSCF ${R1} ${R8}  molA "${input}" +1 2 ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        makeSCF ${R2} ${R9}  molA "${input}"  0 1 ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        makeSCF ${R3} ${R10} molA "${input}"  0 1 ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        
        makeSCF ${R4} ${R11} molB "${input}" -1 2 ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        makeSCF ${R4} ${R12} molB "${input}" +1 2 ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        makeSCF ${R5} ${R13} molB "${input}"  0 1 ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        makeSCF ${R6} ${R14} molB "${input}"  0 1 ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}

        if [ ${extra_SP} == "true" ]; then
            makeSCF ${R1} ${R15}  molA "${input}" 0 1 ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
            makeSCF ${R2} ${R16}  molA "${input}" -1 2 ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
            makeSCF ${R3} ${R17}  molA "${input}" +1 2 ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
            
            makeSCF ${R4} ${R18}  molB "${input}" 0 1 ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
            makeSCF ${R5} ${R19}  molB "${input}" -1 2 ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
            makeSCF ${R6} ${R20}  molB "${input}" +1 2 ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        fi
        
	elif [ "${oniom_ok}" == "true" ]; then
        makeSCF_oniom ${R1} ${R7}  molA "${input}" "-1 2 -1 2 0 1" ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        makeSCF_oniom ${R1} ${R8}  molA "${input}" "+1 2 +1 2 0 1" ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        makeSCF_oniom ${R2} ${R9}  molA "${input}" " 0 1  0 1 0 1" ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        makeSCF_oniom ${R3} ${R10} molA "${input}" " 0 1  0 1 0 1" ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        
        makeSCF_oniom ${R4} ${R11} molB "${input}" "-1 2 -1 2 0 1" ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        makeSCF_oniom ${R4} ${R12} molB "${input}" "+1 2 +1 2 0 1" ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        makeSCF_oniom ${R5} ${R13} molB "${input}" " 0 1  0 1 0 1" ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        makeSCF_oniom ${R6} ${R14} molB "${input}" " 0 1  0 1 0 1" ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        
        if [ ${extra_SP} == "true" ]; then
            makeSCF_oniom ${R1} ${R15}  molA "${input}" " 0 1  0 1 0 1"  ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
            makeSCF_oniom ${R2} ${R16}  molA "${input}" "-1 2 -1 2 0 1"  ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
            makeSCF_oniom ${R3} ${R17}  molA "${input}" "+1 2 +1 2 0 1"  ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
            
            makeSCF_oniom ${R4} ${R18}  molB "${input}" " 0 1  0 1 0 1"  ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
            makeSCF_oniom ${R5} ${R19}  molB "${input}" "-1 2 -1 2 0 1"  ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
            makeSCF_oniom ${R6} ${R20}  molB "${input}" "+1 2 +1 2 0 1"  ${foldername} ${nproc} ${mem} ${CURDIR} ${DFTB_ok} ${zieldatei}
        fi
        
        
    fi
    
	if [[ ${sub_ok} == "true" ]] ; then #Submittierung  
        foldername=${foldername##*/} 	
        if [ "${oniom_ok}" == "true" ]; then
            if [ -e ./Oniom/ ] ; then
                cd ./Oniom/
            else
                echo "Error: Directory ./Oniom does not exist"
            fi # exists ./Oniom
            sub-myg09-lambda_in_moria.sh ${monomerA} ${monomerB} ${foldername} ${zieldatei} oniom
            echo "sub-myg09-lambda_in_moria.sh ${monomerA} ${monomerB} ${foldername} ${zieldatei} oniom "
        else 
            if [ -e ./lambda_in/ ] ; then
                cd ./lambda_in/
                sub-myg09-lambda_in_moria.sh ${monomerA} ${monomerB} ${foldername} ${zieldatei}
                echo "sub-myg09-lambda_in_moria.sh ${monomerA} ${monomerB} ${foldername} ${zieldatei} "
            else
                echo "Error: Directory ./lambda_in/ does not exist"
            fi # exists ./lambda_in/
        fi # "${oniom_ok}" == "true"
        cd ../
    fi

else # Existiert Ordner
    echo "Error: The folder ${foldername} does not exist!"
    echo "EXIT";     exit 1		
fi # End 


 
       
