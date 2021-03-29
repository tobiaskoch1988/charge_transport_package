#!/bin/bash

if [ $# -ne 2 ]; then
	echo "  Programm creates the input CPMD for *.cube files of frontier orbitals  e.g. for Donor /Acceptor in CDFT results"
	echo "  You need to provide a CPMD input and the logfile "
	echo "  Start $0 in CDFT folder with converged data. "
	echo "  $0   CPMD.inp     CPMD.out "
	exit 1
fi
CDFT_inputfile=$1
CDFT_logfile=$2


echo "Use: CDFT_inputfile: ${CDFT_inputfile}        CDFT_logfile: ${CDFT_logfile} "
if [[ ! -e "${CDFT_inputfile}" ]] ; then
        echo "Error the file does not exist in the current folder: ${CDFT_inputfile} "
        exit 1
fi
if [[ ! -e "${CDFT_logfile}" ]] ; then
        echo "Error the file does not exist in the current folder: ${CDFT_logfile} "
        exit 1
fi
if [[ "${CDFT_inputfile: -4}" != ".inp" ]] ; then
        echo "Error; the kmc options file is not a inp file: ${CDFT_inputfile}"
        exit 1
fi
if [[ "${CDFT_logfile: -4}" != ".out" ]] ; then
        echo "Error; the kmc options file is not a .out file: ${CDFT_logfile}"
        exit 1
fi

acc_file="$( basename ${CDFT_inputfile} .inp )_acc.inp"
don_file="$( basename ${CDFT_inputfile} .inp )_don.inp"
if [[ -e "${acc_file}" ]] ; then
        echo "Error the file already exists: ${acc_file} "
        exit 1
fi
if [[ -e "${don_file}" ]] ; then
        echo "Error the file already exists: ${don_file} "
        exit 1
fi


BETA=$( grep "NUMBER OF BETA STATES:" ${CDFT_logfile} | awk '{ print $(5) }' )
##### Check if BETA is an interger
if ! [[ "$BETA" =~ ^[0-9]+$ ]] ; then
	echo "Error: NUMBER OF BETA STATES: not found in  ${CDFT_logfile}  "
        echo "Sorry integers only:     $BETA"
	exit 1
else
	echo "NUMBER OF BETA STATES: $BETA "
fi




cat ${CDFT_inputfile} |  awk -v BETA="${BETA}" 'BEGIN{ print_on=1 ; CPMD_skip=0 } 
{
	if( index($(1),"&CPMD") != 0){  CPMD_skip=1 ;  print_on=0 } 
	
	if(print_on==1) { print $(0) }
	
	
	if( (CPMD_skip==1) && ( index($(1),"&END") != 0) ){
		print("&CPMD")										
		print("                CENTER MOLECULE ON")
		print("                PROPERTIES")
		print("                LSD")
		print("                RESTART WAVEFUNCTION COORDINATES")
		print("                RHOOUT")
		print("&END")
		print(" ")
		print(" ")
		print("&PROP ")
		print("        CUBEFILE ORBITALS ")
		print("        3 ")
		print("        "BETA-2"  "BETA-1"  "BETA" ")
		print("&END ")										
		
		print_on=1
		CPMD_skip=0
		} 

	} END{}' > ${acc_file}


cat ${CDFT_inputfile} |  awk -v BETA="${BETA}" 'BEGIN{ print_on=1 ; CPMD_skip=0 } 
{
	if( index($(1),"&CPMD") != 0){  CPMD_skip=1 ;  print_on=0 } 
	
	if(print_on==1) { print $(0) }
	
	
	if( (CPMD_skip==1) && ( index($(1),"&END") != 0) ){
		print("&CPMD")										
		print("                CENTER MOLECULE ON")
		print("                PROPERTIES")
		print("                LSD")
		print("                RESTART WAVEFUNCTION COORDINATES")
		print("                RHOOUT")
		print("&END")
		print(" ")
		print(" ")
		print("&PROP ")
		print("        CUBEFILE ORBITALS ")
		print("        3 ")
		print("        "BETA-2"  "BETA-1"  "BETA" ")
		print("&END ")										
		
		print_on=1
		CPMD_skip=0
		} 

	} END{}' > ${don_file}


if [[ ! -e "${acc_file}" ]] ; then
        echo "Error the file was not created: ${acc_file} "
        exit 1
else
	echo "New file: ${acc_file}"
	pfadmove ${acc_file}
	folder=$( basename ${acc_file} .inp )
	if [[ -e "RESTART.STATE_1" ]]  ; then
		echo "The file exists: RESTART.STATE_1 "
		if [[ -e "$folder" ]]  ; then
			cp RESTART.STATE_1  $folder/RESTART
			cp *.psp $folder
			echo "--- START ---"
			echo " cd $folder "
			echo " sub-cpmd_CDFT_highmem ${acc_file} 70 1 1 express 2 "
			echo " cd .. "
		fi
	else
		echo "File does not exist: RESTART.STATE_1 "
		exit 1
	fi	
fi ### acc_file


if [[ ! -e "${don_file}" ]] ; then
        echo "Error the file was not created: ${don_file} "
        exit 1
else
	echo "New file: ${don_file}"
	pfadmove ${don_file}
	folder=$( basename ${don_file} .inp )
	if [[ -e "RESTART.STATE_2" ]]  ; then
		echo "The file exists: RESTART.STATE_2 "
		if [[ -e "$folder" ]]  ; then
			cp RESTART.STATE_2  $folder/RESTART
			cp *.psp $folder
			echo " cd $folder "
			echo " sub-cpmd_CDFT_highmem ${don_file} 70 1 1 express 2 "
			echo " cd .. "
		fi
	else
		echo "File does not exist: RESTART.STATE_2 "
		exit 1
	fi
	
fi ### don_file
