#!/bin/bash
charge="T1"
for Geo in ` ls *.xyz `  ; do 
	file=$( basename $Geo .xyz )
	newfilename="Opt_${charge}_${file}_PBE0_SDD_D3BJ_SDD_vac.inp"
	echo "%NProcShared=8"    >> $newfilename
	echo "%Mem=8GB"    >> $newfilename
	##echo "%Chk=SCF_${charge}_PBE0_6-31Gs_CHelpG_$file.chk"  >> $newfilename
	echo "# Opt(tight) PBE1PBE/SDD SCF(XQC,MaxConventionalCycles=400,MaxCycle=800) EmpiricalDispersion=GD3BJ   pop=(minimal,CHelpG) nosymm IOp(6/7=3)" >> $newfilename
	##echo "#pop=minimal B3LYP/6-31G* polar CPHF(RdFreq)"  >> $newfilename
	##echo "#SCF(XQC,MaxConventionalCycles=400,MaxCycle=800) pop=(minimal,CHelpG) nosymm punch=mo IOp(6/7=3)" >> $newfilename
	##### echo " PBE1PBE/6-31G* pop(hirshfeld) density=current iop(7/33=1) iop(6/80=1)  TD(NStates=50)  NoSymm EmpiricalDispersion=GD3BJ  SCRF(PCM,Solvent=Chloroform,Read)  '
    #### tetrahydrofuran
	echo " "                    >> $newfilename
	echo "monomer $Geo ${charge} $file generated with SCF_make_Inputfiles.sh"  >> $newfilename
	echo " "                    >> $newfilename
	if [ "${charge}" == 'n' ]; then
		echo "0 1"     >> $newfilename
	elif [ "${charge}" == 'el' ]; then
		 echo "-1 2"  >> $newfilename 
	elif [ "${charge}" == 'lo' ]; then	
                 echo "+1 2"  >> $newfilename
    elif [ "${charge}" == 'S0' ]; then
			echo "0 1"     >> $newfilename
    elif [ "${charge}" == 'T1' ]; then
			echo "0 3"     >> $newfilename
    elif [ "${charge}" == 'S1' ]; then
			echo "0 1"     >> $newfilename
    elif [ "${charge}" == 'TDDFT' ]; then
			echo "0 1"     >> $newfilename
	fi	
	linenumber=$( head -1 ${Geo} | awk '{ print $(1) }' )
	tail -$linenumber ${Geo}     >> $newfilename
	echo " "                    >> $newfilename 
	echo "" >> $newfilename
	echo " " >> $newfilename

done
