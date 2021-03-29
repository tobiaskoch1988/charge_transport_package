#!/bin/bash
# bash-script creates dimer.inp and monomer.inp files and folder structure for the DIPRO methode
# you need monomerA.xyz and monomerB.xyz files in foldername 
# $1 provides molA.xyz
# $2 provides molB.xyz
# $3 provides dimername
# $4 provides functional/basisset
# $5 provides the foldername
# dimer run; prepare orbitals from monomers

if [ $# -eq 0 ]; then
	echo "Usage: $0 < monomerA.xyz monomerB.xyz dimername functional/basisset foldername >  Optional: < Geo_data_info > "
	exit 1
fi

if [ ! -d "$5/" ]; then
    echo "Folder $5 doesn't exist !"
    echo "Usage: $0 < monomerA.xyz monomerB.xyz dimername functional/basisset foldername >   Optional: < Geo_data_info > "
    exit 1
fi

if [ ! -f "$5/$1" ]; then
	echo "File $1 doesn't exist in $5/ !"
	echo "Usage: $0 < monomerA.xyz monomerB.xyz dimername functional/basisset foldername >  Optional: < Geo_data_info > "
	exit 1
fi

if [ ! -f "$5/$2" ]; then
    echo "File $2 doesn't exist in $5/ !"
    echo "Usage: $0 < monomerA.xyz monomerB.xyz dimername functional/basisset foldername >   Optional: < Geo_data_info > "
    exit 1
fi

monomerA=`basename $1 .xyz`
monomerB=`basename $2 .xyz`
dim0=${3}

calculate_molA_again='True'
calculate_molB_again='True'

nproc=8	
mem=8	
CURDIR=`pwd`	


if [[ ${HOSTNAME} == *"edoras"* ]] ; then
	nproc=1
	mem=1
elif [[ ${HOSTNAME} == *"palma"* ]] ; then
	nproc=12
	mem=24
elif [[ ${HOSTNAME} == *"balrog"* ]] || [[ ${HOSTNAME} == *"kankra"* ]] ; then
	nproc=8
	mem=8
elif [[ ${HOSTNAME} == *"sl"* ]] ; then ### phi
	nproc=24
	mem=48
elif [[ ${HOSTNAME} == *"morgoth"* ]] ; then 
	nproc=40
	mem=64
fi


if [ -n "${4}" ]; then
    calculation_method=${4}
else
    calculation_method="B3LYP/6-311G**"
fi

if [ -n "${5}" ]; then	
    foldername=${5}
else   	
    foldername="ordner"
fi		    

### to provide some extra information in the data files
if [ $# -eq 6 ]; then
    Geo_data_info=${6}
else
    Geo_data_info=" "
fi




if [[ "${calculation_method}" == *"DFTB"* ]]; then 
    DFTB_ok='true'
    echo "DFTB_ok with MIO-1-1 = ${DFTB_ok} "
else
    DFTB_ok='false'
fi

# DFTB Anhang funktion fuer P3MT
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

    echo "@${skf_file_pfad}/O-O.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/H-O.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/O-H.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/C-O.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/O-C.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/O-S.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/S-O.skf/N" >> ${filename_DFTB}

    echo "@${skf_file_pfad}/N-N.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/H-N.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/N-H.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/C-N.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/N-C.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/O-N.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/N-O.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/N-S.skf/N" >> ${filename_DFTB}
    echo "@${skf_file_pfad}/S-N.skf/N" >> ${filename_DFTB}

    echo " " >> ${filename_DFTB}
} # Ende function DFTB_skf_Anhang
	




function DIPRO_neighbourcheck() {
### function to check if the other DIPRO calculations were performed and copy folder to the directory; 
### resid_X is the resid which needs to be checked if it has been calculated before.
### calculate_molA/molB_again= True /False say if there was another folder found.
### function needs to be executed in the DIPRO folder(Bsp: foldername=Dim_DIPBI_DIPBI_0_207_PM3_G0_N2 )
### if foldername/${monomerX}.xyz exists, the geometry is compaired with the second geometry of the neighbours and only if the files are identical the folder is copied ; in order to avoid shifts due to periodic boundary conditions.
resid_X=${1}
mol_X=${2}
monomerX=${3}
foldername=${4}
calculate_molA_again=${5}
calculate_molB_again=${6}
neighbour_folder_list=($( ls -d ../Dim*_${resid_X}_*_N* ) )

echo "vorher, ${1} ${2} ${3} ${4} ${5} ${6}"
#echo "${#neighbour_folder_list[@]} ${neighbour_folder_list[@]}"

if [[ ${#neighbour_folder_list[@]} -lt 2 ]] ; then ## check uf the list is to short to check for suitable neighbours
    echo " calculate_molA_again: ${calculate_molA_again} ,calculate_molB_again: ${calculate_molB_again} "
    return 
fi 
## delete current foldername from the list
i=0
for name in ${neighbour_folder_list[@]} ; do
        if [ "${name}" == "../${foldername}" ]; then
                ### delete one element ( ${foldername})  in array
                neighbour_folder_list=(${neighbour_folder_list[@]:0:$i} ${neighbour_folder_list[@]:$(($i + 1))})
                break
        fi # foldername 
        i=$(($i+1))
done

echo "${neighbour_folder_list[@]}"

for neighbour_folder in ${neighbour_folder_list[@]} ; do
    for mol_C in 'molA' 'molB' ; do # mol_C sucht die moeglichen Nachbarn ab.
     if [ -e ${neighbour_folder}/${mol_C}/${monomerX}'.xyz' ] && [ -e ./${monomerX}'.xyz' ] ; then
        if [ -e ${neighbour_folder}/${mol_C}/'monomer.log' ] && [ -e ${neighbour_folder}/${mol_C}/fort.7 ] && [ -e ./${mol_X} ] ; then
               #echo "${neighbour_folder}/${mol_C}/monomer.log"
               #echo $( ls ${neighbour_folder}/${mol_C}/fort.* )
                normalR1=$(tail -1 ${neighbour_folder}/${mol_C}/monomer.log | cut -c 1-19)
                if [[ ${normalR1} == *"Normal termination"* ]]; then
                        neighbour_xyz=$( ls ${neighbour_folder}/${mol_C}/${monomerX}'.xyz'  )
                        ### check if both input geometries were identical; so that no changes due to the periodic boundary conditions were maid.
                        str=$( diff ./${monomerX}'.xyz'  ${neighbour_xyz}  | tail -1 )
                        if [ -z "${str}" ] ; then # leerer String ?? => keine Geometrieunterschiede?
                                echo " str leer: ${str}"
                                ## Copy folder with data here
                                cp -r ${neighbour_folder}/${mol_C}/* ${mol_X}
                                echo " use/copy old data from ${neighbour_folder}/${mol_C}/* to ./${mol_X} "
                                if [ "${mol_X}" == 'molA' ] ; then
                                        calculate_molA_again='False'
                                        echo " calculate_molA_again: ${calculate_molA_again} ,calculate_molB_again: ${calculate_molB_again} "
                                        return
                                        break 3
                                elif [ "${mol_X}" == 'molB' ] ; then
                                        calculate_molB_again='False'
                                        echo " calculate_molA_again: ${calculate_molA_again} ,calculate_molB_again: ${calculate_molB_again} "
                                        return
                                        break 3
                                fi # set calculate again statements  
                        fi # Geometrien gleich ?
                fi # Normal termination
        fi # existiert *.log
      fi # exiistiert *.xyz     
    done # end loop molX = molA, molB
done # end loop neighbour_folder_list

echo " calculate_molA_again: ${calculate_molA_again} ,calculate_molB_again: ${calculate_molB_again} "
} # End function_DIPRO_neighbourcheck



######Create file structure
    
if [ ! -d ${foldername} ]; then # Create folder, if it does not exist
   mkdir  ${foldername}
fi

cd ${foldername}
if [ ! -d dim/ ]; then # Create folder, if it does not exist
    mkdir -p dim/
fi
if [ ! -d molA/ ]; then # Create folder, if it does not exist
    mkdir -p molA/
fi
if [ ! -d molB/ ]; then # Create folder, if it does not exist
    mkdir -p molB/
fi


###### Start: Check if Resids were already calculated #######

### DIPRO_neighbourcheck Resid1/Resid2 'molA' monomerA Ordner_N1 True/False True/False
output_line=$( DIPRO_neighbourcheck ${monomerA##*_} 'molA' ${monomerA} ${foldername} ${calculate_molA_again} ${calculate_molB_again} )
echo "${output_line}"
output_line=${output_line##*calculate_molA_again:} 
calculate_molA_again="${output_line%%,calculate_molB_again:*}"
calculate_molB_again=${output_line##*calculate_molB_again:}

 
output_line=$( DIPRO_neighbourcheck ${monomerB##*_} 'molB' ${monomerB} ${foldername} ${calculate_molA_again} ${calculate_molB_again} )
echo "${output_line}"
output_line=${output_line##*calculate_molA_again:} 
calculate_molA_again="${output_line%%,calculate_molB_again:*}"
calculate_molB_again=${output_line##*calculate_molB_again:}
echo "${output_line}"

echo " calculate_molA_again: ${calculate_molA_again} "
echo " calculate_molB_again: ${calculate_molB_again} "
###### End check if Resids were already calculated #######


###Create dimer.xyz file

NummerA=$(head -1 ${monomerA}'.xyz') 
NummerB=$(head -1 ${monomerB}'.xyz') 
NummerC=$((${NummerA} + ${NummerB}))   
 echo ${NummerC} > ${dim0}.xyz
 echo ' ' >> ${dim0}.xyz
 tail -$((${NummerA})) ${monomerA}.xyz >> ${dim0}.xyz 
 tail -$((${NummerB})) ${monomerB}.xyz >> ${dim0}.xyz 

# prepare dimer .inp  file
mv ${dim0}.xyz dim/
cd dim/

XYZ=`ls *.xyz`
echo "%nproc="$nproc > dimer.inp
echo "%mem="$mem"GB" >> dimer.inp
echo "# pop=minimal ${calculation_method} nosymm  IOp(3/33=1) punch=mo" >> dimer.inp
echo "" >> dimer.inp
echo "dimer $dim0 run generated from $monomerA and $monomerB by mol_all_setup.sh  with ${Geo_data_info} " >> dimer.inp
echo "" >> dimer.inp
echo "0 1" >> dimer.inp
sed -e "1,2d" $XYZ >> dimer.inp
echo "" >> dimer.inp
if [ "${DFTB_ok}" == 'true' ] && [ "${monomerA}" != "${monomerA%P3*}" ] && [ "${monomerB}" != "${monomerB%P3*}" ] ; then ##DFTB und P3MT/P3HT
    DFTB_skf_Anhang dimer.inp
fi 
    
    
if [[ "${calculate_molA_again}" == *'True'* ]] ; then
    echo " prepare A"   
    #Prepare monomer.inp file for molA
    cd ../molA
    mv ../${monomerA}.xyz .
    
    XYZ=`ls *.xyz` 

    echo "%nproc="$nproc > monomer.inp
    echo "%Mem="$mem"GB" >> monomer.inp
    echo "#pop=minimal ${calculation_method} nosymm punch=mo IOp(6/7=3)" >> monomer.inp
    echo "" >> monomer.inp
    echo "monomer $monomerA generated from mol_all_setup.sh with ${Geo_data_info}  " >> monomer.inp
    echo "" >> monomer.inp
    echo "0 1" >> monomer.inp
    sed -e "1,2d" $XYZ >> monomer.inp
    echo "" >> monomer.inp
    if [ "${DFTB_ok}" == 'true' ] && [ "${monomerA}" != "${monomerA%P3*}" ] ; then ##DFTB und P3MT/P3HT
        DFTB_skf_Anhang monomer.inp
    fi 
fi # calculate_molA_again        	



if [[ "${calculate_molB_again}" == *'True'* ]] ; then
    echo " prepare B"
    #Prepare monomer.inp file for molB
    cd ../molB
    mv ../${monomerB}.xyz .
    
    XYZ=`ls *.xyz` 

    echo "%nproc="$nproc > monomer.inp
    echo "%Mem="$mem"GB" >> monomer.inp
    echo "#pop=minimal ${calculation_method} nosymm punch=mo IOp(6/7=3)" >> monomer.inp
    echo "" >> monomer.inp
    echo "monomer $monomerB generated from mol_all_setup.sh with ${Geo_data_info} " >> monomer.inp
    echo "" >> monomer.inp
    echo "0 1" >> monomer.inp
    sed -e "1,2d" $XYZ >> monomer.inp
    echo "" >> monomer.inp
    if [ "${DFTB_ok}" == 'true' ] && [ "${monomerB}" != "${monomerB%P3*}" ]; then ##DFTB und P3MT/P3HT
        DFTB_skf_Anhang monomer.inp
    fi 
    cd ${CURDIR}
fi # calculate_molB_again
