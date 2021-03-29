#/bin/bash


p_occ_filename="occP_collection_jortner_UxUyUz.dat"
tmp_p_occ_filename="tmp_${p_occ_filename}"
U_sample_max=4
U_sample_min=6
N_max_sim=70 ##70

if [ $# -ge 4 ]; then 
        p_occ_filename=${1}	
	N_max_sim=${2}
	U_sample_min=${3}
	U_sample_max=${4}
else
	echo " $0 calculates the average occupation numbers for   p_occ_trajectory.dat  files"
	echo " $0           p_occ_filename        N_max_sim     U_sample_min   U_sample_max  "
	echo " $0        occP_collection.dat         70              1              1        "
	exit 1
fi


if [[ -e ${p_occ_filename} ]] ; then
	echo "Error: The occupation file exist in current folder: $p_occ_filename"
	exit 1
fi

if [[ -e ${tmp_p_occ_filename} ]] ; then
	echo "Error: The occupation file exist in current folder: $tmp_p_occ_filename"
	exit 1
fi

echo "--- START ---"
        echo " $0           p_occ_filename        N_max_sim     U_sample_min   U_sample_max  "
	echo " $0           $p_occ_filename       $N_max_sim    $U_sample_min  $U_sample_max "

###defaults
n_sim=0
N_lines=-1
N_ref_data=-2
echo "###   Hopping pair jump list   U_sample=${U_sample}   N_max_sim=${N_max_sim} "  > hopping_pair_list.dat
for U in ` seq $U_sample_min $U_sample_max ` ; do
	for N in ` seq 1 $N_max_sim ` ; do
		DIR="kmc_fortran_U"${U}"_N"${N}
		if [ -d "$DIR" ]; then ### Take action if $DIR exists ###
			if [ -s ${DIR}/p_occ_trajectory.dat ] ; then
				N_lines=$( wc -l ${DIR}/p_occ_trajectory.dat | awk '{ print $(1) }')
				if [ "${n_sim}" == 0 ] ; then
					N_ref_data=$N_lines
					cat ${DIR}/p_occ_trajectory.dat > ${p_occ_filename}
					n_sim=$((${n_sim}+1))
					N_lines=-1
				fi
				if [[ "${N_lines}" == "${N_ref_data}" ]] ;then
					paste ${p_occ_filename}   ${DIR}/p_occ_trajectory.dat   | awk '{ if ( index(inputline,"#") == 0 ){ printf("%20.12E\n",$(1) + $(2)) } } ' > ${tmp_p_occ_filename}	
					n_sim=$((${n_sim}+1))
					tmp_N_lines=$( wc -l ${tmp_p_occ_filename} | awk '{ print $(1) }')
					if [[ "${N_lines}" != "${tmp_N_lines}" ]] ;then
					       	echo "Error: Missmatch in the number of lines ${N_lines} != ${tmp_N_lines} "	
						exit 1
					else
						mv ${tmp_p_occ_filename}  ${p_occ_filename} 
					fi
				fi  ### compar line numbers
				
			fi 
		else
			###  Control will jump here if $DIR does NOT exists ###
		echo "Error: ${DIR} not found. Can not continue."
			exit 1
		fi
	done ### loop_N
done ### loop_U

echo "n_sim done: $n_sim "
occP_summe=$( cat ${p_occ_filename} | awk 'BEGIN{ sum=0}{ sum=sum+$(1)  }END{ print sum }' ) 
echo "Sum of all occupation numbers:  $occP_summe    ==  $n_sim ???"
echo "New average occupation file created:  ${p_occ_filename} "
echo "Normalize occupation numbers"
mv ${p_occ_filename} ${tmp_p_occ_filename} 
cat ${tmp_p_occ_filename} | awk -v occP_summe="${occP_summe}" ' { printf( "%20.12E\n", $(1)/occP_summe) } ' > ${p_occ_filename}
rm ${tmp_p_occ_filename}
echo "Normalized occupation numbers:   ${p_occ_filename}"
