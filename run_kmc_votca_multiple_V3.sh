#!/bin/bash
### bash script to calculate multipe kMC-simulations with different Voltages U
### needs 1) options_VOTCA_kMC_file= options.xml
###       2) kMC_VOTCA_startline: Bsp: " kmc_run -e kmcmultiple -o options.xml -f state.sql "
###       3) Modus = [increasing/multidirection/multidirection_scaled]
###       4) |U_ext| = external voltage [V/m]
### Directions can be selected with U_list ; Caution: U_ext is independent from / overrides the voltage in the inputfile: options_VOTCA_kMC_file.xml
### With "N_kMC_random_insertions =  5" you can loop for one voltage direction over (5) start configurations if random_charge_insertion = True in the options_VOTCA_kMC_file
###  
### Use: rate_type = "jortner" or "marcus" or "weiss_dorsey"     or "rate_type_scan" always in combination with use_my_rates_to_votca_sql='True'
### 
### The insertion modus can be tuned, to select only insertion on different P3HT: use_P3HT_charge_insejection="True" or on DIPBI use_DIPBI_charge_insejection='False', 
### Therefore the random charge insertion needs to be switched off: use_random_charge_insejection="False"
### Number of CPUS per Node used for gnu parallel ${NCPUS} ; NCPUS is determined in the script
### - can read options from the options.xml file to bash variables. (New readable options need to be appended to keywordlist)
### 
### Also run_kmc_fortran can be executed in parallel mode if you set use_kmc_fortran='True'; default: use_kmc_fortran='False' =>  use  kmc_run (VOTCA).
### Using:  run_kmc_fortran             KMC_FOR_MULTIPLE_CHARGES_xml   options.xml   statefile.sql
###  or     charge_transport_package    KMC_FOR_MULTIPLE_CHARGES_xml   options.xml   statefile.sql
### 
#set -vx

### read options from the options.xml file
read_options_from_xml="True"
### keywordlist for all readable options
keywordlist=$( echo "N_kMC_random_insertions" "use_random_charge_insejection"  "calc_parallel" "use_my_rates_to_votca_sql"  "use_kmc_fortran"  "use_P3HT_charge_insejection"  "use_DIPBI_charge_insejection"  "use_PPDI_charge_insejection"  "use_PBDT_TS1_charge_insejection" "system"  "modify_rates"  "rate_type"  "dump_rates"  "cutoff_dump_rate"  "N_xy_plane"  "readsystem")



### make sure that options_VOTCA_kMC_filename is correct!
N_kMC_random_insertions=4

### Start simulations at specific injection sites use_P3HT_charge_insejection, use_DIPBI_charge_insejection default use_random_charge_insejection="True"
use_random_charge_insejection="True"
### default: use_P3HT_charge_insejection="False"
use_P3HT_charge_insejection="False"
### default: use_DIPBI_charge_insejection='False'
use_DIPBI_charge_insejection='False'
### default: use_PPDI_charge_insejection='False'
use_PPDI_charge_insejection='False'
### default: use_PBDT_TS1_charge_insejection='False'
use_PBDT_TS1_charge_insejection='False'


### calculate many kMC trajectories in parallel.
calc_parallel="True" 

### use xtp_createStatefile_fromtxt_t_koch08.py to create sql files
use_my_rates_to_votca_sql='False'

### use kmc_fortran program for kmc simulations e.g. for multicharge kmc simulations with VSSM method.
### read (some) setup data from the votca options.xml file format and the rate data from a votca *.sql file. : Default: use_kmc_fortran='False' ==> VOTCA kmc_multiple is used.
use_kmc_fortran='True'

### SYSTEM= DIPBI_P3HT_500K or DIPBI_P3HT_700K, DIPBI_P3HT_900K; or  P3HT_pur ; PPDI_PBDT_TS1
SYSTEM='DIPBI_P3HT_500K' 

### adapt the rate calculation e.g. necessary for different U_ext or temperature T
modify_rates='False'

### rate_type = "jortner" or "marcus" or "weiss_dorsey"     or rate_type_scan
rate_type="rate_type_scan"

### dump rates to 1 for a certain cutoff
dump_rates='True'

### default: cutoff_dump_rate=1.0E-15
cutoff_dump_rate="1.0E-10"

### also for xy, xz and yz-plane
N_xy_plane=36



### default initials 
U_ext="1.0E+7"
NCPUS=72 
U_array_new=" 0.003,0,0 "
U_list=" 0,0,1 " 
NGeo=4
N_frame=0
NPROCESSORS=1
multidirection_modus='single'
temperature=300
T_list=" ${temperature} "
rate_type_list=" ${rate_type} "

### key to read the system from the xml file
readsystem='False'



if [ $# -eq 0 ]; then
	echo "Start $0 with: ${args[@]} "
	echo "   "
    echo "   Modus=[single,singleUx,singleUy,singleUz,UxUyUz,xy-plane,yz-plane,xz-plane,increasing,increasing_small_steps,increasing_small_stepsUx,increasing_small_stepsUy,increasing_small_stepsUz,multidirection,multidirection_scaled,temperature_increasing] " 
    echo "Usage:    $0  options_VOTCA_kMC_file  kMC_VOTCA_startline    Modus        U_ext    "
    echo "Example:  $0        options.xml            state.sql        single       1.0E+7    "
    exit 1
else
	args=("$@")
	echo "Start  $0 "
	echo "       $0    ${args[@]} "
fi

##### function U_xy_plane rotates the vector in the xy-plane for and gives a list of U_list.
##### The argument N_xy_plane gives the number of evaluated segments in the circle.
####   Examples: plane 10 degrees steps 
####	N_xy_plane=36
#####   U_list=$( U_xy_plane "${N_xy_plane}" )
function U_xy_plane(){
        if [ $# -lt 1 ] ; then
            echo "> Argument missing in U_xy_plane. Expected: 1"
            exit 1
        fi

        N_xy_plane=${1}
        U_list=$( echo " " | awk -v N_xy_plane="${N_xy_plane}" 'function abs(v){return v < 0 ? -v : v}
         BEGIN{ PI=3.141592}{
         U_list=""
         for( i =1; i <=N_xy_plane ; i++){
                 x=sin((2*i*PI)/N_xy_plane)
                 y=cos((2*i*PI)/N_xy_plane)
                 if ( abs(x) < 1.0E-5){ x=0.0}
                 if ( abs(y) < 1.0E-5){ y=0.0}   
                 U_step=sprintf("%.4g,%.4g,0.0",x,y)
                 ###print U_step
                 U_list=U_list" "U_step 
                 ### sprintf("%s", printf("%s %.4g,%.4g,0.0",U_list,sin((2*i*PI)/N_xy_plane),cos((2*i*PI)/N_xy_plane)) )
         } #loop        
         } END { print U_list }' )
        echo "${U_list}"
}
##### End Function  U_xy_plane


function U_yz_plane(){
        if [ $# -lt 1 ] ; then
            echo "> Argument missing in U_yz_plane. Expected: 1"
            exit 1
        fi

        N_xy_plane=${1}
        U_list=$( echo " " | awk -v N_yz_plane="${N_xy_plane}" 'function abs(v){return v < 0 ? -v : v}
         BEGIN{ PI=3.141592}{
         U_list=""
         for( i =1; i <=N_yz_plane ; i++){
                 y=sin((2*i*PI)/N_yz_plane)
                 z=cos((2*i*PI)/N_yz_plane)
                 if ( abs(y) < 1.0E-5){ y=0.0}
                 if ( abs(z) < 1.0E-5){ z=0.0}   
                 U_step=sprintf("0.0,%.4g,%.4g",y,z)
                 ###print U_step
                 U_list=U_list" "U_step 
                 ### sprintf("%s", printf("%s 0.0,%.4g,%.4g",U_list,sin((2*i*PI)/N_yz_plane),cos((2*i*PI)/N_yz_plane)) )
         } #loop        
         } END { print U_list }' )
        echo "${U_list}"
}
##### End Function  U_yz_plane


function U_xz_plane(){
        if [ $# -lt 1 ] ; then
            echo "> Argument missing in U_xz_plane. Expected: 1"
            exit 1
        fi

        N_xy_plane=${1}
        U_list=$( echo " " | awk -v N_xz_plane="${N_xy_plane}" 'function abs(v){return v < 0 ? -v : v}
         BEGIN{ PI=3.141592}{
         U_list=""
         for( i =1; i <=N_xz_plane ; i++){
                 x=sin((2*i*PI)/N_xz_plane)
                 z=cos((2*i*PI)/N_xz_plane)
                 if ( abs(x) < 1.0E-5){ x=0.0}
                 if ( abs(z) < 1.0E-5){ z=0.0}   
                 U_step=sprintf("%.4g,0.0,%.4g",x,z)
                 ###print U_step
                 U_list=U_list" "U_step 
                 ### sprintf("%s", printf("%s %.4g,0.0,%.4g",U_list,sin((2*i*PI)/N_xz_plane),cos((2*i*PI)/N_xz_plane)) )
         } #loop        
         } END { print U_list }' )
        echo "${U_list}"
}
##### End Function  U_yz_plane





### function to check if another script terminated with an error and exits this script
function check_exit_1(){
	error_in="$1" 
	error_in_line="$2" 
	if [[ "${error_in}" == "1" ]] ; then
		echo "GLOBAL ERROR: set to 1"
		echo "IN: ${error_in_line}: TERMINATE EXECUTION OF ${0} IN "$( pwd ) 
		exit 1
	fi
}
export -f check_exit_1
export -f U_xy_plane




### Start reading data

if [ $# -ge 1 ]; then   ## Einlesen des Starts des KMC-Programmes mit einer Rechnung oder standard
    kMC_options_file=${1}
else 
    echo "Usage: $0 options_VOTCA_kMC_file kMC_VOTCA_statefile Modus=[single,singleUx,singleUy,singleUz,UxUyUz,xy-plane,yz-plane,xz-plane,increasing,increasing_small_steps,increasing_small_stepsUx,increasing_small_stepsUy,increasing_small_stepsUz,multidirection,multidirection_scaled]  U_ext "
    exit 1
fi

if [ $# -ge 2 ]; then   ## Einlesen des kMC_VOTCA_statefile des KMC-Programmes mit einer Rechnung oder standard
    kMC_statefile=${2}
else 
    echo "Usage: $0 options_VOTCA_kMC_file kMC_VOTCA_statefile Modus=[single,singleUx,singleUy,singleUz,UxUyUz,xy-plane,yz-plane,xz-plane,increasing,increasing_small_steps,increasing_small_steps,multidirection,multidirection_scaled]  U_ext "
    exit 1
fi

if [ $# -ge 3 ]; then   ## Einlesen des Modus des KMC-Programmes mit einer Rechnung
    multidirection_modus=${3}
    if [ "${multidirection_modus}" == "single" ] ; then
	U_list=" 0,0,1 "
    elif [ "${multidirection_modus}" == "singleUx" ] ; then
	U_list=" 1,0,0 "
    elif [ "${multidirection_modus}" == "singleUy" ] ; then
	U_list=" 0,1,0 "
    elif [ "${multidirection_modus}" == "singleUz" ] ; then
	U_list=" 0,0,1 "
    elif [ "${multidirection_modus}" == "UxUyUz" ] ; then
	U_list=" 1,0,0   0,1,0   0,0,1 "
    elif [ "${multidirection_modus}" == "xy-plane" ] ; then
	U_list=$( U_xy_plane "${N_xy_plane}" )
    elif [ "${multidirection_modus}" == "yz-plane" ] ; then
	U_list=$( U_yz_plane "${N_xy_plane}" )
    elif [ "${multidirection_modus}" == "xz-plane" ] ; then
	U_list=$( U_xz_plane "${N_xy_plane}" )
    elif [ "${multidirection_modus}" == "increasing_small_steps" ] ; then
	U_list=" 1,0,0   0,1,0   0,0,1   3,0,0   0,3,0   0,0,3   6,0,0   0,6,0   0,0,6   9,0,0   0,9,0   0,0,9   12,0,0   0,12,0   0,0,12   15,0,0   0,15,0   0,0,15   18,0,0   0,18,0   0,0,18   21,0,0   0,21,0   0,0,21   24,0,0   0,24,0   0,0,24   27,0,0   0,27,0   0,0,27   30,0,0   0,30,0   0,0,30 " 
    elif [ "${multidirection_modus}" == "increasing_small_stepsUx" ] ; then
	U_list=" 1,0,0   3,0,0   6,0,0   9,0,0   12,0,0   15,0,0   18,0,0   21,0,0   24,0,0   27,0,0   30,0,0 " 
    elif [ "${multidirection_modus}" == "increasing_small_stepsUy" ] ; then
	U_list=" 0,1,0   0,3,0   0,6,0   0,9,0   0,12,0   0,15,0   0,18,0   0,21,0   0,24,0   0,27,0   0,30,0 " 
    elif [ "${multidirection_modus}" == "increasing_small_stepsUz" ] ; then
	U_list=" 0,0,1   0,0,3   0,0,6   0,0,9   0,0,12   0,0,15   0,0,18   0,0,21   0,0,24   0,0,27   0,0,30 " 
    elif [ "${multidirection_modus}" == "increasing" ] ; then
        #U_list=" 0.001,0,0   0.006,0,0   0.009,0,0   0.012,0,0   0.015,0,0   0.018,0,0 " 
		U_list=" 0,0,1   0,0,5   0,0,10   0,0,50   0,0,100   0,0,500   0,0,1000   0,0,5000    0,0,10000   0,0,50000 "
       	U_list="${U_list} 0,0,1.0E5    0,0,5.0E5   0,0,1.0E6   0,0,5.0E6   0,0,1.0E7   0,0,5.0E7   0,0,1.0E8   0,0,5.0E8   0,0,1.0E9   0,0,5.0E9 "
      	U_list="${U_list} 0,0,1.0E10   0,0,5.0E10   0,0,1.0E11   0,0,5.0E11   0,0,1.0E12   0,0,5.0E12   0,0,1.0E13   0,0,5.0E13   0,0,1.0E14   0,0,5.0E14 "
		U_list="${U_list} 0,0,1.0E14   0,0,5.0E14   0,0,1.0E15   0,0,5.0E15   0,0,1.0E16   0,0,5.0E16 "	
    elif [ "${multidirection_modus}" == "multidirection" ] ; then
        U_list=" 0.009,0,0   0,0.009,0   0,0,0.009   0.003,0.003,0   0,0.003,0.003   0.003,0,0.003 " 
        ## U_list="${U_list} [-0.003,-0.003,0    -0.003,0.003,0   0.003,-0.003,0   0,0.003,-0.003    0,-0.003,0.003   0,-0.003,-0.003 "
        ## U_list="${U_list} [-0.003,0,0.003   0.003,0,-0.003   -0.003,0,-0.003 "
        
        ###24 directions positiv and negativ 
    elif [ "${multidirection_modus}" == "multidirection_scaled" ] ; then ### 26 directions
		 U_list=" 1,0,0   0,1,0   0,0,1   -1,0,0   0,-1,0   0,0,-1 " ### side planes
		 U_list="${U_list}  -1,1,1   1,-1,-1   1,-1,1   -1,1,-1   1,1,1   -1,-1,-1   1,1,-1   -1,-1,1 "   ### corners
		 U_list="${U_list}   1,1,0   -1,-1,0    1,0,1    -1,0,-1    0,1,1    0,-1,-1    -1,1,0    1,-1,0    -1,0,1    1,0,-1    0,-1,1    0,1,-1]  " ### edges
    elif [ "${multidirection_modus}" == "temperature_increasing" ] ; then 
		U_list=" 1,0,0    0,1,0    0,0,1 "  
		T_list=" 250 300 350 400 450 500 550 600 700 800 900 1000 2000 5000 "
    else
        echo "Error: No_appropriante_Modus_selected: ${multidirection_modus}" 
        echo "Usage: $0 kMC_options_file.xml  state.sql   Modus=[single,singleUx,singleUy,singleUz,increasing/multidirection/multidirection_scaled]  U_ext  "       
        exit 1
    fi
    echo "Selected_Voltages: ${U_list}"
else 
    echo "Usage: $0 kMC_options_file.xml state.sql   Modus=[increasing/multidirection]  U_ext "
    exit 1
fi

if [ $# -ge 4 ]; then
            U_ext="${4}"
            echo "Read external Voltage U_ext:  ${U_ext}"
fi



kMC_options_filename=$( basename ${kMC_options_file} .xml)

#convienience function to change xml option
changeoption(){
    sed -i "s&<${1}.*>.*</${1}>&<${1}>${2}</${1}>&" $3
}
#convienience function to delete xml option
deleteoption(){
 sed -i "s&<${1}.*>.*</${1}>&&" $2
}

### convenience function to get an argument xml option
### e.g. getoption /options/run_kmc_votca_multiple/datafiles/readsystem ${kMC_options_file}
getoption(){
	##xml2 < ${2} | grep "${1}" |  sed 's/.*=//'
	xmllint --xpath "string(//${1})" ${2}
}

echo "Using: ${kMC_options_file}  ${kMC_statefile}   ${multidirection_modus}  ${U_ext} "
echo "For every external voltage N_kMC_random_insertions are selected: ${N_kMC_random_insertions} "



if [[ ! -e ${kMC_options_file} ]] ; then
	echo "Error: kmc options file does not exist in current folder: ${kMC_options_file}"
	exit 1
fi

if [[ "${kMC_options_file: -4}" != ".xml" ]] ; then
	echo "Error; the kmc options file is not a xml file: ${kMC_options_file}"
	exit 1
fi


if [[ ! -e ${kMC_statefile} ]] ; then
	echo "Error: kmc state file does not exist in current folder: ${kMC_statefile}"
	exit 1
fi

if [[ "${kMC_statefile: -4}" != ".sql" ]] ; then
	echo "Error; the kmc state file does not posses a .sql termination: ${kMC_statefile}"
	exit 1
fi

### read the local options from the options.xml file
### keywordlist contains all data from xml file
### keywordlist=$( echo "N_kMC_random_insertions" "use_random_charge_insejection"  "calc_parallel" "use_my_rates_to_votca_sql"  "use_kmc_fortran"  "use_P3HT_charge_insejection"  "use_DIPBI_charge_insejection"  "use_PPDI_charge_insejection"  "use_PBDT_TS1_charge_insejection" "system"  "modify_rates"  "rate_type"  "dump_rates"  "cutoff_dump_rate"  "N_xy_plane" )
### converts key=value to a bash variable.
if [[ "${read_options_from_xml}" == "True" ]] ; then
	
	echo "--- Start reading options from ${kMC_options_file} ---"
	prefix="/options/run_kmc_votca_multiple/datafiles"
	readsystem=$( getoption "${prefix}/readsystem" "${kMC_options_file}" )
	echo "readsystem: $readsystem "
	if [ "${readsystem}" == "1" ] || [[ "${readsystem}" == "True" ]]; then
		### add keys for system datafiles to the keywordlist 
		keywordlist=$( echo " ${keywordlist[@]} gro_filename  neighbourlist_filename   no_box_COM_filename  JAB_filename " )
	fi
	prefix=""
	echo "Keywordlist:  ${keywordlist[@]} "
	for key in ${keywordlist[@]} ; do
			### get value for keyword in xml file
			# xml2 < ${kMC_options_file} | grep /$key | grep -v '!' 
			 ### reduct the tree with  sed 's/.*=//'
			 ### filter comments: grep -v "!" 
			#value=$( xml2 < ${kMC_options_file} | grep /$key |  sed 's/.*=//' | grep -v '!' | head --l 1)
			value=$( xmllint --xpath "string(//$key)" ${kMC_options_file} )
			### key and value not empty
			if [[ "${key}" == *['!'@#\$%^\&*()]* ]] && [[ "${value}" == *['!'@#\$%^\&*()]* ]]  ; then
					echo "Warning, option is not selected in ${kMC_options_file}   skip:  <${key}>  ${value} "
					continue
			fi
			### convert [1,0] to [True,False] 
			### Exclude numeric keywords
			if [[ "${key}" != "N_kMC_random_insertions" ]] && [[ "${key}" != "N_xy_plane" ]] ; then 
					if [[ "${value}" == "1" ]] ; then
							value='True'
					elif [[ "${value}" == "0" ]] ; then
							value='False'
					fi
			fi ### convert [1,0] to [True,False]

			
			if [ ! -z "${key}" ] && [ ! -z "${value}" ] ; then
					### Set keyword to the bash varible
					echo "Use the option: <${key}>   $value"
					 ##declare $(${key})="${value}"
					eval "$( echo ${key##*/})=$value"
			else
				echo " Ignore / Not found in xml: key:   ${key}   value: $value"
			fi #key and value exist?

	done ### read keywords
	echo "-----"
fi ### read_options_from_xml
##############################



if [[ "${use_my_rates_to_votca_sql}" == "True" ]] ; then
	### Data files for:  use_my_rates_to_votca_sql
	
	if [[ "${readsystem}" == "True" ]] ;then
		echo "get system data from grofile"
	
	### DIPBI / P3HT
	elif [[ "${SYSTEM}" == "DIPBI_P3HT_500K" ]] || [[ "${SYSTEM}" == "DIPBI_P3HT" ]]  ; then
		gro_filename='data/equi_2000ps_500k_def2_theta_75_G0.gro'
		neighbourlist_filename='data/new_sorted_neighbours_equi_2000ps_500k_def2_theta_75_G0.ngh'
		no_box_COM_filename='data/no_box_equi_2000ps_500k_def2_theta_75_G0.xyz'
		JAB_filename='data/Datenauswertung_SAB_H_H0_Jab_H0_H0_SAB_L0_L0_Jab_L0_L0_G0.dat'

	### DIPBI / P3HT 700K
	elif [[ "${SYSTEM}" == "DIPBI_P3HT_700K" ]] ; then
		gro_filename='data/equi_1500ps_700k_def2_theta_75_G0.gro'
		neighbourlist_filename='data/new_sorted_neighbours_equi_1500ps_700k_def2_theta_75_G0.ngh'
		no_box_COM_filename='data/no_box_equi_1500ps_700k_def2_theta_75_G0.xyz'
		JAB_filename='data/Datenauswertung_SAB_H_H0_Jab_H0_H0_SAB_L0_L0_Jab_L0_L0_G0_700K_KMC_U0_003_75deg_def2_V1_all.dat'

	### DIPBI / P3HT 900K
	elif [[ "${SYSTEM}" == "DIPBI_P3HT_900K" ]] ; then
		gro_filename='data/equi_2000ps_900k_def2_theta_75_G0.gro'
		neighbourlist_filename='data/new_sorted_neighbours_equi_2000ps_900k_def2_theta_75_G0.ngh'
		no_box_COM_filename='data/no_box_equi_2000ps_900k_def2_theta_75_G0.xyz'
		JAB_filename='data/Datenauswertung_SAB_H_H0_Jab_H0_H0_SAB_L0_L0_Jab_L0_L0_G0_900K_KMC_U0_003_75deg_def2_V1_all.dat'

	### P3HT
	elif [[ "${SYSTEM}" == "P3HT_pur" ]] || [[ "${SYSTEM}" == "P3HT_pur_1000ps_900K_65_G0" ]] ; then	
		gro_filename='data/P3HT_pur_1000ps_900K_65_G0.gro'
		neighbourlist_filename='data/new_sorted_neighbours_P3HT_pur_1000ps_900K_65_G0.ngh'
		no_box_COM_filename='data/no_box_P3HT_pur_1000ps_900K_65_G0.xyz'
		JAB_filename='data/Datenauswertung_SAB_H_H0_Jab_H0_H0_SAB_L0_L0_Jab_L0_L0_G0.dat'
		
	### PPDI PBDT_TS1
	elif [[ "${SYSTEM}" == "PPDI_PBDT_TS1" ]] ; then
		gro_filename='data/out_mischbox_900k_30ns-34ns_60deg_G0.gro'
		neighbourlist_filename='data/new_sorted_neighbours_mischbox_900k_30ns-34ns_60deg_12NN_G0.ngh'
		no_box_COM_filename='data/no_box_mischbox_900k_30ns-34ns_60deg_G0.xyz'
		JAB_filename='data/Datenauswertung_SAB_H_H0_Jab_H0_H0_SAB_L0_L0_Jab_L0_L0_G0.dat'

	### multigeo
	elif [[ "${SYSTEM}" == "DIPBI_P3HT_multigeo" ]] ; then
		###N_frame=1
		#if [ $# -ge 5 ]; then 
		#	N_frame=${5}
		#	### test if N_Geo is integer
		#	re='^[0-9]+$'
		#	if ! [[ ${N_frame} =~ $re ]] ; then
		#		   echo "Error: Reading the frame number N_frame, but it is Not a number:  ${N_frame} " >&2; exit 1
		#	fi		
		#fi

		gro_filename="../indexed_G${N_frame}.gro"
		neighbourlist_filename="../new_sorted_neighbours_G${N_frame}.ngh"
		no_box_COM_filename="../no_box_G${N_frame}.xyz"
		JAB_filename="data/Datenauswertung_SAB_H_H0_Jab_H0_H0_SAB_L0_L0_Jab_L0_L0_G${N_frame}.dat"  ####"state_G${N_frame}.sql"	
		
	else
		echo " Error: System not parametrized, modify ${0} ! "
		exit 1
	fi # System
fi ## use_my_rates_to_votca_sql






if [[ "${rate_type}" == "marcus" ]] || [[ "${rate_type}" == "jortner" ]] || [[ "${rate_type}" == "weiss_dorsey" ]] ; then
	echo "Selected rate_type: ${rate_type} "
	rate_type_list=" ${rate_type} "
elif [[ "${rate_type}" == "rate_type_scan" ]] ; then
	echo "Selected rate_type: ${rate_type} "
	rate_type_list=" marcus  jortner  weiss_dorsey "
else
	echo " The selected rate_type is not supported:  ${rate_type}"
	echo " Select rate_type=[marcus,jortner,weiss_dorsey]"
	exit 1
fi








if [[ "${use_my_rates_to_votca_sql}" == 'True' ]] ; then
	echo "Selected use_my_rates_to_votca_sql Use:   ${gro_filename}  ${neighbourlist_filename} ${no_box_COM_filename}   ${JAB_filename} "  

	for data_file_name in ${gro_filename} ${neighbourlist_filename} ${no_box_COM_filename} ${JAB_filename} 
	do
		if [[ ! -e "${data_file_name}" ]] ; then
			echo "Error: The file does not exist : ${data_file_name}"
			echo "Make sure it is provided for  use_my_rates_to_votca_sql"
			exit 1
		fi
	done

	### Create boxfile.xml for xtp_createStatefile_fromtxt_t_koch08.py
	if [[ ! -e boxfile.xml ]] ; then
						echo '<options>' 		>   boxfile.xml
						echo "	<boxX>13.90228</boxX>"	>>  boxfile.xml
						echo "	<boxY>26.96753</boxY>"	>>  boxfile.xml
						echo "	<boxZ>14.58116</boxZ>"	>>  boxfile.xml
						echo "	<boxXY>0</boxXY>"	>>  boxfile.xml
						echo "	<boxXZ>0</boxXZ>"	>>  boxfile.xml
						echo "	<boxYZ>0</boxYZ>"	>>  boxfile.xml
						echo "</options>"		>>  boxfile.xml
						IFS=' ' read -r -a box_i <<< $( tail --l 1 ${gro_filename})
						changeoption "boxX" ${box_i[0]} boxfile.xml
						changeoption "boxY" ${box_i[1]} boxfile.xml
						changeoption "boxZ" ${box_i[2]} boxfile.xml	
	fi
fi

### add GNU-prallel to PATH
if [[ "${HOSTNAME}" == *"balrog"* ]] || [[ "${HOSTNAME}" == *"kankra"* ]]  ; then
	PATH=/opt/parallel-20180122/bin:$PATH
fi

### get the number of cores available.
NPROCESSORS=$(grep -c ^processor /proc/cpuinfo 2>/dev/null || sysctl -n hw.ncpu || echo "$NUMBER_OF_PROCESSORS")
if ! [[ "${NPROCESSORS}" =~ ^[0-9]+$ ]] ; then ### check if is not an integer.
	echo "Number of processores could not be determined: $NPROCESSORS"
	echo "Use: NCPUS=${NCPUS}"
else
	NCPUS=${NPROCESSORS}
fi



# function zur Multiplikation der Arrayrichtung (wird normiert) mit einer konstanten float Zahl (als Betrag) in Bash ueber ein python skript
# use function u_scale:: result = U_ext*U_array_new/|U_array_new|
# u_scale ${U_ext} ${U_array_new} 
# u_scale 4 [1,2,2]
# result [1.3333, 2.6667, 2.6667] as |{U_array_new| = 9
# Format  U_vec=[0.003,0,0] (with comma seperation)
function u_scale(){  
        U_ext=${1}
        U_arr=${2}
        echo ${U_arr} | awk -F"," -v U_ext="${U_ext}" '{ Ux=$(1); Uy=$(2); Uz=$(3) ; U_norm=1 ; Ux=Ux*U_ext/U_norm ; Uy=Uy*U_ext/U_norm ;  Uz=Uz*U_ext/U_norm ;  printf"%8.3E  %8.3E  %8.3E\n",Ux,Uy,Uz } '
}


function u_scale_norm(){  
        U_ext=${1}
        U_arr=${2}
        echo ${U_arr} | awk -F"," -v U_ext="${U_ext}" '{ Ux=$(1); Uy=$(2); Uz=$(3) ; U_norm=sqrt(Ux*Ux+Uy*Uy+Uz*Uz) ; Ux=Ux*U_ext/U_norm ; Uy=Uy*U_ext/U_norm ;  Uz=Uz*U_ext/U_norm ;  printf"%8.3E  %8.3E  %8.3E\n",Ux,Uy,Uz } '
}


#### ende function  u_scale

##### function kmc parallel for VOTCA 
run_kmc_parallel_line(){
	calcNR=${1} 
	UextNR=${2} 
	kMC_options_file=${3}
	kMC_options_filename=$( basename ${kMC_options_file} .xml)_U${UextNR}_N${calcNR}.xml
	kMC_statefile=${4}  
  	kMC_statefilename=$( basename ${kMC_statefile} .sql )_U${UextNR}.sql  
 #echo  "Use arguments:  $@" 
 #echo  "kmc_run -e kmcmultiple -o kmc_votca_U${UextNR}_N${calcNR}/${kMC_options_filename} -f ${kMC_statefilename}  >> kmc_votca_U${UextNR}_N${calcNR}/kmc_votca_U${UextNR}_N${calcNR}.out" 
  	kmc_run -e kmcmultiple -o kmc_votca_U${UextNR}_N${calcNR}/${kMC_options_filename} -f ${kMC_statefilename}  >> kmc_votca_U${UextNR}_N${calcNR}/kmc_votca_U${UextNR}_N${calcNR}.out 
}
######



##### function kmc_fortran_parallel
run_kmc_fortran_parallel_line(){
	calcNR=${1} 
	UextNR=${2} 
	kMC_options_file=${3}
	kMC_options_filename=$( basename ${kMC_options_file} .xml)_U${UextNR}_N${calcNR}.xml
	kMC_statefile=${4}  
  	kMC_statefilename=$( basename ${kMC_statefile} .sql )_U${UextNR}.sql  
 #echo  "Use arguments:  $@" 
 #echo  "kmc_fortran  KMC_FOR_MULTIPLE_CHARGES_xml   ${kmc_foldername}/${options}  ${kMC_statefile0}  >> ${kmc_foldername}/${kmc_foldername}.out" 
	###### Use kmc_fortran programm for kmc simulations and read the input configurations from the options_votca.xml file. 
	if [[ -e ${kMC_statefilename} ]] ; then 
		cp  ${kMC_statefilename}  kmc_fortran_U${UextNR}_N${calcNR}/
	fi
	
	### read gromacs input filename from options.xml file in: gro_input_filename
	### e.g. <gro_inputfile>out_equi_2000ps_500k_def2_theta_75_G0_V1.gro</gro_inputfile>    <!-- gromacs filename *.gro with morphology data -->
	### gro_input_filename=out_equi_2000ps_500k_def2_theta_75_G0_V1.gro
	gro_input_filename=$( grep "<gro_inputfile>"   kmc_fortran_U${UextNR}_N${calcNR}/${kMC_options_filename}  ) ; gro_input_filename=${gro_input_filename%%</*} ;  gro_input_filename=${gro_input_filename##*>}
	#if [[ -e "${gro_input_filename}" ]] && [ -n ${gro_input_filename} ] ; then
	#	cp  ${gro_input_filename}    kmc_fortran_U${UextNR}_N${calcNR}/
	#fi
	
	cd kmc_fortran_U${UextNR}_N${calcNR}
	if [[ -e ../${gro_input_filename} ]] && [ -n ../${gro_input_filename} ] && [ ! -e ${gro_input_filename} ] ; then
		#cp  ../${gro_input_filename}   .
		ln -s ../${gro_input_filename}  ${gro_input_filename}
	fi
	
	
	###./run_kmc_fortran  KMC_FOR_MULTIPLE_CHARGES_xml     kmc_fortran_U${UextNR}_N${calcNR}/${kMC_options_filename}     ${kMC_statefilename}  >>    kmc_fortran_U${UextNR}_N${calcNR}/kmc_fortran_U${UextNR}_N${calcNR}.out 
	if [[ -e ../run_kmc_fortran ]] ; then
		../run_kmc_fortran            KMC_FOR_MULTIPLE_CHARGES_xml     ${kMC_options_filename}     ${kMC_statefilename}  >>    kmc_fortran_U${UextNR}_N${calcNR}.out 
	elif [[ -e ../charge_transport_package ]] ; then
		../charge_transport_package    KMC_FOR_MULTIPLE_CHARGES_xml    ${kMC_options_filename}     ${kMC_statefilename}  >>    kmc_fortran_U${UextNR}_N${calcNR}.out 
	elif [[ -e charge_transport_package ]] ; then
		charge_transport_package    KMC_FOR_MULTIPLE_CHARGES_xml      ${kMC_options_filename}     ${kMC_statefilename}  >>    kmc_fortran_U${UextNR}_N${calcNR}.out 
	fi
	cd ..
}
######




export -f run_kmc_parallel_line
export -f run_kmc_fortran_parallel_line

declare -i UextNR
UextNR=1

### get initial temperature from optionsfile
temperature=$( grep "<temperature>" ${kMC_options_file} ) ; temperature=${temperature##*<temperature>}  ; temperature=${temperature%%</temperature>*}

### get charge carrier type from option file
carrier=$( grep "<carriertype>"    ${kMC_options_file} ) ;  carrier=${carrier##*<carriertype>} ;  carrier=${carrier%%</carriertype>*} ; carrier=${carrier:0:1} 
if [ "{carrier}" == "h" ] || [ "{carrier}" == "e" ] || [ "{carrier}" == "t" ]|| [ "{carrier}" == "s" ]   ; then
	echo "carrier: ${carrier}"
elif [ -z "${carrier}" ] ; then
	echo "carrier: ${carrier}"
	echo "Error: carrier type in optionfile ${kMC_options_file} was not detemined correcty!"
	error_ctp_run_rates="1" 
	check_exit_1 "${error_ctp_run_rate}" "$((${LINENO}-2))" 
else 
	echo "carrier: ${carrier}"
	echo "Error: carrier type in optionfile ${kMC_options_file} was not detemined correcty!"
	error_ctp_run_rates="1" 
	check_exit_1 "${error_ctp_run_rate}" "$((${LINENO}-2))" 
fi


### Schleife ueber mehrere ratentypen
for rate_type_new in  ${rate_type_list} ; do	
	### Schleife ueber mehrer Spannungen 
	for U_array_new in  ${U_list} ; do
		### Schleife ueber mehrere Temperaturen
		for T_new in  ${T_list} ; do	 
			### rate_types 
			rate_type=${rate_type_new}			
			### function for scaling
			if [[ "${multidirection_modus}" == "multidirection_scaled" ]] ; then
					### use function u_scale U_array_new= U_ext*U_array_new/|U_array_new|
					U_array_new2=$( u_scale_norm ${U_ext} ${U_array_new} )
					if [[ -n "${U_array_new}" ]] && [[ ! -z "${U_array_new}" ]] ; then
						echo "New_Voltage: ${U_array_new} with |U_ext|=${U_ext}"
					else
						echo "New_Voltage: ${U_array_new} with |U_ext|=${U_ext}"
						echo 'Error: u_scale_did_not_work,_check_definition_of_:U_array_new'
						echo "ENDE"
						exit 1
					fi
				
				else 
				U_array_new2=$( u_scale ${U_ext} ${U_array_new}  )
				echo "New_Voltage: ${U_array_new}" 
			fi # ${multidirection_modus} == "multidirection_scaled"
			echo "Use electric field: ${U_array_new} "


			options0=$( basename ${kMC_options_file} .xml)"_U"${UextNR}".xml"
			kMC_statefile0=$( basename ${kMC_statefile} .sql)"_U"${UextNR}".sql"
			cp ${kMC_options_file} ${options0}
			if [[ "${use_my_rates_to_votca_sql}" != "True" ]]; then ### Standard art
				cp ${kMC_statefile} ${kMC_statefile0}
			fi

			### modify external field
			### remove brackets
			U_arr2=$( echo "${U_array_new2}" | sed 's/\[//g'  | sed 's/\]//g' )
				
			### get rid of separating ","
			U_arr3=$( echo "${U_arr2}" | sed 's/,//g' | xargs )

			if [ "${multidirection_modus}" == "temperature_increasing" ] ; then
				temperature="${T_new}"
			fi
				
			### check if the field for in the marcus rate expression was set correctly, otherwise update the rates with the new field in U_arr2
			field_in_rate=$( grep "<field>" ${options0} ) ; field_in_rate=${field_in_rate##*<field>}  ; field_in_rate=${field_in_rate%%</field>*} ; field_in_rate=$( echo ${field_in_rate} | xargs )
			temperature_in_rate=$( grep "<temperature>" ${options0} ) ; temperature_in_rate=${temperature_in_rate##*<temperature>}  ; temperature_in_rate=${temperature_in_rate%%</temperature>*} ; temperature_in_rate=$( echo ${temperature_in_rate}  | xargs )
			echo "tmp: U_arr3: ${U_arr3}    field_in_rate: ${field_in_rate}   temperature_in_rate: ${temperature_in_rate} temperature: ${temperature}   ${rate_type}"
			if [ -z "${U_arr3}" ] || [ -z "${field_in_rate}" ] || [ -z "${temperature_in_rate}" ] || [ -z "${temperature_in_rate}" ] ; then
				echo "Failed: U_arr3: ${U_arr3}    field_in_rate: ${field_in_rate}   temperature_in_rate: ${temperature_in_rate} temperature: ${temperature}"
				echo "Error: At leasts one of the strings is empty. Check calculation input or modules."
				exit 1 
			fi



			if [[ "${modify_rates}" == "True" ]] ; then
				if [[ "${field_in_rate}" != "${U_arr3}" ]] || [[ "${temperature_in_rate}" != "${temperature}" ]] ; then
						echo "Updating the rates for kmc calculation for the external field: ${U_arr3}  [V/m]  and temperature: ${temperature}  [K]"
						### <!-- field x y z -->
						changeoption "field" "${U_arr3}" ${options0}
						### <!-- T [K] -->
						changeoption "temperature" "${temperature}" ${options0}

						if [[ "${use_my_rates_to_votca_sql}" != "True" ]]; then ### Standard art			
							echo "#7) Update Rates"
							ctp_run -e rates -o ${options0} -f ${kMC_statefile0} -t ${NCPUS}
							error_ctp_run_rates="$?" 
							check_exit_1 "${error_ctp_run_rate}" "$((${LINENO}-2))" 
						else ### use_my_rates_to_votca_sql
							echo "#7) Create new sql file"
							### prepare data for statefile
							if [[ -e ./my_rates_to_votca_sql.out ]] ; then
								echo " ./my_rates_to_votca_sql.out ${gro_filename}  ${neighbourlist_filename} ${no_box_COM_filename} data_collection_$( basename ${kMC_statefile0} .sql).dat  ${JAB_filename}  ${U_arr3}   ${temperature}  ${rate_type} "
								./my_rates_to_votca_sql.out        ${gro_filename}  ${neighbourlist_filename} ${no_box_COM_filename} data_collection_$( basename ${kMC_statefile0} .sql).dat  ${JAB_filename}  ${U_arr3}   ${temperature}  ${rate_type}  
								error_my_rates_to_votca_sql="$?" 
								check_exit_1 "${error_my_rates_to_votca_sql}" "$((${LINENO}-2))" 
							elif [[ -e ./charge_transport_package ]] ; then
								echo " ./charge_transport_package  rates_to_VOTCA_sqlfile  ${gro_filename}  ${neighbourlist_filename} ${no_box_COM_filename} data_collection_$( basename ${kMC_statefile0} .sql).dat  ${JAB_filename}  ${U_arr3}   ${temperature}  ${rate_type}  ${kMC_options_file}"
								./charge_transport_package         rates_to_VOTCA_sqlfile  ${gro_filename}  ${neighbourlist_filename} ${no_box_COM_filename} data_collection_$( basename ${kMC_statefile0} .sql).dat  ${JAB_filename}  ${U_arr3}   ${temperature}  ${rate_type}  ${kMC_options_file}
								error_my_rates_to_votca_sql="$?" 
								check_exit_1 "${error_my_rates_to_votca_sql}" "$((${LINENO}-2))" 
							else
								echo " charge_transport_package    rates_to_VOTCA_sqlfile  ${gro_filename}  ${neighbourlist_filename} ${no_box_COM_filename} data_collection_$( basename ${kMC_statefile0} .sql).dat  ${JAB_filename}  ${U_arr3}   ${temperature}  ${rate_type}   ${kMC_options_file}"
								       charge_transport_package    rates_to_VOTCA_sqlfile  ${gro_filename}  ${neighbourlist_filename} ${no_box_COM_filename} data_collection_$( basename ${kMC_statefile0} .sql).dat  ${JAB_filename}  ${U_arr3}   ${temperature}  ${rate_type}   ${kMC_options_file}
								error_my_rates_to_votca_sql="$?" 
								check_exit_1 "${error_my_rates_to_votca_sql}" "$((${LINENO}-2))" 
							fi
								
							if [[ -e votca_pairfile.dat ]] && [[ -e votca_segfile.dat ]]  ; then
								 #### create new statefile  -c create    -t h hole transport
								echo "xtp_createStatefile_fromtxt_t_koch08.py -f ${kMC_statefile0} -b boxfile.xml -p votca_pairfile.dat -s votca_segfile.dat -t ${carrier} -c"
								      xtp_createStatefile_fromtxt_t_koch08.py -f ${kMC_statefile0} -b boxfile.xml -p votca_pairfile.dat -s votca_segfile.dat -t ${carrier} -c
								error_my_rates_to_votca_sql="$?" 
								check_exit_1 "${error_my_rates_to_votca_sql}" "$((${LINENO}-2))" 
								mv votca_pairfile.dat "votca_pairfile_U${UextNR}.dat"
								mv votca_segfile.dat  "votca_segfile_U${UextNR}.dat"
							else
								echo "Error: Can not find new votca_pairfile.dat or votca_segfile.dat! EXIT"
								exit 1
							fi
						fi   ### use_my_rates_to_votca_sql

						if [[ "${dump_rates}" == "True" ]] ; then 
							echo "#8) Set all rates to 1.0 s**-1 for rate <  ${cutoff_dump_rate}   "
							sqlite3 ${kMC_statefile0} " UPDATE pairs SET rate12h=${cutoff_dump_rate}  WHERE rate12h < ${cutoff_dump_rate} "
							sqlite3 ${kMC_statefile0} " UPDATE pairs SET rate21h=${cutoff_dump_rate}  WHERE rate21h < ${cutoff_dump_rate} "
							sqlite3 ${kMC_statefile0} " UPDATE pairs SET rate12e=${cutoff_dump_rate}  WHERE rate12e < ${cutoff_dump_rate} "
							sqlite3 ${kMC_statefile0} " UPDATE pairs SET rate21e=${cutoff_dump_rate}  WHERE rate21e < ${cutoff_dump_rate} "
						fi ### dump ratee

					fi ### check <field> in rates
			fi ## "${modify_rates}" 


			echo "_U${UextNR}  ${U_arr3}  _Nmax ${N_kMC_random_insertions} Start  " $( pwd ) >> Uext_iterations_list.dat
			echo "_U${UextNR}  ${U_arr3}  _Nmax ${N_kMC_random_insertions} Start  " $( pwd ) 
			### erstelle ordnerstruktur
			for i in ` seq 1 ${N_kMC_random_insertions} ` ; do
				if [[ "${use_kmc_fortran}" == 'True' ]] ; then 
					kmc_foldername="kmc_fortran_U"${UextNR}"_N"${i}
				else ### VOTCA
					kmc_foldername="kmc_votca_U"${UextNR}"_N"${i}
				fi ### kmc_fortran or kmc_votca 
				
				options=$( basename ${options0} .xml)"_N"${i}".xml"
				if [[ ! -e ${kmc_foldername} ]] ; then
					mkdir ${kmc_foldername}
				fi
				### create new optoions file
				cp ${options0} ${options} 
				### change random initial seed
				rand_seed=$( shuf -i 2000-65000 -n 1)
				rand_seed=${i}${rand_seed}
				echo "Modify seed in: ${options} ${rand_seed}"
				changeoption "seed" ${rand_seed} ${options}
				echo "tmp: test ${options}"	
				if [ "${use_random_charge_insejection}" == "True" ] ; then
					xml_injectionmethod=$( getoption  "kmcmultiple/injectionmethod"  ${options} ) 
					if [[ "${xml_injectionmethod}" == "injection_segment_type" ]] ; then
						echo "Select injection at injection_segment_type: "$(  getoption  "kmcmultiple/injection_segment_type"  ${options}) 
					else ### standard random insertion!
						changeoption "injectionmethod"  "random"   ${options}
					fi
				elif [ "${use_P3HT_charge_insejection}" == "True" ] || [ "${use_DIPBI_charge_insejection}" == "True" ] || [ "${use_PPDI_charge_insejection}" == "True" ] || [ "${use_PBDT_TS1_charge_insejection}" == "True" ] ; then
				### use P3HT charge injection at PH or PM sites
						### get list of injection node ids: list_injection_ids
						if [ ${i} == 1 ] ; then 
							echo "List of selected injection nodes"
							if [ "${use_P3HT_charge_insejection}" == "True" ] ; then
									sqlite3 ${kMC_statefile0}  " SELECT id,name from segments WHERE name LIKE '%PH%S%' OR '%PM%S%'  "  | head --l ${N_kMC_random_insertions} 
									list_injection_ids=( $( sqlite3 ${kMC_statefile0}  " SELECT id from segments WHERE name LIKE '%PH%S%' OR '%PM%S%'  "  | head --l ${N_kMC_random_insertions}  ) )
							elif  [ "${use_DIPBI_charge_insejection}" == "True" ] ; then
									sqlite3 ${kMC_statefile0}  " SELECT id,name from segments WHERE name LIKE '%DIPBI%' OR '%DPB%'  "  | head --l ${N_kMC_random_insertions} 
									list_injection_ids=( $( sqlite3 ${kMC_statefile0}  " SELECT id from segments WHERE name LIKE '%DIPBI%' OR '%DPB%'  "  | head --l ${N_kMC_random_insertions}  ) )
							elif  [ "${use_PPDI_charge_insejection}" == "True" ] ; then
									sqlite3 ${kMC_statefile0}  " SELECT id,name from segments WHERE name LIKE '%PPDI%%'  "  | head --l ${N_kMC_random_insertions} 
									list_injection_ids=( $( sqlite3 ${kMC_statefile0}  " SELECT id from segments WHERE name LIKE '%PPDI%%'  "  | head --l ${N_kMC_random_insertions}  ) )
							elif  [ "${use_PBDT_TS1_charge_insejection}" == "True" ] ; then
									sqlite3 ${kMC_statefile0}  " SELECT id,name from segments WHERE name LIKE '%A%B%S%'  "  | head --l ${N_kMC_random_insertions} 
									list_injection_ids=( $( sqlite3 ${kMC_statefile0}  " SELECT id from segments WHERE name LIKE '%A%B%S%'  "  | head --l ${N_kMC_random_insertions}  ) )
							else ### create list with sequence 1 to number of charges
									list_injection_ids=( $( seq 1 ${N_kMC_random_insertions} ) )
							fi
							echo " Liste ist_injection_id: ${list_injection_id[@]} "
						fi ### first iteration get list_injection_ids
						
						### switch injectionmethod from random to injectionnode
						changeoption    "injectionmethod"  "injectionnode"            ${options}
						### insert id of selceted start kmc node
						start_injection_id=${list_injection_ids[$((${i}-1))]} 
						start_injection_name="START_${i}"
						sqlite3 ${kMC_statefile0} " UPDATE segments SET name='${start_injection_name}' WHERE _id='${start_injection_id}' "
						
						if [[ "${use_kmc_fortran}" == 'True' ]] ; then ### kmc_fortran injection selection via id 
							changeoption    "injection"    "${start_injection_id}"       ${options}
							echo "changeoption"    "injection"    "${start_injection_id}"       ${options}
						else ### Votca injection selection via name
							changeoption    "injection"    "${start_injection_name}"       ${options}
							echo " changeoption "   "injection"        "START=${list_injection_ids[$((${i}-1))]}"    ${options}
						fi
						 
				fi
				
				### modify external field
				### remove brackets
				U_arr2=$( echo "${U_array_new2}" | sed 's/\[//g'  | sed 's/\]//g' )
				
				### get rid of separating "," and trim empty spaces with xargs
				U_arr3=$( echo "${U_arr2}" | sed 's/,//g' | xargs )

				if [ "${multidirection_modus}" == "temperature_increasing" ] ; then
					temperature="${T_new}"
				fi
				
				### check if the field for in the marcus rate expression was set correctly, otherwise update the rates with the new field in U_arr2
				field_in_rate=$( grep "<field>" ${options} ) ; field_in_rate=${field_in_rate##*<field>}  ; field_in_rate=${field_in_rate%%</field>*} ; field_in_rate=$( echo ${field_in_rate} | xargs )
				temperature_in_rate=$( grep "<temperature>" ${options} ) ; temperature_in_rate=${temperature_in_rate##*<temperature>}  ; temperature_in_rate=${temperature_in_rate%%</temperature>*} ; temperature_in_rate=$( echo ${temperature_in_rate}  | xargs )
				echo "tmp: U_arr3: ${U_arr3}    field_in_rate: ${field_in_rate}   temperature_in_rate: ${temperature_in_rate} temperature: ${temperature}   ${rate_type}"
				if [ -z "${U_arr3}" ] || [ -z "${field_in_rate}" ] || [ -z "${temperature_in_rate}" ] || [ -z "${temperature_in_rate}" ] ; then
					echo "Failed: U_arr3: ${U_arr3}    field_in_rate: ${field_in_rate}   temperature_in_rate: ${temperature_in_rate} temperature: ${temperature}"
					echo "Error: At leasts one of the strings is empty. Check calculation input or modules."
					exit 1 
				fi

				echo "modify_rates: ${modify_rates}" 

				if [[ "${modify_rates}" == "True" ]] ; then
					if [[ "${field_in_rate}" != "${U_arr3}" ]] || [[ "${temperature_in_rate}" != "${temperature}" ]] ; then
						echo "Updating the rates for kmc calculation for the external field: ${U_arr3}  [V/m]  and temperature: ${temperature}  [K]"
						### <!-- field x y z -->
						changeoption "field" "${U_arr3}" ${options}
						### <!-- T [K] -->
						changeoption "temperature" "${temperature}" ${options}
					

						if [[ "${use_my_rates_to_votca_sql}" != "True" ]]; then ### Standard art			
							echo "#7) Update Rates"
							ctp_run -e rates -o ${options} -f ${kMC_statefile0} -t ${NCPUS}
							error_ctp_run_rates="$?" 
							check_exit_1 "${error_ctp_run_rate}" "$((${LINENO}-2))" 
						else ### use_my_rates_to_votca_sql
							echo "#7) Create new sql file"
							### prepare data for statefile
							if [[ -e ./my_rates_to_votca_sql.out ]] ; then
								echo "./my_rates_to_votca_sql.out ${gro_filename}  ${neighbourlist_filename} ${no_box_COM_filename} " data_collection_$( basename ${kMC_statefile0} .sql).dat " ${JAB_filename}  ${U_arr3}   ${temperature}  ${rate_type} "
								./my_rates_to_votca_sql.out       ${gro_filename}  ${neighbourlist_filename} ${no_box_COM_filename}   data_collection_$( basename ${kMC_statefile0} .sql).dat   ${JAB_filename}  ${U_arr3}   ${temperature}  ${rate_type}
								error_my_rates_to_votca_sql="$?" 
								check_exit_1 "${error_my_rates_to_votca_sql}" "$((${LINENO}-2))" 
							elif [[ -e ./charge_transport_package ]] ; then
								echo "./charge_transport_package  rates_to_VOTCA_sqlfile  ${gro_filename}  ${neighbourlist_filename} ${no_box_COM_filename} " data_collection_$( basename ${kMC_statefile0} .sql).dat " ${JAB_filename}  ${U_arr3}   ${temperature}  ${rate_type}  ${kMC_options_file}"
								./charge_transport_package        rates_to_VOTCA_sqlfile  ${gro_filename}  ${neighbourlist_filename} ${no_box_COM_filename}   data_collection_$( basename ${kMC_statefile0} .sql).dat   ${JAB_filename}  ${U_arr3}   ${temperature}  ${rate_type}  ${kMC_options_file}
								error_my_rates_to_votca_sql="$?" 
								check_exit_1 "${error_my_rates_to_votca_sql}" "$((${LINENO}-2))" 
							else
								echo "  charge_transport_package   rates_to_VOTCA_sqlfile  ${gro_filename}  ${neighbourlist_filename} ${no_box_COM_filename} " data_collection_$( basename ${kMC_statefile0} .sql).dat " ${JAB_filename}  ${U_arr3}   ${temperature}  ${rate_type}  ${kMC_options_file}"
								        charge_transport_package   rates_to_VOTCA_sqlfile  ${gro_filename}  ${neighbourlist_filename} ${no_box_COM_filename}   data_collection_$( basename ${kMC_statefile0} .sql).dat   ${JAB_filename}  ${U_arr3}   ${temperature}  ${rate_type}  ${kMC_options_file}
								error_my_rates_to_votca_sql="$?" 
								check_exit_1 "${error_my_rates_to_votca_sql}" "$((${LINENO}-2))" 
							fi ### rates to votca.sql
							
							
							if [[ -e votca_pairfile.dat ]] && [[ -e votca_segfile.dat ]]  ; then
								 #### create new statefile  -c create    -t h hole transport
								echo "xtp_createStatefile_fromtxt_t_koch08.py -f ${kMC_statefile0} -b boxfile.xml -p votca_pairfile.dat -s votca_segfile.dat -t ${carrier} -c"
								      xtp_createStatefile_fromtxt_t_koch08.py -f ${kMC_statefile0} -b boxfile.xml -p votca_pairfile.dat -s votca_segfile.dat -t ${carrier} -c
								mv votca_pairfile.dat "votca_pairfile_U${UextNR}.dat"
								mv votca_segfile.dat "votca_segfile_U${UextNR}.dat"
							else
								echo "Error: Can not find new votca_pairfile.dat or votca_segfile.dat! EXIT"
								exit 1
							fi
						fi   ### use_my_rates_to_votca_sql


						if [[ "${dump_rates}" == "True" ]] ; then
							echo "#8 dump rates:  cutoff: ${cutoff_dump_rate}"
							sqlite3 ${kMC_statefile0} " UPDATE pairs SET rate12h=${cutoff_dump_rate}  WHERE rate12h < ${cutoff_dump_rate} "
							sqlite3 ${kMC_statefile0} " UPDATE pairs SET rate21h=${cutoff_dump_rate}  WHERE rate21h < ${cutoff_dump_rate} "
							sqlite3 ${kMC_statefile0} " UPDATE pairs SET rate12e=${cutoff_dump_rate}  WHERE rate12e < ${cutoff_dump_rate} "
							sqlite3 ${kMC_statefile0} " UPDATE pairs SET rate21e=${cutoff_dump_rate}  WHERE rate21e < ${cutoff_dump_rate} "
						fi ### dump rates

					fi ### check <field> in rates
				fi ### "${modify_rates}"

				### split into components
				IFS=', ' read -r -a U_i <<< "$U_arr2"

				changeoption "fieldX" ${U_i[0]} ${options}
				changeoption "fieldY" ${U_i[1]} ${options}
				changeoption "fieldZ" ${U_i[2]} ${options}
				mv ${options} ${kmc_foldername}

				#### calculation in serial mode
				if [[ "${calc_parallel}" != "True" ]] ; then
					if [[ "${use_kmc_fortran}" == 'True' ]] ; then ### run kmc_fortran programm 
							if [[ -e ./run_kmc_fortran ]] ; then
								echo "./run_kmc_fortran  KMC_FOR_MULTIPLE_CHARGES_xml   ${kmc_foldername}/${options}  ${kMC_statefile0} "
								./run_kmc_fortran        KMC_FOR_MULTIPLE_CHARGES_xml   ${kmc_foldername}/${options}  ${kMC_statefile0}  >> ${kmc_foldername}/${kmc_foldername}.out
							elif [[ -e ./charge_transport_package ]] ; then
								echo "./charge_transport_package  KMC_FOR_MULTIPLE_CHARGES_xml   ${kmc_foldername}/${options}  ${kMC_statefile0} "
								./charge_transport_package        KMC_FOR_MULTIPLE_CHARGES_xml   ${kmc_foldername}/${options}  ${kMC_statefile0}  >> ${kmc_foldername}/${kmc_foldername}.out
							else
								echo " charge_transport_package  KMC_FOR_MULTIPLE_CHARGES_xml   ${kmc_foldername}/${options}  ${kMC_statefile0} "
								 charge_transport_package        KMC_FOR_MULTIPLE_CHARGES_xml   ${kmc_foldername}/${options}  ${kMC_statefile0}  >> ${kmc_foldername}/${kmc_foldername}.out
							fi
					else ### run VOTCA 
						###	--timeout 600 
						echo " Start kmc single trajectory run"
						kmc_run -e kmcmultiple -o ${kmc_foldername}/${options} -f ${kMC_statefile0}  >> ${kmc_foldername}/${kmc_foldername}.out
					fi ### run VOTCA 
				fi ### calc_parallel /= True
			done 

			if [[ "${calc_parallel}" == "True" ]] ; then
				### Prepare kmc_run_parallel_file: run_kmc_prepare.sh
				if [[ "${use_kmc_fortran}" == 'True' ]] ; then ### run kmc_fortran programm 
					parallel  --env ${NCPUS} -j $((${NCPUS}-2)) --delay 10.0  --no-run-if-empty  run_kmc_fortran_parallel_line {1} {2} ${kMC_options_file} ${kMC_statefile}  ::: ` seq 1 ${N_kMC_random_insertions} ` ::: ${UextNR}
					echo " loop done for N Startconfigurations and ${N_kMC_random_insertions} and ${U_ext} "
				
				else ### run VOTCA 
						###  USE VOTCA kmc programm 
						### --timeout 900%  ### timeout is given in seconds --timeout 1000000
						### high delay time, to avoid double access to sql files 
						parallel --env ${NCPUS} -j $((${NCPUS}-2)) --delay 3.0  --no-run-if-empty  run_kmc_parallel_line {1} {2} ${kMC_options_file} ${kMC_statefile}  ::: ` seq 1 ${N_kMC_random_insertions} ` ::: ${UextNR}

						echo " loop done for N Startconfigurations and ${N_kMC_random_insertions} and ${U_ext} "
						### END loop N_Geo
				fi ### kmc_fortran OR run VOTCA 

						### Test if mobility analysis program is available and execute it, if it is possible.
				if [ -x "$(command -v  average_votca_mobilities.sh )" ] ; then
							data_mobility=$( average_votca_mobilities.sh ${UextNR} )
							N_mobilities_done=$( grep  "Overall average mobilit"  *_U${UextNR}_*/*.out | wc -l )
							echo  "${data_mobility}          ${U_arr3}    "  $( pwd ) "T= ${temperature_in_rate}    ${rate_type}    ${N_mobilities_done}  " >> results_mobility_steps.dat
				fi ### average_votca_mobilities.sh
				


				echo "_U${UextNR}  ${U_arr3}  _Nmax ${N_kMC_random_insertions} Done " $( pwd ) "T= ${temperature_in_rate}   ${rate_type}"
				### increment
				UextNR=$((UextNR+1))  
			fi ### calc_parallel
			
			
		
		done ### loop for temperatures 
	done ### loop extermal fields
done ### loop for different rate_types



### 			find "kmc_votca_U"${UextNR}"_"* -type f | parallel --verbose --env $NCPUS -j $((${NUS}-2)) --delay 2.0 --timeout 600 --no-run-if-empty ' cd {//} && bash run_kmc_parallel_line.sh {1} {2} ' ::: ` seq 1 ${N_kMC_random_insertions} ` ::: ${UextNR}	
#find kmc_votca_* -type f | parallel -j2 'cd {//} && kmc_run -e kmcmultiple -o kmc_votca_{}/${kMC_options_filename}_{}.xml -f ${kMC_statefilename} > kmc_votca_{}/kmc_votca_{}.out ' ::: `seq 1 ${NGeo}`
#find kmc_votca_* -type f | parallel -j2 'cd {//} && bash tmp_parallel.sh {} > kmc_votca_{}/kmc_votca_{}.out ' ::: `seq 1 ${NGeo}`
#find kmc_votca_* -type f | parallel -j2 'cd {//} && bash tmp_parallel.sh {} > kmc_votca_{}/kmc_votca_{}.out ' ::: `seq 1 ${NGeo}`
#run_kmc_prepare "${outdir}" "${kMC_VOTCA_startline}" {2} > kmc_votca_{}.out
#parallel -j 2 --results $( run_kmc_prepare "${outdir}" "${kMC_VOTCA_startline}" {1} ) ::: `seq 1 ${NGeo}`
#parallel -j 2 $( run_kmc_prepare ${outdir}_{1}/votca_kmc_calc_{1}  ${kMC_VOTCA_startline} {1} ) ::: `seq 1 ${NGeo}`

#echo "Use:  kmc_run -e kmcmultiple -o kmc_votca_{}/${kMC_options_filename}_{}.xml -f ${kMC_statefilename} " 


U=$((UextNR-1))  

if [[ "${multidirection_modus}" == "multidirection" ]] || [[ "${multidirection_modus}" == "multidirection_scaled" ]] || [[ "${multidirection_modus}" == "xy-plane" ]]  ; then
   filename="Overall_average_mobilities_data_U_all.dat"
   grep "Overall average mobilit"  *_U*/*.out > ${filename}
   N_mobilities_done=$( grep "Overall average mobilit"  *_U*/*.out | wc -l )
else
   filename="Overall_average_mobilities_data_U${U}.dat"
   grep "Overall average mobilit"  *_U${U}_*/*.out > ${filename}  
   N_mobilities_done=$( grep "Overall average mobilit"  *_U${U}_*/*.out | wc -l )
fi
echo "N_mobilities_done:  ${N_mobilities_done}  "

mittelwert_ref=$( cat ${filename} | sed 's/ Overall average mobilities//g' | sed 's/ Overall average mobility//g'  | awk  -F"[=/ ]"  'BEGIN{ mu_x=0.0 ; mu_y=0.0 ; mu_z=0.0 ; var=0.0 ; ix=0 ; iy=0 ; iz=0 ; sigma=0.0 ; j=0}
                {  inputline=$0
		   for( i =1; i <=NF ; i++){
            if ( index($(i),"<muX>") != 0){ 
                           X=$(i+1)
                           mu_x=mu_x+X
                           ix=ix+1
                           #print "X:",X,mu_x,ix,"    ",inputline	  
                        }
   
            if ( index($(i),"<muY>") != 0){
                           Y=$(i+1)
                           mu_y=mu_y+Y
                           iy=iy+1
                           #print "Y:",mu_y,iy            
                     }
            if ( index($(i),"<muZ>") != 0){ 
                           Z=$(i+1)
                           mu_z=mu_z+Z
                           iz=iz+1
                           #print "Z:",z,mu_z 
                        }
            } ### End split inputstring to numbers

            j=j+1 
            mittelwert=(mu_x^2+mu_y^2+mu_z^2)^(0.5) 
            
            ### PRINT OUTPUT
            #if( ix != 0) { printf " %.4E ",mu_x/ix} else {printf " %.4E ",mu_x}
            #if( iy != 0) { printf " %.4E ",mu_y/iy} else {printf " %.4E ",mu_y}
            #if( iz != 0) { printf " %.4E ",mu_z/iz} else {printf " %.4E ",mu_z}
            #if( j  != 0) { printf " %.4E \n",mittelwert/j} else {printf " %.4E ",mittelwert}
            #printf "%.4e  %.4e  %.4e  %.4e \n",mu_x,mu_y,mu_z,mittelwert
            
         } 
		   END{
            if( ix != 0 ){ mu_x=mu_x/ix}
            if( iy != 0 ){ mu_y=mu_y/iy}			   
            if( iz != 0 ){ mu_z=mu_z/iz}
            mittelwert=sqrt(mu_x^2+mu_y^2+mu_z^2)  
            printf "%.4e  %.4e  %.4e  %.4e \n",mu_x,mu_y,mu_z,mittelwert}' )


#echo "${mittelwert_ref}"
read -a mittelwert_arr <<< ${mittelwert_ref}
#echo "ARRAY:  ${mittelwert_arr[0]}  ${mittelwert_arr[1]}  ${mittelwert_arr[2]}  ${mittelwert_arr[3]}"

mu_and_sigma=$( cat ${filename} | sed 's/ Overall average mobilities //g' |  sed 's/ Overall average mobility//g' | awk  -F"[=/ ]"  -v  mu_x_mean="${mittelwert_arr[0]}" -v mu_y_mean="${mittelwert_arr[1]}" -v mu_z_mean="${mittelwert_arr[2]}" -v mu_ges="${mittelwert_arr[3]}"  'BEGIN{ mu_x=0.0 ; mu_y=0.0 ; mu_z=0.0 ; var_x=0.0 ; var_y=0.0 ; var_z=0.0 ; ix=0 ; iy=0 ; iz=0 ; i=0 ; sigma_x=0.0 ; sigma_y=0.0 ; sigma_z=0.0 ; j=0  }
         {  inputline=$0
		   #print inputline 
		   for( i =1; i <=NF ; i++){
            if ( index($(i),"<muX>") != 0){ 
                           X=$(i+1)
                           mu_x=mu_x+X
                           ix=ix+1
                           var_x=var_x+(X-mu_x_mean)^2 
			   if( ix > 1 ){ sigma_x=sqrt(var_x/(ix-1)) }
                           #if( ix > 1 ){ print "X:",ix,mu_x/ix,sigma_x }
		   				}
            if ( index($(i),"<muY>") != 0){
                           Y=$(i+1)
                           mu_y=mu_y+Y
                           iy=iy+1
                           var_y=var_y+(Y-mu_y_mean)^2 
			   if( iy > 1 ){ sigma_y=sqrt(var_y/(iy-1)) }
                           #if(iy > 1 ){ print "Y:",iy,mu_y/iy,sigma_y}             
						  }
                    
            if ( index($(i),"<muZ>") != 0){ 
                           Z=$(i+1)
                           mu_z=mu_z+Z
                           iz=iz+1
                           var_z=var_z+(Z-mu_z_mean)^2
                           if( iz > 1 ){ sigma_z=sqrt(var_z/(iz-1)) }
                           #if( iz > 1 ){ print "Z:",iz,mu_z/iz,sigma_z}
						   }
         } ### End split inputstring to numbers
         
		   j=j+1 
         if (j > 1) {
            if( ix > 1 ){ sigma_x=sqrt(var_x/(ix-1)) }
            if( iy > 1 ){ sigma_y=sqrt(var_y/(iy-1)) }
            if( iz > 1 ){ sigma_z=sqrt(var_z/(iz-1)) }
            }
                  
		   mittelwert=sqrt(mu_x^2+mu_y^2+mu_z^2)    
         	   sigma_ges=sqrt(sigma_x^2+sigma_y^2+sigma_z^2)
		   #if (j > 1) {print "mu_ges:",j,mittelwert,sigma_ges }
         } 
		   END{ 
            sigma_ges=sqrt(sigma_x^2+sigma_y^2+sigma_z^2)
            printf "%.4e %-6s %.4e    %.4e %-6s %.4e    %.4e %-6s %.4e    %.4e %-6s %.4e \n",mu_x_mean,"   ",sigma_x, mu_y_mean,"   ",sigma_y, mu_z_mean,"   ",sigma_z, mu_ges,"   ",sigma_ges 
            }'  )    

echo ${mu_and_sigma}   

### theta, phi line for polar plot
###for i in $U_list ; do echo $i | sed 's/\[//g' | sed 's/\]//g' | awk -F"[,]" '{  r=($1^2+$2^2+$3^2)^(0.5) ;  x=($3/r) ;  theta=atan2((1.-x^2)^0.5,x) ; phi=atan2($2,$1) ; print theta,"     ",phi }' ; done




