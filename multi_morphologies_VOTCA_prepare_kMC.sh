#!/bin/bash
### Pragamm zur Durchfuehrung von kMC-Rechnungen mit veraenderung der zugrundeliegenden Geometrie

 if [ $# -eq 0 ]; then
    echo "Usage: $0 < kMC Basisordner > "
    exit 1
 fi

 if [ ! -d "$1" ]; then
    echo "directory $1 doesn't exist!"
    exit 1
 fi

CURDIR=`pwd`
kmc_programmname='run_kmc_votca_multiple_V3.sh'

### default=15
r_cut_lambda_out_ang=30

kMC_Basisordner='dipbi_p3ht_900K_6ns_to_300K_500ps_G0.gro'
kMC_Ordner_old=${kMC_Basisordner} ## %%_G*}
indexed_gro_old="dipbi_p3ht_900K_6ns_to_300K_500ps_G1_to_G12.gro" 
typeset -i Geo_Start=1  ## Variable um Geometrieschritt in der *.gro file zu zaehlen  G${Geo_Start} e.g: G1
typeset -i N_Geo_Stepsmax=1

### critical_cutoff_angle_theta_SCCS  [in degree]  
theta_cutoff=75  

#### default: start_VOTCA_setup='False'
start_VOTCA_setup='False'
#### evaluate data for reference 32mer     default: calc_data_32mer="False" 
calc_data_32mer="True" 
start_VOTCA_setup_32mer="True" 
### dumps the data from the *.sql file to readable files.
dump_data_32mer="False"

### run_kmc_votca_multiple.sh script.
run_KMC_local="False"
### use sub-KMC-VOTCA   for submission to queue        default:        sub_KMC="False"
sub_KMC="False"
### local wait functions, is not recommended           default: wait_for_queue="False"
wait_for_queue="False"

### use reference data for theta=0.0 deg  default:use_theta0_reference_data='True'
use_theta0_reference_data='True'

### Cutoff for reference data at theta=0.0, so one can find sufficient neighbours in a subsequent search with theta_cutoff=xx in  neighbour_list_with_reference_resid.py
### r_cutoff_lambda_for_reference_data [Ang]
r_cutoff_lambda_for_reference_data=30.0

### Path to excess all data files for DIPRO calculations.
path_to_OUZO='../get_started_OUZO_DIPRO_PM3'
### setup_gro_to_DIPRO with OUZO programm. Default: setup_DIPRO='False'
setup_DIPRO='False'
### g09 sting calculation method; default g09_calculation_method_string=' PM3 ' ; or use B3LYP_6-31Gs ; DO not use "*" or other special characters.
g09_calculation_method_string='PM3'
### default g09_DIPRO_method_line=' PM3 SCF(XQC,MaxCycles=2024,MaxConventionalCycles=1024) '
g09_DIPRO_method_line='   PM3   SCF(XQC,MaxCycles=4024,MaxConventionalCycles=2024) '

### Add a symbolic link to the folder; so QC_FILES and MP_FILES are available.
add_symbolic_link_to_folder="False"


### setup local VOTCA environment
if [ "${HOSTNAME}" == "r06m01" ] ; then
	echo "palma HOSTNAME:  ${HOSTNAME}"
	source /home/t/t_koch08/bin/start_config_VOTCA
	executable_folder='/home/t/t_koch08/bin/def2'
	grom_pfad='/home/t/t_koch08/my_votca_scratch/multi_geo'
else
	echo "HOSTNAME: ${HOSTNAME}"
	source /home/t/t_koch08/bin/startconfig_votca_ctp
	executable_folder='/home/t/t_koch08/bin/def2/'
	grom_pfad='/home/t/t_koch08/votca_my_systems/multi_geo'   # Pfad fuer die gromacs Daten *.trp.xtc ; *.gro ; *.tpr
fi

##########
# Initialisierung
last_kMC_run_ok='True' 

###############################################################################################
#### START  define local functions
###############################################################################################

###############################################################################################
#### Function to check if another script terminated with an error and exits this script
###############################################################################################

	function check_exit_1(){
		error_in="$1" 
		error_in_line="$2" 
		if [[ "${error_in}" == "1" ]] ; then
			echo "GLOBAL ERROR: set to 1"
			echo "IN: ${error_in_line}: TERMINATE EXECUTION OF ${0} IN "$( pwd ) 
			exit 1
		fi
	}
###############################################################################################
	

###############################################################################################	
### Function to dump data from VOTCA *sql data base into fortran readable data files and move them to the data folder, if availanle
### e.g.  dump_data_from_sql_file  statefile_sql   filename          modus[all,l_out,E_out,iZINDO]
### e.g.  dump_data_from_sql_file  state.sql       dipbi_p3ht_system all
function  dump_data_from_sql_file(){
      statefile_sql="${1}"
      filename="${2}"
      dump_modus="${3}"
      
      if [[ "${dump_modus}" == "all" ]] || [[ "${dump_modus}" == "l_out" ]] || [[ "${dump_modus}" == "E_out" ]] || [[ "${dump_modus}" == "iZINDO" ]] ; then
		 echo "Use dump_data_from_sql_file:  ${dump_modus}"
      else
         echo "Error: modus[all,l_out,E_out,iZINDO] is not supported in dump_data_from_sql_file: ${dump_modus} "
         exit 1
      fi       
         
      #echo "${1}  ${2}  ${3}  "
      if [ ! -e "${statefile_sql}" ] ; then
         echo "Error: ${statefile_sql} does not exist."
         exit 1
      else
         if [[ "${statefile_sql: -4}" != ".sql" ]] ; then
            echo "Error: ${statefile_sql} is not a *.sql file."
         else 
               if [[ "${dump_modus}" == "all" ]] || [[ "${dump_modus}" == "l_out" ]] ; then
                     l_out_da=$( sqlite3 ${statefile_sql}  ' SELECT id,seg1,seg2,lOh,lOe FROM pairs ' | head -100 | awk -F "|" 'BEGIN{data_da="False"}{ if( $4 != "0.0" || $5 != "0.0" ){ data_da="True" } }END{ print data_da }' )
                     if [ "${l_out_da}" == "True" ] ; then
                        echo "dump l_out  from ${statefile_sql}     to  l_out_${filename}.dat"
                        echo "### sqlite3 ${statefile_sql}  ' SELECT id,seg1,seg2,lOh,lOe FROM pairs '"          > l_out_${filename}.dat
                        sqlite3 ${statefile_sql}  ' SELECT id,seg1,seg2,lOh,lOe FROM pairs '                    >> l_out_${filename}.dat
                        sed -i 's/|/    /g'  l_out_${filename}.dat
                        if [ -d data ] && [ -e data ] ; then
                           mv l_out_${filename}.dat             data/
                        fi # mv l_out
                     fi # l_out_da
               fi #dump_modus 
               
               
               if [[ "${dump_modus}" == "all" ]] || [[ "${dump_modus}" == "E_out" ]] ; then
                     E_out_da=$( sqlite3 ${statefile_sql}  ' SELECT id,eCation,eNeutral,eAnion FROM segments ' | head -100 | awk -F "|" 'BEGIN{data_da="False"}{ if( $2 != "0.0" || $4 != "0.0" ){ data_da="True" } }END{ print data_da }' )
                     if [ "${E_out_da}" == "True" ] ; then
                        echo "dump E_out  from ${statefile_sql}     to  E_out_emultipoles_${filename}.dat"
                        echo "### sqlite3 ${statefile_sql}  ' SELECT id,eCation,eNeutral,eAnion FROM segments '" > E_out_emultipoles_${filename}.dat
                        sqlite3 ${statefile_sql}  ' SELECT id,eCation,eNeutral,eAnion FROM segments ' >>           E_out_emultipoles_${filename}.dat
                        sed -i 's/|/    /g'  E_out_emultipoles_${filename}.dat
                        if [ -d data ] && [ -e data ] ; then
                           mv E_out_emultipoles_${filename}.dat data/
                        fi # mv E_out
                     fi # E_out_da
               fi #dump_modus    


               if [[ "${dump_modus}" == "all" ]] || [[ "${dump_modus}" == "iZINDO" ]] ; then
				  if [ ! -e J_AB_iZINDO_MOO_${filename}.dat ] ; then
					  iZINDO_out_da=$( sqlite3 ${statefile_sql}  ' SELECT id,seg1,seg2,Jeff2h,Jeff2e FROM pairs ' | head -100 | awk -F "|" 'BEGIN{data_da="False"}{ if( $4 != "0.0" || $5 != "0.0" ){ data_da="True" } }END{ print data_da }' )
					  if [ "${iZINDO_out_da}" == "True" ] ; then
						 echo "dump iZINDO from ${statefile_sql}     to  J_AB_iZINDO_MOO_${filename}.dat"
						 echo " ### sqlite3 ${statefile_sql}  ' SELECT id,seg1,seg2,Jeff2h,Jeff2e ' > J_AB_iZINDO_MOO_${filename}.dat     " > J_AB_iZINDO_MOO_${filename}.dat
						 sqlite3 ${statefile_sql}  ' SELECT id,seg1,seg2,Jeff2h,Jeff2e FROM pairs '    >> J_AB_iZINDO_MOO_${filename}.dat
						 sed -i 's/|/    /g' J_AB_iZINDO_MOO_${filename}.dat
						 if [ -d data ] && [ -e data ] ; then
							mv J_AB_iZINDO_MOO_${filename}.dat   data/
						 fi # mv iZINDO
					  fi # iZINDO_out_da
					fi # file already exists.
               fi #dump_modus 
         fi ### *.sql ?
      fi # statefile_sql exists?

}
### END function  dump_data_from_sql_file    from   VOTCA *.sql file
### END Function to dump data 
###############################################################################################





###############################################################################################    
###			make_gro_to_VOTCA_90					START									###
###############################################################################################  
function make_gro_to_VOTCA_90(){
  Geo_Step=${1}
 if [ ! -e gro_to_VOTCA_G${Geo_Step}_90 ] ; then
	mkdir gro_to_VOTCA_G${Geo_Step}_90
	echo "mkdir gro_to_VOTCA_G${Geo_Step}_90"
 fi 
 
	cd gro_to_VOTCA_G${Geo_Step}_90
	if [ -e ../indexed_G${Geo_Step}_theta_90.gro ] ; then
		cp ../indexed_G${Geo_Step}_theta_90.gro  .
	fi
	
	if [ -e ../../enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_V5 ] && [ ! -e state_indexed_G${Geo_Step}_theta_90.sql ] ; then
			echo "Start ../../enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_V5 ..."
			../../enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_V5    gro_to_VOTCA indexed_G${Geo_Step}_theta_90.gro  out_indexed_G${Geo_Step}_theta_90.gro  1 ../new_sorted_neighbours_G${Geo_Step}.ngh > enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_G${Geo_Step}.log
			error_gro_to_VOTCA_90="$?" 
			check_exit_1 "${error_gro_to_VOTCA_90}" "$((${LINENO}-2))" 
	fi
	
	### ADD a symbolic link for   QC_FILES   and   MP_FILES
	if [ "${add_symbolic_link_to_folder}" == "True" ] ; then
			if [ ! -e QC_FILES ] ; then
				echo "---- Link to reference data "
				echo "ln -s /opt/files/t_koch08/votca_reference_data/dipbi_p3ht/QC_FILES/ QC_FILES"
				ln -s /opt/files/t_koch08/votca_reference_data/dipbi_p3ht/QC_FILES/ QC_FILES
			fi
			
			if [ ! -e MP_FILES ] ; then
				echo "---- Link to reference data "
				echo "ln -s /opt/files/t_koch08/votca_reference_data/dipbi_p3ht/MP_FILES_log2mps_CHelpG_PBE0_6-31Gs_PHxxA_crystal/    MP_FILES"
				ln -s /opt/files/t_koch08/votca_reference_data/dipbi_p3ht/MP_FILES_log2mps_CHelpG_PBE0_6-31Gs_PHxxA_crystal/    MP_FILES
			fi
	fi #### add_symbolic_link to folder
	
	
	if [ ! -e state_indexed_G${Geo_Step}_theta_90.sql ] ; then 
			/opt/gromacs-5.0.4-d-WP/bin/grompp_d  -c out_indexed_G${Geo_Step}_theta_90.gro -p fake_topology_indexed_G${Geo_Step}_theta_90.top -f gromp.mdp
			error_grompp_90="$?" 
			check_exit_1 "${error_grompp_90}" "$((${LINENO}-2))" 
			
			echo "grompp_d done"
			
			  ### 2) GROMACS *.tpr
			/opt/gromacs-5.0.4-d-WP/bin/mdrun_d  -s topol.tpr -x indexed_G${Geo_Step}_theta_90.xtc
			echo "mdrun_d done"
			
			### maping
			ctp_map -t topol.tpr -c indexed_G${Geo_Step}_theta_90.xtc -s map_indexed_G${Geo_Step}_theta_90.xml -f state_indexed_G${Geo_Step}_theta_90.sql
			error_ctp_map_90="$?" 
			check_exit_1 "${error_ctp_map_90}" "$((${LINENO}-2))" 
			echo "ctp_map done"
			
			if [ -e state_indexed_G${Geo_Step}_theta_90.sql ]; then
				echo "File created: state_indexed_G${Geo_Step}_theta_90.sql"
			else
				echo "Error in file production: state_indexed_G${Geo_Step}_theta_90.sql"
				exit 1
			fi
			
			cp  map_indexed_G${Geo_Step}_theta_90.xml system.xml
	fi ### statefile exists?
	
	

	
	pairs_exist=$(sqlite3 state_indexed_G${Geo_Step}_theta_90.sql " SELECT * FROM pairs " | wc -l ) 
	if [ "${pairs_exist}" == "0" ] ; then
			### neighborlist
			changeoption cutoff 2.3 neighbours_constrained_VOTCA__indexed_G${Geo_Step}_theta_90.xml
			ctp_run -e neighborlist -o neighbours_constrained_VOTCA__indexed_G${Geo_Step}_theta_90.xml -f state_indexed_G${Geo_Step}_theta_90.sql
			error_ctp_run_90="$?" 
			check_exit_1 "${error_ctp_run_90}" "$((${LINENO}-2))" 
	else
		echo " neighborlist already exists"
	fi

	if [ -e ../../options_cpt_run_DIPBI_P3HT_inter_and_intra_h.xml ]; then
		cp ../../options_cpt_run_DIPBI_P3HT_inter_and_intra_h.xml .
	else
		echo "Option file is not available: ../../options_cpt_run_DIPBI_P3HT_inter_and_intra_h.xml"
		exit 1
	fi
	
	echo '####### submit these lines !!!  #####################################################################################################################################'
	
	if [ "${HOSTNAME}" == "r06m01" ] ; then
			echo "  sub-CTP-VOTCA-palma2c-scratch  emultipoles_dipbi_p3ht_900K_6ns_to_300K_500ps_G${Geo_Step}_90   "'"   ctp_run -e emultipole   -o  ./options_cpt_run_DIPBI_P3HT_inter_and_intra_h.xml -f state_indexed_G${Geo_Step}_theta_90.sql -t 72 "'" 72 1 moria  120"
			echo "  sub-CTP-VOTCA-palma2c-scratch  eoutersphere_dipbi_p3ht_900K_6ns_to_300K_500ps_G${Geo_Step}_90  "'"   ctp_run -e eoutersphere -o  ./options_cpt_run_DIPBI_P3HT_inter_and_intra_h.xml -f state_indexed_G${Geo_Step}_theta_90.sql -t 72 "'" 72 1 moria  120 "
			#sub-CTP-VOTCA-palma2c-scratch  eoutersphere_dipbi_p3ht_900K_6ns_to_300K_500ps_G${Geo_Step}_90  "   ctp_run -e eoutersphere -o  ./options_cpt_run_DIPBI_P3HT_inter_and_intra_h.xml -f state_indexed_G${Geo_Step}_theta_90.sql -t 72 " 72 1 normal  120
	else
			echo "  sub-CTP-VOTCA_rsync  emultipoles_dipbi_p3ht_900K_6ns_to_300K_500ps_G${Geo_Step}_90   "'"   ctp_run -e emultipole   -o  ./options_cpt_run_DIPBI_P3HT_inter_and_intra_h.xml -f state_indexed_G${Geo_Step}_theta_90.sql -t 8 "'" 8 1 moria"
			echo "  sub-CTP-VOTCA_rsync  eoutersphere_dipbi_p3ht_900K_6ns_to_300K_500ps_G${Geo_Step}_90  "'"   ctp_run -e eoutersphere -o  ./options_cpt_run_DIPBI_P3HT_inter_and_intra_h.xml -f state_indexed_G${Geo_Step}_theta_90.sql -t 8 "'" 8 1 moria"
			#sub-CTP-VOTCA_rsync             eoutersphere_dipbi_p3ht_900K_6ns_to_300K_500ps_G${Geo_Step}_90  "   ctp_run -e eoutersphere -o  ./options_cpt_run_DIPBI_P3HT_inter_and_intra_h.xml -f state_indexed_G${Geo_Step}_theta_90.sql -t 8 " 8 1 batch
	fi ### palma or edoras / else
	
	
	### rename *.sql file
	if [ -e state_indexed_G${Geo_Step}_theta_90.sql ] && [ "${dump_data_32mer}" == "True" ] ; then
									if [ "${calc_data_32mer}" == "True" ] ; then 
											### function dump_data_from_sql_file in *.dat format.
											echo "dump_data_from_sql_file  state_indexed_G${Geo_Step}_theta_90.sql  "
											dump_data_from_sql_file  state_indexed_G${Geo_Step}_theta_90.sql      "${prefix}_32mer_G${Geo_Step}"  "l_out"	
											dump_data_from_sql_file  state_indexed_G${Geo_Step}_theta_90.sql      "${prefix}_32me_G${Geo_Step}"  "E_out"	
											dump_data_from_sql_file  state_indexed_G${Geo_Step}_theta_90.sql      "${prefix}_32me_G${Geo_Step}"  "iZINDO"
									else ### STANDARD 
											### function dump_data_from_sql_file in *.dat format.
											echo "dump_data_from_sql_file  state_indexed_G${Geo_Step}_theta_90.sql     ${prefix}  all"
											dump_data_from_sql_file  state_indexed_G${Geo_Step}_theta_90.sql           ${prefix}  "all"
									fi
	else 
							echo "Error: No *.sql file produced." 
	fi
	
	
 ### go back to G${Geo_Step} folder 
 cd ..
}   #### make_gro_to_VOTCA_90
######################  make_gro_to_VOTCA_90     END ##########################################################################################################################################################################################








###############################################################################################
#### END define local functions
###############################################################################################






Geo_Step=${Geo_Start}
## Hier Schleifenbeginn fuer multigeo_kMC
for i in `seq ${Geo_Start} ${N_Geo_Stepsmax}`; do  # Schleife ueber mehrer Geimetruen


kMC_Ordner_Step=${kMC_Basisordner%%_G*}_G${Geo_Step}
prefix=${kMC_Basisordner%%_G*}
indexed_gro_step="indexed_G${Geo_Step}.gro"

echo ${kMC_Ordner_Step}   ${indexed_gro_step}

current_grofilename=${kMC_Basisordner%%_G*}_G${Geo_Step}.gro
box_size_ang=$( tail --l 1 ${current_grofilename} | awk '{ printf "%10.3f  %10.3f   %10.3f \n", 10.0*$(1),10.0*$(2),10.0*$(3) }' )

echo ${current_grofilename}     ${box_size_ang}

if [ ! -e "${kMC_Ordner_Step}" ] ; then
	mkdir ${kMC_Ordner_Step}
fi # Ordner existiert

if [ ! -e "${kMC_Ordner_Step}/data" ] ; then
		mkdir ${kMC_Ordner_Step}/data
fi ###

if [ ! -e "${kMC_Ordner_Step}/${current_grofilename}" ] ; then
			cp  -p ${current_grofilename} ${kMC_Ordner_Step}/
fi ### skip current morphology as file already exists.

echo " From ${kMC_Ordner_old} to ${kMC_Ordner_Step} "
##cp -r ${kMC_Ordner_old} ${kMC_Ordner_Step}

# copy files
#for file in ${kMC_Ordner_old}/data ; do
#    if [ -e ${file} ] ; then
#        cp -r ${file} ${kMC_Ordner_Step}/
#    else
#        echo "Fehler: Die_Datei_konnte_nicht_kopiert_werden: ${file} "
#        exit 1
#    fi
#done 




if [ -d ${kMC_Ordner_Step} ] ; then
	cd  ${kMC_Ordner_Step}
	
	### generate indexed *.gro file with SCCS-cutoff angle for the Polymer; use reference data with theta=0.0 with "${use_theta0_reference_data}" == 'True' to get 1) resid_cut0_to_resid_r_cut_XX.X.dat file and  2) a more accurate neighbour list. with neighbour_list_with_reference_resid.py.
	if [ ! -e "indexed_G${Geo_Step}.gro" ] ; then
						
						if [ "${use_theta0_reference_data}" == 'True' ] && [ ! -e indexed_G${Geo_Step}_theta_0.gro ] ; then
							###### 1a) Einteilung der P3HT in sites. Erstellung indexed.gro
							echo " 1a) Einteilung in Hopping-sites      0.0   [in degree]"
							${executable_folder}/planar_dieder_def2   ${current_grofilename}   0.0
							if [ -e indexed.gro ] ; then
								mv indexed.gro   indexed_G${Geo_Step}_theta_0.gro
							else
								echo "Fehler: 1a) Die_Datei_konnte_nicht_erzeugt_werden: indexed.gro"
								echo "EXIT";  exit 1
							fi	#### gro?
						else
							echo " 1a) indexed_G${Geo_Step}_theta_0.gro already exists. "
						fi ## ! -e indexed_G${Geo_Step}_theta_0.gro
						
						
						# 1) Einteilung der P3HT in sites. Erstellung indexed.gro
						echo " 1) Einteilung in Hopping-sites      ${theta_cutoff}  [in degree]"
						if [ "${use_theta0_reference_data}" == 'True' ] && [ -e indexed_G${Geo_Step}_theta_0.gro ] ; then ### Beide methoden erzeugen die gleiche *.gro Datei. Nur hat die Referenzdatei, die mit indexed_G${Geo_Step}_theta_0.gro erzeugt wird die passende  resid_cut0_to_resid_r_cut_  datei.
							echo " use:   ${executable_folder}/planar_dieder_def2   indexed_G${Geo_Step}_theta_0.gro     ${theta_cutoff} "
							${executable_folder}/planar_dieder_def2   indexed_G${Geo_Step}_theta_0.gro     ${theta_cutoff} 
							
						elif [ -e "${current_grofilename}" ]; then
							echo " use:   ${executable_folder}/planar_dieder_def2   ${current_grofilename}     ${theta_cutoff} "  
							${executable_folder}/planar_dieder_def2   ${current_grofilename}     ${theta_cutoff} 
						else
							echo "Error: in planar_dieder_def2 "
							echo "EXIT"; exit 1
						fi
						
						
						if [ -e "indexed.gro" ] ; then
							mv indexed.gro   indexed_G${Geo_Step}.gro
						else
							echo "Fehler: Die_Datei_konnte_nicht_erzeugt_werden: indexed.gro"
							echo "EXIT";  exit 1
						fi
	else
		echo " 1) indexed_G${Geo_Step}.gro already exists. "
	fi # -e indexed_G${Geo_Step}.gro 





	###############################################################################################
    #### START   Reference data for theta=0
	### 1b)  data for reference theta_cut 0
	if [ "${use_theta0_reference_data}" == 'True' ] ; then
	
						if [ ! -e "indexed_G${Geo_Step}_theta_0.gro" ] ; then
							###### 1b) Einteilung der P3HT in sites. Erstellung indexed.gro
							echo " 1b) Einteilung in Hopping-sites      0.0   [in degree]"
							${executable_folder}/planar_dieder_def2   ${current_grofilename}   0.0
							if [ -e indexed.gro ] ; then
								mv indexed.gro   indexed_G${Geo_Step}_theta_0.gro
							else
								echo "Fehler: 1b) Die_Datei_konnte_nicht_erzeugt_werden: indexed.gro"
								echo "EXIT";  exit 1
							fi	#### gro?
						else
							echo " 1b) indexed_G${Geo_Step}_theta_0.gro already exists. "
						fi ## ! -e indexed_G${Geo_Step}_theta_0.gro
						
						if [ ! -e sorted_neighbours_G${Geo_Step}_theta_0.ngh ] || [ ! -e sorted_lambda_neighbours_G${Geo_Step}_theta_0.ngh ] ; then
										### create COM_G${Geo_Step}_theta_0.xyz
										if [ ! -e "COM_G${Geo_Step}_theta_0.xyz" ]; then
											${executable_folder}/COM_per_mol_def2   indexed_G${Geo_Step}_theta_0.gro
											echo "${executable_folder}/COM_per_mol_def2   indexed_G${Geo_Step}_theta_0.gro done"
											if [ -e COM.xyz ] ; then
												mv COM.xyz COM_G${Geo_Step}_theta_0.xyz
											else
												echo "Fehler: Die_Datei_konnte_nicht_erzeugt_werden: COM_G${Geo_Step}.xyz"
												echo "EXIT";  exit 1
											fi ####COM.xyz ?
										else 
											echo "File COM_G${Geo_Step}_theta_0.xyz already exists."
										fi
										
										
										### create no_box_G${Geo_Step}_theta_0.xyz
										if [ ! -e "no_box_G${Geo_Step}_theta_0.xyz" ]; then
													${executable_folder}/removebox        COM_G${Geo_Step}_theta_0.xyz  ${box_size_ang}
													echo "${executable_folder}/removebox  COM_G${Geo_Step}_theta_0.xyz  ${box_size_ang} done"
													if [ -e no_box.xyz ] ; then
															mv no_box.xyz no_box_G${Geo_Step}_theta_0.xyz
													else
															echo "Fehler: Die_Datei_konnte_nicht_erzeugt_werden: no_box_G${Geo_Step}.xyz"
															echo "EXIT";  exit 1
													fi ### no_box.xyz?
										else
												echo "File no_box_G${Geo_Step}_theta_0.xyz already exists!"
										fi
										
										### Currently not necessary
										if [ ! -e sorted_neighbours_G${Geo_Step}_theta_0.ngh ] && [ "${use_theta0_reference_data}" != 'True' ]  ; then
												echo "Start neighbour_list_sorted "
												${executable_folder}/neighbour_list_sorted          no_box_G${Geo_Step}_theta_0.xyz   indexed_G${Geo_Step}_theta_0.gro  ${box_size_ang}
												echo "${executable_folder}/neighbour_list_sorted    no_box_G${Geo_Step}_theta_0.xyz   indexed_G${Geo_Step}_theta_0.gro  ${box_size_ang}  done"
												if [ -e neighbours.ngh ] && [ -e sorted_neighbours.ngh ] ; then
													mv        neighbours.ngh          neighbours_G${Geo_Step}_theta_0.ngh 
													mv sorted_neighbours.ngh   sorted_neighbours_G${Geo_Step}_theta_0.ngh 
													if [ -e lambda_neighbours.ngh ]; then
														mv lambda_neighbours.ngh   lambda_neighbours_G${Geo_Step}_theta_0.ngh 
													fi
												else
													echo "Error: file missing."
													echo "EXIT";  exit 1
												fi
												
												if [ -e sorted_neighbours_G${Geo_Step}_theta_0.ngh ] ; then
													echo " 1c) Reference data for theta=0.0 was produced successfully!"
													#rm COM_G${Geo_Step}_theta_0.xyz
													#rm no_box_G${Geo_Step}_theta_0.xyz
													#rm neighbours_G${Geo_Step}_theta_0.ngh
												else
													echo " Error: 1c) No file produced sorted_neighbours_G${Geo_Step}_theta_0.ngh"
													exit 1
												fi
										fi ### ! -e sorted_neighbours_G${Geo_Step}_theta_0.ngh 
										
										if [ ! -e sorted_lambda_neighbours_G${Geo_Step}_theta_0.ngh ] ; then
												echo "Start  lambda_list sorted neighbours with theta=0.0   with NN_cutoff:   ${r_cutoff_lambda_for_reference_data}   "
												#### lambda_list [input.xyz] [input.gro] [cutoff [Angstr]]
														${executable_folder}/lambda_list         no_box_G${Geo_Step}_theta_0.xyz   indexed_G${Geo_Step}_theta_0.gro  ${r_cutoff_lambda_for_reference_data}
												echo   "${executable_folder}/lambda_list         no_box_G${Geo_Step}_theta_0.xyz   indexed_G${Geo_Step}_theta_0.gro  ${r_cutoff_lambda_for_reference_data} done"
												if [ -e sorted_lambda_neighbours.ngh ] ; then
													mv sorted_lambda_neighbours.ngh    sorted_lambda_neighbours_G${Geo_Step}_theta_0.ngh
													
												else
													echo "Error: file missing."
													echo "EXIT";  exit 1
												fi
												
												if [ -e lambda_neighbours.ngh ] ; then
													mv  lambda_neighbours.ngh    lambda_neighbours_G${Geo_Step}_theta_0.ngh
												fi
												
												
												if [ -e sorted_lambda_neighbours_G${Geo_Step}_theta_0.ngh ] ; then
													echo " 1c) Reference data for theta=0.0 was produced successfully!"
													#rm COM_G${Geo_Step}_theta_0.xyz
													#rm no_box_G${Geo_Step}_theta_0.xyz
													#rm neighbours_G${Geo_Step}_theta_0.ngh
												else
													echo " Error: 1c) No file produced sorted_lambda_neighbours_G${Geo_Step}_theta_0.ngh"
													exit 1
												fi
										fi ###  ! -e sorted_lambda_neighbours_G${Geo_Step}_theta_0.ngh 
									
									
									
									
								else 
										echo " 1c) File sorted_neighbours_G${Geo_Step}_theta_0.ngh already exists."
								fi ### ! -e sorted_neighbours_G${Geo_Step}_theta_0.ngh ???
		
	else
			echo "use_theta0_reference_data:  ${use_theta0_reference_data} "
	fi


    #### END  Reference data for theta=0
	######################################################################################################################################################################

	########################################
    #### 1d) calc data with reference 32mer
	   if [ "${calc_data_32mer}" == "True" ]; then
		   if [ ! -e indexed_G${Geo_Step}_theta_90.gro ] ; then
								# 1d) Einteilung der P3HT in sites. Erstellung indexed.gro
								echo " 1d) Einteilung in Hopping-sites      90.0  [in degree]"
								${executable_folder}/planar_dieder_def2   ${current_grofilename}     90.0 
								if [ -e "indexed.gro" ] ; then
									mv indexed.gro   indexed_G${Geo_Step}_theta_90.gro 
								else
									echo "Fehler: 1d) Die_Datei_konnte_nicht_erzeugt_werden: indexed_G${Geo_Step}_theta_90.gro "
									echo "EXIT";  exit 1
								fi
		   else
				echo " 1d) indexed_G${Geo_Step}_theta_90.gro already exists. "
			fi # -e indexed_G${Geo_Step}.gro 
			
			
			if [ "${start_VOTCA_setup_32mer}" == "True" ] ; then
					echo "START make_gro_to_VOTCA_90 for ${Geo_Step}"
					make_gro_to_VOTCA_90 ${Geo_Step}
			fi ###  make_gro_to_VOTCA_90
		   
		fi ## "calc_data_32mer" == "True" 
	
	#### calc data with reference 32mer
	########################################



		if [ ! -e "COM_G${Geo_Step}.xyz" ] ; then    
						# 2) Erstellung der Center-of-Mass.xyz Datei (COM.xyz)
						echo " 2) Erstellung COM_G${Geo_Step}.xyz" 
						${executable_folder}/COM_per_mol_def2   indexed_G${Geo_Step}.gro   
						if [ -e COM.xyz ] ; then
							mv COM.xyz COM_G${Geo_Step}.xyz
						else
							echo "Fehler: Die_Datei_konnte_nicht_erzeugt_werden: COM_G${Geo_Step}.xyz"
							echo "EXIT";  exit 1
						fi
		else
					echo " 2) COM_G${Geo_Step}.xyz already exists. "
		fi ## -e COM_G${Geo_Step}.xyz
	
	
	
		if [ ! -e "no_box_G${Geo_Step}.xyz" ] ; then 
							# 3) Alles in die Box schieben, sodass nichts herausragt, => no_box.xyz
							echo " 3) Erstellung no_box_G${Geo_Step}.xyz "
							echo " ${executable_folder}/removebox  COM_G${Geo_Step}.xyz  ${box_size_ang} "
							${executable_folder}/removebox  COM_G${Geo_Step}.xyz  ${box_size_ang}
							if [ -e no_box.xyz ] ; then
								mv no_box.xyz no_box_G${Geo_Step}.xyz
							else
								echo "Fehler: Die_Datei_konnte_nicht_erzeugt_werden: no_box_G${Geo_Step}.xyz"
								echo "EXIT";  exit 1
							fi
		else	
			echo " 3) no_box_G${Geo_Step}.xyz already exists. "		
		fi # -e  -e no_box_G${Geo_Step}.xyz
	
	
	
		if [ "${use_theta0_reference_data}" != 'True' ] ; then    
			if [ ! -e "lambda_neighbours_G${Geo_Step}.ngh" ] && [ "${use_theta0_reference_data}" != 'True' ]  ; then    
									# 4) Erstellung der lambda_list
									echo " 4) Erstellung der lambda_neighbours_G${Geo_Step} "
									${executable_folder}/lambda_list no_box_G${Geo_Step}.xyz indexed_G${Geo_Step}.gro ${r_cut_lambda_out_ang} ${box_size_ang}
									error_lambda_list="$?" 
									check_exit_1 "${error_lambda_list}" "$((${LINENO}-2))" 
									if [ -e lambda_neighbours.ngh ] ; then
										mv lambda_neighbours.ngh lambda_neighbours_G${Geo_Step}.ngh
									else
										echo "Fehler: Die_Datei_konnte_nicht_erzeugt_werden: lambda_neighbours_G${Geo_Step}.ngh"
										echo "EXIT";  exit 1
									fi
									
									if [ -e sorted_lambda_neighbours.ngh ] ; then
										mv sorted_lambda_neighbours.ngh sorted_lambda_neighbours_G${Geo_Step}.ngh
									else
										echo "Fehler: Die_Datei_konnte_nicht_erzeugt_werden: sorted_lambda_neighbours_G${Geo_Step}.ngh"
										echo "EXIT";  exit 1
									fi
									
			else
						echo " 4) lambda_neighbours_G${Geo_Step}.ngh already exists."
			fi ### lambda_neighbours_G${Geo_Step}.ngh  exists?  
		fi ### use_theta0_reference_data
		
		
		
		if [ "${use_theta0_reference_data}" != 'True' ] ; then    
				if [ -e sorted_neighbours_G${Geo_Step}.ngh ] && [ -e neighbours_G${Geo_Step}.ngh ] && [ "${use_theta0_reference_data}" != 'True' ]  ; then  ###&& [ -e sorted_lambda_neighbours_G${Geo_Step}.ngh ]
						echo " 5) Neighbour files already exist. "
				else
							
							# 5) Erstellung einer sortierten Nachbarschaftsliste
							echo " 5) Erstellung der Nachbarschaftsliste sorted_neighbours_G${Geo_Step}.ngh "
							${executable_folder}/neighbour_list_sorted no_box_G${Geo_Step}.xyz indexed_G${Geo_Step}.gro ${box_size_ang}
							if [ -e neighbours.ngh ] && [ -e sorted_neighbours.ngh ] ; then   ### [ -e sorted_lambda_neighbours.ngh ] &&
								##mv sorted_lambda_neighbours.ngh sorted_lambda_neighbours_G${Geo_Step}.ngh
								mv neighbours.ngh neighbours_G${Geo_Step}.ngh
								mv sorted_neighbours.ngh sorted_neighbours_G${Geo_Step}.ngh 
							else
								echo "Fehler: Eine_der_Datei_konnte_nicht_erzeugt_werden: neighbours.ngh,sorted_neighbours.ngh "
								echo "EXIT";  exit 1
							fi
				fi ### sorted_neighbours
		fi ###  use_theta0_reference_data != True
		
		
		
		### 6) Schwefelcounter / Zaehlt die Anzahl der zusammenhaengenden Resids der Thiophenringe und vergibt fuer DIPBI=33
		if [ -e Resid_Molname_Schwefel_G${Geo_Step}.dat ] ; then
					echo " 6) Resid_Molname_Schwefel_G${Geo_Step}.dat exists."
					#echo "Fehler: 6)_Schwefelcounter_Die_Dateien_existieren_bereits: Resid_Molname_Schwefel_G${Geo_Step}.dat"
					#echo "EXIT";  exit 1
		else 
					echo  "6)   ${executable_folder}/Schwefelcounter ${current_grofilename}   Resid_Molname_Schwefel_G${Geo_Step}.dat "
					${executable_folder}/Schwefelcounter ${current_grofilename}   Resid_Molname_Schwefel_G${Geo_Step}.dat
		fi ## Schwefelcounter






	### 7) new_sorted_neighbours; Erstellung der Nachbarschaftesliste fuer theta_cut anhand der Abstaende, die im 0deg-file fuer theta=0 berechnet wurden 
	### echo ' 7) new_sorted_neighbours wird erstellt '
	if [ -e new_sorted_neighbours_G${Geo_Step}.ngh ] ; then
			echo " 7) new_sorted_neighbours_G${Geo_Step}.ngh exists."
			#echo 'Fehler: Die_Datei_existiert_bereits,_und_soll_nicht_ueberschrieben_werden, ENDE'
			#exit 1
	else 
			if [ -e resid_cut0_to_resid_r_cut_${theta_cutoff}.0.dat ] && [ "${use_theta0_reference_data}" == 'True' ] ; then
					echo ' 7) new_sorted_neighbours wird erstellen '
						
					if [ -e sorted_lambda_neighbours_G${Geo_Step}_theta_0.ngh ] ; then
							### Better choice to calculated for a bigger number of accessible neighbours
						    echo "${executable_folder}/neighbour_list_with_reference_resid.py sorted_lambda_neighbours_G${Geo_Step}_theta_0.ngh  resid_cut0_to_resid_r_cut_${theta_cutoff}.0.dat"	
								  ${executable_folder}/neighbour_list_with_reference_resid.py sorted_lambda_neighbours_G${Geo_Step}_theta_0.ngh  resid_cut0_to_resid_r_cut_${theta_cutoff}.0.dat	
								
					elif [ -e sorted_neighbours_G${Geo_Step}_theta_0.ngh ] ; then	
							echo "${executable_folder}/neighbour_list_with_reference_resid.py       sorted_neighbours_G${Geo_Step}_theta_0.ngh   resid_cut0_to_resid_r_cut_${theta_cutoff}.0.dat"
								  ${executable_folder}/neighbour_list_with_reference_resid.py       sorted_neighbours_G${Geo_Step}_theta_0.ngh   resid_cut0_to_resid_r_cut_${theta_cutoff}.0.dat
								
					fi #### calculate neighbour_list_with_reference_resid.py for  sorted_lambda_neighbours_G${Geo_Step}_theta_0.ngh OR sorted_neighbours_G${Geo_Step}_theta_0.ngh
					
					
					
					if [ -e new_sorted_neighbours.ngh ] ; then
							echo "Datei_erzeugt: new_sorted_neighbours_G${Geo_Step}.ngh"
							mv new_sorted_neighbours.ngh new_sorted_neighbours_G${Geo_Step}.ngh
					else 
							echo "Fehler: Die_Datei_wurde_nicht_erzeugt: new_sorted_neighbours.ngh"
							echo "Fehler: Ende_beim_Erstellen_von: new_sorted_neighbours_G${Geo_Step}.ngh"
							## TO DO check error!!!
							exit 1
					fi 
			else
					echo "Fehler: Die_Datei_fuer_0_deg_existiert_nicht_im_Ordner: resid_cut0_to_resid_r_cut_${theta_cutoff}.0.dat"
					echo "Ende" ; exit 1
			fi # existiert *0.ngh ?
	fi #  new_sorted_neighbours.ngh exists ?
 
  
 
	### 8) run enumerated gro to votca
	if [[ "${start_VOTCA_setup}" == "True" ]] ; then
					if [ -e state_G${Geo_Step}.sql ]; then
							echo "8) state_G${Geo_Step}.sql already exists."
					else
							#if [ ! -e new_sorted_neighbours_G${Geo_Step}.ngh ] ; then
							#	echo "8) Error: new_sorted_neighbours_G${Geo_Step}.ngh does not exist, can not execute enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare "
							#	exit 1
							#else 
									if [ ! -e ../enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_V5 ] ; then
										echo "Error: gro_to_votca_executable is not available:   ../enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_V5"
									else
										if [ -e enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_G${Geo_Step}.log ] ; then
											rm enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_G${Geo_Step}.log
										fi ## rm *.log
										
										if [ "${calc_data_32mer}" == "True" ] && [ -e indexed_G${Geo_Step}_theta_90.gro ] ; then
											echo "8)   Start 32mer case"
											echo "8)  ../enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_V5     gro_to_VOTCA                	     indexed_G${Geo_Step}_theta_90.gro              out_indexed_G${Geo_Step}_theta_90.gro        1     "
											../enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_V5               gro_to_VOTCA                      indexed_G${Geo_Step}_theta_90.gro              out_indexed_G${Geo_Step}_theta_90.gro        1         > enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_G${Geo_Step}.log
											error_gro_to_VOTCA="$?" 
											check_exit_1 "${error_gro_to_VOTCA}" "$((${LINENO}-2))" 									
										else ### STANDARTD CASE
											#### enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_V5         :        method     enumerated_inputfilename.gro    outputfilename.gro  I_shift       [ neighbourlist_filename.ngh ] 
											echo "8)  ../enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_V5     gro_to_VOTCA                	     indexed_G${Geo_Step}.gro              out_indexed_G${Geo_Step}.gro        1            new_sorted_neighbours_G${Geo_Step}.ngh "
											../enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_V5               gro_to_VOTCA                      indexed_G${Geo_Step}.gro              out_indexed_G${Geo_Step}.gro        1            new_sorted_neighbours_G${Geo_Step}.ngh  > enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_G${Geo_Step}.log
											error_gro_to_VOTCA="$?" 
											check_exit_1 "${error_gro_to_VOTCA}" "$((${LINENO}-2))" 
										fi ## 32mer?
									fi ### -e ../../enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_V5
							#fi ### -e new_sorted_neighbours_G${Geo_Step}.ngh 
					fi ## state_G${Geo_Step}.sql ?
 
				###Next Steps TO DO:
				   if [[ -e enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_G${Geo_Step}.log ]] ; then 				
						grep "Next Steps TO DO:" -A 500 enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_G${Geo_Step}.log > bash_setup_votca_G${Geo_Step}.sh
						chmod +rwx bash_setup_votca_G${Geo_Step}.sh
						
						
						### ADD a symbolic link for   QC_FILES   and   MP_FILES    / or copy from lower layer 
						if [ "${add_symbolic_link_to_folder}" == "True" ] ; then
								if [ ! -e QC_FILES ] ; then
									echo "---- Link to reference data "
									echo "ln -s /opt/files/t_koch08/votca_reference_data/dipbi_p3ht/QC_FILES/ QC_FILES"
									ln -s /opt/files/t_koch08/votca_reference_data/dipbi_p3ht/QC_FILES/ QC_FILES
								fi
								
								if [ ! -e MP_FILES ] ; then
									echo "---- Link to reference data "
									echo "ln -s /opt/files/t_koch08/votca_reference_data/dipbi_p3ht/MP_FILES_log2mps_CHelpG_PBE0_6-31Gs_PHxxA_crystal/    MP_FILES"
									ln -s /opt/files/t_koch08/votca_reference_data/dipbi_p3ht/MP_FILES_log2mps_CHelpG_PBE0_6-31Gs_PHxxA_crystal/    MP_FILES
								fi
						else   ### copy from lower layer 
												MP_files_copied='False'
												### start copy folder with MP_FILES
												if [ ! -d ./MP_FILES ] && [ ! -e ./MP_FILES ] ; then ## check if MP_FILES are available.
													if [ -d ../MP_FILES ] && [ -e ../MP_FILES ] ; then
														echo "mv ../MP_FILES  . "
														mv ../MP_FILES  .
														MP_files_copied='True'
													else 
														echo "Error: No MP_FILES supplied."
														exit 1
													fi 
												fi ## MP_FILES ??
												### start copy folder with MP_FILES
												
												QC_FILES_copied='False'
												### start copy folder with QC_FILES
												if [ ! -d ./QC_FILES ] && [ ! -e ./QC_FILES ] ; then ## check if QC_FILES are available.
													if [ -d ../QC_FILES ] && [ -e ../QC_FILES ] ; then
														echo "mv ../QC_FILES  . "
														mv ../QC_FILES  .
														QC_FILES_copied='True'
													else 
														echo "Error: No QC_FILES supplied."
														exit 1
													fi 
												fi ## QC_FILES ??
												### start copy folder with QC_FILES				
						fi #### add_symbolic_link to folder   / or copy from lower layer 
						
						
						
						
						if [ -e bash_setup_votca_G${Geo_Step}.sh ] ; then
							./bash_setup_votca_G${Geo_Step}.sh 
						else
							echo "Error: File does not exist:  bash_setup_votca_G${Geo_Step}.sh "
						fi ### bash_setup_votca_G${Geo_Step}.sh ???
						 
						 
						 ### rename *.sql file
						 if [ -e state.sql ] ; then
									if [ "${calc_data_32mer}" == "True" ] ; then 
											echo "mv state.sql  state_G${Geo_Step}_32mer.sql"
											mv state.sql  state_G${Geo_Step}_32mer.sql
											### function dump_data_from_sql_file in *.dat format.
											dump_data_from_sql_file  state_G${Geo_Step}_32mer.sql      "${prefix}_32mer"  "l_out"	
											dump_data_from_sql_file  state_G${Geo_Step}_32mer.sql      "${prefix}_32mer"  "E_out"	
									else ### STANDARD 
											echo "mv state.sql  state_G${Geo_Step}.sql"
											mv state.sql  state_G${Geo_Step}.sql
											### function dump_data_from_sql_file in *.dat format.
											dump_data_from_sql_file  state_G${Geo_Step}.sql      ${prefix}  "all"
									fi
						 else 
							echo "Error: No *.sql file produced." 
						 fi
						 
						 						 
						### clean up 
						if [ "${add_symbolic_link_to_folder}" == "True" ] ; then
							continue
						else
							if [ "${MP_files_copied}" == 'True' ]; then
								### copy backwards
								echo "mv ../MP_FILES  .."
								mv ../MP_FILES  ..	
							fi ### MP_files_copied
							
							if [ "${QC_files_copied}" == 'True' ]; then
								### copy backwards
								echo "mv ../QC_FILES  .."
								mv ../QC_FILES  ..	
							fi ### QC_files_copied			
						fi ### clean up 
						
				   fi ### enumerated_gro_to_VOTCA_DIPBI_P3HT_prepare_G${Geo_Step}.log ???
   fi ### "${start_VOTCA_setup}" == "True" 
   ### END 8) run enumerated gro to votca 
   ################################################################################################################################################# 



	if [ "${setup_DIPRO}" == 'True' ] ; then
				echo " 10) Run setup DIPRO "
				if [ -e no_box_G${Geo_Step}.xyz ] && [ -e indexed_G${Geo_Step}.gro ] && [ -e new_sorted_neighbours_G${Geo_Step}.ngh ] ; then
					if [ -e subroutines/ ] && [ -e subroutines/DIPRO ] ;then
						echo "Folder  subroutines/ already exists."
						
					###path_to_OUZO='../get_started_OUZO_DIPRO_PM3'
					elif [ -e "${path_to_OUZO}" ]; then
						if [ -e "${path_to_OUZO}/subroutines" ] ; then
							echo "copy DIPRO subroutines folder:    cp -r "${path_to_OUZO}"/subroutines  "
							cp -r "${path_to_OUZO}"/subroutines     .
							
							if [ -e ${path_to_OUZO}/kmc_multigeo_V4 ] && [ ! -e kmc_multigeo_V4 ]; then
								echo "cp  ${path_to_OUZO}/kmc_multigeo_V4  ."
								cp  ${path_to_OUZO}/kmc_multigeo_V4  .
							fi
							
						else
							echo "Error: Dipro subroutines folder not found:  ${path_to_OUZO}/subroutines "; exit 1
						fi
					else 
						echo "Error: Dipro folder not found:  ${path_to_OUZO} "; exit 1
					fi
					
					if [ -e subroutines/aufruf.out ] && [ -e subroutines/DIPRO ] && [ -e subroutines/DIPRO/mol_all_setup_kankra.sh ] && [ -e subroutines/DIPRO/DIPRO-nMO-V10 ] ; then 
						if [ ! -e subroutines/Datenauswertung_SAB_H_H0_Jab_H0_H0_SAB_L0_L0_Jab_L0_L0_G${Geo_Step}.dat ] ; then
							echo "### DIPRO data ${kMC_Ordner_Step} with new_sorted_neighbours_G${Geo_Step}.ngh with method: ${g09_DIPRO_method_line} " >>  subroutines/Datenauswertung_SAB_H_H0_Jab_H0_H0_SAB_L0_L0_Jab_L0_L0_G${Geo_Step}.dat 
						fi
					else
							echo "DIPRO data missing"
							exit 1
						
					fi
					
					
					
					
					
					first_Resid=$( head --l 3 indexed_G${Geo_Step}.gro | tail --l 1 | cut -c 1-5 | awk '{ print $(1) }' )
					last_Resid=$(  tail --l 2 indexed_G${Geo_Step}.gro | head --l 1 | cut -c 1-5 | awk '{ print $(1) }' )
					### modify and replace pattern in   "kmc_multigeo_V4"
					#calc_complete_list=True 
					### this is the index from which the complete list calculation starts
					#start_index_for_complete_list_calculation= 700
					#end_index_for_complete_list_calculation= 1283
					declare -i ind ;ind=0
					for grep_pattern in "calc_complete_list=" "start_index_for_complete_list_calculation="  "end_index_for_complete_list_calculation=" "g09_DIPRO_method_line="  "calculation_method=" "frame_number=" ; do
						#ind=$((${ind}+1))
						ind+=1
						pattern=$( grep ${grep_pattern} kmc_multigeo_V4 | head --l 1 )
						#echo "${pattern}"
						if [[ ${ind} == 1 ]] ; then
								replace_pattern="${grep_pattern} True"  
						elif [[ ${ind} == 2 ]] ; then
								replace_pattern="${grep_pattern} ${first_Resid}" 
						elif [[ ${ind} == 3 ]] ; then
								replace_pattern="${grep_pattern} ${last_Resid}" 
						elif [[ ${ind} == 4 ]] ; then
								replace_pattern="${grep_pattern} '"${g09_DIPRO_method_line}"' " 
						elif [[ ${ind} == 5 ]] ; then
								replace_pattern="${grep_pattern} ' "${g09_calculation_method_string}" ' " 
						elif [[ ${ind} == 6 ]] ; then
								replace_pattern="${grep_pattern} ${Geo_Step} " 
						fi

						#### Replace first matching pattern with replace_pattern in kmc_file
						sed -i '0,/'"${pattern}"'/ s/'"${pattern}"'/'"${replace_pattern}"'/'  kmc_multigeo_V4
						### check successful replacement
						
						pattern=$( grep ${grep_pattern} kmc_multigeo_V4 | head --l 1 )
						
						if [ "${pattern}" != "${replace_pattern}" ] ; then
									echo "Error: substitution of the ${replace_pattern} was not successful in kmc_multigeo_V4 !"
									echo "Tried to replace ${pattern}  by  ${replace_pattern} "
									echo "Check the data "; exit 1
						else 
									echo "Using: ${pattern}  in kmc_multigeo_V4"
						fi
						
						
						
					done 
					########## modified kmc_multigeo_V4 for current _G${Geo_Step}.gro file
					
					
					
					### parameterexpected: kmc_multigeo_V4  expected: [distances.ngh] [COM.xyz] [morphology.gro] [debug = True/False] [lambda_out_neighbour_list] optional: [old matrix elem List] [old hole matrix elem List]
					  echo "START_DIPRO_LINE:       /usr/bin/python ./kmc_multigeo_V4   new_sorted_neighbours_G${Geo_Step}.ngh    no_box_G${Geo_Step}.xyz   indexed_G${Geo_Step}.gro  False   new_sorted_neighbours_G${Geo_Step}.ngh "
					  
					  echo "SUBMISSION_INPUTLINE:  sub-KMC   ${kMC_Ordner_Step}  "  "' /usr/bin/python ./kmc_multigeo_V4   new_sorted_neighbours_G${Geo_Step}.ngh    no_box_G${Geo_Step}.xyz   indexed_G${Geo_Step}.gro  False   new_sorted_neighbours_G${Geo_Step}.ngh '"  True  False   "${executable_folder}" 
					   
				else
					echo "Error: Not all needed input files for DIPRO calcualations are available."
					exit 1
				fi ### all input files available?
				
									echo 'My_EXIT'
					exit 1
	else
		echo "No setup_DIPRO selected: ${setup_DIPRO}"
	fi ### setup_DIPRO

    
    
    ## 12) Run kmc local 
     if [ "${run_KMC_local}" == "True" ] ; then  
             if [ -e ../run_kmc_votca_multiple_V3.sh ] ; then # Abgrage ob neue multigeo_kmc.dat vorhanden ist
	
				if [ ! -e ${prefix}_G${Geo_Step}.sql ]; then
					echo " " >> ${prefix}_G${Geo_Step}.sql
				fi		
				###  Modus=[single,singleUx,singleUy,singleUz,UxUyUz,xy-plane,yz-plane,xz-plane,increasing,increasing_small_steps,increasing_small_stepsUx,increasing_small_stepsUy,increasing_small_stepsUz,multidirection,multidirection_scaled,temperature_increasing] 
				#### Usage:    ./run_kmc_votca_multiple_V3.sh  options_VOTCA_kMC_file  kMC_VOTCA_startline    Modus        U_ext    
				#### Example:  ./run_kmc_votca_multiple_V3.sh        options.xml            state.sql        single       1.0E+7   
	
				echo "../run_kmc_votca_multiple_V3.sh  ../options_cpt_run_DIPBI_P3HT_inter_and_intra_h.xml     ${prefix}_G${Geo_Step}.sql   	UxUyUz  	1.0E+7"
				error_run_kmc_votca_multiple="$?" 
				check_exit_1 "${error_run_kmc_votca_multiple}" "$((${LINENO}-2))" 				
				
			else
				echo "Error: File is not available: ../run_kmc_votca_multiple_V3.sh "
				exit 1
			fi
	fi # run_KMC_local
    
    

    
    
    
    ## 13) Submit kMC-Rechnung        
    echo " ## 13) Submit kMC Step ${Geo_Step} "
    echo " ## 13) sub-KMC-VOTCA ../options_cpt_run_DIPBI_P3HT_inter_and_intra_h.xml  ${prefix}_G${Geo_Step}.sql   	UxUyUz  	1.0E+7   72  1 normal 48 "
    
	if [ "${sub_KMC}" == "True" ] ; then  
		if [ -e ./run_kmc_votca_multiple_V3.sh ] ; then # Abfrage ob  kmc mulitple da ist.
			
				###  Modus=[single,singleUx,singleUy,singleUz,UxUyUz,xy-plane,yz-plane,xz-plane,increasing,increasing_small_steps,increasing_small_stepsUx,increasing_small_stepsUy,increasing_small_stepsUz,multidirection,multidirection_scaled,temperature_increasing] 
				#### Usage:    ./run_kmc_votca_multiple_V3.sh  options_VOTCA_kMC_file  kMC_VOTCA_startline    Modus        U_ext    
				#### Example:  ./run_kmc_votca_multiple_V3.sh        options.xml            state.sql        single       1.0E+7   
				if [ ! -e ${prefix}_G${Geo_Step}.sql ]; then
					echo " " >> ${prefix}_G${Geo_Step}.sql
				fi
				echo "Start 13) sub-KMC-VOTCA ../options_cpt_run_DIPBI_P3HT_inter_and_intra_h.xml  ${prefix}_G${Geo_Step}.sql   	UxUyUz  	1.0E+7   72  1 normal 48 "
				sub-KMC-VOTCA ../options_cpt_run_DIPBI_P3HT_inter_and_intra_h.xml  ${prefix}_G${Geo_Step}.sql   	UxUyUz  	1.0E+7       72  1 normal 48 
				error_sub-KMC-VOTCA ="$?" 
				check_exit_1 "${error_sub-KMC-VOTCA}" "$((${LINENO}-2))" 
				
					            
				if [ wait_for_queue == "True" ]; then
						last_kMC_run_ok="false"
						## Hole Job ID um zu schauen, ob der alte Job in der Queue beendet wurde 
						### Check ob der JOB schon fertig ist
						sleep 20 # Zeit fuer die Initialisierung der Submittierung
						Job_queue_id=$(qstat -n | grep t_koch08 | grep ${kMC_Ordner_Step::12} | cut -c 1-15 ) ; echo ${Job_queue_id}
						status_symbol=$( qstat | grep ${Job_queue_id} | cut -c 70-70)   # Vorsicht hier beim cut
						echo ${status_symbol}
						#status_symbol="R"
						if [ "${status_symbol}" == "R" -o "${status_symbol}" == "Q" ] ; then # Job Running or waiting in the Queue 
								echo " Job bei der Rechnung: ${Job_queue_id}"
								echo " Warten ..."
								while [ `qstat | grep ${Job_queue_id} | cut -c 70-70` != "C" ]; do
										#echo "sleep" 
										sleep 30
								done
								sleep 30 ## Zeit um das Zurueckkopieren zu ermoeglichen
								if [ `qstat | grep ${Job_queue_id} | cut -c 70-70` == "C" ] ; then
									echo "Job finised"
									if [ -e $(ls ${kMC_Ordner_Step}/kmc_traj*.xyz) ] ; then 
										finishing_time=$( grep "finishing_time:" $( ls ${kMC_Ordner_Step}/kmc_info* | tail -1 ) )
										echo "${finishing_time} ${kMC_Ordner_Step}"
										if [ "${finishing_time%% *}" == "finishing_time:" ] ; then 
											tail -50 $( ls ${kMC_Ordner_Step}/kmc_info* | tail -1 ) ## Anzeige der letzten Output Daten
											last_kMC_run_ok='True'
										else
											echo "Fehler: Rechnung_nicht_richtig_beendet: ${kMC_Ordner_Step}/kmc_info"
										fi # finishing_time in kmc_info
									else
										echo "Fehler: Das_kmc_traj.xyz_file_existiert_nicht_${kMC_Ordner_Step}_Exit" ; exit 1 
									fi # Existiert kmc_traj file
								fi # Job compleated
						else
							echo "Fehler: The Job is not in the Queue after submission, check input and submission file. Exit ";  exit 1
						fi 
						### Ende Check ob der JOB schon fertig ist
					fi ## wait_for_queue == "True" 

			else # Datendatei fuer neuen Start nicht vorhanden: multigeo_kmc.dat? 
				echo "Fehler: Die_Datei_existiert_nicht: ../run_kmc_votca_multiple_V3.sh  "
				echo "EXIT"
				exit 1
			fi # Existiert neue multigeo_kmc.dat Datei ?
		fi # submit neuen kMC_multigeo Step
		
		
		
	cd ${CURDIR} # zurueck in kMC_Ordner_Step
	echo "TESTS ALL OK"
   
    
else
    echo "Fehler: Der Ordner ${kMC_Ordner_Step} existiert nicht. "
    echo "EXIT"
    exit 1
fi


## Abschnitt fuer Analyse der Daten nach der kMC-Rechnung
if [ "${last_kMC_run_ok}" == 'True' ] ; then
	echo "Start Analyse ${kMC_Ordner_Step}"
	echo " ....TO DO ADD ANALYSIS scripts here"
	echo "Ende Analyer ${kMC_Ordner_Step}"
else
	echo "Fehler: KMC_run_nicht_ok!"
	exit 1
fi ## Ende Analyse der kMC-Rechnungs-Daten


### Reset der Daten fuer den naechsten Geometrie-Schritt
  Geo_Step=$(($Geo_Step + 1))              # mit jeder Iteration, Geo_Step hochzaehlen
  kMC_Ordner_old=${kMC_Ordner_Step}
  indexed_gro_old=${indexed_gro_step}
if [ "${sub_KMC}" == "True" ] ; then
	### reset variable to restart.
	last_kMC_run_ok='False' 
fi
  
done # Ende Schleife ueber die Geometrien


 echo " Ende nach ${Geo_Step} von ${N_Geo_Stepsmax} "
 
