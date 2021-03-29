#!/bin/bash
### The programm evalutates the mobility and corrects the error in the VOTCA implementation concerning the implemente projection of the displacement vector on the external field vektor.
### Use VOTCA only if the external field vector is parrallel to the box sides, or use this skript.
### Execute this skript  in the folder with the kmc simulations created by run_kmc_votca_multiple_V3.sh
### The script will get the applied external field and the mobility from the options.xml files and the kmc_votca_outputfile.
### Maybe one has tto adapt the names below : optionsfilexml and kmc_votca_outputfile
### CAUTION!! Change the filenames TWO times !!!
### Set the number of calculations using:
### Ni_max=64
### Ui_max=16 
### 

### The script is needed for correct polar mobility plots using VOTCA!  
function check_file_exists(){     # Schaut ob Datei existiert 
        filename=${1}
        if [ -e ${filename} ] ; then
            echo "True"
        else
            echo "False"
            echo "Fehler: ${filename} ist nicht im Startordner verfuegbar." 
            echo ${filename}
            echo "EXIT"
            exit 1
        fi
    } # file_exists

### Change the parameters here !!!  
Ni_max=64
Ui_max=16 


### dafault values all zero:
average_mobilityX=0.0
average_mobilityY=0.0
average_mobilityZ=0.0
average_mobility=0.0

sigma_mobilityX=0.0
sigma_mobilityY=0.0
sigma_mobilityZ=0.0
sigma_mobility=0.0

N_mobilities=0
N_mobilities_total=0


for Ui in ` seq 1 ${Ui_max} `; do

average_mobilityX=0.0
average_mobilityY=0.0
average_mobilityZ=0.0
average_mobility=0.0

sigma_mobilityX=0.0
sigma_mobilityY=0.0
sigma_mobilityZ=0.0
sigma_mobility=0.0

N_mobilities=0


  for Ni in ` seq 1 ${Ni_max} ` ; do
   
   ### CAUTION!! Change the filenames here and BELOW !!!
   optionsfilexml=" kmc_votca_U${Ui}_N${Ni}/options_cpt_run_DIPBI_P3HT_inter_and_intra_PM3_U${Ui}_N${Ni}.xml "
   kmc_votca_outputfile="kmc_votca_U${Ui}_N${Ni}/kmc_votca_U${Ui}_N${Ni}.out"

      if [ ! -e ${optionsfilexml} ] ; then
         echo "Skip  ${optionsfilexml}"
         continue      # Skip rest of this particular loop iteration.
      fi
      if [ ! -e ${kmc_votca_outputfile} ] ; then
         continue      # Skip rest of this particular loop iteration.
      fi

      
      
   fieldX=$( grep "fieldX" ${optionsfilexml} ) ; fieldX=${fieldX#*>} ; fieldX=${fieldX%<*} ##; echo ${fieldX}
   fieldY=$( grep "fieldY" ${optionsfilexml} ) ; fieldY=${fieldY#*>} ; fieldY=${fieldY%<*} ##; echo ${fieldY}
   fieldZ=$( grep "fieldZ" ${optionsfilexml} ) ; fieldZ=${fieldZ#*>} ; fieldZ=${fieldZ%<*} ##; echo ${fieldZ} 

   #echo "External field in [V/m]:    ${fieldX}   ${fieldY}   ${fieldZ} "
   result_line=$( grep "Average velocities (m/s):" -A 1 ${kmc_votca_outputfile}   | sed 's/,/ /g' | sed 's&]& &g'  | sed 's/charge 1: / /g' | sed 's/\[/ /g'  | tail --l 1 | awk  -v fieldX=${fieldX} -v fieldY=${fieldY} -v  fieldZ=${fieldZ} -v kmc_votca_outputfile=${kmc_votca_outputfile} '
      function abs(v){return v < 0 ? -v : v}
      BEGIN{}
      {

      vx=$(1)
      vy=$(2)
      vz=$(3)
      ### print $(1),$(2),$(3),fieldX,fieldY,fieldZ 
      absolute_field = sqrt(fieldX*fieldX + fieldY*fieldY + fieldZ*fieldZ)	
      fieldfactorsX=(1.0E4*fieldX)/absolute_field/absolute_field
      fieldfactorsY=(1.0E4*fieldY)/absolute_field/absolute_field
      fieldfactorsZ=(1.0E4*fieldZ)/absolute_field/absolute_field

      mobility_new= vx*fieldfactorsX + vy*fieldfactorsY + vz*fieldfactorsZ 
      
      mobilityX=fieldX * mobility_new / absolute_field 
      mobilityY=fieldY * mobility_new / absolute_field 
      mobilityZ=fieldZ * mobility_new / absolute_field 

      fieldfactorsX=0
      fieldfactorsY=0
      fieldfactorsZ=0
      if ( abs(fieldX) != 0.0 ) { 	fieldfactorsX=(1.0E4/fieldX) }
      if ( abs(fieldY) != 0.0 ) {	fieldfactorsY=(1.0E4/fieldY) }
      if ( abs(fieldZ) != 0.0 ) {	fieldfactorsZ=(1.0E4/fieldZ) }

      mobilityX_old=vx*fieldfactorsX
      mobilityY_old=vy*fieldfactorsY
      mobilityZ_old=vz*fieldfactorsZ
      mobility_old=sqrt( mobilityX_old^2  +  mobilityY_old^2 +  mobilityZ_old^2 )
      
     ### print "mobility_new: ",mobility_new,"mobilityX: ",mobilityX," mobilityY: ",mobilityY," mobilityZ: ",mobilityZ, "mobility_old: ",mobility_old, "mobilityX: ",mobilityX_old," mobilityY: ",mobilityY_old," mobilityZ: ",mobilityZ_old,"   ",fieldX," ",fieldY," ",fieldZ,"  ",kmc_votca_outputfile 
         printf " mobility_new: %.10E  mobilityX:  %.10E  mobilityY:  %.10E  mobilityZ: %.10E  mobility_old:   %.10E  mobilityX: %.10E  mobilityY:  %.10  mobilityZ: %.10E           %.4E   %.4E %.4E  \n",mobility_new,mobilityX,mobilityY,mobilityZ, mobility_old, mobilityX_old, mobilityY_old, mobilityZ_old, fieldX, fieldY, fieldZ,kmc_votca_outputfile 
   }'  )

      N_data=$( echo ${result_line} | wc | awk '{ print $(3) }' )
      if [[ ${N_data} -ge 10 ]] ; then
            N_mobilities=$(  echo ${N_mobilities} | awk '{ print $(1)+1 }' )
            
            average_mobilityX=$( echo ${average_mobilityX}  ${result_line} | awk '{ printf "%.10e", $(1) + $(5) }' )
            average_mobilityY=$( echo ${average_mobilityY}  ${result_line} | awk '{ printf "%.10e", $(1) + $(7) }' )
            average_mobilityZ=$( echo ${average_mobilityZ}  ${result_line} | awk '{ printf "%.10e", $(1) + $(9) }' )

           ### echo "SUM: ${average_mobilityX}       ${average_mobilityY}    ${average_mobilityZ}     ${average_mobility}   ${N_mobilities} "
      fi 
      
      if [ ${Ni} == ${Ni_max} ] ; then
            average_mobilityX_ref=$( echo ${average_mobilityX}  ${N_mobilities} | awk '{ var=$(1); i=$(2) ; if( i > 1){ printf "%.10e", var/i } else{ printf "%.10e", var} }' )
            average_mobilityY_ref=$( echo ${average_mobilityY}  ${N_mobilities} | awk '{ var=$(1); i=$(2) ; if( i > 1){ printf "%.10e", var/i } else{ printf "%.10e", var} }' )
            average_mobilityZ_ref=$( echo ${average_mobilityZ}  ${N_mobilities} | awk '{ var=$(1); i=$(2) ; if( i > 1){ printf "%.10e", var/i } else{ printf "%.10e", var} }' )
                      

            average_mobility_ref=$( echo  ${average_mobilityX_ref}   ${average_mobilityY_ref}  ${average_mobilityZ_ref}    | awk '{ printf "%.10e", sqrt($(1)^2 + $(2)^2 + $(3)^2) }' )
            
            #echo "Averaged reference data:"  ${average_mobilityX_ref}   ${average_mobilityY_ref}  ${average_mobilityZ_ref} ${average_mobility_ref}
      fi
      
   done ### loop Ni kmc files 

   average_mobilityX=0.0
   average_mobilityY=0.0
   average_mobilityZ=0.0
   average_mobility=0.0
   
   sigma_mobilityX=0.0
   sigma_mobilityY=0.0
   sigma_mobilityZ=0.0
   sigma_mobility=0.0
   
   
   #### calculate sigma with same function
   for Ni in ` seq 1  ${Ni_max}   ` ; do
	
      ### CAUTION!! Change the filenames here and ABOVE !!!  
      optionsfilexml=" kmc_votca_U${Ui}_N${Ni}/options_cpt_run_DIPBI_P3HT_inter_and_intra_PM3_U${Ui}_N${Ni}.xml "
      kmc_votca_outputfile="kmc_votca_U${Ui}_N${Ni}/kmc_votca_U${Ui}_N${Ni}.out"

      if [ ! -e ${optionsfilexml} ] ; then
         continue      # Skip rest of this particular loop iteration.
      fi
      if [ ! -e ${kmc_votca_outputfile} ] ; then
         continue      # Skip rest of this particular loop iteration.
      fi
         
         
      fieldX=$( grep "fieldX" ${optionsfilexml} ) ; fieldX=${fieldX#*>} ; fieldX=${fieldX%<*} ##; echo ${fieldX}
      fieldY=$( grep "fieldY" ${optionsfilexml} ) ; fieldY=${fieldY#*>} ; fieldY=${fieldY%<*} ##; echo ${fieldY}
      fieldZ=$( grep "fieldZ" ${optionsfilexml} ) ; fieldZ=${fieldZ#*>} ; fieldZ=${fieldZ%<*} ##; echo ${fieldZ} 

      ###echo "External field in [V/m]:    ${fieldX}   ${fieldY}   ${fieldZ} "
      ### grep "Average velocities (m/s):" -A 1 ${kmc_votca_outputfile}  
      result_line=$( grep "Average velocities (m/s):" -A 1 ${kmc_votca_outputfile}   | sed 's/,/ /g' | sed 's&]& &g'  | sed 's/charge 1: / /g' | sed 's/\[/ /g'  | tail --l 1 | awk  -v fieldX=${fieldX} -v fieldY=${fieldY} -v  fieldZ=${fieldZ} -v kmc_votca_outputfile=${kmc_votca_outputfile} '
         function abs(v){return v < 0 ? -v : v}
         BEGIN{}
         {

         vx=$(1)
         vy=$(2)
         vz=$(3)
         	
         ### print $(1),$(2),$(3),fieldX,fieldY,fieldZ ,absolute_field
         absolute_field = sqrt(fieldX*fieldX + fieldY*fieldY + fieldZ*fieldZ)	
         fieldfactorsX=(1.0E4*fieldX)/absolute_field/absolute_field
         fieldfactorsY=(1.0E4*fieldY)/absolute_field/absolute_field
         fieldfactorsZ=(1.0E4*fieldZ)/absolute_field/absolute_field

         mobility_new= vx*fieldfactorsX + vy*fieldfactorsY + vz*fieldfactorsZ 
      
         mobilityX=fieldX * mobility_new / absolute_field 
         mobilityY=fieldY * mobility_new / absolute_field 
         mobilityZ=fieldZ * mobility_new / absolute_field 

         ### old implementation BuG in VOTCA
         fieldfactorsX=0
         fieldfactorsY=0
         fieldfactorsZ=0
         if ( abs(fieldX) != 0.0 ) { 	fieldfactorsX=(1.0E4/fieldX) }
         if ( abs(fieldY) != 0.0 ) {	fieldfactorsY=(1.0E4/fieldY) }
         if ( abs(fieldZ) != 0.0 ) {	fieldfactorsZ=(1.0E4/fieldZ) }

         mobilityX_old=vx*fieldfactorsX
         mobilityY_old=vy*fieldfactorsY
         mobilityZ_old=vz*fieldfactorsZ
         mobility_old=sqrt( mobilityX_old^2  +  mobilityY_old^2 +  mobilityZ_old^2 )
         
         ###print "mobility_new: ",mobility_new,"mobilityX: ",mobilityX," mobilityY: ",mobilityY," mobilityZ: ",mobilityZ, "mobility_old: ",mobility_old, "mobilityX: ",mobilityX_old," mobilityY: ",mobilityY_old," mobilityZ: ",mobilityZ_old,"   ",fieldX," ",fieldY," ",fieldZ,"  ",kmc_votca_outputfile 
          printf " mobility_new: %.10E  mobilityX:  %.10E  mobilityY:  %.10E  mobilityZ: %.10E  mobility_old:   %.10E  mobilityX: %.10E  mobilityY:  %.10E  mobilityZ: %.10E           %.4E   %.4E %.4E  \n",mobility_new,mobilityX,mobilityY,mobilityZ, mobility_old, mobilityX_old, mobilityY_old, mobilityZ_old, fieldX, fieldY, fieldZ,kmc_votca_outputfile 
      }'  )


         
         N_data=$( echo ${result_line} | wc | awk '{ print $(3) }' )
         if [[ ${N_data} -ge 10 ]] ; then
               ##N_mobilities=$(  echo ${N_mobilities} | awk '{ print $(1)+1 }' )
               N_mobilities_total=$(  echo ${N_mobilities_total} | awk '{ print $(1)+1 }' )
               
               #echo "RESULT: ${result_line} "

               MuX=$( echo  ${result_line} | awk '{ printf "%.10e",$(4) }' )
               MuY=$( echo  ${result_line} | awk '{ printf "%.10e",$(6) }' )
               MuZ=$( echo  ${result_line} | awk '{ printf "%.10e",$(8) }' )
               #### echo "mu:  $MuX, $MuY, $MuZ"
               sigma_mobilityX=$( echo ${average_mobilityX_ref}   ${MuX}     ${sigma_mobilityX} | awk '{ printf "%.10e",$(3) + ($(1) - $(2))^2 }'  )
               sigma_mobilityY=$( echo ${average_mobilityY_ref}   ${MuY}     ${sigma_mobilityY} | awk '{ printf "%.10e",$(3) + ($(1) - $(2))^2 }'  )
               sigma_mobilityZ=$( echo ${average_mobilityZ_ref}   ${MuZ}     ${sigma_mobilityZ} | awk '{ printf "%.10e",$(3) + ($(1) - $(2))^2 }'  )
                       
               ###echo  "VAR:" ${sigma_mobilityX}   ${sigma_mobilityY}  ${sigma_mobilityZ}

         fi 
         
         if [ "${Ni}" == "${Ni_max}" ] ; then

            sigma_mobilityX=$( echo ${sigma_mobilityX}   ${N_mobilities}      | awk   '{ var=$(1) ; i=$(2) ;  if(i > 1) { sigma=(var/(i-1))^0.5 } else{ sigma=0 } printf "%.10e", sigma }' )
            sigma_mobilityY=$( echo ${sigma_mobilityY}   ${N_mobilities}      | awk   '{ var=$(1) ; i=$(2) ;  if(i > 1) { sigma=(var/(i-1))^0.5 } else{ sigma=0 } printf "%.10e", sigma }' )
            sigma_mobilityZ=$( echo ${sigma_mobilityZ}   ${N_mobilities}      | awk   '{ var=$(1) ; i=$(2) ;  if(i > 1) { sigma=(var/(i-1))^0.5 } else{ sigma=0 } printf "%.10e", sigma }' )
            sigma_mobility=$(  echo ${sigma_mobilityX}   ${sigma_mobilityY}   ${sigma_mobilityZ}    | awk '{ printf "%.10e", sqrt( $(1)^2 + $(2)^2 + $(3)^2) }' )
            #### echo "---------RESULT-----------"
            data_line=$(  echo " ${average_mobilityX_ref}  ${sigma_mobilityX}   ${average_mobilityY_ref}  ${sigma_mobilityY}  ${average_mobilityZ_ref}  ${sigma_mobilityZ}  ${average_mobility_ref}  ${sigma_mobility}    ${N_mobilities} " | awk '{ printf "%.4E  %.4E    %.4E  %.4E    %.4E  %.4E    %.4E  %.4E               \n",$(1),$(2),$(3),$(4),$(5),$(6),$(7),$(8) }' )
            echo   "${data_line}            ${fieldX}   ${fieldY}   ${fieldZ}    ${N_mobilities}"
         fi
      done ### loop Ni kmc files 

done ### loop output files Ui
