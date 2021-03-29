#!/bin/bash
# Pragamm zur Erzeugung der Startdaten einer kMC-Rechnung
 if [ $# -eq 0 ]; then
    echo "Usage: $0 < kMC Basisordner > "
    exit 1
 fi

 if [ ! -d "$1" ]; then
    echo "directory $1 doesn't exist!"
    exit 1
 fi

skript_path='..'
gro_filename_start='equi_1500ps_700k_def2_theta' ###'final_416xTHP32_pur' ### Verkuerzter String! ohne Endung
box_size_ang=' 139.4119  270.4301  146.2196  '  ##$( ${skript_path}/boxsize | head -1 )
r_cut_lambda_out_ang=20
typeset -i Geo_Step=0  ## Variable um Geometrieschritt in der *.gro file zu zählen
for crit_angle_theta_SCCS in 70 75 80 85  ; do   ###45 50 55 60 65 70 75 80 85

    echo " 1) Einteilung in Hopping-sites"
    if [ -e ${gro_filename_start}.gro ] ; then
    	${skript_path}/planar_dieder_def2 ${gro_filename_start}.gro  ${crit_angle_theta_SCCS} 
    else
		echo "Fehler: Die_gro_Startdatei_existiert_nicht: ${gro_filename_start}.gro"
	        echo "EXIT";  exit 1
    fi	

    if [ -e indexed.gro ] ; then
        gro_filename="${gro_filename_start}_${crit_angle_theta_SCCS}"
        mv indexed.gro ${gro_filename}_G${Geo_Step}.gro
    else
        echo "Fehler: Die_Datei_konnte_nicht_erzeugt_werden: indexed.gro"
        echo "EXIT";  exit 1
    fi
    
    
    # 2) Erstellung der Center-of-Mass.xyz Datei (COM.xyz)
    echo " 2) Erstellung COM_G${Geo_Step}.xyz" 
    ${skript_path}/COM_per_mol_def2 ${gro_filename}_G${Geo_Step}.gro   
    if [ -e COM.xyz ] ; then
        mv COM.xyz COM_${gro_filename}_G${Geo_Step}.xyz
    else
        echo "Fehler: Die_Datei_konnte_nicht_erzeugt_werden: COM_${gro_filename}_G${Geo_Step}.xyz"
        echo "EXIT";  exit 1
    fi
    
    
    # 3) Alles in die Box schieben, sodass nichts herausragt, => no_box.xyz
    echo " 3) Erstellung no_box_G${Geo_Step}.xyz "

    ${skript_path}/removebox  COM_${gro_filename}_G${Geo_Step}.xyz  ${box_size_ang}
    if [ -e no_box.xyz ] ; then
        mv no_box.xyz no_box_${gro_filename}_G${Geo_Step}.xyz
    else
        echo "Fehler: Die_Datei_konnte_nicht_erzeugt_werden: no_box_${gro_filename}_G${Geo_Step}.xyz"
        echo "EXIT";  exit 1
    fi


    
if [ -e lambda_neighbours_${gro_filename}_G${Geo_Step}.ngh ]; then    
    echo " 4) Die Datei existiert bereits lambda_neighbours_${gro_filename}_G${Geo_Step}.ngh"
else    
    # 4) Erstellung der lambda_list
    echo " 4) Erstellung der lambda_neighbours_G${Geo_Step} "
    ${skript_path}/lambda_list no_box_${gro_filename}_G${Geo_Step}.xyz ${gro_filename}_G${Geo_Step}.gro ${r_cut_lambda_out_ang} ${box_size_ang}
    if [ -e lambda_neighbours.ngh ] ; then
        mv lambda_neighbours.ngh lambda_neighbours_${gro_filename}_G${Geo_Step}.ngh
    else
        echo "Fehler: Die_Datei_konnte_nicht_erzeugt_werden: lambda_neighbours_${gro_filename}_G${Geo_Step}.ngh"
        echo "EXIT";  exit 1
    fi
fi 
 
   
if [ -e sorted_neighbours_${gro_filename}_G${Geo_Step}.ngh ] && [ -e neighbours_${gro_filename}_G${Geo_Step}.ngh ] && [ -e  sorted_lambda_neighbours_${gro_filename}_G${Geo_Step}.ngh ]  ; then
    echo " 5) Die_Dateien_existieren_bereits"
else
    # 5) Erstellung einer sortierten Nachbarschaftsliste
    echo " 5) Erstellung der Nachbarschaftsliste sorted_lambda_neighbours_G${Geo_Step}.ngh "
    ${skript_path}/neighbour_list_sorted no_box_${gro_filename}_G${Geo_Step}.xyz ${gro_filename}_G${Geo_Step}.gro ${box_size_ang}
    if [ -e sorted_lambda_neighbours.ngh ] && [ -e neighbours.ngh ] && [ -e sorted_neighbours.ngh ] ; then
        mv sorted_lambda_neighbours.ngh sorted_lambda_neighbours_${gro_filename}_G${Geo_Step}.ngh
        mv neighbours.ngh neighbours_${gro_filename}_G${Geo_Step}.ngh
        mv sorted_neighbours.ngh sorted_neighbours_${gro_filename}_G${Geo_Step}.ngh 
    else
        echo "Fehler: Eine_der_Datei_konnte_nicht_erzeugt_werden: sorted_lambda_neighbours.ngh,neighbours.ngh,sorted_neighbours.ngh "
        echo "EXIT";  exit 1
    fi
fi    


### 6) Schwefelcounter / Zaehlt die Anzahl der zusammenhaengenden Resids der Thiophenringe und vergibt fuer DIPBI=33
if [ -e Resid_Molname_Schwefel_${gro_filename}_G${Geo_Step}.dat ] ; then
    echo "Fehler: 6)_Schwefelcounter_Die_Dateien_existieren_bereits: Resid_Molname_Schwefel_${gro_filename}_G${Geo_Step}.dat"
    echo "EXIT";  exit 1
else 
	${skript_path}/Schwefelcounter ${gro_filename}_G${Geo_Step}.gro Resid_Molname_Schwefel_${gro_filename}_G${Geo_Step}.dat
fi ## Schwefelcounter


### 7) new_sorted_neighbours; Erstellung der Nachbarschaftesliste fuer theta_cut anhand der Abstaende, die im 0deg-file fuer theta=0 berechnet wurden 
echo ' 7) new_sorted_neighbours wird erstellt '
if [ -e 'new_sorted_neighbours.ngh' ] ; then
	echo 'Fehler: Die_Datei_existiert_bereits,_und_soll_nicht_ueberschrieben_werden, ENDE'
	exit 1
else 
	if [ -e "sorted_neighbours_${gro_filename_start}_0.ngh" ] && [ -e resid_cut0_to_resid_r_cut_${crit_angle_theta_SCCS}.0.dat ] ; then
		echo "${skript_path}/neighbour_list_with_reference_resid.py ${gro_filename_start}_0.ngh resid_cut0_to_resid_r_cut_${crit_angle_theta_SCCS}.0.dat"
		${skript_path}/neighbour_list_with_reference_resid.py sorted_neighbours_${gro_filename_start}_0.ngh resid_cut0_to_resid_r_cut_${crit_angle_theta_SCCS}.0.dat
		if [ -e new_sorted_neighbours.ngh ] ; then
			echo "Datei_erzeugt: new_sorted_neighbours_${gro_filename}_G${Geo_Step}.ngh"
			mv new_sorted_neighbours.ngh new_sorted_neighbours_${gro_filename}_G${Geo_Step}.ngh
		else 
			echo "Fehler: Die_Datei_wurde_nicht_erzeugt: new_sorted_neighbours.ngh"
			echo "Fehler: Ende_beim_Erstellen_von: new_sorted_neighbours_${gro_filename}_G${Geo_Step}.ngh"
			exit 1
		fi 
	else
		echo "Fehler: Die_Datei_fuer_0_deg_existiert_nicht_im_Ordner: sorted_neighbours_${gro_filename_start}_0.ngh"
		echo "Ende" ; exit 1
	fi # existiert *0.ngh ?
fi #  new_sorted_neighbours.ngh exists ?


### 8) Auswertung der intramolekularen Daten fuer die Verteilungen der energetic disorder \sigma
echo "8) Berechnung der energetic disorder "
if [ -e Daten_energetic_disorder_${gro_filename}_G${Geo_Step}.dat ] ; then
	echo "Fehler: Die_Datei_existiert_bereits,_und_soll_nicht_ueberschrieben_werden  Daten_energetic_disorder_${gro_filename}_G${Geo_Step}.dat  , ENDE"
	exit 1
else 
    if [ -e ${gro_filename}_G${Geo_Step}.gro ] && [ -e new_sorted_neighbours_${gro_filename}_G${Geo_Step}.ngh ]  && [ -e no_box_${gro_filename}_G${Geo_Step}.xyz ] ; then
        ${skript_path}/calc_energetic_disorder_V2 ${gro_filename}_G${Geo_Step}.gro  new_sorted_neighbours_${gro_filename}_G${Geo_Step}.ngh no_box_${gro_filename}_G${Geo_Step}.xyz Daten_energetic_disorder_${gro_filename}_G${Geo_Step}.dat
    else
		echo "Fehler: calc_energetic_disorder:_Eine_benoetigte_Datei_existiert_nicht_im_Ordner:  ${gro_filename}_G${Geo_Step}.gro  new_sorted_neighbours_${gro_filename}_G${Geo_Step}.ngh no_box_${gro_filename}_G${Geo_Step}.xyz "
		echo "Ende" ; exit 1
	fi # existiert *0.ngh ?
fi 


done # Schleife crit_angle_theta_SCCS
    
