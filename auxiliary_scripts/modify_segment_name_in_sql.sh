#!/bin/bash

if [ $# -ne 2 ] ; then
	echo " Programm to read a list of names for the segments into the *.sql file "
	echo " Make sure the lists have equal length!"
	echo " Usinge:  $0         readfile            sqlfile   "
	echo "          $0         Resid_list.dat      state.sql "
	exit 1 
fi

readfile=${1}
sqlfile=${2}

N_lines_readfile=$( cat ${readfile} | wc -l ) 
N_lines_sqlfile=$( sqlite3 ${sqlfile} " SELECT name FROM segments " | wc -l )

echo " ${N_lines_readfile}   ${N_lines_sqlfile} "
if [[ "${N_lines_readfile}" == "${N_lines_sqlfile}" ]]; then
	echo " Using: ${readfile}     ${sqlfile} "
else
	echo "Warning the number of lines in the files do not match!"
	echo " ${N_lines_readfile} ${readfile} "
	echo " ${N_lines_sqlfile}  ${sqlfile} "
	exit 1
fi
i=0
while IFS= read -r cmd; do
	i=$((${i}+1))	
    	printf '%s %s \n' "$i" "$cmd"
	newname="$cmd"
	sqlite3 ${sqlfile} " UPDATE segments SET name='$newname'  WHERE id='${i}' "

     
done < "$readfile"
