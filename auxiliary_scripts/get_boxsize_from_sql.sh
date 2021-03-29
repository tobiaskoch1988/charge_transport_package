#!/bin/bash
statefile='state.sql'
if [[ $# -eq 1 ]]; then   ## Einlesen des Names des statefiles.
	statefile=${1}
else 
	echo " Gets the box dimensions for a state.sql   IN: state.sql    OUT: xbox, ybox zbox "
	echo " Usage:   ${0}   state.sql " 
	exit 1
fi


if [ -e "${statefile}" ] ; then
	echo "Use:    statefile: ${statefile} "
else 
	echo "Error; the statefile does not exist:  ${statefile}"
	exit 1
fi


if [[ "${statefile: -4}" != ".sql" ]] ; then
	echo "Error; the kmc state file does not posses a .sql termination: ${statefile}"
	exit 1
fi



### boxsize in nm
xbox=$( sqlite3  ${statefile}  " SELECT box11 FROM frames " )
ybox=$( sqlite3  ${statefile}  " SELECT box22 FROM frames " )
zbox=$( sqlite3  ${statefile}  " SELECT box33 FROM frames " )


echo "Use boxsize   [nm]:  ${xbox}  ${ybox} ${zbox} "

### nm2Ang = 10.0
xbox=$( echo ${xbox} | awk '{ print $(1)*10.0 }' )
ybox=$( echo ${ybox} | awk '{ print $(1)*10.0 }' )
zbox=$( echo ${zbox} | awk '{ print $(1)*10.0 }' )
echo "Use boxsize  [Ang]:  ${xbox}  ${ybox} ${zbox} "
