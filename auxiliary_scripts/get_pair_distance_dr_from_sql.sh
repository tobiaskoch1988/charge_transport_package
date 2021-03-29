#!/bin/bash
if [ $# -ne 1 ]; then
	echo " Programm writes distances from sql file to a data list "
	echo " Usage:   $0      sqlfilename  "
	echo "          $0         state.sql "
	exit 1
fi


sqlfilename=${1}
echo "${sqlfilename}"
if [[ "${sqlfilename: -4}" != ".sql" ]] ; then
	echo "Error; the kmc state file does not posses a .sql termination: ${sqlfilename}"
	exit 1
fi
outputfilename=dr_Ang_pairs_${sqlfilename}.dat
if [ -e "${outputfilename}" ] ; then
	echo "File already exits. Make sure it is not overwritten.! ${outputfilename} "
	exit 1
fi


xbox=$( sqlite3 ${sqlfilename} " SELECT box11 FROM frames " )
ybox=$( sqlite3 ${sqlfilename} " SELECT box22 FROM frames " )
zbox=$( sqlite3 ${sqlfilename} " SELECT box33 FROM frames " )


echo "BOX ${xbox} ${ybox} ${zbox}  [nm] in sqlfile:     ${sqlfilename}"
sqlite3  ${sqlfilename}  " SELECT drX,drY,drZ FROM pairs " | sed "s/|/ /g"  |  awk -v xbox="${xbox}" -v ybox="${ybox}" -v zbox="${zbox}"  'BEGIN { FS = " " ; nm2Ang=10.0 
xbox=xbox*nm2Ang
ybox=ybox*nm2Ang
zbox=zbox*nm2Ang
}{ 

drX=$(1)*nm2Ang
drY=$(2)*nm2Ang
drZ=$(3)*nm2Ang
## print drX,drY,drZ
##### correct periodic boundary conditions r > rbox
if( drX > xbox) { drX=drX-xbox }
if( drY > ybox) { drY=drY-ybox }
if( drZ > zbox) { drZ=drZ-zbox }

### correct  periodic boundary conditions r < rbox
#if( (-1.0*drX) > xbox) { drX=drX+xbox }
#if( (-1.0*drY) > ybox) { drY=drY+ybox }
#if( (-1.0*drZ) > zbox) { drZ=drZ+zbox }

dr=sqrt(drX*drX + drY*drY + drZ*drZ)
print drX, drY, drZ, dr
 }' > ${outputfilename}
