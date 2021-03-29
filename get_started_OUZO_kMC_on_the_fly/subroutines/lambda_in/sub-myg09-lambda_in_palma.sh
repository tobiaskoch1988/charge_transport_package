#!/bin/bash

 if [ $# -le 3 ]; then
    echo "Usage: $0  monomerA monomerB foldername zieldatei [waitfor_jobid]"
    exit 1
 fi

 if [ ! -d "$3" ]; then
    echo "Error: directory $1 doesn't exist!"
    exit 1
 fi


 KEY=".ssh/stahs"
 monomerA=${1}
 monomerB=${2}
 JOB=${3}
 NAM=${3}
 zieldatei=${4}

 DEPEND="#" 
 if [ $# -ge 5 ]
 then
      DEPEND="#PBS -W depend=afterany:$3"
 fi

 CURDIR=`pwd`
 WAITFOR='#'
 QUEUE='p0doltsi'  #####'default'
 NODE=1
 NCPUS=72
 MEM=72
 TIME=160
 echo PPN = ${NCPUS}
 echo MEM = ${MEM}

 qsub -V << EOF
#!/bin/bash
#PBS -l nodes=${NODE}:westmere:ppn=${NCPUS}
#PBS -l walltime=${TIME}:00:00
##PBS -l mem=${MEM}
#PBS -A p0doltsi
#PBS -q ${QUEUE}
#PBS -N ${NAM}
#PBS -o ${JOB}.out
#PBS -e ${JOB}.err
#PBS -j oe
##PBS -m ae
##PBS -M t_koch08@uni-muenster.de
${WAITFOR}


# activate verbose shell
 set -vx

 HOME=/home/t/t_koch08
 DATDIR=`pwd`  
 SCRATCH="/scratch_gss/tmp/t_koch08/g09.\${PBS_JOBID}"
 mkdir -p \${SCRATCH}

# Betrachtung des tail-Prozess
#echo "ssh palma001 \"ssh \${HOSTNAME} 'tail \${SCRATCH}/${JOB}.log -n 100 -f'\"" > tail.\${PBS_JOBID}.sh
#chmod 755 tail.\${PBS_JOBID}.sh
#scp tail.\${PBS_JOBID}.sh palma001:\${DATDIR}

#
# set g09 environment palma
#
  source /Applic.PALMA/gruppen/q0grimme/bin/set_pgi-9.0.sh

  g09root=/Applic.PALMA/gruppen/q0grimme/g09-D01/g09/ 
 #g09root=/data/shared/apps/chemie
 #g09root=/Applic.ZIV/gaussian        #Zivsmp-Pfad
  export GAUSS_EXEDIR=\$g09root
  export PATH=\$g09root:\$PATH
  export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$g09root/g09:/share/Applic.PALMA/compiler/intel/12.1/composerxe/lib/intel64



#shorten host name
S_HOST=\${PBS_O_HOST%.uni-muenster.de}


# ! need to use scp/ssh via gondor !
# create temporary input/output directory
#TMPINOUT="/exports/queuehome/m/mboec_02/queue_inout/g09.\${PBS_JOBID}"
#TMPINOUT="queue_inout/g09.\${PBS_JOBID}"
# ssh gondor mkdir -p \${TMPINOUT}

# create SCRATCH directory
# SCRATCH="/scratch/g09.\${PBS_JOBID}"
# mkdir -p \${SCRATCH}


 echo "Starting ${JOB} from ${CURDIR} on \${PARENT}"

# run JOB in SCRATCH directory on PARENT
# cd \${SCRATCH}
# cp -r \${DATDIR}/${JOB} .
cd \${DATDIR}

# fetch input data 
#cp -p ${CURDIR}/${JOB} .
# ssh gondor " scp -i ${KEY} -pr \${S_HOST}:${CURDIR}/${JOB} \${TMPINOUT}"
# ssh gondor " scp -i ${KEY} -p \${S_HOST}:${CURDIR}/${JOB%.com}.chk \${TMPINOUT}"

# scp -pr gondor:\${TMPINOUT}/* .

# set GAUSSIAN environment
# rev. D01 @gondor
#   export g09root=/opt/intel
#   export GAUSS_EXEDIR=/opt/intel/g09   
#   export PATH=/opt/intel/g09:$PATH

# set INTEL environment
#currently @gondor
#   export LD_LIBRARY_PATH=/opt/intel/lib/intel64:/opt/intel/mkl/lib/intel64
   

### Start run_lambda_in
cd ${JOB}
./run_lambda_in.sh ${monomerA} ${monomerB} ${JOB} ${zieldatei} "false" "true" "sub_g09"
cd ../
 
# put output data & clean up
# scp -pr \${SCRATCH}/* gondor:\${TMPINOUT}
# ssh gondor " scp -i ${KEY} -pr \${TMPINOUT}/* \${S_HOST}:${CURDIR}"
# ssh gondor " ssh -i ${KEY} \${S_HOST} \"cat ${CURDIR}/${JOB}/${zieldatei} >> ${CURDIR}/${zieldatei} \" "

#Zurueckkopieren der Daten in das Ausgangsverzeichnis
 cp -pr ${JOB} \${DATDIR}


##rm -r   \${SCRATCH}/*
cd ../
rmdir \${SCRATCH}

# Loeschen des tail Prozesses nach der Rechnung
rm \${DATDIR}/tail.\${PBS_JOBID}.sh

# retrieve PBS output
# scp -p /var/spool/torque/spool/\${PBS_JOBID}.OU gondor:\${TMPINOUT}
# scp -p /var/spool/torque/spool/\${PBS_JOBID}.ER gondor:\${TMPINOUT}

# ssh  gondor " scp -i ${KEY} \${TMPINOUT}/\${PBS_JOBID}.OU \${S_HOST}:${CURDIR}/${JOB} "
# ssh  gondor " scp -i ${KEY} \${TMPINOUT}/\${PBS_JOBID}.ER \${S_HOST}:${CURDIR}/${JOB} "

# ##ssh gondor " rm -r \${TMPINOUT}/* "
# ssh gondor " rmdir \${TMPINOUT} "


EOF

