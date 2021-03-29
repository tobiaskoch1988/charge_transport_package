#!/bin/bash

 if [ $# -le 3 ]; then
    echo "Usage: $0  monomerA monomerB foldername zieldatei oniom [waitfor_jobid]"
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
 oniom=${5}
  if [ $# -ge 5 ]; then
    oniom=${5}
  else
    oniom='false'
  fi
  
 DEPEND="#" 
 if [ $# -ge 6 ]
 then
      DEPEND="#PBS -W depend=afterany:$6"
 fi

 CURDIR=`pwd`


 PPN=8
 MEM=8
 echo PPN = ${PPN}
 echo MEM = ${MEM}

 qsub -V << EOF
#!/bin/bash
#PBS -l nodes=1:ppn=${PPN}
#PBS -l walltime=00:48:00:00
##PBS -l mem=${MEM}
#PBS -q batch
#PBS -N ${NAM}
##PBS -e ${HOME}/QSUBIO
##PBS -o ${HOME}/QSUBIO
##PBS -m e
${DEPEND}
#

#shorten host name
S_HOST=\${PBS_O_HOST%.uni-muenster.de}

# activate verbose shell
# set -vx

# ! need to use scp/ssh via gondor !
# create temporary input/output directory
#TMPINOUT="/exports/queuehome/m/mboec_02/queue_inout/g09.\${PBS_JOBID}"
 TMPINOUT="queue_inout/g09.\${PBS_JOBID}"
 ssh gondor mkdir -p \${TMPINOUT}

# create SCRATCH directory
 SCRATCH="/scratch/g09.\${PBS_JOBID}"
 mkdir -p \${SCRATCH}


 echo "Starting ${JOB} from ${CURDIR} on \${PARENT}"

# run JOB in SCRATCH directory on PARENT
 cd \${SCRATCH}
 #touch _START_


# fetch input data 
#cp -p ${CURDIR}/${JOB} .
 ssh gondor " scp -i ${KEY} -pr \${S_HOST}:${CURDIR}/${JOB} \${TMPINOUT}"
# ssh gondor " scp -i ${KEY} -p \${S_HOST}:${CURDIR}/${JOB%.com}.chk \${TMPINOUT}"

 scp -pr gondor:\${TMPINOUT}/* .

# set GAUSSIAN environment
 # rev. D01 @gondor
   export g09root=/opt/intel
   export GAUSS_EXEDIR=/opt/intel/g09   
   export PATH=/opt/intel/g09:$PATH

# set INTEL environment
 #currently @gondor
   export LD_LIBRARY_PATH=/opt/intel/lib/intel64:/opt/intel/mkl/lib/intel64
   

### Start run_lambda_in
cd ${JOB}
./run_lambda_in.sh ${monomerA} ${monomerB} ${JOB} ${zieldatei} ${oniom}
cd ../
 
# put output data & clean up
 scp -pr \${SCRATCH}/* gondor:\${TMPINOUT}
 ssh gondor " scp -i ${KEY} -pr \${TMPINOUT}/* \${S_HOST}:${CURDIR}"
 ssh gondor " ssh -i ${KEY} \${S_HOST} \"cat ${CURDIR}/${JOB}/${zieldatei} >> ${CURDIR}/${zieldatei} \" "


##rm -r   \${SCRATCH}/*
rmdir \${SCRATCH}


# retrieve PBS output
 scp -p /var/spool/torque/spool/\${PBS_JOBID}.OU gondor:\${TMPINOUT}
 scp -p /var/spool/torque/spool/\${PBS_JOBID}.ER gondor:\${TMPINOUT}

 ssh  gondor " scp -i ${KEY} \${TMPINOUT}/\${PBS_JOBID}.OU \${S_HOST}:${CURDIR}/${JOB} "
 ssh  gondor " scp -i ${KEY} \${TMPINOUT}/\${PBS_JOBID}.ER \${S_HOST}:${CURDIR}/${JOB} "

 ##ssh gondor " rm -r \${TMPINOUT}/* "
 ssh gondor " rmdir \${TMPINOUT} "


EOF

