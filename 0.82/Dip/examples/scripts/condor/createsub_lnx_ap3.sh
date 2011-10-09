#!/bin/bash

NAME=ap3p

#this stuff might change from run to run
EXECUTABLE=${HOME}/bin/decomp/alps_ap3
PARAM_FILE=${HOME}/bin/decomp/${NAME}.parm
INIT_DIR=${HOME}/running/decomp/ap3/${NAME}
ARGS="--param ${NAME}.parm"
INST_DIR=$HOME/COIN/coin-Decomp/Decomp/data/AP3/Grundel
instance_sets=ap3
OUT_FILE=condor.${NAME}.submit

#this stuff will not change often from run to run
UNIVERSE=vanilla
OUT_SUFF=".condor.out"
ERR_SUFF=".condor.err.\$(Process)"
HOLD=False
EXEC_STATUS=`ls -lt ${EXECUTABLE}`
MAX_SECONDS=3600
#TIME_LIMIT=600
LOG=${INIT_DIR}/condor.log
EXTEN=""



echo "
universe     = $UNIVERSE
notification = Error
log          = ${LOG}
notify_user  = matthew.galati@sas.com
executable   = ${EXECUTABLE}
initialdir   = ${INIT_DIR}
MAX_SECONDS  = ${MAX_SECONDS}
periodic_remove = (JobStatus == 2) && (CurrentTime - EnteredCurrentStatus > ${MAX_SECONDS})
requirements = (Arch == \"X86_64\") && (OpSys == \"LINUX\")
# executable: ${EXEC_STATUS}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
">>$OUT_FILE
let count=0
for i in $instance_sets
do
   echo "Generating queue commands for $i"
   for file in `cat $INST_DIR/$i.list`
   do     
     filename=${INST_DIR}/${file}.txt
     filesol=${INST_DIR}/${file}.sol
     fileparam=${PARAM_FILE}
     echo "
     output      = ${file}${OUT_SUFF}
     error       = ${file}${ERR_SUFF}
     arguments   = "${ARGS} --AP3:Instance ${file}"
     should_transfer_files = YES
     transfer_input_files = $filename, $filesol, $fileparam
     WhenToTransferOutput = ON_EXIT_OR_EVICT
     hold        = $HOLD

     queue 1

# -----------------------------------------------------------------------------
      ">>$OUT_FILE
      echo -n "."
      let count=$count+1
   done
   echo ""
done
echo
echo "$OUT_FILE has $count jobs. Use condor_submit to submit it."

