#!/bin/bash

NAME=$1

#this stuff might change from run to run
EXECUTABLE=${HOME}/bin/decomp/alps_tsp
PARAM_FILE=${HOME}/bin/decomp/${NAME}.parm
INIT_DIR=${HOME}/running/decomp/tsp/${NAME}
ARGS="--param ${NAME}.parm"
INST_DIR=$HOME/COIN/coin-Decomp/Decomp/data/TSP/TSPLIB
OUT_FILE=condor.${NAME}.submit

#this stuff will not change often from run to run
UNIVERSE=vanilla
OUT_SUFF=".condor.out"
ERR_SUFF=".condor.err.\$(Process)"
HOLD=False
EXEC_STATUS=`ls -lt ${EXECUTABLE}`
MAX_SECONDS=600
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
requirements = (OpSys == \"LINUX\")
# executable: ${EXEC_STATUS}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
">>$OUT_FILE
let count=0
echo "Generating queue commands for TSPLIB"
for file in `cat $INST_DIR/tsplib.list`
   do     
     filename=${INST_DIR}/${file}.tsp
     filesol=${INST_DIR}/TSPLIB_opt
     fileparam=${PARAM_FILE}
     echo "
     output      = ${file}${OUT_SUFF}
     error       = ${file}${ERR_SUFF}
     arguments   = "${ARGS} --TSP:Instance ${file}"
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
echo
echo "$OUT_FILE has $count jobs. Use condor_submit to submit it."

