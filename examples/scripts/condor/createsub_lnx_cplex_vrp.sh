#!/bin/bash

NAME=$1

#this stuff might change from run to run
EXECUTABLE=/usr/bin/perl
WRAPPER=${HOME}/bin/decomp/license_wrapper_vrp.pl
PARAM_FILE=${HOME}/bin/decomp/${NAME}.parm
INIT_DIR=${HOME}/running/decomp/vrp/${NAME}
ARGS="${WRAPPER} --param ${NAME}.parm"
INST_DIR=$HOME/COIN/coin-Decomp/Decomp/data/VRP/vrplib
OUT_FILE=condor.${NAME}.submit

#this stuff will not change often from run to run
UNIVERSE=vanilla
OUT_SUFF=".condor.out"
ERR_SUFF=".condor.err.\$(Process)"
HOLD=False
EXEC_STATUS=`ls -lt ${EXECUTABLE}`
MAX_SECONDS=1200
HOLD_SECONDS=600
LOG=${INIT_DIR}/condor.log
EXTEN=""



echo "
getenv       = TRUE
universe     = $UNIVERSE
notification = Error
log          = ${LOG}
notify_user  = matthew.galati@sas.com
executable   = ${EXECUTABLE}
initialdir   = ${INIT_DIR}
MAX_SECONDS  = ${MAX_SECONDS}
HOLD_SECONDS  = ${HOLD_SECONDS}
periodic_remove = (JobStatus == 2) && (CurrentTime - EnteredCurrentStatus > ${MAX_SECONDS})
periodic_release = (JobStatus == 5) && (CurrentTime - EnteredCurrentStatus > ${HOLD_SECONDS})
requirements = (OpSys == \"LINUX\")
# executable: ${EXEC_STATUS}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
">>$OUT_FILE
let count=0
echo "Generating queue commands for VRPLIB"
while read line;
   do     
     instance=`echo ${line} | awk '{print $1}'`
     nRoutes=`echo ${line} | awk '{print $3}'`
     echo "instance=" $instance
     echo "nRoutes =" $nRoutes

     filename=${INST_DIR}/${instance}.vrp
     filesol=${INST_DIR}/vrplib.opt
     fileparam=${PARAM_FILE}
     echo "
     output      = ${instance}${OUT_SUFF}
     error       = ${instance}${ERR_SUFF}
     arguments   = "${ARGS} --VRP:Instance ${instance} --VRP:NumRoutes ${nRoutes}"
     should_transfer_files = YES
     transfer_input_files = $filename, $filesol, $fileparam
     WhenToTransferOutput = ON_EXIT_OR_EVICT
     hold        = $HOLD
     transfer_executable = False
     queue 1

# -----------------------------------------------------------------------------
      ">>$OUT_FILE
      echo -n "."
      let count=$count+1
   done < $INST_DIR/vrplib.opt
   echo ""
echo
echo "$OUT_FILE has $count jobs. Use condor_submit to submit it."

