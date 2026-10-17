#!/bin/bash

NAME=$1

#this stuff might change from run to run
#EXECUTABLE=${HOME}/bin/decomp/alps_mad
EXECUTABLE=/usr/bin/perl
WRAPPER=${HOME}/bin/decomp/license_wrapper.pl
PARAM_FILE=${HOME}/bin/decomp/${NAME}.parm
INIT_DIR=${HOME}/running/decomp/mad/${NAME}
ARGS="${WRAPPER} --param ${NAME}.parm"
INST_DIR=$HOME/COIN/coin-Decomp/Decomp/data/MAD/
instance_sets=miplib
OUT_FILE=condor.${NAME}.submit

#this stuff will not change often from run to run
UNIVERSE=vanilla
OUT_SUFF=".condor.out"
ERR_SUFF=".condor.err.\$(Process)"
HOLD=False
EXEC_STATUS=`ls -lt ${EXECUTABLE}`
MAX_SECONDS=3600
HOLD_SECONDS=600
#TIME_LIMIT=600
LOG=${INIT_DIR}/condor.log
EXTEN=""



echo "
getenv       = TRUE
universe     = $UNIVERSE
notification = Error
log          = ${LOG}
notify_user  = magh@lehigh.edu
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
for i in $instance_sets
do
   echo "Generating queue commands for $i"
   for file in `cat $INST_DIR/$i/madlib.$i.list`
   do     
     filename=${INST_DIR}/$i/${file}.p.lp
     filesol=${INST_DIR}/$i/madlib.$i.opt
     fileparam=${PARAM_FILE}
     echo "
     output      = ${file}${OUT_SUFF}
     error       = ${file}${ERR_SUFF}
     arguments   = "${ARGS} --MAD:Instance ${file}"
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
   done
   echo ""
done
echo
echo "$OUT_FILE has $count jobs. Use condor_submit to submit it."

