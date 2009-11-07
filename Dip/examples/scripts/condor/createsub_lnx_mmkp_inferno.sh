#!/bin/bash

NAME=$1 #example=mmkpp
OPT=$2  #example=-g or -O
VERSION=$3 #example=   or -10 (the latter uses cpx10.2)

#inferno
EXECUTABLE=/usr/local/bin/perl

PARAM_FILE=${HOME}/bin${VERSION}/decomp/${NAME}.parm
INIT_DIR=${HOME}/running${VERSION}/decomp${OPT}/mmkp/${NAME}
WRAP_DIR=${HOME}/running${VERSION}
ARGS="--param ${NAME}.parm"


DATA_DIR=${HOME}/COIN/coin-Dip/Dip/data/MMKP/hifi
OUT_FILE=condor${OPT}.${NAME}.submit

rm ${OUT_FILE}

#this stuff will not change often from run to run
UNIVERSE=vanilla
OUT_SUFF=".condor.out"
ERR_SUFF=".condor.err.\$(Process)"
HOLD=False
EXEC_STATUS=`ls -lt ${EXECUTABLE}`
MAX_SECONDS=3600
HOLD_SECONDS=300
LOG=${INIT_DIR}condor.log
EXTEN=""

#orclus
#requirements = (Arch == \"X86_64\") && (OpSys == \"LINUX\")
#notify_user  = matthew.galati@sas.com
#inferno
#requirements = (Subnet == "192.168.3")
#notify_user  = magh@lehigh.edu

#release the job on hold after so many minutes
#periodic_release = (JobStatus == 5) && (CurrentTime - EnteredCurrentStatus > ${HOLD_SECONDS})


echo "
#get your environment
getenv       = TRUE
universe     = $UNIVERSE
notification = Error
log          = ${LOG}
notify_user  = magh@lehigh.edu
executable   = ${EXECUTABLE}
initialdir   = ${INIT_DIR}

requirements = (Arch == \"X86_64\") && (subnet == \"192.168.3\") && (Memory >= 0)

#remove job after so many minutes anyway.
periodic_remove = (JobStatus == 2) && (CurrentTime - EnteredCurrentStatus > ${MAX_SECONDS})

#release the job on hold after so many minutes
periodic_release = (JobStatus == 5) && (CurrentTime - EnteredCurrentStatus > ${HOLD_SECONDS})


# executable: ${EXEC_STATUS}
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
">>$OUT_FILE

let count=0


for i in `seq 13`
  do
   echo "Generating queue commands for I${i}"
   filename=I${i}
   fileparam=${PARAM_FILE}
   echo "
     output      = I${i}${OUT_SUFF}
     error       = I${i}${ERR_SUFF}
     arguments   = ${WRAP_DIR}/license_wrapper_mmkp.pl ${ARGS} --MMKP:Instance I${i}
     on_exit_remove = (ExitBySignal==FALSE) && (ExitCode == 0)
     on_exit_hold   = (ExitBySignal==TRUE) || (ExitCode != 0)
     should_transfer_files = YES
     transfer_input_files = ${DATA_DIR}/$filename, $fileparam
     WhenToTransferOutput = ON_EXIT_OR_EVICT
     hold        = $HOLD
     transfer_executable = False

   queue 1

# -----------------------------------------------------------------------------
   ">>$OUT_FILE
   echo -n "."
   let count=$count+1
done


for i in `seq 9`
  do
   echo "Generating queue commands for INST0${i}"
   filename=INST0${i}
   fileparam=${PARAM_FILE}
   echo "
     output      = INST0${i}${OUT_SUFF}
     error       = INST0${i}${ERR_SUFF}
     arguments   = ${WRAP_DIR}/license_wrapper_mmkp.pl ${ARGS} --MMKP:Instance INST0${i}
     on_exit_remove = (ExitBySignal==FALSE) && (ExitCode == 0)
     on_exit_hold   = (ExitBySignal==TRUE) || (ExitCode != 0)
     should_transfer_files = YES
     transfer_input_files = ${DATA_DIR}/$filename, $fileparam
     WhenToTransferOutput = ON_EXIT_OR_EVICT
     hold        = $HOLD
     transfer_executable = False

   queue 1

# -----------------------------------------------------------------------------
   ">>$OUT_FILE
   echo -n "."
   let count=$count+1
done

for i in `seq 11 20`
  do
   echo "Generating queue commands for INST${i}"
   filename=INST${i}
   fileparam=${PARAM_FILE}
   echo "
     output      = INST${i}${OUT_SUFF}
     error       = INST${i}${ERR_SUFF}
     arguments   = ${WRAP_DIR}/license_wrapper_mmkp.pl ${ARGS} --MMKP:Instance INST${i}
     on_exit_remove = (ExitBySignal==FALSE) && (ExitCode == 0)
     on_exit_hold   = (ExitBySignal==TRUE) || (ExitCode != 0)
     should_transfer_files = YES
     transfer_input_files = ${DATA_DIR}/$filename, $fileparam
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
echo ""
echo "$OUT_FILE has $count jobs. Use condor_submit to submit it."

