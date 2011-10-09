#!/bin/bash

NAME=$1 #example=mmkpp
OPT=$2  #example=-g or -O

EXECUTABLE=${HOME}/bin/decomp/decomp_mmkp${OPT}
PARAM_FILE=${HOME}/bin/decomp/${NAME}.parm
INIT_DIR=${HOME}/running/decomp${OPT}/mmkp/${NAME}
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
MAX_SECONDS=1000
LOG=${INIT_DIR}condor.log
EXTEN=""

#orclus
#requirements = (Arch == \"X86_64\") && (OpSys == \"LINUX\")
#inferno
#requirements = (Subnet == "192.168.3")

echo "
universe     = $UNIVERSE
notification = Error
log          = ${LOG}
notify_user  = matthew.galati@sas.com
executable   = ${EXECUTABLE}
initialdir   = ${INIT_DIR}
requirements = (Arch == \"X86_64\") && (OpSys == \"LINUX\")
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
     arguments   = "${ARGS} --MMKP:Instance I${i}"
     should_transfer_files = YES
     transfer_input_files = ${DATA_DIR}/$filename, $fileparam
     WhenToTransferOutput = ON_EXIT_OR_EVICT
     hold        = $HOLD

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
     arguments   = "${ARGS} --MMKP:Instance INST0${i}"
     should_transfer_files = YES
     transfer_input_files = ${DATA_DIR}/$filename, $fileparam
     WhenToTransferOutput = ON_EXIT_OR_EVICT
     hold        = $HOLD

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
     arguments   = "${ARGS} --MMKP:Instance INST${i}"
     should_transfer_files = YES
     transfer_input_files = ${DATA_DIR}/$filename, $fileparam
     WhenToTransferOutput = ON_EXIT_OR_EVICT
     hold        = $HOLD

   queue 1

# -----------------------------------------------------------------------------
   ">>$OUT_FILE
   echo -n "."
   let count=$count+1
done

echo ""
echo ""
echo "$OUT_FILE has $count jobs. Use condor_submit to submit it."

