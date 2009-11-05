#!/bin/bash

NAME=$1 #example=gapp
OPT=$2  #example=-g or -O

EXECUTABLE=${HOME}/bin/decomp/decomp_gap${OPT}
PARAM_FILE=${HOME}/bin/decomp/${NAME}.parm
INIT_DIR=${HOME}/running/decomp${OPT}/gap/${NAME}
ARGS="--param ${NAME}.parm"


DATA_DIR=${HOME}/COIN/coin-Dip/Dip/data/GAP
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


for i in a05100  b10100  c10400 c40400 e10400 e40400 gap0515-5 a05200  b10200  c15900 c60900 e15900 e60900 gap0520-1 a10100  b20100  c20100 c801600 e20100 e801600 a10200  b20200  c201600 e201600 a20100  c05100  c20200 e05100 e20200 gap0515-1 a20200 c05200 c20400 e05200 e20400 gap0515-2 b05100 c10100 c30900  e10100 e30900 gap0515-3 b05200  c10200  c401600  e10200 e401600 gap0515-4
do
   echo "Generating queue commands for ${i}"
   filename=${i}
   fileparam=${PARAM_FILE}
   echo "
     output      = ${i}${OUT_SUFF}
     error       = ${i}${ERR_SUFF}
     arguments   = "${ARGS} --GAP:Instance ${i}"
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

