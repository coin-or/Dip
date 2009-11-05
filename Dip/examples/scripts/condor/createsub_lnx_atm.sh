#!/bin/bash

NAME=$1 #example=atmp
OPT=$2  #example=-g or -O

#this stuff might change from run to run
EXECUTABLE=${HOME}/bin/decomp/decomp_atm${OPT}
PARAM_FILE=${HOME}/bin/decomp/${NAME}.parm
INIT_DIR=${HOME}/running/decomp${OPT}/atm/${NAME}
ARGS="--param ${NAME}.parm"
INST_DIR=$HOME/COIN/coin-Dip/Dip/data/ATM
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

#MAX_SECONDS  = ${MAX_SECONDS}
#periodic_remove = (JobStatus == 2) && (CurrentTime - EnteredCurrentStatus > ${MAX_SECONDS})

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
for ad in 5_25 5_50 10_50 10_100 20_100 20_200
do
  for i in {1..5}
  do
   echo "Generating queue commands for ${ad}_${i}"
   filenameA=atm_randA_${ad}_${i}.txt
   filenameD=atm_randD_${ad}_${i}.txt
   filenameAD=atm_randAD_${ad}_${i}.txt
   fileparam=${PARAM_FILE}
   echo "
     output      = ${ad}_${i}${OUT_SUFF}
     error       = ${ad}_${i}${ERR_SUFF}
     arguments   = "${ARGS} --ATM:DataAtm ${filenameA} --ATM:DataDate ${filenameD} --ATM:DataAtmDate ${filenameAD}"
     should_transfer_files = YES
     transfer_input_files = ${INST_DIR}/$filenameA, ${INST_DIR}/$filenameD, ${INST_DIR}/$filenameAD, $fileparam
     WhenToTransferOutput = ON_EXIT_OR_EVICT
     hold        = $HOLD

   queue 1

# -----------------------------------------------------------------------------
   ">>$OUT_FILE
   echo -n "."
   let count=$count+1
   done
done
echo ""
echo ""
echo "$OUT_FILE has $count jobs. Use condor_submit to submit it."

