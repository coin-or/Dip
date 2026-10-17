#!/bin/bash

NAME=$1 #example=milpblockp
OPT=$2  #example=-g or -O
VERSION=$3 #example=   or -10 (the latter uses cpx10.2)


#inferno
EXECUTABLE=/usr/local/bin/perl

PARAM_FILE=${HOME}/bin${VERSION}/decomp/${NAME}.parm
INIT_DIR=${HOME}/running${VERSION}/decomp${OPT}/milpblock/${NAME}
WRAP_DIR=${HOME}/running${VERSION}
ARGS="--param ${NAME}.parm"


DATA_DIR=${HOME}/COIN/coin-Dip/Dip/data
DRUG_DIR=${DATA_DIR}/MILP/block/drugport
RETAIL_DIR=${DATA_DIR}/MILP/block/retail
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
#inferno
#requirements = (Subnet == "192.168.3")


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

#Retail
for i in 3 4 6 20 22 27 31 33
  do
   echo "Generating queue commands for ${i}"
   filename=retail${i}.mps
   filenameB=retail${i}.block
   fileparam=${PARAM_FILE}
   echo "
     output      = retail${i}${OUT_SUFF}
     error       = retail${i}${ERR_SUFF}
     on_exit_remove = (ExitBySignal==FALSE) && (ExitCode == 0)
     on_exit_hold   = (ExitBySignal==TRUE) || (ExitCode != 0)
     arguments   = ${WRAP_DIR}/license_wrapper_milpblock.pl ${ARGS} --MILPBlock:Instance retail${i} --MILPBlock:BlockFile retail${i}.block
     should_transfer_files = YES
     transfer_input_files = ${RETAIL_DIR}/$filename, ${RETAIL_DIR}/$filenameB, $fileparam
     WhenToTransferOutput = ON_EXIT_OR_EVICT
     hold        = $HOLD

   queue 1

# -----------------------------------------------------------------------------
   ">>$OUT_FILE
   echo -n "."
   let count=$count+1
done

# #Drug Portfolio
# for i in 20 100 500
#   do
#    echo "Generating queue commands for ${i}"
#    filename=drug_port_${i}.mps
#    filenameB=drug_port_${i}.block
#    fileparam=${PARAM_FILE}
#    echo "
#      output      = drug_port_${i}${OUT_SUFF}
#      error       = drug_port_${i}${ERR_SUFF}
#      arguments   = "${ARGS} --MILPBlock:Instance drug_port_${i} --MILPBlock:BlockFile drug_port_${i}.block"
#      should_transfer_files = YES
#      transfer_input_files = ${DRUG_DIR}/$filename, ${DRUG_DIR}/$filenameB, $fileparam
#      WhenToTransferOutput = ON_EXIT_OR_EVICT
#      hold        = $HOLD

#    queue 1

# # -----------------------------------------------------------------------------
#    ">>$OUT_FILE
#    echo -n "."
#    let count=$count+1
# done
echo ""
echo ""
echo "$OUT_FILE has $count jobs. Use condor_submit to submit it."

