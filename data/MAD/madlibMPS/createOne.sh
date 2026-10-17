#!/bin/bash

EXE=~/COIN/coin-Decomp/build/Decomp/examples/MAD/alps_mad
DATADIR=..
DATASUBDIR=`echo $1 | awk -F/ '{print $1}'`
INSTANCE=`echo $1 | awk -F/ '{print $2}'`
NUMBLOCKS=$2

$EXE --MAD:DataDir $DATADIR --MAD:DataSubDir $DATASUBDIR --MAD:Instance $INSTANCE --MAD:NumBlocks $NUMBLOCKS



