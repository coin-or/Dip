EXE=/home/magh/AFS/coin-Decomp/build/Decomp/examples/MILPBlock/decomp_milpblock
PARM=/home/magh/AFS/coin-Decomp/Decomp/data/MILP/block/retail.direct.parm

#for i in 3 4 6 20 22 27 31 33
for i in 33
do
  echo "Running " ${i}
${EXE} --param ${PARM} --MILPBlock:Instance retail${i} --MILPBlock:BlockFile retail${i}.block > retail${i}.direct.txt 2> retail${i}.direct.err
done
