EXE=/home/magh/AFS/COIN/coin-Dip/build-O/Dip/examples/MILPBlock/decomp_milpblock
PARM=/home/magh/AFS/COIN/coin-Dip/Dip/data/MILP/block/retail.pc.parm

for i in 3 4 6 20 22 27 31 33
do
  echo "Running " ${i}
${EXE} --param ${PARM} --MILPBlock:Instance retail${i} --MILPBlock:BlockFile retail${i}.block > log/retail${i}.pc.txt 2> log/retail${i}.pc.err
done
