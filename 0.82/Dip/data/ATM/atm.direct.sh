EXE=/home/magh/AFS/coin-Decomp/build-g/Decomp/examples/ATM/decomp_atm
PARM=/home/magh/AFS/coin-Decomp/Decomp/data/ATM/atm.direct.parm

for ad in 5_25 5_50 10_50 10_100 20_100 20_200
do
  for i in {1..5}
  do
  echo "Running " ${ad}_${i}
${EXE} --param ${PARM} --ATM:DataAtm atm_randA_${ad}_${i}.txt --ATM:DataDate atm_randD_${ad}_${i}.txt --ATM:DataAtmDate atm_randAD_${ad}_${i}.txt > atm_${ad}_${i}.direct.txt 2> atm_${ad}_${i}.direct.err
  done
done
