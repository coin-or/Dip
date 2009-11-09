EXE=/home/magh/AFS/COIN/coin-Dip/build-O/Dip/examples/ATM/decomp_atm
PARM=/home/magh/AFS/COIN/coin-Dip/Dip/data/ATM/atm.pc.parm

#20_200
for ad in 5_25 5_50 10_50 10_100 20_100
do
  for i in {1..5}
  do
  echo "Running " ${ad}_${i}
${EXE} --param ${PARM} --ATM:DataAtm atm_randA_${ad}_${i}.txt --ATM:DataDate atm_randD_${ad}_${i}.txt --ATM:DataAtmDate atm_randAD_${ad}_${i}.txt > log/atm_${ad}_${i}.pc.txt 2> log/atm_${ad}_${i}.pc.err
  done
done
