#$1 = pc or cpm
EXE=/users/magala/COIN/coin-Dip/build-g/Dip/examples/MMKP/decomp_mmkp
PARM=/users/magala/COIN/coin-Dip/Dip/data/MMKP/gsimon/mmkp.$1.parm

PREFIX=10-5-5
for t in G-CL-DS G-CL-DU G-CL-DW G-CL-S G-CL-U G-CL-W G-CU-DS G-CU-DU G-CU-DW G-CU-S G-CU-U G-CU-W G-L-DSI G-L-DSUI G-L-DSU G-L-DS G-L-DU G-L-DW G-L-S G-L-U G-L-W G-R-DSI G-R-DSU G-R-S G-R-U G-R-W G-U-DS G-U-DU G-U-DW G-U-S G-U-U G-U-W
do
  for d in `seq 25`
  do
  echo "Running " ${t}.${d}
  ${EXE} --param ${PARM} --MMKP:DataDir ${t} --MMKP:Instance ${t}.${d} > gsimon_${PREFIX}_${t}.${d}.$1.txt 2> gsimon_${PREFIX}_${t}.${d}.$1.err 
  done
done
