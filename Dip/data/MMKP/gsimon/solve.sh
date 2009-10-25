EXE=/home/magh/AFS/coin-Decomp/build-g/Dip/examples/MMKP/decomp_mmkp
PARM=/home/magh/AFS/coin-Decomp/Dip/data/MMKP/gsimon/mmkp.pc.parm

PREFIX=10-5-5
for t in G-CL-DS G-CL-DU G-CL-DW G-CL-S G-CL-U G-CL-W G-CU-DS G-CU-DU G-CU-DW G-CU-S G-CU-U G-CU-W G-L-DSI G-L-DSUI G-L-DSU G-L-DS G-L-DU G-L-DW G-L-S G-L-U G-L-W G-R-DSI G-R-DSU G-R-S G-R-U G-R-W G-U-DS G-U-DU G-U-DW G-U-S G-U-U G-U-W
do
  echo "Running " ${t}
${EXE} --param ${PARM} --MMKP:DataDir ${t} --MMKP:Instance 

> gsimon_${PREFIX}_${t}.pc.txt 2> gsimon_${PREFIX}_${t}.pc.err 
done
