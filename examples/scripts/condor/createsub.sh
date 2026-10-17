#c  = cpm
#c2 = cpm+dc
#p  = pc
#p2 = pc+solveMasterAsIp 
#p3 = pc+solverMasterAsIp+InitVarsWithIP(10s)
#d  = direct

EXT=_inferno

OPT=$1       #example=-g or -O
VERSION=$2   #example=   or -10 (the latter uses CPA10.2)


for i in atm milpblock gap
do
  ./createsub_lnx_${i}${EXT}.sh ${i}p ${OPT} ${VERSION}
  ./createsub_lnx_${i}${EXT}.sh ${i}c ${OPT} ${VERSION}
  ./createsub_lnx_${i}${EXT}.sh ${i}c2 ${OPT} ${VERSION}
  ./createsub_lnx_${i}${EXT}.sh ${i}d ${OPT} ${VERSION}
  ./createsub_lnx_${i}${EXT}.sh ${i}p2 ${OPT} ${VERSION}
done
for i in mmkp
do
  ./createsub_lnx_${i}${EXT}.sh ${i}p ${OPT} ${VERSION}
  ./createsub_lnx_${i}${EXT}.sh ${i}c ${OPT} ${VERSION}
  ./createsub_lnx_${i}${EXT}.sh ${i}c2 ${OPT} ${VERSION}
  ./createsub_lnx_${i}${EXT}.sh ${i}d ${OPT} ${VERSION}
done
for i in atm milpblock gap mmkp
do
  ./createsub_lnx_${i}${EXT}.sh ${i}p3 ${OPT} ${VERSION}
done


# for i in #vrp
# do
#   ./createsub_lnx_${i}.sh ${i}p ${OPT} ${VERSION}
#   ./createsub_lnx_${i}.sh ${i}p2 ${OPT} ${VERSION}
#   ./createsub_lnx_${i}.sh ${i}c ${OPT} ${VERSION}
#   ./createsub_lnx_${i}.sh ${i}c2 ${OPT} ${VERSION}
# done
