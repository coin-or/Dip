#c  = cpm
#c2 = cpm+dc
#p  = pc
#p  = pc+solveMasterAsIp 
#d  = direct

EXT=_inferno

for i in atm milpblock gap
do
  ./createsub_lnx_${i}${EXT}.sh ${i}p $1
  ./createsub_lnx_${i}${EXT}.sh ${i}c $1
  ./createsub_lnx_${i}${EXT}.sh ${i}c2 $1
  ./createsub_lnx_${i}${EXT}.sh ${i}d $1
  ./createsub_lnx_${i}${EXT}.sh ${i}p2 $1
done
for i in mmkp
do
  ./createsub_lnx_${i}${EXT}.sh ${i}p $1
  ./createsub_lnx_${i}${EXT}.sh ${i}c $1
  ./createsub_lnx_${i}${EXT}.sh ${i}c2 $1
  ./createsub_lnx_${i}${EXT}.sh ${i}d $1
done
# for i in #vrp
# do
#   ./createsub_lnx_${i}.sh ${i}p $1
#   ./createsub_lnx_${i}.sh ${i}p2 $1
#   ./createsub_lnx_${i}.sh ${i}c $1
#   ./createsub_lnx_${i}.sh ${i}c2 $1
# done
