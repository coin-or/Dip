for i in atm gap milpblock mmkp
do
  ./createsub_lnx_${i}.sh ${i}p $1
  ./createsub_lnx_${i}.sh ${i}p2 $1
  ./createsub_lnx_${i}.sh ${i}c $1
  ./createsub_lnx_${i}.sh ${i}c2 $1
done
