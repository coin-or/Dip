#!/bin/bash
PREFIX=10-5-5
for t in G-CL-DS G-CL-DU G-CL-DW G-CL-S G-CL-U G-CL-W G-CU-DS G-CU-DU G-CU-DW G-CU-S G-CU-U G-CU-W G-L-DSI G-L-DSUI G-L-DSU G-L-DS G-L-DU G-L-DW G-L-S G-L-U G-L-W G-R-DSI G-R-DSU G-R-S G-R-U G-R-W G-U-DS G-U-DU G-U-DW G-U-S G-U-U G-U-W
do   
  v=1
  echo "Running " ${t} " generator"  
# to get 25 instances each
#   if random elemebts, instances=5, con_levels=5
#   else                             con_level=25
  `sed 's/instances="20"/instances="5"/' ${PREFIX}-${t}.xml > tmp.xml`
  `sed 's/constraint_levels="100"/constraint_levels="5"/' tmp.xml > tmp2.xml`
  python gen_mmkp.py --def_file tmp2.xml
  mkdir ${t}
  for k in *data
  do
    mv ${k} ${t}/${t}.${v}
    v=$((v+1))
  done
done
rm tmp.xml
rm tmp2.xml