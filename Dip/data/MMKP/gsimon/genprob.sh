PREFIX=10-5-5
for t in G-CL-DS G-CL-DU G-CL-DW G-CL-S G-CL-U G-CL-W G-CU-DS G-CU-DU G-CU-DW G-CU-S G-CU-U G-CU-W G-L-DSI G-L-DSUI G-L-DSU G-L-DS G-L-DU G-L-DW G-L-S G-L-U G-L-W G-R-DSI G-R-DSU G-R-S G-R-U G-R-W G-U-DS G-U-DU G-U-DW G-U-S G-U-U G-U-W
do
  echo "Running " ${t} " generator"
  `sed 's/instances="20"/instances="5"/' ${PREFIX}-${t}.xml > tmp.xml`
  `sed 's/constraint_levels="100"/constraint_levels="5"/' tmp.xml > tmp2.xml`
  python gen_mmkp.py --def_file tmp2.xml
  mkdir ${t}
  mv *.data ${t}
done
