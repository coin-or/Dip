./parse.awk.retail.sh ../data/MILP/block/results/retail_pc_1.3600.txt 3600 
mv pc.3600 pc1.3600
./parse.awk.retail.sh ../data/MILP/block/results/retail_pc_2.3600.txt 3600
mv pc.3600 pc2.3600
./parse2.awk.retail.sh ../data/MILP/block/results/retail_cpx11.3600.txt 3600


#awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\\\\ \\hline"}' pc2.3600 > pc2.3600.2
#awk '{print $2 "\t" $3 "\t" $4 "\t" $5 "\t" $6 "\t" $7 "\\\\ \\hline"}' cpx.3600 > cpx.3600.2

#paste pc1.3600 pc2.3600.2 > t
#paste t  cpx.3600.2 > results.3600

#rm pc1.3600 pc2.3600 t cpx.3600 cpx.3600.2
