./parse.awk.sh ../data/ATM/results/atm_pc_1.3600.txt 3600 1    #price+branch
./parse.awk.sh ../data/ATM/results/atm_pc_2.3600.txt 3600 2    #pc
./parse2.awk.sh ../data/ATM/results/atm_cpx11.3600.txt 3600    


awk '{print $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12}' pc.3600.2 > tmp.2
awk '{print $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\\\\ \\hline"}' pc.3600.1 > tmp.1

paste cpx.3600 tmp.2 tmp.1 > results.3600

