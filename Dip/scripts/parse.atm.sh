./parse.awk.sh ../data/ATM/results/atm_pc1.3600.txt 3600
./parse.awk.sh ../data/ATM/results/atm_pc1.300.txt 300
./parse2.awk.sh ../data/ATM/results/atm_cpx11.3600.txt 3600
./parse2.awk.sh ../data/ATM/results/atm_cpx11.300.txt 300

awk '{print $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\\\\ \\hline"}' cpx.3600 > cpx.3600.2
awk '{print $6 "\t" $7 "\t" $8 "\t" $9 "\t" $10 "\t" $11 "\t" $12 "\\\\ \\hline"}' cpx.300 > cpx.300.2

paste pc.300 cpx.300.2 > results.300
paste pc.3600 cpx.3600.2 > results.3600


rm pc.3600 pc.300
rm cpx.3600 cpx.300
rm cpx.3600.2 cpx.300.2
