OPT=$1          #example=-g or -O
VERSION=$2      #example=   or -10 (the latter uses CPA10.2)


DIR=${HOME}/running${VERSION}/decomp${OPT}
DIRM=${DIR}/mmkp/mmkp

#direct vs cpm
join -1 1 -2 1 ${DIRM}d/mmkpd.600.direct ${DIRM}c/mmkpc.600 > tmp1
#pc vs cpm+dc
join -1 1 -2 1 ${DIRM}p/mmkpp.600 ${DIRM}c2/mmkpc2.600      > tmp2
#hybrid vs pc-nested(m2kp0)
join -1 1 -2 1 ${DIRM}p3/mmkpp3.600 ${DIRM}p4/mmkpp4.600      > tmp3
#pc-nested(mmkp)
cat ${DIRM}p5/mmkpp5.600 > tmp4


#11 columns
#instance  1                        2               
#          time, lb, ub, gap, nodes time, lb, ub, gap, nodes 
awk '{ print $0 " \\\\ \\hline";}' tmp  > mmkp1.tex
awk '{ print $0 " \\\\ \\hline";}' tmp2 > mmkp2.tex
awk '{ print $0 " \\\\ \\hline";}' tmp3 > mmkp3.tex
awk '{ print $0 " \\\\ \\hline";}' tmp4 > mmkp4.tex
