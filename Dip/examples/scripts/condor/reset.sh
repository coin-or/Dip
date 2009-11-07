OPT=$1       #example=-g or -O
VERSION=$2   #example=   or -10 (the latter uses cpx10.2)

rm -f ${HOME}/running${VERSION}/decomp${OPT}
mkdir ${HOME}/running${VERSION}/decomp${OPT}
for i in atm gap milpblock mmkp 
do
  mkdir ${HOME}/running${VERSION}/decomp${OPT}/${i}
  mkdir ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}p
  mkdir ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}p2
  mkdir ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}p3
  mkdir ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}c
  mkdir ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}c2
  mkdir ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}d
done
