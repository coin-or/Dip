OPT=$1       #example=-g or -O
VERSION=$2   #example=   or -10 (the latter uses CPA10.2)

cp report.awk ${HOME}/running${VERSION}/decomp${OPT}
cp report.root.awk ${HOME}/running${VERSION}/decomp${OPT}
cp report.direct.awk ${HOME}/running${VERSION}/decomp${OPT}

for i in mmkp
do
  for t in p p3 p4 p5 p6 c c2
  do
    echo ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}${t}        
    cd ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}${t}    
    ../../report.awk report.txt ${i}${t} 600
    ../../report.root.awk report.root.txt ${i}${t} 600
  done
  echo ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}d
  cd ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}d
  ../../report.direct.awk report.direct.txt ${i}d 600
done

