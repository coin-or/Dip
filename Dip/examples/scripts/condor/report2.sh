for i in atm gap milpblock mmkp
do
  for t in p p2 c c2
  do
    echo ${HOME}/running/decomp${1}/${i}/${i}${t}        
    cp report.awk ${HOME}/running/decomp${1}/${i}/${i}${t}/
    cd ${HOME}/running/decomp${1}/${i}/${i}${t}    
    ./report.awk report.out ${i}${t} 600
  done
done
