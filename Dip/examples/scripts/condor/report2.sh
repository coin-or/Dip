cp report.awk ${HOME}/running/decomp${1}
cp report.root.awk ${HOME}/running/decomp${1}
cp report.direct.awk ${HOME}/running/decomp${1}

for i in atm milpblock mmkp gap
do
  for t in p p2 c c2
  do
    echo ${HOME}/running/decomp${1}/${i}/${i}${t}        
    cd ${HOME}/running/decomp${1}/${i}/${i}${t}    
    ../../report.awk report.txt ${i}${t} 600
    ../../report.root.awk report.root.txt ${i}${t} 600
  done
  echo ${HOME}/running/decomp${1}/${i}/${i}d
  cd ${HOME}/running/decomp${1}/${i}/${i}d
  ../../report.direct.awk report.direct.txt ${i}d 600
done


# for i in #vrp
# do
#   for t in p p2 c c2
#   do
#     echo ${HOME}/running/decomp${1}/${i}/${i}${t}        
#     cp report.awk ${HOME}/running/decomp${1}/${i}/${i}${t}/
#     cp report.root.awk ${HOME}/running/decomp${1}/${i}/${i}${t}/
#     cd ${HOME}/running/decomp${1}/${i}/${i}${t}    
#     ./report.awk report.txt ${i}${t} 600
#     ./report.root.awk report.root.txt ${i}${t} 600
#   done
# done
