OPT=$1       #example=-g or -O
VERSION=$2   #example=   or -10 (the latter uses CPA10.2)

for i in atm milpblock mmkp gap
do
  echo ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}d    
  cd ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}d    
  for f in *.err*
  do
    if [ -s $f ]
      then
        echo "======================== ERROR $f"
        cat $f | grep -v "licensed to" | grep -v "my_name_in"
    fi
  done
  grep "DIRECT SOLVE" *.out > report.direct.txt
done


for i in atm milpblock gap
do
  for t in p p2 p3 c c2
  do
    echo ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}${t}    
    cd ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}${t}    
    for f in *.err*
    do
      if [ -s $f ]
        then
          echo "======================== ERROR $f"
          cat $f | grep -v "licensed to" | grep -v "my_name_in"
      fi
    done
    grep "TotalCPU" *.out > report.txt
    grep "Process Node 1 " *.out > report.root.txt
  done
done

for i in mmkp
do
  for t in p p3 p4 c c2
  do
    echo ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}${t}    
    cd ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}${t}    
    for f in *.err*
    do
      if [ -s $f ]
        then
          echo "======================== ERROR $f"
          cat $f | grep -v "licensed to" | grep -v "my_name_in"
      fi
    done
    grep "TotalCPU" *.out > report.txt
    grep "Process Node 1 " *.out > report.root.txt
  done
done




# for i in vrp
# do
#   for t in p p2 c c2
#   do
#     echo ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}${t}    
#     cd ${HOME}/running${VERSION}/decomp${OPT}/${i}/${i}${t}    
#     for f in *.err*
#     do
#       if [ -s $f ]
#         then
#           echo "======================== ERROR $f"
#           cat $f
#       fi
#     done
#     grep "TotalCPU" *.out > report.txt
#     grep "Process Node 1 " *.out > report.root.txt
#   done
# done
