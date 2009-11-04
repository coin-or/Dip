for i in atm milpblock mmkp #gap
do
  for t in p p2 c c2
  do
    echo ${HOME}/running/decomp${1}/${i}/${i}${t}    
    cd ${HOME}/running/decomp${1}/${i}/${i}${t}    
    for f in *.err*
    do
      if [ -s $f ]
        then
          echo "======================== ERROR $f"
          cat $f
      fi
    done
    grep "TotalCPU" *.out > report.txt
    grep "Process Node 1 " *.out > report.root.txt
  done

  echo ${HOME}/running/decomp${1}/${i}/${i}d    
  cd ${HOME}/running/decomp${1}/${i}/${i}d    
  for f in *.err*
  do
    if [ -s $f ]
      then
        echo "======================== ERROR $f"
        cat $f
    fi
  done
  grep "DIRECT SOLVE" *.out > report.direct.txt
done
# for i in vrp
# do
#   for t in p p2 c c2
#   do
#     echo ${HOME}/running/decomp${1}/${i}/${i}${t}    
#     cd ${HOME}/running/decomp${1}/${i}/${i}${t}    
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
