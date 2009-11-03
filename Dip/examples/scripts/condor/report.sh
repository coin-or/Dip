for i in atm gap milpblock mmkp
do
  for t in p p2 c c2 d
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
done
for i in vrp
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
done
