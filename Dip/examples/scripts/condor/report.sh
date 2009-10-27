for i in atm gap milpblock mmkp
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
    grep "TotalCPU" *.out > report.out
  done
done
