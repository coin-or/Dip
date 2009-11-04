for i in atm milpblock mmkp #gap
do
  echo "condor_submit condor${1}.${i}p.submit"
  condor_submit condor${1}.${i}p.submit
  echo "condor_submit condor${1}.${i}p2.submit"
  condor_submit condor${1}.${i}p2.submit
  echo "condor_submit condor${1}.${i}c.submit"
  condor_submit condor${1}.${i}c.submit
  echo "condor_submit condor${1}.${i}c2.submit"
  condor_submit condor${1}.${i}c2.submit
  echo "condor_submit condor${1}.${i}d.submit"
  condor_submit condor${1}.${i}d.submit
done
# for i in vrp
# do
#   echo "condor_submit condor${1}.${i}p.submit"
#   condor_submit condor${1}.${i}p.submit
#   echo "condor_submit condor${1}.${i}p2.submit"
#   condor_submit condor${1}.${i}p2.submit
#   echo "condor_submit condor${1}.${i}c.submit"
#   condor_submit condor${1}.${i}c.submit
#   echo "condor_submit condor${1}.${i}c2.submit"
#   condor_submit condor${1}.${i}c2.submit
# done
