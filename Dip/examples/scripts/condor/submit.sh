for i in atm gap milpblock mmkp
do
  condor_submit condor${1}.${i}p.submit
  condor_submit condor${1}.${i}p2.submit
  condor_submit condor${1}.${i}c.submit
  condor_submit condor${1}.${i}c2.submit
done
