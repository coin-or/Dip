for i in milpblock gap atm
do
  echo "condor_submit condor${1}.${i}p.submit"
  condor_submit condor${1}.${i}p.submit
sleep 20
  echo "condor_submit condor${1}.${i}p2.submit"
  condor_submit condor${1}.${i}p2.submit
sleep 20
  echo "condor_submit condor${1}.${i}c.submit"
  condor_submit condor${1}.${i}c.submit
sleep 20
  echo "condor_submit condor${1}.${i}c2.submit"
  condor_submit condor${1}.${i}c2.submit
sleep 20
  echo "condor_submit condor${1}.${i}d.submit"
  condor_submit condor${1}.${i}d.submit
sleep 20
done

for i in mmkp
do
  echo "condor_submit condor${1}.${i}p.submit"
  condor_submit condor${1}.${i}p.submit
sleep 20
  echo "condor_submit condor${1}.${i}c.submit"
  condor_submit condor${1}.${i}c.submit
sleep 20
  echo "condor_submit condor${1}.${i}c2.submit"
  condor_submit condor${1}.${i}c2.submit
sleep 20
  echo "condor_submit condor${1}.${i}d.submit"
  condor_submit condor${1}.${i}d.submit
sleep 20
done
for i in milpblock gap atm mmkp
do
  echo "condor_submit condor${1}.${i}p3.submit"
  condor_submit condor${1}.${i}p3.submit
sleep 20
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
