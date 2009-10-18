for i in *_r*
do
  base=`basename $i .lp`
  base2=`basename $base .p`
  base3=`basename $base2 _r`
  echo $base3
  mv $i $base3.p.lp
done
