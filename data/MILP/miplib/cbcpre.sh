CBC=/users/magala/COIN/coin-Decomp/build-g/bin/cbc
DATA=/users/magala/COIN/coin-Decomp/Decomp/data/MILP
while read i
do
  base=`basename $i .mps`
  echo $base
  ${CBC} ${DATA}/miplib3/$i -preprocess save -heuristic off -maxnode -1 -solve > ${DATA}/miplib3_pre/$base.prelog
  mv presolved.mps ${DATA}/miplib3_pre/${base}_pre.mps
done < $1