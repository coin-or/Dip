DIR=/users/magala/COIN/coin-Dip/build-O/Dip/examples/MILPBlock

echo "Running 100-dip1"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm100 > rand100.dip1.txt
echo "Running 100-dip2"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm100 --PRICE_AND_CUT:CompressColumns 0 > rand100.dip2.txt
echo "Running 100-dip3"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm100 --PRICE_AND_CUT:CompressColumns 0 --PRICE_AND_CUT:SolveMasterAsIpFreqPass 2 > rand100.dip3.txt
echo "Running 100-dip4"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm100 --PRICE_AND_CUT:CompressColumns 0 --PRICE_AND_CUT:SolveMasterAsIpFreqPass 2 --PRICE_AND_CUT:BreakOutPartial 1 > rand100.dip4.txt
echo "Running 100-dip5"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm100 --PRICE_AND_CUT:CompressColumns 0 --PRICE_AND_CUT:BreakOutPartial 1 > rand100.dip5.txt
