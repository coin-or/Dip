DIR=/users/magala/COIN/coin-Dip/build-O/Dip/examples/MILPBlock

echo "Running 200-dip1"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm200 > rand200.dip1.txt
echo "Running 200-dip2"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm200 --PRICE_AND_CUT:CompressColumns 0 > rand200.dip2.txt
echo "Running 200-dip3"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm200 --PRICE_AND_CUT:CompressColumns 0 --PRICE_AND_CUT:SolveMasterAsIpFreqPass 2 > rand200.dip3.txt
echo "Running 200-dip4"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm200 --PRICE_AND_CUT:CompressColumns 0 --PRICE_AND_CUT:SolveMasterAsIpFreqPass 2 --PRICE_AND_CUT:BreakOutPartial 1 > rand200.dip4.txt
echo "Running 200-dip5"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm200 --PRICE_AND_CUT:CompressColumns 0 --PRICE_AND_CUT:BreakOutPartial 1 > rand200.dip5.txt
