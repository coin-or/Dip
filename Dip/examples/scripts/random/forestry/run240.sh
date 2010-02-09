DIR=/users/magala/COIN/coin-Dip/build-O/Dip/examples/MILPBlock

echo "Running 240-dip1"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm240 > rand240.dip1.txt
echo "Running 240-dip2"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm240 --PRICE_AND_CUT:CompressColumns 0 > rand240.dip2.txt
echo "Running 240-dip3"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm240 --PRICE_AND_CUT:CompressColumns 0 --PRICE_AND_CUT:SolveMasterAsIpFreqPass 2 > rand240.dip3.txt
echo "Running 240-dip4"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm240 --PRICE_AND_CUT:CompressColumns 0 --PRICE_AND_CUT:SolveMasterAsIpFreqPass 2 --PRICE_AND_CUT:BreakOutPartial 1 > rand240.dip4.txt
echo "Running 240-dip5"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm240 --PRICE_AND_CUT:CompressColumns 0 --PRICE_AND_CUT:BreakOutPartial 1 > rand240.dip5.txt
