DIR=/users/magala/COIN/coin-Dip/build-O/Dip/examples/MILPBlock

#echo "Running 40-dip1"
#$DIR/decomp_milpblock --param $DIR/milpblockR.parm40 > rand40.dip1.txt
#echo "Running 40-dip2"
#$DIR/decomp_milpblock --param $DIR/milpblockR.parm40 --PRICE_AND_CUT:CompressColumns 0 > rand40.dip2.txt
#echo "Running 40-dip3"
#$DIR/decomp_milpblock --param $DIR/milpblockR.parm40 --PRICE_AND_CUT:CompressColumns 0 --PRICE_AND_CUT:SolveMasterAsIpFreqPass 2 > rand40.dip3.txt
#echo "Running 40-dip4"
#$DIR/decomp_milpblock --param $DIR/milpblockR.parm40 --PRICE_AND_CUT:CompressColumns 0 --PRICE_AND_CUT:SolveMasterAsIpFreqPass 2 --PRICE_AND_CUT:BreakOutPartial 1 > rand40.dip4.txt
echo "Running 40-dip5"
$DIR/decomp_milpblock --param $DIR/milpblockR.parm40 --PRICE_AND_CUT:CompressColumns 0 --PRICE_AND_CUT:BreakOutPartial 1 > rand40.dip5.txt
