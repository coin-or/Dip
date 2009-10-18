# read edges from file and call addedge.awk to "solidify" in .dot file
let count=0
BN=`basename $1 .dot`
cp $BN.dot $BN.0.dot
while read i j; do
    echo $i
    echo $j
    let count1=$count+1
    awk -f ../../../scripts/addedge.awk $i $j $BN.$count.dot > $BN.$count1.dot
    let count=$count1
done < $2
cp $BN.$count1.dot tmp.dot
