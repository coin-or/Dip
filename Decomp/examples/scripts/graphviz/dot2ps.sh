base=`basename $1 .dot`
neato -Gsep=.1 -Gsplines=true -Gmaxiter=500 -Goverlap=false -Tdot $1 -o tmp
neato -Gsep=.1 -Gsplines=true -Gmaxiter=500 -Goverlap=false -Tjpg  tmp -o $base.jpg
convert -size 500x500 $base.jpg $base.jpg
rm tmp
