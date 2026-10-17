../../../scripts/addedges.sh sg$1.dot d$1.$2.list
sed 's/dotted/invis/' tmp.dot > tmp2.dot
neato -Gsep=.1 -Gsplines=true -Goverlap=false/scale -Tps tmp2.dot -o sgD$1.$2.ps