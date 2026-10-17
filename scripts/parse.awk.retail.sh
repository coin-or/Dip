#parser for decomp output
TIME=$2

awk -F'.' '{print $1 " ->" $0'} $1 > tmp
awk '
{printf "Instance= %15s LB= %12.3f UB= %12.3f Gap= %10.4f Nodes= %10d Cpu= %10.2f Real= %10.2f\n",$1,$5,$7,($7-$5)/$7,$9,$15,$21;
}' tmp > tmp2

#awk -F'_' '{print $2 " " $3 " " $4}' tmp2 > tmp3

awk -v time="$TIME" '
{
timeLim = time-(time/20);
if($6 == 999999){
  if($12 >= timeLim){
    printf "%10s & T       & $\\infty$  & %8d\n", 
    $2,$10;
  }
  else{
    printf "%10s & %8.2f & $\\infty$  & %8d\n",   
    $2,$12,$10;
  }
}
else{
  if($12 >= timeLim){
    printf "%10s & T        & %8.2f\\% & %8d\n", 
    $2,($8 *100),$10; 
  }  
  else{
    if($8 <= 0.00){
      printf "%10s & %8.2f & OPT & %8d\n", 
      $2,$12,$10;
    }
    else{
      printf "%10s & %8.2f & %8.2f\\% & %8d\n", 
      $2,$12,($8 *100),$10;
    }
 }
}
}' tmp2 > tmp3
#sort -n +0 -1 +2 -3 tmp3 > pc.${TIME}
sort tmp3 > pc.${TIME}

#rm tmp
#rm tmp2
#rm tmp3
#rm tmp4


