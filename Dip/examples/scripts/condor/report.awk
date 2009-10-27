#parser for decomp output
REPORT=$1
NAME=$2
TIME=$3

awk -F'.' '{print $1 " ->" $0'} $1 > tmp
awk '
{
if($7 > 1.0e15){
   printf "Instance= %15s LB= %12.3f UB= %12s Gap= %10s Nodes= %10d Cpu= %10.2f Real= %10.2f\n",$1,$5,999999,999999,$9,$15,$21;
}
else{
   if($7 > 0){
      printf "Instance= %15s LB= %12.3f UB= %12.3f Gap= %10.4f Nodes= %10d Cpu= %10.2f Real= %10.2f\n",$1,$5,$7,($7-$5)/$7,$9,$15,$21;
   }
   else{
      printf "Instance= %15s LB= %12.3f UB= %12.3f Gap= %10.4f Nodes= %10d Cpu= %10.2f Real= %10.2f\n",$1,$5,$7,($7-$5)/(-$7),$9,$15,$21;
   }
 }
}' tmp > tmp2
 
#awk -F'_' '{print $2 " " $3 " " $4}' tmp2 > tmp3

# $2     , $12,  $4, $6, $8,  $10
#instance, time, LB, UB, gap, nodes
awk -v time="$TIME" '
{
timeLim = time-(time/20);
if($6 == 999999){ #no ub
  if($12 >= timeLim){
    printf "%15s & T        & %10.2f & $\\infty$ & $\\infty$  & %8d\n", 
       $2,$4,$10;
  }
  else{
    printf "%15s & %8.2f & & %10.2f & $\\infty$ & $\\infty$  & %8d\n",   
       $2,$12,$4,$10;
  }
}
else{
  if($12 >= timeLim){ #exceed time
    printf "%15s & T        & %10.2f & %10.2f & %8.2f\\% & %8d\n", 
       $2,$4,$6,($8 *100),$10; 
  }  
  else{
    if($8 <= 0.00){
      printf "%15s & %8.2f & %10.2f & %10.2f & OPT & %8d\n", 
         $2,$12,$4,$6,$10;
    }
    else{
      printf "%15s & %8.2f & %10.2f & %10.2f & %8.2f\\% & %8d\n", 
         $2,$12,$4,$6,($8 *100),$10;
    }
 }
}
}' tmp2 > tmp3
   sort tmp3 > ${NAME}.${TIME}
#   sort -n +0 -1 tmp3 > ${NAME}.${TIME}
#   sort -n +0 -1 +2 -3 tmp3 > ${NAME}.${TIME}

#rm tmp
#rm tmp2
#rm tmp3
#rm tmp4


