#parser for direct milp output
TIME=$2

awk -F'.' '{print $1 " ->" $0'} $1 > tmp
#$1 Instance
#$11 LB
#$13 UB
#$9  Nodes
#$5  RealTime
#$7  CpuTime

#retail20 ->retail20.direct.txt:DIRECT SOLVE RealTime= 3629.44399 CpuTime= 3600.01000 Nodes=    4646270 LB=        80.02 UB=        80.18
awk '
{
if($13>1.0e15){
printf "Instance= %15s LB= %12.3f UB= %12s Gap= %10s Nodes= %10d Cpu= %10.2f Real= %10.2f\n",$1,$11,999999,999999,$9,$7,$5;
}
else{
printf "Instance= %15s LB= %12.3f UB= %12.3f Gap= %10.4f Nodes= %10d Cpu= %10.2f Real= %10.2f\n",$1,$11,$13,($13-$11)/$13,$9,$7,$5;
}
}' tmp > tmp2


#$1 Instance
#$4 LB
#$6 UB
#$8 Gap
#$10 Nodes
#$12 Cpu
#$14 Real
#Instance=        retail20 LB=       80.020 UB=       80.180 Gap=     0.0020 Nodes=    4646270 Cpu=    3600.01 Real=    3629.44

awk -v time="$TIME" '
{
timeLim = time-(time/20);
if($6 == 999999){
  if($12 >= timeLim){
    printf "%10s  & T       & $\\infty$  & %8d\n", 
    $2,$10;
  }
  else{
    printf "%10s  & %8.2f & $\\infty$  & %8d\n", 
    $2,$12,$10;
  }
}
else{
  if($12 >= timeLim){
    printf "%10s  & T        & %8.2f\\% & %8d\n", 
    $2,($8 *100),$10;
  }
  else{
    if($8 <= 0.01){
      printf "%10s  & %8.2f & OPT & %8d\n", 
      $2,$12,$10;
    }
    else{
      printf "%10s  & %8.2f & %8.2f\\% & %8d\n", 
      $2,$12,($8 *100),$10;
    }
  }
}
}' tmp2 > tmp3
#sort -n +0 -1 +2 -3 tmp3 > cpx.${TIME}
sort tmp3 > cpx.${TIME}

#rm tmp
#rm tmp2
#rm tmp3
#rm tmp4
