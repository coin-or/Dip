#parser for decomp output
REPORT=$1
NAME=$2
TIME=$3

#break up name of file to give instance name in first column
awk -F'.' '{print $1 " ->" $0'} $1 > tmp

#change INF to 999999
sed 's/INF/999999/' tmp > tmp2

# $1     , $15,  $5, $7,      $9
#instance, time, LB, UB, gap, nodes

#Instance=1
#Statu=$3
#LB=$5
#UB=$7
#Nodes=$9
#CPU=$15
#Real=$21

awk -v time="$TIME" '
{
   timeLim = time-(time/20);
   if($5 == -999999){ # no lb
      if($15 >= timeLim){ #exceed time
         printf "%15s & T        & $\\infty$ & $\\infty$ & $\\infty$  & %8d\n", 
            $1,$9;
      }
      else{
         #no lb, did not exceed time limit <-- error?
         printf "%15s & %8.2f & & $\\infty$ & $\\infty$ & $\\infty$  & %8d\n",   
            $1,$15,$9;
      }
   }   
   else if($7 == 999999){ #no ub
      if($15 >= timeLim){ #exceed time
         printf "%15s & T        & %10.2f & $\\infty$ & $\\infty$  & %8d\n", 
            $1,$5,$9;
      }
      else{
         #no ub, did not exceed time limit <-- error?
         printf "%15s & %8.2f & & %10.2f & $\\infty$ & $\\infty$  & %8d\n",   
            $1,$15,$5,$9;
      }
   }
   else{
      #if ub is negative, will need to flip in gap calc for absolute value
      if($7 >= 0){ mult=1 } else {mult=-1};
      gap = 100*($7-$5)/(mult*$7);
      if($15 >= timeLim){ #exceed time
         printf "%15s & T        & %10.2f & %10.2f & %8.2f\\% & %8d\n", 
            $1,$5,$7,gap,$9; 
      }  
      else{
         if(gap <= 0.0000001){
            printf "%15s & %8.2f & %10.2f & %10.2f & OPT & %8d\n", 
               $1,$15,$5,$7,$9;
         }
         else{
            printf "%15s & %8.2f & %10.2f & %10.2f & %8.2f\\% & %8d\n", 
               $1,$13,$5,$7,gap,$9;
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


