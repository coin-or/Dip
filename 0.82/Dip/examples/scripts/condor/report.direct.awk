#parser for decomp output
# $2                    $3    $4             $5    $6         $7      $8       $9     $10       $11      $12
#I10.condor.out:DIRECT SOLVE Real=624.92614  Cpu= 600.28874  Nodes= 144236   objLB= -61485.119 objUB= -61426.000
#10_100_1.condor.out:DIRECT SOLVE Real=603.64554  Cpu= 600.21175  Nodes= 5081     objLB= 1137632.843 objUB=  INF


REPORT=$1
NAME=$2
TIME=$3

#break up name of file to give instance name in first column
awk -F'.' '{print $1 " ->" $0'} $1 > tmp.direct

#change INF to 999999
sed 's/INF/999999/' tmp.direct > tmp2.direct

# $1     , $6,   $10, $12,     , $8
#instance, time, LB,  UB,  gap, nodes

awk -v time="$TIME" '
{
   timeLim = time-(time/20);
   if($12 == 999999){ #no ub
      if($6 >= timeLim){ #exceed time
         printf "%15s & T        & %10.2f & $\\infty$ & $\\infty$ & %8d\n", 
	   $1,$10,$8;
      }
      else{
         #no ub, did not exceed time limit <-- error?
         printf "%15s & %8.2f & & %10.2f & $\\infty$ & $\\infty$ & %8d\n",   
	   $1,$6,$10,$8;
      }
   }
   else{
      #if ub is negative, will need to flip in gap calc for absolute value
      if($12 >= 0){ mult=1 } else {mult=-1};
      gap = 100*($12-$10)/(mult*$12);
      if($6 >= timeLim){ #exceed time
         printf "%15s & T        & %10.2f & %10.2f & %8.2f\\% & %8d\n", 
	   $1,$10,$12,gap,$8; 
      }  
      else{
         if(gap <= 0.0000001){
            printf "%15s & %8.2f & %10.2f & %10.2f & OPT & %8d\n", 
	      $1,$6,$10,$12,$8;
         }
         else{
            printf "%15s & %8.2f & %10.2f & %10.2f & %8.2f\\% & %8d\n", 
	      $1,$6,$10,$12,gap,$8;
         }
      }
   }
}' tmp2.direct > tmp3.direct
sort tmp3.direct > ${NAME}.${TIME}.direct


