# Sum and average a column of numbers
# usage: if numbers are in column 1:
# gawk -f average.awk 1 data.in

BEGIN {
  n = ARGV[1];
  ARGV[1] = "";
  max = -1e38;
  min = 1e38;
 
}

{
  if($n > max) max = $n;
  if($n < min) min = $n;
  sum += $n;
  sumq += $n * $n;
}
  
END {
  if(NR > 0){
    average = sum/NR;
    sigq = sumq / NR - average * average;
    if(sigq > 0){
      sigma = sqrt(sigq);
    }
    else sigma = 0;
    
    printf("%10.4f\t%10.4f\t%10.4f\t%10.4f\t%10.4f %d\n", sum, average, sigma, 
	   min, max, NR);
  }
   
}
