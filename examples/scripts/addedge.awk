BEGIN {
  i = ARGV[1]
  j = ARGV[2]
  ARGV[1] = ""
  ARGV[2] = ""
  regex1 = sprintf("\t%d -- %d ",i,j)
  regex2 = sprintf("\t%d -- %d ",j,i)
  #printf regex1
  #printf regex2
}

{
#  printf NF "\n"
  where1 = match($0, regex1)
  where2 = match($0, regex2)
  if(where1 != 0 || where2 != 0){
      string = substr($(NF), 0, length($(NF))-2)
      printf "\t" $1 " " $2 " " $3 " "
      for (i = 4; i <= NF-1; i=i+1) 
           printf $i " "
      printf string ", style=solid, label=\"\"];\n"
  }
  else print $0;
}

END{
}
