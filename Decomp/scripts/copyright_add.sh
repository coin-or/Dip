#!/bin/bash

cat Copyright  > tmp
cat $1        >> tmp
mv tmp $1

# for i in ../src/*.cpp 
# do
#   cat Copyright  > tmp
#   cat $i        >> tmp
#   mv tmp $i
# done

# for i in ../src/*.h
# do
#   cat Copyright  > tmp
#   cat $i        >> tmp
#   mv tmp $i
# done
