#!/bin/bash

LIST_FILE=$1

while read instance numBlocks
do
  echo $instance " " $numBlocks
  ./createOne.sh $instance $numBlocks
done < $LIST_FILE
