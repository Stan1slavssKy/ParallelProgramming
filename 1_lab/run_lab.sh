#! /usr/bin/bash

for (( n = 1; n < 9; n++))
do
  echo "n: $n"
  echo "X: 50"
  echo "T: 1"
  for (( i = 0; i < 10; i++))
  do
    mpirun -np $n ./cross
  done
done > ./build/statistic.txt
