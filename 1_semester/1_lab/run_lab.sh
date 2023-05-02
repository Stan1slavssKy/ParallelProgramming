#! /usr/bin/bash

for (( n = 1; n < 9; n++))
do
  echo "n: $n"
  echo "X: 50"
  echo "T: 100"
  for (( i = 0; i < 10; i++))
  do
    mpirun -np $n ./build/cross
  done
done > ./build/statistic.txt
