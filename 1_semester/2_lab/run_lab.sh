#! /usr/bin/bash

for (( n = 1; n < 17; n++))
do
    echo "n: $n"
    for (( i = 0; i < 50; i++))
    do
        ./build/integration $n 1e-11
    done
done > ./build/statistic.txt
