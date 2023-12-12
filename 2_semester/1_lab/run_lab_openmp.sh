#! /usr/bin/bash

THREAD_NUMBER=$1
RUN_OBJ_NAME=$2

for (( n = 1; n < ${THREAD_NUMBER} + 1; n++))
do
    echo "n: $n"
    for (( i = 0; i < 15; i++))
    do
        ./build/${RUN_OBJ_NAME} ${n}
    done
done > ./build/statistic_${RUN_OBJ_NAME}.txt
