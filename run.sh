#!/bin/bash

for p in 2 4 8
do
    make clean
    make np=$p -kj
    for i in 1000 2000 5000 7000 10000
    do
        if [ ! -f $i.inp ]; then
            python3 testgen.py $i
        fi

        ./main $i.inp $i.out
        python3 checkSolution.py $i.inp $i.out
    done
done
