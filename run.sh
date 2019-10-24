#!/bin/bash

for p in 4
do
    make clean
    make np=$p -kj
    for i in 4941 10680 23166 #34761 77360 154908
    #for i in 2000 5000 7000 10000 20000 30000
    do
        cp tests/$i.inp .
        cp ans/$i.ans .
        #if [ ! -f $i.inp ]; then
        #    python3 testgen.py $i
        #fi

        ./main $i.inp $i.out
        if [ $i -lt 15000 ]; then
            python3 checkSolution.py $i.inp $i.out
            mv $i.ans ans/
        fi
        rm ${i}*
    done
done
