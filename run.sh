#!/bin/bash

make -j -k -s
for i in 100 200 400 800
do
    if [ ! -f $i ]; then
        ./testgen.py $i
    fi

    ./main $i.inp $i.out 2>$i.err;
    ./checkSolution.py $i.inp $i.out
done
