#!/bin/bash

for i in 10 25 50 100 250 500 1000 2500 5000 10000
do
    ./testgen.py $i

    make -j -k -s && \
    ./main $i.inp $i.out;
    ./checkSolution.py $i.inp $i.out
done
