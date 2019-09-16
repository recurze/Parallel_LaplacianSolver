#!/bin/bash

for i in 100 200 400 800
do
    ./testgen.py $i

    make -j -k -s && \
    ./main $i.inp $i.out;
    ./checkSolution.py $i.inp $i.out
done
