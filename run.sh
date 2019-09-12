#!/bin/bash

if [ ! -f $1.inp ]; then
    ./testgen.py $1
fi

make -j -k -s && \
./main $1.inp $1.out && \
./checkSolution.py $1.inp $1.out
