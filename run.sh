#!/bin/bash

if [ "$1" == "d" ];then
    mpirun -gdb -np $2 ./main
else
    mpirun -np $1 ./main data/mat data/vec
fi
