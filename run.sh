#!/bin/bash

if [ "$1" == "d" ];then
    gdb ./main
else
    ./main data/mat data/vec
fi
