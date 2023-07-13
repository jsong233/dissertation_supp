#!/bin/bash

or_loop=5
init_loop=5
OR[0]=1
OR[1]=2
OR[2]=4
OR[3]=6
OR[4]=8

s=0
while test $s -lt $or_loop
do
    #select or
    echo "or loop = $s"
    rm or_value.dat
    echo ${OR[$s]} >>or_value.dat
    #average over different initialization
    h=0
    while test $h -lt $init_loop
    do
        echo "init loop = $h"
        ./motion_diff_or >tmp.dat
        mv tmp.dat tmp"$h"_"$s".dat
        mv motion.dat motion"$h"_or"$s".dat
        h=`expr $h + 1`
    done
    s=`expr $s + 1`
done