#!/bin/bash

nl_loop=6
init_loop=5
NL[0]=1
NL[1]=1.5
NL[2]=2
NL[3]=2.5
NL[4]=3
NL[5]=3.5

s=0
while test $s -lt $nl_loop
do
    #select nl
    echo "nl loop = $s"
    rm nl_value.dat
    echo ${NL[$s]} >>nl_value.dat
    #average over different initialization
    h=0
    while test $h -lt $init_loop
    do
        echo "init loop = $h"
        ./motion_diff_nl >tmp.dat
        mv tmp.dat tmp"$h"_"$s".dat
        mv motion.dat motion"$h"_nl"$s".dat
        h=`expr $h + 1`
    done
    s=`expr $s + 1`
done