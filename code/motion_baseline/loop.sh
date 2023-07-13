#!/bin/bash

init_loop=5
#average over different initialization
h=0
while test $h -lt $init_loop
do
    echo "init loop = $h"
    ./motion_baseline >tmp.dat
    mv tmp.dat tmp"$h".dat
    mv motion.dat motion"$h".dat
    h=`expr $h + 1`
done
s=`expr $s + 1`
