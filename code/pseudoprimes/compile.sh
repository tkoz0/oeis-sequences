#!/bin/bash
g++ -O3 -Wall -Werror pp.c -o pp
g++ -O3 -Wall -Werror sorted_diff.c -o sorted_diff
mkdir fpp_tmp
parallel g++ -O3 -Wall -Werror fpp_tmp.cpp -o ./fpp_tmp/fpp_{} \
    -ftemplate-depth=1100 -DBASE={} ::: $(seq 2 1023)
