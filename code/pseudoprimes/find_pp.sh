#!/bin/bash
echo TYPE=FPP
echo BASE=$1
echo LO_BOUND=$2
echo HI_BOUND=$3
./sorted_diff <(./fpp_tmp/fpp_$1 $2 $3) \
              <(primesieve --threads=1 --print=1 $2 $3 && echo done)
