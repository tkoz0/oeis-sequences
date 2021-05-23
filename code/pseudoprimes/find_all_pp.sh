#!/bin/bash
echo BASE=$1
mkdir out_$1
python3 job_generator.py \
    | parallel --slf .. \
        --basefile fpp_tmp/fpp_$1 \
        --basefile sorted_diff \
        --basefile find_pp.sh \
        --cleanup \
        --header : \
        --colsep , \
        --tag \
        --files \
        --tmpdir out_$1 \
        --joblog jobs_$1.log \
        --bar \
        --resume-failed \
        ./find_pp.sh $1 {lo} {hi} \
    >> parallel_$1.stdout
