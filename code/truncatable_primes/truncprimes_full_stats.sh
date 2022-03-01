#!/bin/bash

PRIME_TYPE=$1
BASE=$2
SPLIT_LEN=$3

echo PRIME_TYPE=$PRIME_TYPE
echo BASE=$BASE
echo SPLIT_LEN=$SPLIT_LEN

if ! [ $PRIME_TYPE == 'r' -o $PRIME_TYPE == 'l' -o $PRIME_TYPE == 'lor' -o $PRIME_TYPE == 'lar' ]
then
    echo invalid prime type
    exit
fi

if [ $BASE -lt 2 -o $BASE -gt 255 ]
then
    echo invalid base
    exit
fi

if [ $SPLIT_LEN -lt 2 ]
then
    echo split length must be \>= 2
    exit
fi

# prompt before running
read -e -p "continue?(y/n) " YN
echo YN=$YN
if [[ ! $YN =~ [yY] ]]
then
    echo exiting
    exit
fi

DIR=${PRIME_TYPE}_$BASE
echo DIR=$DIR

mkdir ./$DIR
mkdir ./$DIR/tmp

# generate tree up to given limit
./tp_tree -p $PRIME_TYPE -b $BASE -l $SPLIT_LEN > ./$DIR/root.bin
cat ./$DIR/root.bin | ./tree_convert -p $PRIME_TYPE -i $BASE > ./$DIR/root.txt
echo computed $(cat ./$DIR/root.txt | wc -l) for recursion splitting

# generate correct stats for the partial tree
if [[ $PRIME_TYPE =~ 'lar' ]]
then
    ./tp_stats -p $PRIME_TYPE -b $BASE -l $(expr $SPLIT_LEN + 2) > ./$DIR/root.csv
else
    ./tp_stats -p $PRIME_TYPE -b $BASE -l $(expr $SPLIT_LEN + 1) > ./$DIR/root.csv
fi
echo computed stats for the nodes in the partial tree

# filter to get roots of proper length
if [[ $PRIME_TYPE =~ 'lar' ]]
then
    TMP=$(expr $SPLIT_LEN - 1)
    cat ./$DIR/root.txt | python3 num_len_filter.py $BASE $SPLIT_LEN $TMP > ./$DIR/job_roots.txt
else
    cat ./$DIR/root.txt | python3 num_len_filter.py $BASE $SPLIT_LEN > ./$DIR/job_roots.txt
fi

# remove duplicate roots
if [[ $PRIME_TYPE =~ 'lor' ]]
then
    cat ./$DIR/job_roots.txt | sort | uniq > ./$DIR/job_roots_2.txt
    mv ./$DIR/job_roots.txt ./$DIR/job_roots_all.txt
    mv ./$DIR/job_roots_2.txt ./$DIR/job_roots.txt
else
    cp ./$DIR/job_roots.txt ./$DIR/job_roots_all.txt
fi

echo splitting into $(cat ./$DIR/job_roots.txt | wc -l) jobs

# run a separate process for each root
cat ./$DIR/job_roots.txt \
    | parallel \
        --tag \
        --files \
        --group \
        --tmpdir ./$DIR/tmp \
        --joblog ./$DIR/jobs.log \
        --resume-failed \
        --bar \
        ./tp_stats -p $PRIME_TYPE -b $BASE -r {} '>' ./$DIR/root_{}.csv \
    >> ./$DIR/parallel.stdout
#        ./truncprimes -p $PRIME_TYPE -b $BASE -r {} '|' xz '>' ./$DIR/root_{}.bin.xz \

# combine subtrees together when all jobs complete successfully
if [ $? -eq 0 ]
then
    echo merging job files into stats.csv
    python3 stats_combine.py ./$DIR $PRIME_TYPE $BASE $SPLIT_LEN > ./$DIR/stats.csv
else
    echo \>= 1 failed job, not merging
fi
cat ./$DIR/stats.csv | md5sum
echo END
exit
