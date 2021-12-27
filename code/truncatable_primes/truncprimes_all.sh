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
./truncprimes -p $PRIME_TYPE -b $BASE -l $SPLIT_LEN > ./$DIR/root.bin
cat ./$DIR/root.bin | ./tree_convert -p $PRIME_TYPE -i $BASE > ./$DIR/root.txt
echo computed $(cat ./$DIR/root.txt | wc -l) for recursion splitting

# filter to get roots of proper length
if [[ $PRIME_TYPE =~ 'lar' ]]
then
    TMP=$(expr $SPLIT_LEN - 1)
    cat ./$DIR/root.txt | python3 filter_numbers.py $BASE $SPLIT_LEN $TMP > ./$DIR/job_roots.txt
else
    cat ./$DIR/root.txt | python3 filter_numbers.py $BASE $SPLIT_LEN > ./$DIR/job_roots.txt
fi

# remove duplicate roots
if [[ $PRIME_TYPE =~ 'lor' ]]
then
    cat ./$DIR/job_roots.txt | sort | uniq > ./$DIR/job_roots_2.txt
    rm ./$DIR/job_roots.txt
    mv ./$DIR/job_roots_2.txt ./$DIR/job_roots.txt
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
        ./truncprimes -p $PRIME_TYPE -b $BASE -r {} '>' ./$DIR/root_{}.bin \
    >> ./$DIR/parallel.stdout
#        ./truncprimes -p $PRIME_TYPE -b $BASE -r {} '|' xz '>' ./$DIR/root_{}.bin.xz \

# combine job subtrees together
echo merging job files into tree.bin
python3 tree_combine.py ./$DIR $PRIME_TYPE $BASE

echo END
exit

# combine all results into a list
cp ./$DIR/root.txt ./$DIR/list.txt
for r in $(cat ./$DIR/job_roots.txt)
do
    #echo root $r has $(cat ./$DIR/root_$r.bin | ./tree_convert -p $PRIME_TYPE -i $BASE -r $r | wc -l) numbers
    cat ./$DIR/root_$r.bin | ./tree_convert -p $PRIME_TYPE -i $BASE -r $r >> ./$DIR/list.txt
done
echo total of $(cat ./$DIR/list.txt | wc -l) numbers

echo END
