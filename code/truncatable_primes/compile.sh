#!/bin/bash
echo --- Compiling tp_tree \(-DWRITE_TREE\)
gcc -Wall -Werror -g -O3 truncprimes.c tp_util.c -lgmp -o tp_tree -DWRITE_TREE
echo --- Compiling tp_stats \(-DWRITE_STATS\)
gcc -Wall -Werror -g -O3 truncprimes.c tp_util.c -lgmp -o tp_stats -DWRITE_STATS
echo --- Compiling tree_convert
gcc -Wall -Werror -g -O3 tree_convert.c -lgmp -o tree_convert
