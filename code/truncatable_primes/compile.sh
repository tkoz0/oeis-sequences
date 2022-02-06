#!/bin/bash
gcc -Wall -Werror -g -O3 truncprimes.c -lgmp -o tp_tree -DWRITE_TREE
gcc -Wall -Werror -g -O3 truncprimes.c -lgmp -o tp_stats -DWRITE_STATS
