# num_len_filter.py <base> <length1> [length2...]
# read integers (as base 10) from stdin (1 per line) and using given base
# output (in base 10) the numbers with length1 (or length2...) digits
import sys
base = int(sys.argv[1])
assert base >= 2
def has_n_digits(num,n):
    return base**(n-1) <= num < base**n
nvals = [int(a) for a in sys.argv[2:]]
for nval in nvals:
    assert nval > 0
for line in sys.stdin:
    num = int(line)
    if any(has_n_digits(num,nval) for nval in nvals):
        print(line,end='')
