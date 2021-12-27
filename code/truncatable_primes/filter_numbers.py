# filter numbers by their length in a base
import sys
base = int(sys.argv[1])
assert base >= 2
def has_n_digits(num,n):
    return base**(n-1) <= num < base**n
nvals = [int(a) for a in sys.argv[2:]]
for line in sys.stdin:
    num = int(line)
    assert 0 <= num < 2**64 # root values must fit in uint64_t
    if any(has_n_digits(num,nval) for nval in nvals):
        print(line,end='')
