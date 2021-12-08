import sys
from typing import Generator

# PRP test with base a, n > 1, usually 1 < a < n-1
# should ensure gcd(a,n)=1 but not required
def prp(n: int, a: int = 2) -> bool:
    # assert n > 1 and a >= 1
    return pow(a,n-1,n) == 1

# SPRP test with base a, n > 2 odd, usually 1 < a < n-1
def sprp(n: int, a: int = 2) -> bool:
    # assert n > 2 and n % 2 == 1 and a >= 1
    s = 0
    d = n-1
    while d % 2 == 0: # n-1 == d * 2**s
        s += 1
        d //= 2
    res = pow(a,d,n)
    if res == 1 or res == n-1:
        return True
    return any((res := (res*res) % n) == n-1 for _ in range(s-1))

# probable prime test to use when computing numbers
# subtrees can be pruned by a proper primality test afterward
# all primes will be output, but some pseudoprimes may be present
def prob_prime_test(n: int) -> bool:
    return n == 2 or (n % 2 == 1 and n > 2 and prp(n))
    #return n == 2 or (n % 2 == 1 and n > 2 and sprp(n))

# right truncatable primes (A024770)
# b >= 2 is the base
# v is the recursion value, initially 0
# digit is appended right by incrementing from b*v
def right_trunc(b: int = 10, v: int = 0) -> Generator[int,None,None]:
    v *= b
    for _ in range(1,b):
        v += 1
        if prob_prime_test(v):
            yield v
            yield from right_trunc(b,v)

# left truncatable primes (A024785)
# b >= 2 is the base
# d is the place value to add for appending a left digit, initially 1
# v is the recursion value, initially 0
# digit is appended left by adding the place value d
def left_trunc(b: int = 10, d: int = 1, v: int = 0) -> Generator[int,None,None]:
    for _ in range(1,b):
        v += d
        if prob_prime_test(v):
            yield v
            yield from left_trunc(b,d*b,v)

# TODO have a helper function that eliminates single digit primes being formed
#   by both left and right digit appends when v=0
# TODO handle multiple recursion paths to same numbers (ex: 7 -> 73, 3 -> 73)
#   this may be better done with BFS to avoid duplicating numbers
# left-or-right truncatable primes (A137812)
# b >= 2 is the base
# d is the place value to add for appending a left digit, initially 1
# v is the recursion value, initially 0
def lor_trunc(b: int = 10, d: int = 1, v: int = 0) -> Generator[int,None,None]:
    d2 = d*b
    for _ in range(1,b): # append left
        v += d
        if prob_prime_test(v):
            yield v
            yield from lor_trunc(b,d2,v)
    v -= (b-1)*d
    v *= b
    for _ in range(1,b): # append right
        v += 1
        if prob_prime_test(v):
            yield v
            yield from lor_trunc(b,d2,v)

# helper function to do the recursion for left-and-right truncatable primes
# b >= 2 is the base
# d is the place value of the leftmost digit of v
# v is the recursion value, initially a 1 or 2 digit prime
def lar_trunc_h(b: int, d: int, v: int) -> Generator[int,None,None]:
    v *= b
    d *= b*b
    for _ in range(1,b):
        v += d
        for _ in range(1,b):
            v += 1
            if prob_prime_test(v):
                yield v
                yield from lar_trunc_h(b,d,v)
        v -= b-1

# left-and-right truncatable primes (A077390)
# initialize the recursion with 1 or 2 digit primes
def lar_trunc(b: int = 10) -> Generator[int,None,None]:
    for v in range(1,b):
        if prob_prime_test(v):
            yield v
            yield from lar_trunc_h(b,1,v)
    for v in range(b,b*b):
        if prob_prime_test(v):
            yield v
            yield from lar_trunc_h(b,b,v)

ptype = sys.argv[1]
base = int(sys.argv[2])
assert ptype in ['l','r','lor','lar']
assert base >= 2

if ptype == 'l':
    primes = left_trunc(base)
elif ptype == 'r':
    primes = right_trunc(base)
elif ptype == 'lor':
    primes = lor_trunc(base)
else:
    primes = lar_trunc(base)

for p in primes:
    print(p)
