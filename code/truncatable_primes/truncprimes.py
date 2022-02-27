import argparse
import gmpy2
from typing import Dict, List, Tuple

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
    #return n == 2 or (n % 2 == 1 and n > 2 and prp(n))
    #return n == 2 or (n % 2 == 1 and n > 2 and sprp(n))
    return gmpy2.is_bpsw_prp(n)

# counts digits of n in base b, returns 0 if n == 0
def count_digits(n: int, b: int) -> int:
    assert n >= 0
    assert b > 1
    d = 0
    while n:
        d += 1
        n //= b
    return d

# based on lower 64 bits of the prime
def hash_init(num: int) -> int:
    return (num%2**64)//2

# rotate 32 bits
def swap_halves(num: int) -> int:
    return ((num >> 32) | (num << 32)) & (2**64-1)

# mix in the new values to scramble the hash value
def hash_update(h: int, d: int, c: int) -> int:
    return h ^ swap_halves((8191*(127*h-d)+c)%2**64)

'''
Recursive functions return tuple (tree,bytes,hash)
'''

# node is (length,number) and list of children
TPNode = Tuple[int,int,int] # len,num,hash(subtree)
TPTree = Tuple[TPNode,Dict[int,'TPTree']]
TPRet = Tuple[TPTree,bytes,int]

def tp_r(base: int, len: int, num: int, root: bytes, maxlen: int) -> TPRet:
    children = dict()
    output = root
    hash = hash_init(num)
    for d in range(1,base) if len+1 <= maxlen else []:
        num2 = num*base+d
        if prob_prime_test(num2):
            ctree,cbytes,chash = tp_r(base,len+1,num2,bytes([d]),maxlen)
            children[d] = ctree
            output += cbytes
            hash = hash_update(hash,d,chash)
    output += b'\xff' # end
    return (((len,num,hash),children),output,hash)

def tp_l(base: int, len: int, num: int, root: bytes, maxlen: int) -> TPRet:
    children = dict()
    output = root
    hash = hash_init(num)
    for d in range(1,base) if len+1 <= maxlen else []:
        num2 = d*(base**len)+num
        if prob_prime_test(num2):
            ctree,cbytes,chash = tp_l(base,len+1,num2,bytes([d]),maxlen)
            children[d] = ctree
            output += cbytes
            hash = hash_update(hash,d,chash)
    output += b'\xff' # end
    return (((len,num,hash),children),output,hash)

def tp_lor(base: int, len: int, num: int, root: bytes, maxlen: int) -> TPRet:
    children = dict()
    output = root
    hash = hash_init(num)
    for d in range(1,base) if len+1 <= maxlen else []: # left
        num2 = d*(base**len)+num
        if prob_prime_test(num2):
            ctree,cbytes,chash = tp_lor(base,len+1,num2,bytes([0,d]),maxlen)
            children[d] = ctree
            output += cbytes
            hash = hash_update(hash,d,chash)
    for d in range(1,base) if num and len+1 <= maxlen else []: # right
        num2 = num*base+d
        if prob_prime_test(num2):
            ctree,cbytes,chash = tp_lor(base,len+1,num2,bytes([1,d]),maxlen)
            children[base+d] = ctree
            output += cbytes
            hash = hash_update(hash,base+d,chash)
    output += b'\xff' # end
    return (((len,num,hash),children),output,hash)

def tp_lar(base: int, len: int, num: int, root: bytes, maxlen: int) -> TPRet:
    children = dict()
    output = root
    hash = hash_init(num)
    for dl in range(1 if num else 0, base): # allow 0 for 0 root
        for dr in range(1 if num or dl == 0 else 0, base):
            num2 = dl*(base**(len+1))+base*num+dr
            len2 = len + (2 if num or dl else 1)
            if len2 <= maxlen and prob_prime_test(num2):
                ctree,cbytes,chash = tp_lar(base,len2,num2,
                                            bytes([dl,dr]),maxlen)
                children[dl*base+dr] = ctree
                output += cbytes
                hash = hash_update(hash,dl*base+dr,chash)
    output += b'\xff' # end
    return (((len,num,hash),children),output,hash)

class TPStats:
    def __init__(self):
        self.hashes = []
        self.pmin = dict() # length -> (children -> min prime)
        self.pmax = dict() # length -> (children -> max prime)
        self.counts = dict() # length -> (children -> count)
    def insert_hash(self, hash: int):
        self.hashes.append(hash)
    def update(self, num: int, len: int, children: int):
        if len not in self.counts:
            self.counts[len] = dict()
            self.pmin[len] = dict()
            self.pmax[len] = dict()
        if children not in self.counts[len]:
            self.counts[len][children] = 0
            self.pmin[len][children] = num
            self.pmax[len][children] = num
        self.counts[len][children] += 1
        self.pmin[len][children] = min(self.pmin[len][children],num)
        self.pmax[len][children] = max(self.pmax[len][children],num)
    def get_hash_modular_distribution(self, mod: int) -> List[int]:
        assert mod > 0
        ret = [0 for _ in range(mod)]
        for hash in self.hashes:
            ret[hash%mod] += 1
        return ret
    def get_hash_bit_distribution(self) -> List[int]:
        ret = [0 for _ in range(64)]
        for hash in self.hashes:
            for i in range(64):
                ret[i] += hash%2
                hash //= 2
        return ret

# adds information to stats object
def tp_tree_stats(root: TPTree, stats: TPStats):
    (length,number,hash),subtrees = root
    if len(subtrees) > 0: # leaf hashes bloat lower bit counts
        stats.insert_hash(hash)
    stats.update(number,length,len(subtrees))
    for subtree in subtrees.values():
        tp_tree_stats(subtree,stats)

# helper function for print_stats
def get_values(max_children: int, obj: Dict[int,int]) -> List[int]:
    return [obj[i] if i in obj else 0 for i in range(max_children)]

# print stats like the C program
def print_stats(stats: TPStats, args: argparse.Namespace):
    print(f'# prime_type = {args.prime_type}')
    print(f'# base = {args.base}')
    print(f'# root = {args.root}')
    print(f'# max_length = {args.max_length}')
    max_children = args.base
    if args.prime_type == 'lor':
        print(f'# NOTE: counts are not applicable')
        max_children *= 2
    if args.prime_type == 'lar':
        max_children *= args.base
        len_order = sorted(filter(lambda x : x % 2 == 1, stats.counts.keys())) \
                + sorted(filter(lambda x : x % 2 == 0, stats.counts.keys()))
    else:
        len_order = sorted(stats.counts.keys())
    print('digits,all'+''.join(f',{k}' for k in range(max_children)))
    for l in len_order:
        all_min = min(stats.pmin[l].values())
        all_max = max(stats.pmax[l].values())
        all_count = sum(stats.counts[l].values())
        if l == 0 or all_count == 0:
            continue # skip unnecessary row
        print(f'{l},{all_count}'+''.join(f',{i}'
            for i in get_values(max_children,stats.counts[l])))
        print(f',{all_min}'+''.join(f',{i}'
            for i in get_values(max_children,stats.pmin[l])))
        print(f',{all_max}'+''.join(f',{i}'
            for i in get_values(max_children,stats.pmax[l])))
    print(f'# hash = {stats.hashes[0]}')

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-b','--base',default='10',
        help='primes base (2 <= base <= 255) (default 10)')
    parser.add_argument('-l','--max_length',default='4294967295',
        help='max length in base digits (default unlimited)')
    parser.add_argument('-p','--prime_type',
        help='type of truncatable primes (r, l, lor, lar)')
    parser.add_argument('-r','--root',default='0',
        help='root of recursion tree (>= 0, 0 for full tree)')
    args = parser.parse_args()
    #print('args =',args)
    args.base = int(args.base)
    args.max_length = int(args.max_length)
    args.root = int(args.root)
    funcs = \
    {
        'r': tp_r,
        'l': tp_l,
        'lor': tp_lor,
        'lar': tp_lar
    }
    tptree,tpbytes,tphash = funcs[args.prime_type](args.base,count_digits(
        args.root,args.base),args.root,
        b'\xff'*((len(args.prime_type)+1)//2),args.max_length)
    #print('bytes length =',len(tpbytes))
    #print('bytes hash (md5) =',hashlib.md5(tpbytes).hexdigest())
    #print('tree hash (custom) =',tphash)
    stats = TPStats()
    tp_tree_stats(tptree,stats)
    print_stats(stats,args)
    #print('hash bit distribution =',stats.get_hash_bit_distribution())
    #while True:
    #    m = int(input())
    #    print(stats.get_hash_modular_distribution(m))


'''
=== OLD CODE BELOW ===

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

def lor_trunc_init(b: int = 10) -> Generator[int,None,None]:
    for v in range(1,b):
        if prob_prime_test(v):
            yield v
            yield from lor_trunc(b,b,v)

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
    primes = lor_trunc_init(base)
else:
    primes = lar_trunc(base)

for p in primes:
    print(p)

'''
