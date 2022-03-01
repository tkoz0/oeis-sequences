# combine stats from subtree files (root_#.csv) into 1
# root.csv has correct stats up to the split length
# root_#.csv files have correct stats above split length, must be combined
# python3 stats_combine.py <input_dir> <prime_type> <base> <split_len>

# output is written to stdout
# up to split length, stats come from root.csv
# beyond that, combine stats from all root_#.csv files
# use job_roots_all.txt to find the root_#.csv files to add
# (job_roots_all.txt includes duplicates for left-or-right primes)

import os
import re
import sys
from typing import Callable, Dict, List

job_re = re.compile(r'root_(\d+).csv')
input_dir = os.path.normpath(sys.argv[1])
prime_type = sys.argv[2]
base = int(sys.argv[3])
split_len = int(sys.argv[4])

assert prime_type in ['r','l','lor','lar']
assert 2 <= base <= 255
assert split_len >= 2

# determine max children count
max_children = base
if prime_type == 'lor':
    max_children *= 2
elif prime_type == 'lar':
    max_children *= base

# for each root_#.csv file, take digit lengths > split_len
job_roots = open(input_dir+'/job_roots_all.txt','r').read().splitlines()

class StatsFile:
    def __init__(self,file:str):
        lines = open(file,'r').read().splitlines()
        table_lines = [line for line in lines if not line.startswith('#')]
        comment_lines = [line for line in lines if line.startswith('#')]
        # set properties from comments
        props = dict()
        for comment in comment_lines:
            comment = comment.split()
            if len(comment) == 4: # lines with "# prop = value"
                props[comment[1]] = comment[3]
        self.prime_type = props['prime_type']
        self.base = int(props['base'])
        self.root = int(props['root'])
        self.max_length = int(props['max_length'])
        self.hash = int(props['hash'])
        # parse table stuff
        self.pmin = dict() # len -> (all,[0,1,.. (num children)])
        self.pmax = dict()
        self.counts = dict()
        assert len(table_lines) % 3 == 1
        assert table_lines[0].startswith('digits,all')
        table_lines = table_lines[1:]
        for i in range(0,len(table_lines),3):
            row = table_lines[i].split(',')
            num_len = int(row[0])
            self.counts[num_len] = (int(row[1]),[int(k) for k in row[2:]])
            row = table_lines[i+1].split(',')
            self.pmin[num_len] = (int(row[1]),[int(k) for k in row[2:]])
            row = table_lines[i+2].split(',')
            self.pmax[num_len] = (int(row[1]),[int(k) for k in row[2:]])
            assert len(self.counts[num_len][1]) == max_children
            assert len(self.pmin[num_len][1]) == max_children
            assert len(self.pmax[num_len][1]) == max_children

# accumulator functions for updating the values
umin = lambda acc,x: min(acc,x) if acc != 0 and x != 0 else max(acc,x)
umax = lambda acc,x: max(acc,x)

# stats objects
pmin = dict() # length -> (children -> min prime)
pmax = dict() # length -> (children -> max prime)
counts = dict() # length -> (children -> count)

# update stats objects with new information
def apply_updates(sf: StatsFile, len_filter: Callable[[int],bool]):
    for length in filter(len_filter,sf.counts.keys()):
        if length not in counts:
            counts[length] = (0,[0]*max_children)
            pmin[length] = (0,[0]*max_children)
            pmax[length] = (0,[0]*max_children)
        a,l = counts[length]
        a2,l2 = sf.counts[length]
        counts[length] = (a+a2,[l[i]+l2[i] for i in range(len(l))])
        a,l = pmin[length]
        a2,l2 = sf.pmin[length]
        pmin[length] = (umin(a,a2),[umin(l[i],l2[i]) for i in range(len(l))])
        a,l = pmax[length]
        a2,l2 = sf.pmax[length]
        pmax[length] = (umax(a,a2),[umax(l[i],l2[i]) for i in range(len(l))])

sf = StatsFile(input_dir+'/root.csv')
apply_updates(sf, lambda x: x <= split_len)

root2hash = dict()

for job_root in job_roots:
    sf = StatsFile(input_dir+'/root_'+job_root+'.csv')
    apply_updates(sf, lambda x: x > split_len)
    root2hash[int(job_root)] = sf.hash

# output result
print(f'# prime_type = {prime_type}')
print(f'# base = {base}')
print(f'# root = 0')
print(f'# max_length = 4294967295')
if prime_type == 'lor':
    print(f'# NOTE: counts are not applicable')
if prime_type == 'lar':
    len_order = sorted(filter(lambda x: x % 2 == 1, counts.keys())) \
            + sorted(filter(lambda x: x % 2 == 0, counts.keys()))
else:
    len_order = sorted(counts.keys())
print('digits,all'+''.join(f',{k}' for k in range(max_children)))
for l in len_order:
    print(f'{l},{counts[l][0]}'+''.join(f',{s}' for s in counts[l][1]))
    print(f',{pmin[l][0]}'+''.join(f',{s}' for s in pmin[l][1]))
    print(f',{pmax[l][0]}'+''.join(f',{s}' for s in pmax[l][1]))

# compute hash
import truncprimes
depth = split_len + (2 if prime_type == 'lar' else 1)
funcs = \
{
    'r': truncprimes.tp_r,
    'l': truncprimes.tp_l,
    'lor': truncprimes.tp_lor,
    'lar': truncprimes.tp_lar
}
rootbytes = b'\xff'*((len(prime_type)+1)//2)
tptree,tpbytes,tphash = funcs[prime_type](base,0,0,rootbytes,depth)

def cust_hash(root: truncprimes.TPTree, hashes: Dict[int,int]) -> int:
    (_,number,_),children = root
    if number in hashes:
        return hashes[number]
    hash = truncprimes.hash_init(number)
    for d in children:
        hash = truncprimes.hash_update(hash,d,cust_hash(children[d],hashes))
    return hash

print(f'# hash = {cust_hash(tptree,root2hash)}')
