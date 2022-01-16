# combines the job subtree files into 1 (original files are not deleted)
# python3 tree_combine.py <input_dir> <prime_type> <base>
# expects "root.bin" to have the tree evaluated partially
# (up to the partial depth used to split the recursion into jobs)
# this program writes a new file "tree.bin" with the contents of "root.bin"
# every time an empty subtree is found, it looks for a file "root_#.bin"
# (where # is the number of the current node of the tree)
# the subtree content is written from "root_#.bin"
# essentially, jobs are determined with a partial depth tree
# this program fills in the uncomputed trees left to the jobs
# for left-or-right truncatable primes, duplicate subtrees are filled

import os
import re
import sys

job_re = re.compile(r'root_(\d+).bin')
input_dir = os.path.normpath(sys.argv[1])
prime_type = sys.argv[2]
base = int(sys.argv[3])
job_files = dict() # root (int) -> path (str)
used_jobs = set() # set of integers (keys in job_files)

# collect mapping of root to job file
for f in os.listdir(input_dir):
    m = job_re.fullmatch(f)
    if m:
        job_files[int(m.groups()[0])] = input_dir+'/'+f

BUFFER_SIZE = 2**24
root_bin = open(input_dir+'/root.bin','rb')
tree_bin = open(input_dir+'/tree.bin','wb',buffering=BUFFER_SIZE)

# read a byte from root_bin
def read_byte() -> int:
    return root_bin.read(1)[0]

# write a byte to tree_bin
def write_byte(b: int):
    tree_bin.write(bytes([b]))

# copy entire contents to tree_bin skipping skip bytes (the root value)
# otherwise if file does not exist, just write the end byte
# after writing the value byte, this function completes the subtree
def write_job_subtree(val: int, skip: int):
    if val in job_files:
        job_file = open(job_files[val],'rb',buffering=BUFFER_SIZE)
        job_file.read(skip)
        # TODO perhaps write chunks instead of reading entire file
        data = job_file.read()
        #print(f'writing {len(data)} bytes from {job_files[val]}')
        tree_bin.write(data)
        job_file.close()
        used_jobs.add(val)
    else:
        write_byte(255)

# val = current node number
# base = number base
# power = left append digit value
# each reads subtrees and end byte, writes subtrees and end byte
# the caller reads the value to calculate the next number for recursion
# the subtrees and end byte read are rewritten, but the job files
# are used to fill in the appropriate leaf recursion nodes

def primes_r(val: int, base: int):
    d = read_byte() # subtree value
    if d == 255: # leaf node, copy appropriate file if it exists
        write_job_subtree(val,1)
        return
    while d != 255:
        assert 0 < d < base
        write_byte(d)
        primes_r(base*val+d,base)
        d = read_byte()
    write_byte(255) # end

def primes_l(val: int, power: int, base: int):
    d = read_byte()
    if d == 255: # leaf node
        write_job_subtree(val,1)
        return
    while d != 255:
        assert 0 < d < base
        write_byte(d)
        primes_l(d*power+val,power*base,base)
        d = read_byte()
    write_byte(255) # end

def primes_lor(val: int, power: int, base: int):
    side = read_byte()
    if side == 255: # leaf node
        write_job_subtree(val,2)
        return
    while side != 255:
        digit = read_byte()
        assert 0 < digit < base
        write_byte(side)
        write_byte(digit)
        if side == 0:
            val2 = digit*power+val
        elif side == 1:
            val2 = base*val+digit
        else:
            val2 = 0 # to suppress error
            assert 0
        primes_lor(val2,power*base,base)
        side = read_byte()
    write_byte(255) # end

def primes_lar(val: int, power: int, base: int):
    ld = read_byte() # left digit
    if ld == 255: # leaf node
        write_job_subtree(val,2)
        return
    while ld != 255:
        rd = read_byte() # right digit
        if val == 0: # root of entire tree allows zeroes
            assert 0 <= ld < base
            assert 0 <= rd < base
        else: # otherwise digits cannot be zero
            assert 0 < ld < base
            assert 0 < rd < base
        write_byte(ld)
        write_byte(rd)
        if ld != 0: # double digit, multiply power by base^2
            primes_lar(ld*power+base*val+rd,power*base*base,base)
        else: # single digit, multiply power by base^1
            primes_lar(ld*power+base*val+rd,power*base,base)
        ld = read_byte()
    write_byte(255) # end

# extract root value bytes and write root bytes
assert read_byte() == 255
write_byte(255)
if prime_type == 'lor' or prime_type == 'lar':
    assert read_byte() == 255
    write_byte(255)

# call routine based on prime type
if prime_type == 'r':
    primes_r(0,base)
elif prime_type == 'l':
    primes_l(0,1,base)
elif prime_type == 'lor':
    primes_lor(0,1,base)
elif prime_type == 'lar':
    primes_lar(0,base,base)
else:
    assert 0

root_bin.close()
tree_bin.close()

unused = set(job_files.keys()) - used_jobs
if len(unused) > 0:
    print(f'warning: used job files for roots {unused}')
