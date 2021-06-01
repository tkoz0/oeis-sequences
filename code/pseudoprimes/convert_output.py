import sys
import zipfile

def factors(n):
    assert type(n) == int and n >= 2
    result = []
    while n % 2 == 0:
        result.append(2)
        n //= 2
    d = 3
    while d*d <= n:
        while n % d == 0:
            result.append(d)
            n //= d
        d += 2
    if n != 1:
        result.append(n)
    return result

assert sys.argv[1].startswith('fpp_base_')
i = 9
while sys.argv[1][i].isdigit(): i += 1
base = int(sys.argv[1][9:i])

input = zipfile.ZipFile(sys.argv[1])
contents = input.namelist() # list of files and dirs
#contents_info = input.infolist() # ZipInfo object for each

skipped_primes = [] # not in probable prime output -> factors of base
job_pseudoprimes = [None for _ in range(4096)]

sys.stderr.write('processing ".par" files in archive\n')

for arcname in contents:
    if not arcname.endswith('.par'): continue
    job_output = input.read(arcname).decode().splitlines()
    if len(job_output) < 5 or job_output[-1] != 'done':
        sys.stderr.write('warning: "%s" does not end with "done"\n'%arcname)
    else:
        assert job_output[0] == 'TYPE=FPP'
        assert job_output[1].startswith('BASE=')
        assert job_output[2].startswith('LO_BOUND=')
        assert job_output[3].startswith('HI_BOUND=')
        out_base = int(job_output[1][5:])
        assert out_base == base
        out_lo = int(job_output[2][9:])
        out_hi = int(job_output[3][9:])
        index = out_lo // 2**30 # job number based on pseudoprime range
        assert 0 <= index < 4096
        if index == 0:
            assert out_lo == 2 and out_hi == 2**30-1
        else:
            assert out_lo % 2**30 == 0 and out_hi == out_lo + 2**30-1
        job_nums = []
        prev_num = 0
        for i in range(4,len(job_output)-1):
            side = job_output[i][0]
            num = int(job_output[i][1:])
            assert num > prev_num
            prev_num = num
            assert out_lo <= num <= out_hi
            if side == '>': skipped_primes.append(num)
            elif side == '<':
                assert pow(base,num-1,num) == 1 # verify probable primality
                job_nums.append(num)
            else: assert 0
        assert job_pseudoprimes[index] == None
        job_pseudoprimes[index] = job_nums

input.close()
sys.stderr.write('finalizing\n')
assert set(skipped_primes) == set(factors(base))
assert all(job_nums != None for job_nums in job_pseudoprimes)
out_pp = []
for job_nums in job_pseudoprimes: out_pp += job_nums
#out_pp = sum(job_pseudoprimes,[])
#assert all(out_pp[i] < out_pp[i+1] for i in range(len(out_pp)-1))
sys.stderr.write('writing\n')
sys.stdout.write('\n'.join(map(str,out_pp)) + '\n')
sys.stderr.write('done\n')
