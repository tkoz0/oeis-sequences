import csv
import subprocess

# for each sequence run get_oeis.py -if ID

input = csv.reader(open('oeis_pseudoprimes.csv','r'))

header = next(input)

assert header[0] == 'base'

print('header =',header)

for row in input:
    base = row[0]
    for i in range(1,len(header)):
        if not row[i]: continue
        print('downloading OEIS A%s to %s_base_%s.txt'%(row[i],header[i],base))
        subprocess.run(['python3','get_oeis.py','-if',row[i]],
                       stdout=open('%s_base_%s.txt'%(header[i],base),'w'))
