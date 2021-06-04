import argparse
import requests
import sys

def get_oeis_b_file(id):
    req = requests.get('https://oeis.org/A%s/b%s.txt'%(id,id))
    assert req.ok
    return req.text

ap = argparse.ArgumentParser(description='download OEIS B file')
ap.add_argument('id',help='the 6 digit OEIS sequence number')
ap.add_argument('-c','--keep-comments',
    help='keeps comments (starting with #, removed by default)',
    action='store_true')
ap.add_argument('-b','--keep-blanks',
    help='keeps blank lines (remove them by default)',
    action='store_true')
ap.add_argument('-i','--keep-indexes',
    help='keeps sequence indexes (removed by default)',
    action='store_true')
ap.add_argument('-f','--check-format',
    help='runs extra checks on B file format',
    action='store_true')
args = ap.parse_args()

assert len(args.id) == 6 and args.id.isdigit()

lines = get_oeis_b_file(args.id).splitlines()

def check_line(line):
    if line == '': return True # blank
    if line.startswith('#'): return True # comment
    if line != line.strip(): return False # extra whitespace not allowed
    line2 = line.split()
    if len(line2) != 2: return False # 2 integers
    if not line2[0].isdigit(): return False
    if not line2[1].isdigit(): return False
    return '%d %d'%(int(line2[0]),int(line2[1])) == line # proper format

if args.check_format:
    for line in lines:
        assert check_line(line)

if not args.keep_comments:
    lines = list(filter(lambda l : not l.startswith('#'), lines))

if not args.keep_blanks:
    lines = list(filter(lambda l : len(l) > 0, lines))

def remove_index(l):
    if l == '' or l.startswith('#'): return l
    return l[l.find(' ')+1:]

if not args.keep_indexes:
    lines = list(map(remove_index,lines))

print('\n'.join(lines))
