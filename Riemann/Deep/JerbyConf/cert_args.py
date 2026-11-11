import re, sys
from decimal import Decimal, getcontext
getcontext().prec = 40
L = open('logs/zeros_low.log').read().split('\n')
want = [int(x) for x in sys.argv[1:]]
args = []
for i, line in enumerate(L):
    if not line.startswith('ZG'): continue
    n = int(re.search(r'n=(\d+)', line).group(1))
    if n not in want: continue
    tc = re.search(r'tc=([\d.e+]+)', line).group(1); W = int(re.search(r'W=(\d+)', line).group(1))
    gz = [x.split(':') for x in line.split()[5:]]; dg = [float(a) for a, b in gz]
    zz = [float(x) for x in L[i+1].split()[2:]]
    dn = dg[W]; left = [z for z in zz if dg[W-1] < z < dn]; right = [z for z in zz if dn < z < dg[W+1]]
    pair = left if len(left) == 2 else right
    mid = Decimal(tc) + (Decimal(repr(pair[0])) + Decimal(repr(pair[1])))/2
    args.append(f'{n}:{mid:.25f}')
print(' '.join(args))
