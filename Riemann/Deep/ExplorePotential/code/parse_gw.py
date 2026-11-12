# parse_gw.py: gram_win -zeros output -> zero ordinates (tc + offsets), contiguity check
import sys, re, numpy as np
tc = {}; zs = []; gram = {}
for line in open(sys.argv[1]):
    if line.startswith('ZG'):
        m = re.match(r'ZG n=(\d+) tc=(\S+) W=(\d+)', line); n = int(m.group(1)); tc[n] = float(m.group(2)); W = int(m.group(3))
        parts = line.split()[5:]
        d = [float(p.split(':')[0]) for p in parts]
        gram[n] = (tc[n] + d[0], tc[n] + d[-1])
    elif line.startswith('ZZ'):
        p = line.split(); n = int(p[1][2:])
        zs.extend(tc[n] + float(x) for x in p[2:])
z = np.array(sorted(zs)); ns = sorted(gram)
gaps = [gram[ns[i+1]][0] - gram[ns[i]][1] for i in range(len(ns) - 1)]
print('windows', len(ns), 'zeros', len(z), 'Gram intervals', 62*len(ns), 'max endpoint mismatch', max(abs(g) for g in gaps))
print('range', gram[ns[0]][0], gram[ns[-1]][1])
with open(sys.argv[2], 'w') as f:
    f.write('# zeta zeros from gram_win -zeros, range %.9f %.9f\n' % (gram[ns[0]][0], gram[ns[-1]][1]))
    for x in z: f.write(repr(float(x)) + '\n')
