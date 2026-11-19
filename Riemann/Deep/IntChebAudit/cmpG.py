import sys, re
mine = {}; 
for l in open(sys.argv[1]):
    m = re.match(r'G (\d+) (\S+) resid (\S+)', l)
    if m: mine[int(m.group(1))] = (float(m.group(2)), m.group(2), float(m.group(3)))
gc = {}
for l in open(sys.argv[2]):
    p = l.split()
    if len(p) == 3 and p[0].isdigit(): gc[int(p[0])] = (float(p[1]), float(p[2]))
common = sorted(set(mine) & set(gc))
from decimal import Decimal as D
diffs = [(D(mine[i][1]) - D(gc[i][0]), i) for i in common]   # D(float) = exact binary value of gcert's floored double
print('compared', len(common), 'of', len(gc), '; mine - gcert: min %.3e max %.3e ; negative (gcert above mine): %d' % (float(min(diffs)[0]), float(max(diffs)[0]), sum(1 for d, i in diffs if d < 0)))
print('max gcert radius %.3e ; max root residual %.3e' % (max(v[1] for v in gc.values()), max(v[2] for v in mine.values())))
lo = sorted(common, key=lambda i: mine[i][0])[:4]
for i in lo: print('lowest G: i=%d mine %s gcert %r' % (i, mine[i][1], gc[i][0]))
