# perhcheck.py: writes perhcheck.gp: for 200 random h <= 10000 (seed 1) the largest d of the per-h table must be a
# fundamental discriminant with class number h (PARI qfbclassno; a consistency check, not part of the proof).
import random, sys
perh = sys.argv[1]
rows = {}
for l in open(perh):
    if l.startswith('#'): continue
    h, c, m = map(int, l.split()); rows[h] = (c, m)
random.seed(1)
hs = sorted(random.sample(range(1, 10001), 200))
with open('perhcheck.gp', 'w') as g:
    g.write('bad = 0; n = 0;\n')
    for h in hs:
        c, m = rows[h]
        g.write('n++; if(!isfundamental(-%d) || qfbclassno(-%d) != %d, bad++; print("BAD h=%d d=%d"));\n' % (m, m, h, h, m))
    g.write('print("perhcheck: ", n, " random h, bad ", bad);\n')
print('wrote perhcheck.gp')
