# rowcheck.py: (1) rows N <= 5000 of our list table must equal ImqClass/lists_sha256.txt; (2) writes rowcheck.gp, which
# checks with PARI that the largest d of every row is a fundamental discriminant with class number h (h taken from
# the per-h table: the h whose largest d is that d), by qfbclassno and quadclassunit.
# Usage: python3 rowcheck.py final_sha.txt final_perh.txt ImqClass/lists_sha256.txt
import sys
sha, perh, old = sys.argv[1], sys.argv[2], sys.argv[3]
def rows(f):
    r = {}
    for l in open(f):
        if l.startswith('#') or not l.strip(): continue
        a = l.split(); r[int(a[0])] = (int(a[1]), int(a[2]), a[3])
    return r
new, prev = rows(sha), rows(old)
bad = 0
for n in sorted(prev):
    ok = new.get(n) == prev[n]; bad += not ok
    print('N=%d ours %s ImqClass %s %s' % (n, new.get(n), prev[n], 'EQUAL' if ok else 'DIFFERENT'))
print('rows N <= 5000 compared: %d, different: %d' % (len(prev), bad))
hof = {}
for l in open(perh):
    if l.startswith('#'): continue
    h, c, m = map(int, l.split()); hof.setdefault(m, []).append(h)
with open('rowcheck.gp', 'w') as g:
    for n in sorted(new):
        cnt, d, _ = new[n]
        hs = [h for h in hof.get(d, []) if h <= n]
        assert len(hs) == 1, (n, d, hs)
        g.write('{my(D=-%d,h=%d); print("N=%d d=%d h=",h," fundamental ",isfundamental(D)," qfbclassno ",qfbclassno(D)," quadclassunit ",quadclassunit(D).no, if(isfundamental(D)&&qfbclassno(D)==h&&quadclassunit(D).no==h," OK"," FAIL"))}\n' % (d, hs[0], n, d))
print('wrote rowcheck.gp (%d rows)' % len(new))
