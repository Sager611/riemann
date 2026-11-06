# t4check.py: tier-4 results of the low sieve chunks against the tables, and the list of sieve-found fields above them.
# T4: isieve2 tier-4 entries "d M" (fundamental d with exact M(d) <= N) of the chunks [LO, HI];
# FOUND: t4.gp output lines "FOUND d h" (classno_upto: h <= N, exact); DUMP: all table fields with h <= N and d >= 1.2e10
# (mjx dump.cfg), complete for d < KTOP = K*2^28 (the tables processed). Checks on the overlap [LO, min(HI, KTOP-1)]:
# the table fields with h <= N there are exactly the FOUND entries there, with the same h, all in T4 with M <= h.
# Writes APPEND: the FOUND fields with d >= KTOP ("d h", increasing), for mjx --append.
# Usage: python3 t4check.py T4 FOUND DUMP N LO HI K APPEND
import sys
t4f, fdf, dumpf, N, LO, HI, K, outf = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]), int(sys.argv[7]), sys.argv[8]
KTOP = K * 2**28
t4 = {}
for l in open(t4f):
    d, M = l.split(); t4[int(d)] = int(M)
found = {}
for l in open(fdf):
    a = l.split()
    if len(a) == 3 and a[0] == 'FOUND': found[int(a[1])] = int(a[2])
tab = {}
for l in open(dumpf):
    a = l.split()
    if len(a) != 2: continue
    d, h = int(a[0]), int(a[1])
    if LO <= d <= HI and d < KTOP and h <= N: tab[d] = h
ovhi = min(HI, KTOP - 1)
fo = {d: h for d, h in found.items() if d <= ovhi}
miss = [d for d in tab if d not in t4]
mbad = [d for d in tab if d in t4 and t4[d] > tab[d]]
diff = sorted(set(tab.items()) ^ set(fo.items()))
notin4 = [d for d in found if d not in t4]
above = sorted((d, h) for d, h in found.items() if d >= KTOP)
print('chunks [%d, %d], tables to %d (k < %d); overlap [%d, %d]' % (LO, HI, KTOP, K, LO, ovhi))
print('t4 entries %d; FOUND (h <= %d) %d, of them in the overlap %d, above the tables %d' % (len(t4), N, len(found), len(fo), len(above)))
print('table fields with h <= N in the overlap: %d; not in t4 (must be 0): %d; M > h (must be 0): %d' % (len(tab), len(miss), len(mbad)))
print('overlap: table fields vs FOUND, symmetric difference (must be empty): %d %s' % (len(diff), diff[:5]))
print('FOUND not in t4 (must be 0): %d' % len(notin4))
with open(outf, 'w') as g:
    for d, h in above: g.write('%d %d\n' % (d, h))
print('wrote %s (%d fields above the tables)' % (outf, len(above)))
print('PASS' if not miss and not mbad and not diff and not notin4 else 'FAIL')
