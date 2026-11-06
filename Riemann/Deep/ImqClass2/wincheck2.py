# wincheck2.py: positive-case check of the sieve against the Mosunov-Jacobson tables on a window [LO, HI].
# DUMP: every (d, h) of the tables with d in a range containing [LO, HI] and h <= HMAX (mjx dump.cfg);
# T4: isieve2 tier-4 output (fundamental d with exact M(d) <= N). Checks: (1) every field with h <= N in [LO, HI] is
# in T4 (the sieve did not certify it) with M(d) <= h; (2) every T4 entry is in [LO, HI]; T4 entries absent from the
# dump have h > HMAX >= N (the dump lists all fields with h <= HMAX; T4 entries are fundamental).
# Usage: python3 wincheck2.py DUMP T4 N LO HI HMAX
import sys
dump, t4f, N, LO, HI, HMAX = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
assert HMAX >= N
tab = {}
for l in open(dump):
    d, h = l.split(); d = int(d)
    if LO <= d <= HI: tab[d] = int(h)
t4 = {}
for l in open(t4f):
    d, M = l.split(); t4[int(d)] = int(M)
small = [d for d, h in tab.items() if h <= N]
missing = [d for d in small if d not in t4]
mbad = [d for d in small if d in t4 and t4[d] > tab[d]]
outside = [d for d in t4 if not (LO <= d <= HI)]
in_tab_big = [d for d in t4 if d in tab and tab[d] > N]
absent = [d for d in t4 if d not in tab]
print('window [%d, %d]: dump fields %d (h <= %d), with h <= N=%d: %d; t4 entries %d' % (LO, HI, len(tab), HMAX, N, len(small), len(t4)))
print('fields with h <= N certified by the sieve (must be 0): %d %s' % (len(missing), missing[:5]))
print('M(d) > h (must be 0): %d %s' % (len(mbad), mbad[:5]))
print('t4 entries outside the window (must be 0): %d' % len(outside))
print('t4 entries that are fields with h <= N: %d; with N < h <= HMAX: %d; absent from the dump (h > HMAX): %d' % (len(small) - len(missing), len(in_tab_big), len(absent)))
print('PASS' if not missing and not mbad and not outside else 'FAIL')
