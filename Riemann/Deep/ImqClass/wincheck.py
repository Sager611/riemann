# wincheck.py: positive-case regression of the sieve against the Mosunov-Jacobson tables on a window.
# dump: all (d, h) of the tables with lo <= d < hi (mjx dump.cfg); t4: isieve tier-4 output (fundamental d with M(d) <= N).
# Checks: every field with h <= N is uncertified (in t4) with M(d) <= h; every other t4 entry has h > N.
# Usage: python3 wincheck.py dump.txt run.t4 N
import sys
tab = {}
for l in open(sys.argv[1]):
    d, h = l.split(); tab[int(d)] = int(h)
t4 = {}
for l in open(sys.argv[2]):
    d, M = l.split(); t4[int(d)] = int(M)
N = int(sys.argv[3])
small = [d for d, h in tab.items() if h <= N]
missing = [d for d in small if d not in t4]
mbad = [d for d in small if d in t4 and t4[d] > tab[d]]
notab = [d for d in t4 if d not in tab]
extra_bad = [d for d in t4 if d in tab and tab[d] <= N and False]
extra = [d for d in t4 if d in tab and tab[d] > N]
print('window fields %d (range %d..%d), with h <= %d: %d; t4 entries %d' % (len(tab), min(tab), max(tab), N, len(small), len(t4)))
print('fields with h <= N certified by the sieve (must be 0): %d %s' % (len(missing), missing[:5]))
print('M(d) > h (must be 0): %d' % len(mbad))
print('t4 entries absent from the table (must be 0): %d %s' % (len(notab), notab[:5]))
print('t4 entries with h > N (correctly kept for exact check): %d' % len(extra))
print('PASS' if not missing and not mbad and not notab else 'FAIL')
