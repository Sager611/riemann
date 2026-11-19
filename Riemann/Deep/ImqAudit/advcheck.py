# advcheck.py (ImqAudit): checks hcount output "d fund h M M2" on the sieve-range samples:
# M = M2 (reduced-form count with a <= X = multiplicative sum), M <= h, h > 10^4 except d = 28154696707,
# and M >= 10^4 + 1 except the tier-3 lanes of the sieve (only d < 3.8e10 can be tier 3; the 33 tier-4 d are known).
import sys
t4 = {int(l.split()[0]) for l in open(sys.argv[2]) if l[0].isdigit()} if len(sys.argv) > 2 else set()
rows = [l.split() for l in open(sys.argv[1])]
bad = []; low = []; mins = []; ratio = []
for d, f, h, M, M2 in rows:
    d = int(d)
    if f != '1': bad.append((d, 'nonfund')); continue
    h, M, M2 = int(h), int(M), int(M2)
    if M != M2: bad.append((d, 'M!=M2', M, M2))
    if M > h: bad.append((d, 'M>h', M, h))
    if h <= 10000 and d != 28154696707: bad.append((d, 'h<=1e4', h))
    if M <= 10000: low.append((d, M, h, d in t4))
    mins.append((h, d, M)); ratio.append(M / h)
mins.sort()
print('rows', len(rows), 'bad', len(bad), bad[:5])
print('M <= 10^4 (must be tier 3, d < 3.8e10):', len(low), [x for x in low if x[0] >= 38_000_000_000 or not x[3]][:5])
print('smallest h:', mins[:5])
print('M/h min %.4f max %.4f' % (min(ratio), max(ratio)))
