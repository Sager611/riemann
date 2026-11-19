# perhsample.py (ImqAudit): (1) perh_h10000.txt consistency: every h <= 10^4 present with count >= 1, total, and the
# REPORT rows (count and largest d for h <= N) recomputed from the per-h table; (2) writes a random sample of 120 h
# (seed 1119) plus h = 1..5 and the 18 rows' h with their largest d to perh_sample_d.txt for hcount.
import random, re
P = 'Riemann/Deep/ImqClass2/'
rows = {}
for l in open(P + 'perh_h10000.txt'):
    if l.startswith('#') or not l.strip(): continue
    h, c, d = map(int, l.split()); rows[h] = (c, d)
missing = [h for h in range(1, 10001) if h not in rows or rows[h][0] < 1]
print('h present', len(rows), 'missing/zero', len(missing), 'total', sum(c for c, d in rows.values()))
rep = {}
for l in open(P + 'lists_sha256.txt'):
    a = l.split()
    if len(a) >= 3 and a[0].isdigit(): rep[int(a[0])] = (int(a[1]), int(a[2]))
bad = []
for N, (cnt, mx) in sorted(rep.items()):
    c2 = sum(rows[h][0] for h in range(1, N + 1)); m2 = max(rows[h][1] for h in range(1, N + 1))
    if (c2, m2) != (cnt, mx): bad.append((N, cnt, c2, mx, m2))
print('sha rows', len(rep), 'row mismatches vs perh', bad[:3])
random.seed(1119)
S = sorted(set(random.sample(range(1, 10001), 120)) | {1, 2, 3, 4, 5, 9819, 10000})
with open('.scratch/imqaudit/perh_sample.txt', 'w') as f:
    for h in S: f.write('%d %d\n' % (rows[h][1], h))
print('sample', len(S))
