#!/usr/bin/env python3
"""chebaudit: compare ChebExp/logs/direct_range.txt (PARI idealprimedec/idealfrobenius/rnfidealprimedec) with
(a) own Kronecker-symbol least primes for the 471 C2 rows, (b) own factormod/coset witnesses logs/direct_check.txt (92 rows)."""
import collections
F = '/Users/adriansager/playground/riemann/Riemann/Deep/ChebExp/logs/direct_range.txt'
rows = [l.split() for l in open(F) if l.strip()]
def kron(d, p):
    if p == 2: return 0 if d % 2 == 0 else (1 if d % 8 in (1, 7) else -1)
    if d % p == 0: return 0
    return 1 if pow(d % p, (p - 1) // 2, p) == 1 else -1
P = [p for p in range(2, 5000) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
def sqfree(m): return all(m % (q * q) for q in P if q * q <= m)
def fund(d): return d not in (0, 1) and ((d % 4 == 1 and sqfree(abs(d))) or (d % 4 == 0 and (d // 4) % 4 in (2, 3) and sqfree(abs(d // 4))))
mine = collections.defaultdict(list)
for d in range(-257, 258):
    if abs(d) >= 3 and fund(d):
        sp = next(p for p in P if kron(d, p) == 1); ine = next(p for p in P if kron(d, p) == -1); ra = next(p for p in P if d % p == 0)
        mine[abs(d)] += [('L', min(sp, ra)), ('Q1', sp), ('Qs', ine)]
fil = collections.defaultdict(list)
for r in rows:
    if r[0] != 'C2': continue
    desc = ' '.join(r[3:-1]); p = int(r[-1])
    fil[int(r[2])].append(('L' if desc.startswith('K=L') else ('Q1' if 'elt[1, 2]' in desc else 'Qs'), p))
eq = all(sorted(fil[D]) == sorted(mine[D]) for D in set(fil) | set(mine))
print('(a) C2: %d file rows, %d discriminants; own Kronecker recomputation identical (as multisets per |d|): %s' % (sum(len(v) for v in fil.values()), len(fil), eq))
# (b)
mc = collections.defaultdict(list)
for l in open('logs/direct_check.txt'):
    t = l.split()
    if len(t) < 5 or t[0] not in ('C3', 'S3', 'C4', 'V4', 'D4', 'A4', 'S4'): continue
    G, D = t[0], int(t[2]); kind = 'L' if t[3] == 'K=L' else ('Q' if t[3] == 'K=Q' else t[4]); mc[(G, D)].append((kind, int(t[-2]), t[-1] == '1'))
fc = collections.defaultdict(list)
for r in rows:
    if r[0] == 'C2': continue
    G, D = r[0], int(r[2]); desc = ' '.join(r[3:-1])
    kind = 'L' if desc.startswith('K=L') else ('Q' if desc.startswith('K=Q') else desc.split()[-1])
    fc[(G, D)].append((kind, int(r[-1])))
nrow = sum(len(v) for v in fc.values()); nmine = sum(len(v) for v in mc.values())
allok = all(ok for v in mc.values() for _, _, ok in v)
ge = []; same = 0; tot = 0
for key in fc:
    for kind in set(k for k, _ in fc[key]):
        a = sorted(p for k, p in fc[key] if k == kind); b = sorted(p for k, p, _ in mc[key] if k == kind)
        tot += len(a)
        if len(a) != len(b): ge.append(('COUNT', key, kind, a, b)); continue
        if a == b: same += len(a)
        elif all(x <= y for x, y in zip(a, b)): ge.append(('file<mine', key, kind, a, b))
        else: ge.append(('FILE>MINE', key, kind, a, b))
print('(b) non-C2: file %d rows, mine %d rows, same field keys: %s; all own witnesses <= (0.62 log 3D)^(1/d2): %s' % (nrow, nmine, set(fc) == set(mc), allok))
print('    rows with identical least prime: %d of %d; differences:' % (same, tot))
for g in ge: print('     ', g)
