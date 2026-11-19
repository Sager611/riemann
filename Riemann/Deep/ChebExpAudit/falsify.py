#!/usr/bin/env python3
"""chebaudit falsification: least primes with prescribed Frobenius vs bounds (A), (B) of chebotarev-numeric.
(1) every quadratic field, 3 <= |d| <= 1e5: K = Q (split = class {1}, inert = class {sigma}) and K = L (least degree-one prime);
(2) Q(zeta_q)/Q, 3 <= q <= QMAX, q != 2 mod 4, every class a mod q; (3) class count sum phi(q), 3 <= q <= 2000;
(4) max theta(p)/p for p <= 3e6 (Rosser-Schoenfeld 1962 Thm 9 constant 1.01624); (5) the C2 rows of ChebExp/logs/direct_range.txt."""
import math, sys, collections
import mpmath as mp
mp.mp.prec = 256
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 600
NP = 3_000_000
b = bytearray([1]) * (NP + 1); b[0] = b[1] = 0
for i in range(2, int(NP ** 0.5) + 1):
    if b[i]: b[i * i::i] = bytearray(len(range(i * i, NP + 1, i)))
P = [i for i in range(NP + 1) if b[i]]
d2 = 1 - (777 - mp.sqrt(1081)) / 852
dA = mp.mpf(1) / 32
def kron(d, p):
    if p == 2: return 0 if d % 2 == 0 else (1 if d % 8 in (1, 7) else -1)
    if d % p == 0: return 0
    return 1 if pow(d % p, (p - 1) // 2, p) == 1 else -1
def sqfree(m):
    m = abs(m)
    for p in P:
        if p * p > m: return True
        if m % (p * p) == 0: return False
    return True
def fund(d):
    if d in (0, 1): return False
    if d % 4 == 1: return sqfree(d)
    if d % 4 == 0: return (d // 4) % 4 in (2, 3) and sqfree(d // 4)
    return False
# bounds: (name, delta, c, exponent, min log Delta)
BOUNDS = [('A cert 0.484^32', dA, '0.484', None, 0), ('A page 0.49^32', dA, '0.49', None, 0), ('A 0.29 (logD>=10)', dA, '0.29', None, 10),
          ('B cert 0.62^(1/d2)', d2, '0.62', None, 0), ('B page 0.73^7.8978', d2, '0.73', mp.mpf('7.8978'), 0), ('B 0.595 (logD>=10)', d2, '0.595', None, 10)]
def check(N, Delta, tag, stats):
    l3 = mp.log(3 * Delta)
    for (nm, de, c, ex, lmin) in BOUNDS:
        if mp.log(Delta) < lmin: continue
        e = ex if ex is not None else 1 / de
        if not (N <= (mp.mpf(c) * l3) ** e): stats['viol'].append((nm, tag, Delta, N))
    for key, de in (('A', dA), ('B', d2)):
        v = mp.exp(de * mp.log(N)) / l3
        if v > stats['max' + key][0]: stats['max' + key] = (v, tag, Delta, N)
stats = dict(viol=[], maxA=(0, None, None, None), maxB=(0, None, None, None))
quad = {}; nq = 0
for d in range(-100000, 100001):
    if abs(d) < 3 or not fund(d): continue
    nq += 1
    sp = next(p for p in P if kron(d, p) == 1); ine = next(p for p in P if kron(d, p) == -1); ra = next(p for p in P if d % p == 0)
    deg1 = min(sp, ra); D = abs(d)
    quad[d] = (sp, ine, deg1)
    check(sp, D, 'quad d=%d split' % d, stats); check(ine, D, 'quad d=%d inert' % d, stats); check(deg1, D, 'quad d=%d K=L' % d, stats)
print('(1) quadratic fields 3<=|d|<=1e5:', nq, 'fields,', 3 * nq, 'cases; violations:', stats['viol'][:5])
print('    max N^(1/32)/log 3D = %s at %s (D=%s, N=%s)' % (mp.nstr(stats['maxA'][0], 10), *stats['maxA'][1:]))
print('    max N^(d2)/log 3D   = %s at %s (D=%s, N=%s)' % (mp.nstr(stats['maxB'][0], 10), *stats['maxB'][1:]))
print('    exact sharpness constants: 7^(1/32)/log 9 =', mp.nstr(mp.mpf(7) ** dA / mp.log(9), 12), ' 7^d2/log 9 =', mp.nstr(mp.mpf(7) ** d2 / mp.log(9), 12),
      ' 3^d2/log 9 =', mp.nstr(mp.mpf(3) ** d2 / mp.log(9), 6), ' 11^d2/log 15 =', mp.nstr(mp.mpf(11) ** d2 / mp.log(15), 6))
print('    (0.484 log 9)^32 =', mp.nstr((mp.mpf('0.484') * mp.log(9)) ** 32, 8), ' (0.49 log 9)^32 =', mp.nstr((mp.mpf('0.49') * mp.log(9)) ** 32, 8),
      ' (0.62 log 9)^(1/d2) =', mp.nstr((mp.mpf('0.62') * mp.log(9)) ** (1 / d2), 8), ' (0.73 log 9)^7.8978 =', mp.nstr((mp.mpf('0.73') * mp.log(9)) ** mp.mpf('7.8978'), 8))
# second-worst quadratic cases for A
worstA = sorted(((mp.mpf(N) ** dA / mp.log(3 * abs(d)), d, k, N) for d, t in quad.items() for k, N in zip(('split', 'inert', 'K=L'), t)), reverse=True)[:5]
print('    top-5 A constants:', [(mp.nstr(v, 6), d, k, N) for v, d, k, N in worstA])
# (2) cyclotomic
st2 = dict(viol=[], maxA=(0, None, None, None), maxB=(0, None, None, None)); ncl = 0
def phi(n):
    r = n; m = n
    for p in P:
        if p * p > m: break
        if m % p == 0:
            r -= r // p
            while m % p == 0: m //= p
    if m > 1: r -= r // m
    return r
for q in range(3, QMAX + 1):
    if q % 4 == 2: continue
    f = phi(q); ps = [p for p in P[:200] if q % p == 0]
    logD = f * math.log(q) - sum(f * math.log(p) / (p - 1) for p in ps)
    Delta = mp.exp(mp.mpf(logD))
    seen = {}; 
    for p in P:
        if q % p == 0: continue
        a = p % q
        if a not in seen: seen[a] = p
        if len(seen) == f: break
    assert len(seen) == f, q
    ncl += f; worst = max(seen.values())
    check(worst, Delta, 'Q(zeta_%d) worst class' % q, st2)
print('(2) Q(zeta_q)/Q, 3<=q<=%d, q != 2 mod 4: %d classes; violations: %s; max A const %s at %s; max B const %s at %s' % (QMAX, ncl, st2['viol'][:3],
      mp.nstr(st2['maxA'][0], 6), st2['maxA'][1], mp.nstr(st2['maxB'][0], 6), st2['maxB'][1]))
s_all = sum(phi(q) for q in range(3, 2001)); s_no2 = sum(phi(q) for q in range(3, 2001) if q % 4 != 2)
print('(3) sum phi(q), 3<=q<=2000: %d; excluding q = 2 mod 4: %d   (page: 1,013,926)' % (s_all, s_no2))
th = 0.0; best = (0, 0)
for p in P:
    th += math.log(p)
    if th / p > best[0]: best = (th / p, p)
thp = mp.fsum(mp.log(p) for p in P if p <= best[1])
print('(4) max theta(p)/p over p <= 3e6: %s at p = %d  (RS62 Thm 9: theta(x) < 1.01624 x)' % (mp.nstr(thp / best[1], 10), best[1]))
rows = [l.split() for l in open('/Users/adriansager/playground/riemann/Riemann/Deep/ChebExp/logs/direct_range.txt') if l.startswith('C2 ')]
file_ms = collections.defaultdict(list)
for r in rows: file_ms[int(r[2])].append((r[3] + (' ' + r[4] if r[3] == 'K=Q' else ''), int(r[-1])))
mine = collections.defaultdict(list)
for d, (sp, ine, deg1) in quad.items():
    if abs(d) <= 257: mine[abs(d)] += [('K=L', deg1), ('K=Q elt[1, 2]', sp), ('K=Q elt[2, 1]', ine)]
# PARI labels: identity is elt[1, 2]
fm = {D: sorted((k.replace('K=L id', 'K=L'), p) for k, p in v) for D, v in file_ms.items()}
mm = {D: sorted(v) for D, v in mine.items()}
print('(5) C2 rows of direct_range.txt: %d rows, %d discriminants; my recomputation equal as multisets: %s' % (len(rows), len(fm), fm == mm))
if fm != mm: print('    first diffs:', [(D, fm.get(D), mm.get(D)) for D in sorted(set(fm) | set(mm)) if fm.get(D) != mm.get(D)][:4])
