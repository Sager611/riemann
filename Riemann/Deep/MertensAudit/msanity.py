# msanity.py -- sanity checks: exact M(x) (own sieve to 1e7; OEIS A084237 M(10^k), k <= 23; Hurst 1610.08551 TeX l.479-512 extremes),
# reported-constant rounding, and the B' rounding of theta141 (exact p/q of Riemann/Tree/R466.lean l.10-12).
import json, math, re
from fractions import Fraction
from decimal import Decimal, getcontext
import numpy as np
N = 10**7
sp = np.ones(N + 1, dtype=bool); sp[:2] = False
for p in range(2, int(N**0.5) + 1):
    if sp[p]: sp[p*p::p] = False
mu = np.ones(N + 1, dtype=np.int8); mu[0] = 0
for p in np.nonzero(sp)[0]:
    mu[p::p] *= -1
    if p*p <= N: mu[p*p::p*p] = 0
M = np.cumsum(mu, dtype=np.int32)
oe = {int(a): int(b) for a, b in (l.split() for l in open('/Users/adriansager/playground/riemann/.scratch/mertensaudit/dl/a084237.txt') if l[0].isdigit())}
assert all(M[10**k] == oe[k] for k in range(8)), [(k, M[10**k], oe[k]) for k in range(8)]
del sp
recs = [('9.69e4', 0.9382), ('1.52e4', 0.9409), ('1.84e6', 0.9382), ('7.91e5', 0.9409)]
qmax = (0.0, 0); rmax = {c: (0.0, 0) for c, _ in recs}
for lo in range(2, N + 1, 10**6):                     # chunks of 1e6: every array <= 8 MB (memory cap)
    x = np.arange(lo, min(lo + 10**6, N + 1), dtype=np.float64); Mm = np.abs(M[lo:lo + len(x)]).astype(np.float64)
    q = Mm/np.sqrt(x); i = int(np.argmax(q))
    if q[i] > qmax[0]: qmax = (float(q[i]), lo + i)
    for c, th in recs:
        r = Mm/(float(c)*x**th); i = int(np.argmax(r))
        if r[i] > rmax[c][0]: rmax[c] = (float(r[i]), lo + i)
print('own sieve: M(10^k) = OEIS for k <= 7; max |M(x)|/sqrt(x) on [2, 1e7] =', round(qmax[0], 5), 'at x =', qmax[1])
for c, th in recs:
    print('bound', c, 'x^%s' % th, ': max |M(x)|/bound on [2, 1e7] =', '%.3e' % rmax[c][0], 'at x =', rmax[c][1], '; at x = 1:', '%.3e' % (1/float(c)))
pts = [(10**k, oe[k]) for k in range(8, 24)] + [(7766842813, 50286), (71578936427177, -4440015), (6287915599821430, -35629003), (7332940231978758, 40371499)]
for c, th in recs[:2]:
    worst = max(abs(m)/(float(c)*n**th) for n, m in pts)
    print('bound', c, 'x^%s' % th, ': max |M|/bound over OEIS 10^8..10^23 and Hurst extremes =', '%.3e' % worst)
# rounding of the reported constants
R = json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/Mertens78/results.json'))
def sci(l): e = math.floor(l); return '%.6fe%d' % (10**(l - e), e)
for r in R: print('M78', r['name'], 'theta', r['theta'], "C'(all x) =", sci(r['allx']['Cprime_log10']), 'C(x>=xa) =', sci(r['best']['log10C']), 'xa = 1e%d' % r['best']['log10xa'])
H = json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/HP141/MertensCH/logs/cert_ch_results_bp.json'))
for r in H: print('HP141', r['theta'], r['c'], [p[0] + ' ' + p[-1][:22] for p in r['pieces']])
# theta141
L = open('/Users/adriansager/playground/riemann/Riemann/Tree/R466.lean').read()
m = re.search(r'theta_exact = (\d+)\s*/\s*(\d+)', L)
TH = Fraction(int(m.group(1)), int(m.group(2)))
getcontext().prec = 40
print('theta141 =', Decimal(TH.numerator)/Decimal(TH.denominator))
for s in ['0.87334585610671858', '0.8733458561067185768', '0.8733458561067185768011848']:
    print(s, '>= theta141 ?', Fraction(s) >= TH, ' diff', float(Fraction(s) - TH))
print('D = 7/8 - Bprime =', Fraction(7, 8) - Fraction('0.87334585610671858'))
