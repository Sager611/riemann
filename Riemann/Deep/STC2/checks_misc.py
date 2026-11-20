# checks_misc.py -- (i) kernels.json vs a fresh extraction of Appendix D (sol.json blocks) from the PDF text;
# (ii) G(N0) = g(N0) >= gt(400) (so G(N0) = max(g(N0), gt(400)) = g(N0) in cz_arb.py); (iii) the sign of nu near L = 653 for l150.
import re, json, math, sys
from fractions import Fraction as Fr
from math import factorial
from flint import arb, ctx
import mpmath as mp
ctx.prec = 80
K = json.load(open('Riemann/Deep/STC2/kernels.json'))
L = open('.scratch/stc2/cg_raw.txt').read().split('\n')
found = {}
i = 0
while i < len(L):
    m = re.match(r'^line_([A-Z]\d+)/sol\.json', L[i].strip())
    if m:
        buf = []; j = i + 1
        while j < len(L) and L[j].strip() != '}': buf.append(L[j]); j += 1
        txt = '\n'.join(ln for ln in buf if not re.search(r'EXPLICIT AND ASYMPTOTIC|CASA|^\s*\d+\s*$', ln)) + '\n}'
        try: found[m.group(1)] = json.loads(txt)
        except Exception as e: print(m.group(1), 'parse fail', e)
        i = j
    i += 1
for nm, v in K.items():
    f = found.get(nm)
    if f is None: print(nm, 'NOT FOUND in fresh extraction'); continue
    same = f['c'] == v['c'] and f['B'] == v['B'] and f['d'] == v['d']
    print(nm, 'kernels.json == fresh Appendix D extraction:', same, ' (n = %d, d = %s)' % (len(f['c']), f['d']))
for nm, k in K.items():
    dq = Fr(k['d']).limit_denominator(10**6); d = arb(dq.numerator)/dq.denominator; s1 = d + arb(1)/2
    B = [(X, int(j)) for X, j in k['B']]; C = [arb(c) for c in k['c']]; DEG = max(j for _, j in B) - 1
    pa = [arb(0)]*(DEG+1); pb = [arb(0)]*(DEG+1)
    for (X, j), c in zip(B, C):
        if X == 'R': pa[j-1] -= c/factorial(j-1)
        else: pb[j-1] -= c/factorial(j-1)
    L0 = arb(10**8).log(); ed = (-d*L0).exp()
    A = sum((pa[i]*L0**i for i in range(DEG+1)), arb(0)); Bv = (ed - 2)/L0 + sum((pb[i]*L0**i for i in range(DEG+1)), arb(0))
    g = (A*A + Bv*Bv).sqrt()/arb.pi()*(-s1*L0).exp()
    LM = arb(400); gt = (sum(((abs(pa[i]) + abs(pb[i]))*LM**i for i in range(DEG+1)), arb(0)) + 3)/arb.pi()*(-s1*LM).exp()
    print(nm, 'g(1e8) = %s, gt(400) = %s, g(1e8) > gt(400): %s' % (g.str(5), gt.str(3), bool(g > gt)))
# (iii) nu(L) = s1 (A^2+B^2) - (A A' + B B') for l150 near L = 653 (mpmath, 50 digits)
mp.mp.dps = 50
k = K['L150']; d = mp.mpf(17)/20; s1 = d + mp.mpf(1)/2
pa = [mp.mpf(0)]*5; pb = [mp.mpf(0)]*5
for (X, j), c in zip(k['B'], k['c']):
    if X == 'R': pa[j-1] -= mp.mpf(c)/mp.factorial(j-1)
    else: pb[j-1] -= mp.mpf(c)/mp.factorial(j-1)
def nu(L):
    A = mp.polyval(pa[::-1], L); A1 = mp.polyval([i*pa[i] for i in range(4, 0, -1)], L)
    ed = mp.e**(-d*L); Bv = (ed - 2)/L + mp.polyval(pb[::-1], L); B1 = -d*ed/L - (ed - 2)/L**2 + mp.polyval([i*pb[i] for i in range(4, 0, -1)], L)
    return s1*(A*A + Bv*Bv) - (A*A1 + Bv*B1), A, Bv
neg = [(L, float(nu(L)[0])) for L in [mp.mpf(x)/10 for x in range(6400, 6700)] if nu(L)[0] < 0]
print('l150: nu < 0 at %d of 300 grid points in [640, 670]; first/last: %s' % (len(neg), (neg[0][0], neg[-1][0]) if neg else None))
Lm = min([mp.mpf(x)/10 for x in range(6400, 6700)], key=lambda L: abs(nu(L)[1]))
print('l150: |pi a(L)| smallest at L = %s: pi a = %.3e, pi b = %.3e, d log g/dL = %.3f' % (Lm, float(nu(Lm)[1]), float(nu(Lm)[2]), float(-(nu(Lm)[0])/(nu(Lm)[1]**2 + nu(Lm)[2]**2))))

