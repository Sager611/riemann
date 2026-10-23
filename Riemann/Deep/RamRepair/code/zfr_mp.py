# zfr_mp.py -- certified lower bounds m_p <= min_x G_p(1, x), G_p(1, x) = sum_{m<=15} p^{-m} P(m x), P(y) = |sum c_k e^{iky}|^2/14912370
# (own computation; compare BTY Table 1 and Padhi Remark 8.2).  Output zfr_mp.json {p: m_p (string)}.
import math, json, sys
import numpy as np
from flint import arb, ctx
from fractions import Fraction as Fr
from rig import fr, ball as rball
ctx.prec = 96
CC = [4,-8,2,20,-9,-34,27,91,-27,-201,32,895,1949,2389,1896,949,239]
NN = 14912370
kv = np.arange(17)
def Pf(y):
    y = np.asarray(y)[..., None]
    C = (CC*np.cos(kv*y)).sum(-1); S = (CC*np.sin(kv*y)).sum(-1); return (C*C + S*S)/NN
def Gf(p, x): return sum(p**(-m)*Pf(m*x) for m in range(1, 16))
def Ga(p, X):         # arb ball X: value via direct evaluation (used at centres)
    acc = arb(0)
    for m in range(1, 16):
        y = m*X; C = sum((CC[k]*(k*y).cos() for k in range(17)), arb(0)); S = sum((CC[k]*(k*y).sin() for k in range(17)), arb(0))
        acc += arb(p)**(-m)*(C*C + S*S)
    return acc/NN
def Gd(p, X):         # derivative enclosure on the ball X
    acc = arb(0)
    for m in range(1, 16):
        y = m*X
        C = sum((CC[k]*(k*y).cos() for k in range(17)), arb(0)); S = sum((CC[k]*(k*y).sin() for k in range(17)), arb(0))
        Cd = sum((-k*CC[k]*(k*y).sin() for k in range(17)), arb(0)); Sd = sum((k*CC[k]*(k*y).cos() for k in range(17)), arb(0))
        acc += arb(p)**(-m)*m*2*(C*Cd + S*Sd)
    return acc/NN
def certify(p, mp_s):
    """G_p(1, x) > every point of the ball arb(mp_s) for all x = pi*u, u in [0, 1] (G is even and 2 pi-periodic); boxes are exact
    rational multiples of pi, mean value enclosure G(c) + G'(X)(X - c) on each box."""
    mp_ = arb(mp_s); PIa = arb.pi(); stack = [(Fr(0), Fr(1))]; n = 0
    while stack:
        lo, hi = stack.pop(); n += 1
        c = PIa*fr((lo + hi)/2); X = PIa*rball(lo, hi)
        if Ga(p, c) + Gd(p, X)*(X - c) - mp_ > 0: continue
        if hi - lo < Fr(1, 10**13): return False, n
        mid = (lo + hi)/2; stack += [(lo, mid), (mid, hi)]
    return True, n
primes = [p for p in range(2, 60) if all(p % q for q in range(2, int(p**0.5) + 1))]
out = {}
for p in primes:
    xs = np.linspace(0, math.pi, 200001); g = Gf(p, xs); i = int(np.argmin(g))
    # refine by golden section around the grid minimum
    lo, hi = xs[max(i-1, 0)], xs[min(i+1, len(xs)-1)]
    for _ in range(200):
        m1, m2 = lo + (hi - lo)*0.382, lo + (hi - lo)*0.618
        if Gf(p, m1) < Gf(p, m2): hi = m2
        else: lo = m1
    gmin = float(Gf(p, (lo + hi)/2))
    mp_ = '%.7e' % (gmin*(1 - 2e-6))               # 8 significant digits, below the minimum; the decimal itself is certified
    ok, n = certify(p, mp_)
    print(p, 'min %.10e at x=%.7f' % (gmin, (lo + hi)/2), 'm_p', mp_, 'certified', ok, n, flush=True)
    assert ok
    out[p] = mp_
json.dump(out, open('zfr_mp.json', 'w'), indent=0)
s = sum(math.log(p)*float(v) for p, v in out.items())
print('sum log p m_p (p < 60) =', s)
