
# KLN Lemma 3.2 (3.5): max_{|t|<=T}|zeta(1/2+it)| <= a1 T^{1/6} log T + a2 for all T > 0, a1 = 0.63, a2 = 2.851.
# For t >= 3 it follows from (3.4) (HPY 0.618 <= 0.63). On [0,3] we check pointwise |zeta(1/2+it)| <= R(t) := min_{T>=t} (a1 T^{1/6} log T + a2)
# = a1 max(t,e^-6)^{1/6} log max(t,e^-6) + a2 (the RHS is increasing for T >= e^-6), by adaptive ball bisection (Arb, 128 bits).
import time
from flint import arb, acb, ctx
ctx.prec = 128
a1, a2 = arb('0.63'), arb('2.851')
e6 = arb(-6).exp()
def Rlow(tl):   # lower bound of R on [tl, th] = R(tl) (R non-decreasing)
    u = tl.max(e6) if tl > e6 else e6
    return a1*u**(arb(1)/6)*u.log() + a2
t0 = time.time()
stack = [(arb(0), arb(3))]; nb = 0; worst = arb(-1); minw = arb(3)
while stack:
    a, b = stack.pop()
    tb = a.union(b)
    v = arb(abs(acb(arb('0.5'), tb).zeta()).upper())
    r = arb(Rlow(a).lower())
    if v < r:
        nb += 1
        if (r - v) < (worst if worst > 0 else arb(10)): worst = r - v
        continue
    w = b - a
    assert w > arb('1e-12'), (a, b, v, r)
    if w < minw: minw = w
    m = (a + b)/2
    stack.append((a, m)); stack.append((m, b))
print('KLN (3.5) with a1=0.63, a2=2.851 on [0,3]: verified on', nb, 'balls; smallest margin', worst.str(5), '; finest width', minw.str(3), ' %.1fs' % (time.time()-t0))
