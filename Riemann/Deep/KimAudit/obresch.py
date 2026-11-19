#!/usr/bin/env python3
"""kimaudit T8: random test of the cited input (I3) = Kim-Lee Thm 3 (Obreschkoff): Z(P) in S(delta) = {|Im z| <= delta|z|},
Q hyperbolic, deg Q <= delta^-2  =>  P(D)Q hyperbolic. Control: deg Q = floor(delta^-2) + 1 (may fail).
Run: bash .scratch/memcap.sh 290 600 uv run --no-project --with mpmath python3 obresch.py"""
import random
import mpmath as mp
mp.mp.dps = 50
random.seed(20261119)
def pdq(P, Q):  # P, Q coefficient lists, highest degree first
    Pc = P[::-1]; out = [mp.mpf(0)]*len(Q); Qk = Q[:]
    for k, pk in enumerate(Pc):
        if len(Qk) == 0: break
        f = pk          # P^{(k)}(0)/k! = k-th coefficient of P
        for i in range(len(Qk)): out[len(Q)-len(Qk)+i] += f*Qk[i]
        Qk = [Qk[i]*(len(Qk)-1-i) for i in range(len(Qk)-1)]
    while len(out) > 1 and out[0] == 0: out = out[1:]
    return out
def hyperbolic(c):
    if len(c) <= 2: return True
    r = mp.polyroots(c, maxsteps=500, extraprec=200)
    return all(abs(mp.im(x)) <= mp.mpf(10)**-12*(1 + abs(x)) for x in r)
res = {'in': [0, 0], 'out': [0, 0]}
for trial in range(800):
    delta = mp.mpf(random.uniform(0.3, 0.75)); N = int(mp.floor(delta**-2))
    roots = []
    for _ in range(random.randint(1, 3)):
        th = mp.asin(delta*mp.mpf(random.choice([1, 1, random.random()])))   # boundary-heavy
        r = mp.mpf(random.uniform(0.2, 3)); sgn = random.choice([1, -1])
        z = sgn*r*mp.expj(th); roots += [z, mp.conj(z)]
    for _ in range(random.randint(0, 2)): roots.append(mp.mpf(random.uniform(-3, 3)))
    co = [mp.mpc(1)]
    for z in roots:
        co = [co[0]] + [co[i] - z*co[i-1] for i in range(1, len(co))] + [-z*co[-1]]
    P = [mp.re(x) for x in co]
    for key, dq in (('in', N), ('out', N+1)):
        qr = [mp.mpf(random.uniform(-2, 2)) for _ in range(dq)]
        Q = [mp.mpf(1)]
        for z in qr: Q = [Q[0]] + [Q[i] - z*Q[i-1] for i in range(1, len(Q))] + [-z*Q[-1]]
        h = hyperbolic(pdq(P, Q)); res[key][0] += 1; res[key][1] += (0 if h else 1)
print('deg Q <= delta^-2 (theorem): %d trials, %d non-hyperbolic' % tuple(res['in']))
print('deg Q = floor(delta^-2)+1 (control, outside): %d trials, %d non-hyperbolic' % tuple(res['out']))
P = [mp.mpf(1), -2*mp.cos(mp.mpf(1)), mp.mpf(1)]; d = int(mp.floor(1/mp.sin(1)**2))
print('sharpness P = z^2 - 2z cos 1 + 1, Q = z^d: d = %d hyperbolic %s, d = %d hyperbolic %s' % (d, hyperbolic(pdq(P, [mp.mpf(1)] + [mp.mpf(0)]*d)), d+1, hyperbolic(pdq(P, [mp.mpf(1)] + [mp.mpf(0)]*(d+1)))))
