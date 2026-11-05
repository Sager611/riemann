# t2_twist.py -- explore3 transfer 2 (N7): Euclidean (Baladi-Vallee) Dirichlet series with arithmetic twists.
# Identities checked by direct summation over rationals u/v in (0,1], gcd(u,v)=1, v <= V:
#  (A) sum v^{-2s}                     = zeta(2s-1)/zeta(2s)                   (poles: s=1; s=rho/2, Re<=1/2)
#  (B) sum chi4(v) v^{-2s}             = L(2s-1,chi4)/L(2s,chi4)              (no pole at s=1)
#  (C) sum chi4(u) v^{-2s}             = D(2s)/L(2s,chi4), D(w)=sum_{m=1,2 mod 4} m^{-w} = 4^{-w}(zeta(w,1/4)+zeta(w,1/2))
# In every case the L-function attached to the arithmetic twist sits at the doubled argument 2s, so its zeros give
# singularities at Re s <= 1/2: half a unit behind the leading pole s = 1 of the untwisted operator.
import sys, math, time
import numpy as np
import mpmath as mp
mp.mp.dps = 30

def chi4(n):
    r = n % 4
    return 0 if r % 2 == 0 else (1 if r == 1 else -1)

def direct(s, V):
    # phi-type sums via a sieve: for each v, S0 = #{u<=v coprime}, S1 = sum chi4(u) over u<=v coprime
    A = B = C = mp.mpf(0)
    # precompute mu
    mu = np.ones(V + 1, dtype=np.int64); isp = np.ones(V + 1, dtype=bool); isp[:2] = False
    for p in range(2, V + 1):
        if isp[p]:
            isp[2*p::p] = False
            mu[p::p] *= -1
            mu[p*p::p*p] = 0
    S = [0, 1, 1, 0]  # partial sums of chi4 by residue of m mod 4 : S(m)=1 if m%4 in (1,2)
    tot0 = np.zeros(V + 1); tot1 = np.zeros(V + 1)
    for d in range(1, V + 1):
        if mu[d] == 0: continue
        cd = chi4(d)
        for v in range(d, V + 1, d):
            m = v // d
            tot0[v] += mu[d] * m
            if cd: tot1[v] += mu[d] * cd * (1 if m % 4 in (1, 2) else 0)
    s = mp.mpc(s)
    for v in range(1, V + 1):
        w = mp.power(v, -2*s)
        A += tot0[v] * w
        B += chi4(v) * tot0[v] * w
        C += tot1[v] * w
    return A, B, C

def L4(w):
    return mp.power(4, -w) * (mp.zeta(w, mp.mpf(1)/4) - mp.zeta(w, mp.mpf(3)/4))

def closed(s):
    s = mp.mpc(s)
    A = mp.zeta(2*s - 1) / mp.zeta(2*s)
    B = L4(2*s - 1) / L4(2*s)
    D = mp.power(4, -2*s) * (mp.zeta(2*s, mp.mpf(1)/4) + mp.zeta(2*s, mp.mpf(1)/2))
    C = D / L4(2*s)
    return A, B, C

if __name__ == '__main__':
    t0 = time.time()
    for s in (complex(1.6, 0.0), complex(1.3, 2.0)):
        for V in (2000, 8000):
            d = direct(s, V); c = closed(s)
            print('s=%s V=%d' % (s, V), ' '.join('%s: direct %s closed %s' % (n, mp.nstr(x, 8), mp.nstr(y, 8)) for n, x, y in zip('ABC', d, c)))
    # singularities of (C): zeros of L(2s,chi4) -> s = (1/2 + i*6.0209489)/2 etc; check |closed C| blows up there
    g1 = mp.findroot(lambda t: L4(mp.mpf(0.5) + 1j*t), 6.02)
    print('first zero of L(.,chi4) at 1/2+%si -> pole of (C) at s=%s' % (mp.nstr(g1, 10), mp.nstr(mp.mpc(0.25, g1/2), 8)))
    for eps in (1e-2, 1e-4):
        print('  |C(s)| at distance %.0e: %s' % (eps, mp.nstr(abs(closed(mp.mpc(0.25, g1/2 + eps))[2]), 6)))
    print('(C) at s=1+1e-9 (no pole):', mp.nstr(closed(mp.mpf(1) + mp.mpf(1e-9))[2], 10), '; |(A)| at s=1+1e-6:', mp.nstr(abs(closed(mp.mpf(1) + mp.mpf(1e-6))[0]), 6))
    print('time %.1fs' % (time.time() - t0))
