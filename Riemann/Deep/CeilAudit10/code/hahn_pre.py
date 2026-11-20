# hahn_pre.py — ceilaudit10 T2c: Hahn pre-image of Xi under the (empirical) lemma z^k -> kappa_k g_k: a_2n = (-1)^n d_2n / kappa_2n.
# Newton k=1 ratio, log-convexity of d_2n, and Arb real-rootedness of the truncated pre-images q_N(w) = sum_{n<=N} |a_2n| (-w)^n.
import mpmath as mp
from flint import arb, acb, acb_poly, ctx
mp.mp.dps = 60; ctx.prec = 700
d = {}
for line in open('Riemann/Deep/ExplorePositivity/logs/d2n_zeta.txt'):
    k, v = line.split()[:2]; d[int(k)] = mp.mpf(v)
M = max(d)//2
print('d_2n available n <=', M, '; d0,d2,d4 =', [mp.nstr(d[2*i], 12) for i in range(3)])
lc = [d[2*n]**2/(d[2*n-2]*d[2*n+2]) for n in range(1, M)]
print('d_2n^2/(d_2n-2 d_2n+2), n=1..8:', [mp.nstr(x, 6) for x in lc[:8]], '; n=20,50,90:', [mp.nstr(lc[i-1], 6) for i in (20, 50, 90) if i < M])
print('log-convex (ratio <= 1) for all n?', all(x <= 1 for x in lc), '; max ratio', mp.nstr(max(lc), 6), 'at n =', 1+lc.index(max(lc)))
kaps = {'romik (kappa=1)': lambda k: mp.mpf(1), 'k!/(3/2)_k': lambda k: mp.factorial(k)/mp.rf(1.5, k),
        '4^k k!/(3/2)_k^2': lambda k: 4**k*mp.factorial(k)/mp.rf(1.5, k)**2, '1/(3/2)_k': lambda k: 1/mp.rf(1.5, k)}
for name, kf in kaps.items():
    a = [d[2*n]/kf(2*n) for n in range(M+1)]
    R = a[1]**2/(a[0]*a[2])
    res = []
    for N in (2, 3, 4, 6, 10, 20, 40):
        co = [a[n]*(-1)**n for n in range(N+1)]
        rts = acb_poly([acb(arb(mp.nstr(v, 55))) for v in co]).roots(maxprec=4000)
        nr = sum(1 for r in rts if not r.imag.contains(0)); npos = sum(1 for r in rts if r.imag.contains(0) and r.real > 0)
        res.append('N=%d:%d nonreal' % (N, nr))
    print('%-18s Newton k=1: a1^2/(a0 a2) = %s (needs >= 2N/(N-1) -> 2) | %s' % (name, mp.nstr(R, 6), ', '.join(res)), flush=True)
