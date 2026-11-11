# Farr-Pauli-Saidak conjecture (Funct. Approx. 64 (2021), restated in arXiv 2301.11209 l. 523-525):
#   |C_alpha(a)| <= 2 |exp(alpha (log w - 1/w))|,  w = W_0(i alpha/(2 pi)),  C_alpha(a) = gamma_alpha(a) - log^alpha(a)/a = gamma_alpha(1+a).
# Integer alpha = n, dyadic a in (0,1]; gamma_n(1+a) in Arb (python-flint), no cancellation.
import sys, mpmath
from flint import acb, arb, ctx
ctx.prec = 192; mpmath.mp.dps = 50
def bound(n):
    w = mpmath.lambertw(1j*n/(2*mpmath.pi), 0)
    return 2*abs(mpmath.exp(n*(mpmath.log(w) - 1/w)))
nmax = int(sys.argv[1]); J = int(sys.argv[2])
ns = list(range(1, 301)) + list(range(301, nmax + 1, 7))
B = {n: bound(n) for n in ns}
best = (0, None, None)
for j in range(1, J + 1):
    b = arb(1) + arb(j)/J
    wr = (0, None)
    for n in ns:
        g = acb.stieltjes(n, acb(b))
        assert abs(g).rad() < abs(g).mid()*arb('1e-6') or abs(g).mid() < 1e-30, (n, j)
        r = mpmath.mpf(abs(g).mid().str(25, radius=False)) / B[n]
        if r > wr[0]: wr = (r, n)
    if wr[0] > best[0]: best = (wr[0], wr[1], j)
    print('a=%d/%d: max ratio %s at n=%d' % (j, J, mpmath.nstr(wr[0], 6), wr[1]), flush=True)
print('OVERALL max |C_n(a)|/FPS bound = %s at n=%d, a=%d/%d' % (mpmath.nstr(best[0], 8), best[1], best[2], J))
