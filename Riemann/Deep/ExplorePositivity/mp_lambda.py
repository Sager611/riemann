# mp_lambda.py — the whole Meixner–Pollaczek family P_n^{(lam)}(x; pi/2):  (n+1)P_{n+1} = 2x P_n - (n+2lam-1) P_{n-1}.
# (a) Polya-type lemma  sum a_k z^k real-rooted => sum k!/(2lam)_k a_k P_k real-rooted  (random tests; Romik states lam = 3/4);
# (b) Xi(t) = sum (-1)^n c_2n^(lam) P_2n^(lam)(t/2),  c_2n^(lam) ∝ int_0^inf omega(e^{2u}) e^{u/2} cosh(u)^(-2lam) tanh(u)^(2n) du;
# (c) Newton k=1 obstruction: needs c2^2/(c0 c4) >= [(2l+2)(2l+3)/(6*2l*(2l+1))] * 2N/(N-1); Cauchy-Schwarz gives <= 1;
# (d) pre-image partial sums p_N(z) = sum_{n<=N} (-1)^n c_2n (2lam)_{2n}/(2n)! z^{2n}: Arb root isolation, N <= NMAX.
import sys, random, mpmath as mp
from mpmath.calculus.quadrature import GaussLegendre
from flint import arb, acb, acb_poly, ctx
from mp_core import omega
mp.mp.dps = 50; ctx.prec = 600
lams = [mp.mpf(v) for v in sys.argv[1].split(',')]; NMAX = int(sys.argv[2])
def Pcoeffs(lam, nmax):
    P = [[mp.mpf(1)], [mp.mpf(0), mp.mpf(2)]]
    for n in range(1, nmax):
        c = [mp.mpf(0)] * (n + 2)
        for k, v in enumerate(P[n]): c[k+1] += 2*v
        for k, v in enumerate(P[n-1]): c[k] -= (n + 2*lam - 1)*v
        P.append([v/(n+1) for v in c])
    return P
EDGES = [0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1, 1.125, 1.25, 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3.0, 3.3]
base = GaussLegendre(mp.mp).calc_nodes(6, mp.mp.prec)
def cl(lam):
    acc = [mp.mpf(0)] * (NMAX + 1)
    for a, b in zip(EDGES[:-1], EDGES[1:]):
        a, b = mp.mpf(a), mp.mpf(b)
        for x, w in base:
            u = (b - a)/2*x + (a + b)/2
            p = omega(mp.exp(2*u)) * mp.exp(u/2) * mp.cosh(u)**(-2*lam) * w*(b - a)/2
            th2 = mp.tanh(u)**2
            for n in range(NMAX + 1):
                acc[n] += p; p *= th2
    return acc
random.seed(7)
for lam in lams:
    P = Pcoeffs(lam, 10)
    worst = 0
    for trial in range(100):
        d = random.randint(2, 8); p = [mp.mpf(1)]
        for r in [random.uniform(-3, 3) for _ in range(d)]:
            q = [mp.mpf(0)]*(len(p)+1)
            for k, v in enumerate(p): q[k+1] += v; q[k] -= r*v
            p = q
        T = [mp.mpf(0)]*(d+1)
        for k, a in enumerate(p):
            wk = mp.factorial(k)/mp.rf(2*lam, k)*a
            for j, c in enumerate(P[k]): T[j] += wk*c
        worst = max(worst, max(abs(mp.im(z)) for z in mp.polyroots(T[::-1], maxsteps=300, extraprec=300)))
    c = cl(lam)
    need = (2*lam+2)*(2*lam+3)/(6*2*lam*(2*lam+1))
    ratio = c[1]**2/(c[0]*c[2])
    fails = []
    for N in range(2, NMAX + 1):
        co = [mp.mpf(0)]*(2*N+1)
        for n in range(N+1): co[2*n] = (-1)**n * c[n] * mp.rf(2*lam, 2*n)/mp.factorial(2*n)
        rts = acb_poly([acb(arb(mp.nstr(v, 45))) for v in co]).roots(maxprec=2400)
        nr = sum(1 for r in rts if not r.imag.contains(0))
        fails.append(nr)
    print('lam=%s lemma: max|Im| over 100 random p = %s | c2^2/(c0c4) = %s, Newton needs >= %s*2N/(N-1) (N->inf: %s) | nonreal zeros of p_N, N=2..%d: %s' % (
        mp.nstr(lam, 4), mp.nstr(worst, 2), mp.nstr(ratio, 5), mp.nstr(need, 5), mp.nstr(2*need, 5), NMAX, fails), flush=True)
