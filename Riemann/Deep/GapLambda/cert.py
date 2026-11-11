"""Certificate: lower bounds for Lambda = limsup (t_{n+1}-t_n)/(2 pi/log t_n) over consecutive critical zeros
from Bui-Hall-Subira Jorge (arXiv 2511.14415): Thm 1.1 (main term c, l.75-100), Thm Wi (l.507-516),
keyineq2 (l.665-668): Lambda > kappa  if  3 lam0(v) c(0,0) > 16 kappa^4 c1(u) + 24 v kappa^2 c2(u).
All main-term integrals are evaluated from the exact rational polynomial integrand (core.py); the large linear
algebra runs in Arb ball arithmetic (outward rounding, prec 320 bits) from exact rational inputs.
The optimiser that proposed P, u, v is not trusted: only the inputs below enter.
Run: uv run --no-project --with python-flint python3 cert.py"""
import sys, time, resource
sys.path.insert(0, '.')
from flint import fmpq, arb, arb_mat, ctx
import core

ctx.prec = 320
A = lambda q: arb(q.p) / arb(q.q)


def shleg_mono(n):
    from math import comb
    return [fmpq((-1) ** (n + k) * comb(n, k) * comb(n + k, k)) for k in range(n + 1)]


def main_terms(theta, r, P):
    """Arb balls: c (scalar), c1 and c2 as lists of u-polynomial coefficients."""
    parts = core.poly_parts(theta)
    terms = {k: core.split_terms(p) for k, p in parts.items()}
    gAs = sorted({k[0] for t in terms.values() for k in t})
    gBs = sorted({k[1] for t in terms.values() for k in t})
    VA = {g: core.V(g, r, P) for g in gAs}
    VB = {g: core.V(g, r, P) for g in gBs}
    colsA = sorted({n for d in VA.values() for n in d})
    colsB = sorted({n for d in VB.values() for n in d})
    nmin = min(min(min(c) for c in colsA), min(min(c) for c in colsB))
    nmax = max(max(max(c) for c in colsA), max(max(c) for c in colsB))
    assert nmin >= 2 * r
    Gt = core.G_exact_table(r, 2 * r, nmax - 2 * r + 1)            # exact rationals, then exact -> ball
    iA = {n: i for i, n in enumerate(colsA)}; iB = {n: i for i, n in enumerate(colsB)}
    MA = arb_mat(len(gAs), len(colsA)); MB = arb_mat(len(colsB), len(gBs)); GM = arb_mat(len(colsA), len(colsB))
    for i, g in enumerate(gAs):
        for n, cf in VA[g].items(): MA[i, iA[n]] = A(cf)
    for j, g in enumerate(gBs):
        for n, cf in VB[g].items(): MB[iB[n], j] = A(cf)
    for i, (n1, n2) in enumerate(colsA):
        for j, (n3, n4) in enumerate(colsB):
            GM[i, j] = A(Gt[(n1, n2, n3, n4)])
    del Gt
    M = MA * GM * MB
    rA = {g: i for i, g in enumerate(gAs)}; rB = {g: j for j, g in enumerate(gBs)}
    out = {}
    for name, tm in terms.items():
        deg = max(k for (_, _, k) in tm)
        cs = [arb(0)] * (deg + 1)
        for (gA, gB, k), cf in tm.items():
            cs[k] += A(cf) * M[rA[gA], rB[gB]]
        out[name] = cs
    return out


def bound(theta, r, P, u, v, kappa=None, label=''):
    t0 = time.time()
    mt = main_terms(theta, r, P)
    ua, va = A(u), A(v)
    c = mt['c'][0]
    c1 = sum((cf * ua ** k for k, cf in enumerate(mt['c1'])), arb(0))
    c2 = sum((cf * ua ** k for k, cf in enumerate(mt['c2'])), arb(0))
    lam0 = (1 + 4 * va + (1 + 8 * va).sqrt()) / 8                    # Thm Wi, l.514
    k2 = ((144 * va ** 2 * c2 ** 2 + 48 * lam0 * c * c1).sqrt() - 12 * va * c2) / (16 * c1)   # keyineq3, l.674-676
    print('[%s] theta=%s r=%d deg P=%d u=%s v=%s' % (label, theta, r, len(P) - 1, u, v))
    print('   c=%s  c1(u)=%s  c2(u)=%s' % (c.str(10), c1.str(10), c2.str(10)))
    print('   kappa_max^2 in %s ; Lambda >= kappa_max in %s' % (k2.str(15), k2.sqrt().str(15)))
    assert c > 0 and c1 > 0
    if kappa is not None:
        ka = A(kappa)
        H = 3 * lam0 * c - 16 * ka ** 4 * c1 - 24 * va * ka ** 2 * c2       # keyineq2, l.666-668
        ok = H > 0                                                          # Arb: true only if the whole ball is > 0
        print('   keyineq2 at kappa=%s: H in %s  ->  %s' % (kappa, H.str(8), 'CERTIFIED Lambda > %s' % kappa if ok else 'NOT certified'))
        assert ok
    print('   (%.1fs, peak RSS %.0f MB)' % (time.time() - t0, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20), flush=True)
    return k2


if __name__ == '__main__':
    # 1) Hall's value: theta=0, r=1, P=1, u=1/2, v=22/49 must give kappa^2 = 11/2 (TeX l.678-682)
    k2 = bound(fmpq(0), 1, [fmpq(1)], fmpq(1, 2), fmpq(22, 49), label='Hall check')
    assert abs(k2 - arb(11) / 2) < arb(10) ** -60
    # 2) BHSJ's value: theta=0.1249, r=1, P=1; their (u,v)=(0.6,0.4) (l.684) and our (u,v)
    bound(fmpq(1249, 10000), 1, [fmpq(1)], fmpq(3, 5), fmpq(2, 5), label='BHSJ params u=.6 v=.4')
    bound(fmpq(1249, 10000), 1, [fmpq(1)], fmpq(5942, 10000), fmpq(4038, 10000), kappa=fmpq(264, 100), label='P=1, tuned u,v')
    # 3) optimised P (degree 11), given by exact shifted-Legendre coefficients (4 decimals), converted exactly
    pl = ['1', '3.2279', '4.6423', '9.7348', '5.6966', '21.5327', '1.8217', '36.1817', '-4.141', '42.4179', '-5.2206', '28.2819']
    def Q(s):
        neg = s.startswith('-'); s = s.lstrip('-'); a, b = (s.split('.') + [''])[:2]
        q = fmpq(int(a + b), 10 ** len(b)); return -q if neg else q
    P = [fmpq(0)] * len(pl)
    for m, s in enumerate(pl):
        for k, cf in enumerate(shleg_mono(m)):
            P[k] += Q(s) * cf
    print('optimised P (monomial, exact):', [str(x) for x in P])
    bound(fmpq(1249, 10000), 1, P, fmpq(6744, 10000), fmpq(40741, 100000), kappa=fmpq(2907, 1000), label='optimised P')
