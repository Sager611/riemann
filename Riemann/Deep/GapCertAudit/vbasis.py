"""gapcertaudit: independent check of the per-gap inequality (V) by a basis different from orch_vnum.py's x(g-x)P_m:
span{ sin(k pi x/g), k <= K } + { x(g-x) P_m(2x/g - 1), m < MP } (MP low-degree polynomials absorb the boundary jets
f''(0), f''(g), f''''(0), ...; with MP = 2 (two cubics) the sine tail converges only like 1/K).
Gram matrices of f^(j), j <= 4, and of the symmetrised brackets [f'f'']_0^g, [f'''f'']_0^g by Gauss-Legendre (192 nodes,
exact for polynomial degree 383) at 50 digits; lambda_min = min Q(f)/|f|^2 on the span (an upper bound for the true inf,
decreasing in K).  Forms: 'cert' (cert4.py's rationals), 'noshift' (a0 - 2e-5), 'sdp' (a0 - 2e-5, phi12 - 2e-6).
usage: vbasis.py g_expr K1,K2,... [MP]"""
import sys
from fractions import Fraction as Fr
import mpmath as mp
mp.mp.dps = 50

A = [Fr('-5899179108889603/6134092019850000'), Fr(1), Fr('14789027537/948575124080'), Fr('-4037122085/176350902122'),
     Fr('1982999731/771179571628')]
PH12 = Fr('3691572375143311/104390384155500000'); PH32 = Fr('4037122085/176350902122')
assert PH32 == -A[3]
FORMS = {'cert': (A, PH12), 'noshift': ([A[0] - Fr(2, 100000)] + A[1:], PH12),
         'sdp': ([A[0] - Fr(2, 100000)] + A[1:], PH12 - Fr(2, 10 ** 6))}


def polys(g, MP):
    """coefficient lists (powers of x) of x(g-x)P_m(2x/g-1), m < MP, and their derivatives 0..4."""
    t = [mp.mpf(-1), 2 / g]                         # 2x/g - 1
    def mul(a, b):
        out = [mp.mpf(0)] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                out[i + j] += x * y
        return out
    P = [[mp.mpf(1)], t]
    while len(P) < MP:
        m = len(P) - 1                              # (m+1)P_{m+1} = (2m+1) t P_m - m P_{m-1}
        a = mul(t, P[m]); b = P[m - 1] + [mp.mpf(0)] * (len(a) - len(P[m - 1]))
        P.append([((2 * m + 1) * x - m * y) / (m + 1) for x, y in zip(a, b)])
    out = []
    for m in range(MP):
        c = mul([mp.mpf(0), g, mp.mpf(-1)], P[m])
        ds = [c]
        for _ in range(4):
            c = [i * c[i] for i in range(1, len(c))] or [mp.mpf(0)]
            ds.append(c)
        out.append(ds)
    return out


def peval(c, x):
    v = mp.mpf(0)
    for a in reversed(c):
        v = v * x + a
    return v


def run(g, K, MP=2):
    gl = mp.calculus.quadrature.GaussLegendre(mp.mp)
    nodes = gl.calc_nodes(7, mp.mp.prec)            # 3*2^6 = 192 nodes on [-1, 1]
    X = [(g / 2) * (1 + t) for t, w in nodes]; Wt = [(g / 2) * w for t, w in nodes]
    nb = K + MP
    PL = polys(g, MP)
    def ders(b, x):                                 # values of the b-th basis function's derivatives 0..4 at x
        if b < K:
            om = (b + 1) * mp.pi / g; s, c = mp.sin(om * x), mp.cos(om * x)
            return [s, om * c, -om ** 2 * s, -om ** 3 * c, om ** 4 * s]
        return [peval(PL[b - K][j], x) for j in range(5)]
    V = [[ders(b, x) for x in X] for b in range(nb)]
    G = [[[mp.fsum(Wt[q] * V[a][q][j] * V[b][q][j] for q in range(len(X))) for b in range(nb)] for a in range(nb)] for j in range(5)]
    E0 = [ders(b, mp.mpf(0)) for b in range(nb)]; Eg = [ders(b, g) for b in range(nb)]
    def br(p, q):                                   # symmetrised [f^(p) f^(q)]_0^g as a bilinear form
        return [[(Eg[a][p] * Eg[b][q] + Eg[b][p] * Eg[a][q] - E0[a][p] * E0[b][q] - E0[b][p] * E0[a][q]) / 2
                 for b in range(nb)] for a in range(nb)]
    B12, B32 = br(1, 2), br(3, 2)
    Lc = mp.cholesky(mp.matrix(G[0])); Li = mp.inverse(Lc)
    res = {}
    for name, (a, p12) in FORMS.items():
        Qm = mp.matrix(nb, nb)
        for i in range(nb):
            for k in range(nb):
                Qm[i, k] = sum(mp.mpf(a[j].numerator) / a[j].denominator * G[j][i][k] for j in range(5)) \
                    + mp.mpf(p12.numerator) / p12.denominator * B12[i][k] + mp.mpf(PH32.numerator) / PH32.denominator * B32[i][k]
        M = Li * Qm * Li.T
        M = (M + M.T) / 2
        ev = mp.eigsy(M, eigvals_only=True)
        res[name] = min(ev)
    return res


if __name__ == '__main__':
    g = mp.mpf(eval(sys.argv[1], {'pi': mp.pi}))
    MP = int(sys.argv[3]) if len(sys.argv) > 3 else 2
    for K in [int(s) for s in sys.argv[2].split(',')]:
        r = run(g, K, MP)
        print('g = %s K = %d MP = %d: ' % (mp.nstr(g, 8), K, MP) + ', '.join('%s %s' % (k, mp.nstr(v, 10)) for k, v in r.items()), flush=True)
