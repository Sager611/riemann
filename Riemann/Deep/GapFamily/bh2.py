"""gapfamily: INDEPENDENT second implementation of Bui-Hall's h_k (arXiv 2304.05181v1), literally from the TeX:
general formula l.418-433 (operators Q^4 R^{2j} applied to c(alpha,beta), l.322-333), evaluated at the point of l.434
(alpha1 = beta1 = 0, alpha2 = -beta2 = pi i lam/L, alpha3 = -beta3 = 2 pi i eta/L); Q(x) = (1/2+x)^k (l.371), R(x) = v+1+x (l.396);
h_k = c_{k,0}/(lam^2 c_{k,1}) (l.441-444).  Original 9 variables x, x1..x4, t1..t4 (no u-substitution, no closed forms):
simplex map x = s0, x1 = (1-x)s1, x2 = (1-x-x1)s2, x3 = (1-x)s3, x4 = (1-x-x3)s4; tensor Gauss-Legendre in all 9 variables.
mode 'general' = exponent of l.418-433 at that point (with eta); mode 'literal' = exponent printed in Prop main1, l.345 (no eta).
usage: python3 bh2.py k lam eta mode n_x,n_t12,n_t34  -> forms M[n][i][l] (n <= 2, basis P = b0 + b1*w) and h at BH's / best (v, P)."""
import sys, math, json
import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.optimize import minimize_scalar
import scipy.linalg as sl

k = int(sys.argv[1]); lam = float(sys.argv[2]); eta = float(sys.argv[3]); mode = sys.argv[4]
nx, n12, n34 = (int(z) for z in sys.argv[5].split(','))
th = 0.2499


def gl(n):
    x, w = leggauss(n)
    return (x + 1) / 2, w / 2


sx, wx = gl(nx); s12, w12 = gl(n12); s34, w34 = gl(n34)
# t-grid (t1, t2, t3, t4)
T = np.stack(np.meshgrid(s12, s12, s34, s34, indexing='ij'), -1).reshape(-1, 4)
WT = np.einsum('a,b,c,d->abcd', w12, w12, w34, w34).ravel()
t1, t2, t3, t4 = (T[:, i][None, :] for i in range(4))
O = np.stack(np.meshgrid(*[np.arange(nx)] * 5, indexing='ij'), -1).reshape(-1, 5)
v0 = 1.2   # centre of the S-moments (BH: v = 0.2, R = v+1)
M = np.zeros((3, 2, 2), dtype=complex)
pl = math.pi * lam
for st in range(0, len(O), 40):
    o = O[st:st + 40]
    x = sx[o[:, 0]]; W = wx[o[:, 0]]
    x1 = (1 - x) * sx[o[:, 1]]; W = W * (1 - x) * wx[o[:, 1]]
    x2 = (1 - x - x1) * sx[o[:, 2]]; W = W * (1 - x - x1) * wx[o[:, 2]]
    x3 = (1 - x) * sx[o[:, 3]]; W = W * (1 - x) * wx[o[:, 3]]
    x4 = (1 - x - x3) * sx[o[:, 4]]; W = W * (1 - x - x3) * wx[o[:, 4]]
    c = lambda z: z[:, None]
    x, x1, x2, x3, x4, W = map(c, (x, x1, x2, x3, x4, W))
    a = 1 - th * (x1 + x3); b = 1 - th * (x2 + x4)
    if mode == 'literal':      # l.345 as printed
        E = th * pl * (x2 - x4 - (x3 - x4) * t3 + (x1 - x2) * t4) - pl * (a * t1 - b * t2) * (t3 - t4)
    else:                      # l.418-433 at the point of l.434; alpha's in units of 1/log T, log y = th log T
        A2, B2, A3, B3 = 1j * pl, -1j * pl, 2j * math.pi * eta, -2j * math.pi * eta
        Ey = th * (-(A3 + B3) * x - A3 * (x1 + x2) - B3 * (x3 + x4) - B2 * x2 - A2 * x4 - A2 * (x3 - x4) * t3 - B2 * (x1 - x2) * t4)
        ET = a * (-A2 * t1 * t3 - B2 * t1 * t4) + b * (-(A2 + B2) * t2 + A2 * t2 * t3 + B2 * t2 * t4)
        E = (Ey + ET) / 1j
    f = np.exp(1j * E) * a * b
    f = f * (th * (x1 - x2) + a * t1 - b * t2) * (th * (x3 - x4) + a * t1 - b * t2)
    q = lambda z: (0.5 + z) ** k   # Q(x) = (1/2 + x)^k, arguments as printed l.427-430
    f = f * q(-th * x3 + th * (x3 - x4) * t3 - a * t1 * (1 - t3) - b * t2 * t3)
    f = f * q(-th * x4 - th * (x3 - x4) * t3 - a * t1 * t3 - b * t2 * (1 - t3))
    f = f * q(-th * x1 + th * (x1 - x2) * t4 - a * t1 * (1 - t4) - b * t2 * t4)
    f = f * q(-th * x2 - th * (x1 - x2) * t4 - a * t1 * t4 - b * t2 * (1 - t4))
    f = f * W * WT[None, :]
    S = th * (x + x1 + x2 + x3 + x4) + a * t1 + b * t2      # R(-S) = v + 1 - S  (l.431)
    w1 = 1 - x - x1 - x2; w2 = 1 - x - x3 - x4               # P(1-x-x1-x2) P(1-x-x3-x4)
    D = S - v0
    for n in range(3):
        g = f * D ** n
        for i in range(2):
            for l in range(2):
                M[n, i, l] += np.sum(g * w1 ** i * w2 ** l)
Ms = 0.5 * (M + np.transpose(M, (0, 2, 1)))
imag = float(np.max(np.abs(Ms.imag)) / np.max(np.abs(Ms.real))); Mr = Ms.real


def h(V, bb=None):
    tt = V - v0
    C2 = Mr[2] - 2 * tt * Mr[1] + tt * tt * Mr[0]
    if bb is not None:
        return (bb @ Mr[0] @ bb) / (lam ** 2 * (bb @ C2 @ bb))
    return sl.eigh(Mr[0], lam ** 2 * C2, eigvals_only=True)[-1]


aBH = {1: 2.5, 2: 2.0}.get(k)
out = {'k': k, 'lam': lam, 'eta': eta, 'mode': mode, 'nodes': [nx, n12, n34], 'imag': imag}
if aBH is not None:
    out['h_BH'] = h(1.2, np.array([1.0, -aBH]))          # BH's printed v = 0.2 (R = v+1), P = 1 - a x
    out['h_BH_Rv'] = h(0.2, np.array([1.0, -aBH]))       # alternative reading R = v + x
Vs = np.linspace(0.0, 2.0, 81); hs = [h(V) for V in Vs]; i0 = int(np.argmax(hs))
r = minimize_scalar(lambda V: -h(V), bounds=(Vs[max(i0 - 1, 0)], Vs[min(i0 + 1, 80)]), method='bounded', options={'xatol': 1e-8})
out['h_best_linearP_v'] = -r.fun; out['v_best'] = r.x - 1
print(json.dumps(out))
