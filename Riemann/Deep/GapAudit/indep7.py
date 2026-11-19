"""gapaudit: independent float evaluation of Bui-Milinovich's c_0, c_1 (arXiv 1410.3635v3, Lemma main l.370-384 with the
Q-operator of l.385-434) for M(s + i*mu*L^{-1}) -- mollifier shifts alpha3 = i*mu/L, beta3 = -i*mu/L, so the extra factor is
y^{-alpha3(x1+x2)-beta3(x3+x4)} = exp(-i*mu*theta*(x1+x2-x3-x4)).  mu = 0: BM's own kernel; mu = kappa*pi/2: the midpoint.
t3, t4 integrated exactly (own derivation): W*(t3,t4-integral) = (1-e^{-ia})(1-e^{ib})/(kappa pi)^2,
a = kappa pi (theta(x3-x4) + D), b = kappa pi (theta(x1-x2) + D), D = (1-theta(x1+x3))t1 - (1-theta(x2+x4))t2.
Then tau_i = A_i t_i (A1 = 1-theta(x1+x3), A2 = 1-theta(x2+x4)) absorbs BM's factor A1*A2; tau moments by Gauss.
Written from BM's TeX only (no GapBar/GapCeil code).  usage: indep7.py r theta v kappa mu_over_c P_spec [n]
P_spec: 'mono:b0,b1,...' (P(x) = sum b_j x^j) or a json file with shifted-Legendre A_n (P(w) = sum (-1)^n A_n P_n(2w-1))."""
import sys, json, math
import numpy as np
from numpy.polynomial.legendre import leggauss, legval

def gl01(n):
    x, w = leggauss(n)
    return (x + 1) / 2, w / 2

def main():
    r = int(sys.argv[1]); th = eval(sys.argv[2]); v = float(sys.argv[3]); kap = eval(sys.argv[4])
    muc = float(sys.argv[5]); spec = sys.argv[6]; n = int(sys.argv[7]) if len(sys.argv) > 7 else 20
    if spec.startswith('mono:'):
        b = [float(s) for s in spec[5:].split(',')]
        P = lambda w: sum(bj * w ** j for j, bj in enumerate(b))
    else:
        js = json.load(open(spec)); A = np.array([float(s) for s in js['A']])
        sg = np.array([(-1) ** k for k in range(len(A))])
        P = lambda w: legval(2 * w - 1, A * sg)
    c = kap * math.pi; mu = muc * c
    xg, xw = gl01(n); tg, tw = gl01(24)
    # collapsed coordinates: x; (x1,x2) = (1-x)(s, (1-s)u); (x3,x4) = (1-x)(s', (1-s')u')
    S1, U1 = np.meshgrid(xg, xg, indexing='ij'); W1 = np.outer(xw, xw) * (1 - S1)
    s, u, wsu = S1.ravel(), U1.ravel(), W1.ravel()
    acc = np.zeros(2, dtype=complex)
    for ix in range(n):
        x = xg[ix]; wx = xw[ix] * (1 - x) ** 4
        x1 = (1 - x) * s; x2 = (1 - x) * (1 - s) * u
        X1 = x1[:, None]; X2 = x2[:, None]; WL = wsu[:, None]
        X3 = x1[None, :]; X4 = x2[None, :]; WR = wsu[None, :]
        wt = wx * WL * WR * x ** (r * r - 1) * (X1 * X2 * X3 * X4) ** (r - 1) * P(1 - x - X1 - X2) * P(1 - x - X3 - X4)
        A1 = 1 - th * (X1 + X3); A2 = 1 - th * (X2 + X4)
        V = v - th * (x + X1 + X2 + X3 + X4)
        # (1-e^{-ia})(1-e^{ib}) = 1 - e^{-ia} - e^{ib} + e^{i(b-a)}; a = c(th(x3-x4) + tau1 - tau2), b = c(th(x1-x2) + tau1 - tau2)
        ph0 = np.exp(1j * c * th * (X2 - X4)) * np.exp(-1j * mu * th * (X1 + X2 - X3 - X4))
        ea = np.exp(-1j * c * th * (X3 - X4)); eb = np.exp(1j * c * th * (X1 - X2))
        def mom(Aa, z, m):  # int_0^Aa tau^m e^{z tau} d tau, Gauss 24
            out = np.zeros(Aa.shape, dtype=complex)
            for k in range(24):
                t = Aa * tg[k]
                out += tw[k] * Aa * t ** m * np.exp(z * t)
            return out
        pm = lambda Aa, m: Aa ** (m + 1) / (m + 1)
        # terms: const 1 + e^{i(b-a)} (tau-free), -e^{-ia} (e^{-ic tau1} e^{+ic tau2}), -e^{ib} (e^{ic tau1} e^{-ic tau2})
        cf = 1 + eb * np.conj(np.exp(1j * c * th * (X3 - X4))) * 1.0  # 1 + e^{i c th (x1-x2-x3+x4)}
        cf = 1 + np.exp(1j * c * th * (X1 - X2 - X3 + X4))
        m1p = [mom(A1, -1j * c, k) for k in range(3)]; m2p = [mom(A2, 1j * c, k) for k in range(3)]
        m1m = [mom(A1, 1j * c, k) for k in range(3)]; m2m = [mom(A2, -1j * c, k) for k in range(3)]
        q1 = [pm(A1, k) for k in range(3)]; q2 = [pm(A2, k) for k in range(3)]
        def box(f1, f2, j):
            if j == 0:
                return f1[0] * f2[0]
            return V * V * f1[0] * f2[0] - 2 * V * (f1[1] * f2[0] + f1[0] * f2[1]) + f1[2] * f2[0] + 2 * f1[1] * f2[1] + f1[0] * f2[2]
        for j in range(2):
            k = cf * box(q1, q2, j) - ea * box(m1p, m2p, j) - eb * box(m1m, m2m, j)
            acc[j] += np.sum(wt * ph0 * k)
    c0, c1 = acc
    print('r=%d theta=%s v=%s kappa=%s mu/c=%s n=%d: c0 = %.12e %+.2ei, c1 = %.12e %+.2ei, h1 = c0/(k^2 c1) = %.10f'
          % (r, sys.argv[2], sys.argv[3], sys.argv[4], muc, n, c0.real, c0.imag, c1.real, c1.imag, (c0 / (kap ** 2 * c1)).real))

if __name__ == '__main__':
    main()

