"""gapcertaudit: independent float evaluation of Bui-Milinovich's main-term moments <S^m>/<1>, m <= 8, for the midpoint
mollifier (alpha3 = -beta3 = i mu/L), from BM's TeX only (arXiv 1410.3635v3, Lemma main l.370-384, operator l.422, l.431):
extends GapAudit/indep7.py (audited, which has m <= 2 through (v-S)^2) to the full S^m weight, S = theta(x+x1+..+x4)+tau1+tau2.
t3, t4 integrated exactly: BM's two linear factors x the t3,t4-integral = (1-e^{-ia})(1-e^{ib})/c^2 (own derivation, FULLPROOF (d)).
usage: indep8.py r theta kappa mu_over_c amp.json [n]   (amp.json: shifted-Legendre A_n, P(w) = sum (-1)^n A_n P_n(2w-1))"""
import sys, json, math, os
from math import factorial
import numpy as np
from numpy.polynomial.legendre import leggauss, legval

MMAX = int(os.environ.get('MMAX', '8'))   # MMAX=12 for the J = 6 moments


def gl01(n):
    x, w = leggauss(n)
    return (x + 1) / 2, w / 2


def main():
    r = int(sys.argv[1]); th = eval(sys.argv[2]); kap = eval(sys.argv[3]); muc = float(sys.argv[4])
    js = json.load(open(sys.argv[5])); n = int(sys.argv[6]) if len(sys.argv) > 6 else 16
    A = np.array([float(s) for s in js['A']]); sg = np.array([(-1) ** k for k in range(len(A))])
    P = lambda w: legval(2 * w - 1, A * sg)
    c = kap * math.pi; mu = muc * c
    xg, xw = gl01(n); tg, tw = gl01(28)
    S1, U1 = np.meshgrid(xg, xg, indexing='ij'); W1 = np.outer(xw, xw) * (1 - S1)
    s, u, wsu = S1.ravel(), U1.ravel(), W1.ravel()
    acc = np.zeros(MMAX + 1, dtype=complex)
    for ix in range(n):
        x = xg[ix]; wx = xw[ix] * (1 - x) ** 4
        x1 = (1 - x) * s; x2 = (1 - x) * (1 - s) * u
        X1 = x1[:, None]; X2 = x2[:, None]; WL = wsu[:, None]
        X3 = x1[None, :]; X4 = x2[None, :]; WR = wsu[None, :]
        wt = wx * WL * WR * x ** (r * r - 1) * (X1 * X2 * X3 * X4) ** (r - 1) * P(1 - x - X1 - X2) * P(1 - x - X3 - X4)
        A1 = 1 - th * (X1 + X3); A2 = 1 - th * (X2 + X4)
        TX = th * (x + X1 + X2 + X3 + X4)
        base = wt * np.exp(1j * c * th * (X2 - X4)) * np.exp(-1j * mu * th * (X1 + X2 - X3 - X4))
        cf = 1 + np.exp(1j * c * th * (X1 - X2 - X3 + X4))
        ea = np.exp(-1j * c * th * (X3 - X4)); eb = np.exp(1j * c * th * (X1 - X2))
        # tau-moments int_0^A tau^k e^{z tau} d tau, k <= MMAX (Gauss 28); e^{+ic} moments are conjugates of e^{-ic} ones
        m1 = [np.zeros(A1.shape, dtype=complex) for _ in range(MMAX + 1)]   # tau1, e^{-ic tau1}
        m2 = [np.zeros(A2.shape, dtype=complex) for _ in range(MMAX + 1)]   # tau2, e^{+ic tau2}
        for k in range(len(tg)):
            t1 = A1 * tg[k]; e1 = tw[k] * A1 * np.exp(-1j * c * t1)
            t2 = A2 * tg[k]; e2 = tw[k] * A2 * np.exp(1j * c * t2)
            p1 = np.ones(A1.shape); p2 = np.ones(A2.shape)
            for m in range(MMAX + 1):
                m1[m] += e1 * p1; m2[m] += e2 * p2
                p1 = p1 * t1; p2 = p2 * t2
        TXp = [np.ones(TX.shape)]
        for m in range(MMAX):
            TXp.append(TXp[-1] * TX)
        for k in range(MMAX + 1):
            q1 = A1 ** (k + 1) / (k + 1)
            for l in range(MMAX + 1 - k):
                q2 = A2 ** (l + 1) / (l + 1)
                # -e^{-ia}: e^{-ic tau1} e^{+ic tau2};  -e^{ib}: e^{+ic tau1} e^{-ic tau2} (conjugate moments)
                Pkl = base * (cf * q1 * q2 - ea * m1[k] * m2[l] - eb * np.conj(m1[k]) * np.conj(m2[l]))
                for p in range(MMAX + 1 - k - l):
                    m = k + l + p
                    coef = factorial(m) / (factorial(k) * factorial(l) * factorial(p))
                    acc[m] += coef * np.sum(Pkl * TXp[p])
    out = {'n': n, 'r': r, 'theta': sys.argv[2], 'kappa': sys.argv[3], 'mu_over_c': muc,
           'S0': [acc[0].real, acc[0].imag], 'ratios': [(acc[m] / acc[0]).real for m in range(MMAX + 1)],
           'imag_rel': [abs(acc[m].imag) / abs(acc[m].real) for m in range(MMAX + 1)]}
    print(json.dumps(out))


if __name__ == '__main__':
    main()
