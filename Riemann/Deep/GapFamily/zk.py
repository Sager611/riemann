"""gapfamily: main-term moment forms for Bui-Hall's f_k = e^{2iv theta} Z^(k)(t) Z^(k)(t+pi lam/L) A(s+2 pi i eta/L).
c_{k,j}(V) = int K (S-V)^{2j}, V = v+1 (BH's R(x) = v+1+x), K = P(w1)P(w2) e^{iE0} G3 G4 (REPORT working notes).
Returns real symmetric forms M[n] (basis w^i, i <= d) of int K w1^i w2^l (S-v0)^n, n <= K, by tensor Gauss-Legendre in 7-D
(X, r1, th1, r2, th2, tau1, tau2) + an inner GL for G3, G4.  Float only.
usage: python3 zk.py k lam eta V0 d Kmax [nodes]  -> prints JSON"""
import sys, json, math
import numpy as np
from numpy.polynomial.legendre import leggauss


def gl(n, a=0.0, b=1.0):
    x, w = leggauss(n)
    return (a + b) / 2 + (b - a) / 2 * x, (b - a) / 2 * w


def forms(k, lam, eta, v0, d=1, Kmax=12, th=0.2499, nodes=(11, 12, 10, 12, 10, 16, 16), nu=None, chunk=600):
    c = math.pi * lam
    n0, n1, nt1, n2, nt2, n3, n4 = nodes
    if nu is None:
        nu = 14 + k
    s0, w0 = gl(n0); s1, w1g = gl(n1); t1, wt1 = gl(nt1); s2, w2g = gl(n2); t2, wt2 = gl(nt2)
    s3, w3 = gl(n3); s4, w4 = gl(n4); uu, wu = gl(nu)
    # outer grid (X, r1, th1, r2, th2)
    O = np.stack(np.meshgrid(np.arange(n0), np.arange(n1), np.arange(nt1), np.arange(n2), np.arange(nt2), indexing='ij'), -1).reshape(-1, 5)
    I3, I4 = np.meshgrid(np.arange(n3), np.arange(n4), indexing='ij'); I3 = I3.ravel(); I4 = I4.ravel()
    S3 = s3[I3][None, :]; S4 = s4[I4][None, :]; W34 = (w3[I3] * w4[I4])[None, :]
    M = np.zeros((Kmax + 1, d + 1, d + 1), dtype=complex)
    for st in range(0, len(O), chunk):
        o = O[st:st + chunk]
        X = th * s0[o[:, 0]]; Wo = th * w0[o[:, 0]]
        R1 = (th - X) * s1[o[:, 1]]; Wo = Wo * (th - X) * w1g[o[:, 1]] * R1
        T1 = t1[o[:, 2]]; Wo = Wo * wt1[o[:, 2]]
        R2 = (th - X) * s2[o[:, 3]]; Wo = Wo * (th - X) * w2g[o[:, 3]] * R2
        T2 = t2[o[:, 4]]; Wo = Wo * wt2[o[:, 4]]
        X1 = R1 * T1; X2 = R1 * (1 - T1); X3 = R2 * T2; X4 = R2 * (1 - T2)
        A1 = 1 - X1 - X3; A2 = 1 - X2 - X4
        col = lambda a: a[:, None]
        tau1 = col(A1) * S3; tau2 = col(A2) * S4
        W = col(Wo * A1 * A2) * W34
        sig = tau1 - tau2
        D3 = col(X3 - X4) + sig; D4 = col(X1 - X2) + sig
        p3 = 0.5 - col(X3) - tau1; q3 = 0.5 - col(X4) - tau2
        p4 = 0.5 - col(X1) - tau1; q4 = 0.5 - col(X2) - tau2
        G3 = np.zeros(D3.shape, dtype=complex); G4 = np.zeros(D4.shape, dtype=complex)
        for ui, wi in zip(uu, wu):
            a3 = ui * D3; a4 = ui * D4
            G3 += wi * ((p3 + a3) * (q3 - a3)) ** k * np.exp(-1j * c * a3)
            G4 += wi * ((p4 + a4) * (q4 - a4)) ** k * np.exp(1j * c * a4)
        G3 *= D3; G4 *= D4
        E0 = c * (X2 - X4) - 2 * math.pi * eta * (X1 + X2 - X3 - X4)
        Kr = W * np.exp(1j * col(E0)) * G3 * G4
        S = col(X + R1 + R2) + tau1 + tau2 - v0
        wv1 = 1 - (X + R1) / th; wv2 = 1 - (X + R2) / th
        Sp = Kr.copy()
        for n in range(Kmax + 1):
            row = Sp.sum(axis=1)
            for i in range(d + 1):
                for l in range(d + 1):
                    M[n, i, l] += np.sum(row * wv1 ** i * wv2 ** l)
            Sp = Sp * S
    herm = float(np.max(np.abs(M - np.conj(np.transpose(M, (0, 2, 1))))) / np.max(np.abs(M)))
    Ms = 0.5 * (M + np.transpose(M, (0, 2, 1)))
    imag = float(np.max(np.abs(Ms.imag)) / np.max(np.abs(Ms.real)))
    return Ms.real, herm, imag


if __name__ == '__main__':
    k = int(sys.argv[1]); lam = float(sys.argv[2]); eta = float(sys.argv[3]); v0 = float(sys.argv[4])
    d = int(sys.argv[5]); Kmax = int(sys.argv[6])
    nodes = tuple(int(x) for x in sys.argv[7].split(',')) if len(sys.argv) > 7 else (11, 12, 10, 12, 10, 16, 16)
    M, herm, imag = forms(k, lam, eta, v0, d, Kmax, nodes=nodes)
    print(json.dumps({'k': k, 'lam': lam, 'eta': eta, 'v0': v0, 'd': d, 'nodes': nodes, 'herm': herm, 'imag': imag,
                      'M': M.tolist()}))
