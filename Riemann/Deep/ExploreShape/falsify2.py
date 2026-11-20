"""explore53 falsification 2: F = H(x) sin(phi(x)), phi' = 1 + beta(1 + cos(x/R)) >= 1 (zero gaps pi/phi' <= pi),
H any trig polynomial of period 2 pi R (|n| <= nmax).  min over H of sum_w p(w)|F^(w)|^2 / sum |F^|^2 (generalized eigenproblem);
also a twist e^{i u x} is free in H.  Must be >= 0 for a valid form."""
import sys, json
import numpy as np


def test(a, R, beta, nmax=40, n=2 ** 14):
    m = R * (1 + beta)
    assert abs(m - round(m)) < 1e-12
    X = 2 * np.pi * R
    x = np.arange(n) * X / n
    phi = (1 + beta) * x + beta * R * np.sin(x / R)
    s = np.sin(phi)
    w = np.fft.fftfreq(n, d=X / n) * 2 * np.pi
    p = sum(a[j] * w ** (2 * j) for j in range(len(a)))
    cols = []
    for k in range(-nmax, nmax + 1):
        cols.append(np.fft.fft(np.exp(1j * k * x / R) * s) / n)
    C = np.array(cols).T
    Q = C.conj().T @ (p[:, None] * C)
    N = C.conj().T @ C
    L = np.linalg.cholesky(N + 1e-13 * np.eye(N.shape[0]))
    Li = np.linalg.inv(L)
    return np.linalg.eigvalsh(Li @ Q @ Li.conj().T)[0]


if __name__ == '__main__':
    forms = json.load(open(sys.argv[1]))
    for name, a in forms.items():
        res = []
        for R, beta in [(2, 0.5), (4, 0.25), (5, 1.0), (8, 0.125), (10, 0.5)]:
            res.append('R%g,b%g:%.4g' % (R, beta, test(a, R, beta)))
        print(name, ' '.join(res), flush=True)
