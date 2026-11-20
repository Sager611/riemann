"""explore53 falsification: every F = sin(x) H(x) (H smooth, period 2 pi M) vanishes at k pi, so all its zero gaps are <= pi.
A valid per-gap form with symbol p(w) = sum_j a_j w^{2j} must give  sum_w p(w) |F^(w)|^2 >= 0 (gaps sum + telescoping).
F^(w) = (h(w-1) - h(w+1)) / 2i  on the frequency grid w in Z/M.  min over h of Q/N = smallest generalized eigenvalue.
Negative => the form is NOT valid (explicit counterexample H).  Also scans |sin|-like F = sin x tanh(sin x / eps)."""
import sys, json
import numpy as np


def mineig(a, M=8, W=12.0):
    nu = np.arange(-int(W * M), int(W * M) + 1) / M          # support of h
    om = np.arange(-int((W + 1) * M), int((W + 1) * M) + 1) / M  # support of F^
    D = np.zeros((len(om), len(nu)), complex)
    idx = {round(x * M): i for i, x in enumerate(om)}
    for j, v in enumerate(nu):
        D[idx[round((v + 1) * M)], j] += 1 / 2j      # h(nu) contributes to w = nu+1 with +1/(2i)
        D[idx[round((v - 1) * M)], j] -= 1 / 2j      # and to w = nu-1 with -1/(2i)
    p = sum(a[j] * om ** (2 * j) for j in range(len(a)))
    Q = D.conj().T @ (p[:, None] * D)
    Nn = D.conj().T @ D
    L = np.linalg.cholesky(Nn + 1e-14 * np.eye(len(nu)))
    Li = np.linalg.inv(L)
    ev, U = np.linalg.eigh(Li @ Q @ Li.conj().T)
    return ev[0]


def tanhtest(a, eps, n=4096):
    x = np.arange(n) * 2 * np.pi / n
    F = np.sin(x) * np.tanh(np.sin(x) / eps)
    Fh = np.fft.fft(F) / n
    w = np.fft.fftfreq(n, d=1.0 / n)
    p = sum(a[j] * w ** (2 * j) for j in range(len(a)))
    return float((p * abs(Fh) ** 2).sum() / (abs(Fh) ** 2).sum())


if __name__ == '__main__':
    forms = json.load(open(sys.argv[1]))
    for name, a in forms.items():
        out = [name]
        for M, W in [(4, 8.0), (8, 12.0)]:
            out.append('mineig(M=%d,W=%g)=%.6g' % (M, W, mineig(a, M, W)))
        out.append('tanh: ' + ' '.join('%.3g:%.4g' % (e, tanhtest(a, e)) for e in (0.02, 0.05, 0.1, 0.3, 1.0, 3.0)))
        print('  '.join(out), flush=True)
