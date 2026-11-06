#!/usr/bin/env python3
"""odlyzko31 experiment 2: nonnegative autocorrelation kernels. q = sum_{j<J} a_j cos((j+1/2) pi t) on [-1,1], h = q^2 >= 0,
g = (h*h)/(h*h)(0) >= 0 with nonnegative Fourier transform; F(x) = g(x/b) sech(w x). J = 1 gives Odlyzko's H.
Float only: optimise (a, b) by Nelder-Mead for given (n, r1, w)."""
import math, sys
import numpy as np
from kern import sech, EG, L8P, consts

def h_coeffs(a):
    """h = q^2 = sum_m c_m cos(m pi t), m = 0..2J-1."""
    J = len(a); c = np.zeros(2 * J)
    for j in range(J):
        for k in range(J):
            c[j + k + 1] += 0.5 * a[j] * a[k]
            c[abs(j - k)] += 0.5 * a[j] * a[k]
    return c

def Cmat_int(u, Mx):
    """C_ml(u) = int_{-1}^{1-u} cos(m pi t) cos(l pi (t+u)) dt, m,l < Mx, u array in [0,2]."""
    om = np.arange(Mx) * np.pi
    u = np.asarray(u)[:, None, None]
    om_m = om[None, :, None]; om_l = om[None, None, :]
    lo, hi = -1.0, 1.0 - u
    def part(freq, phase_coef):
        # int_lo^hi cos(freq t + phase) dt, phase = phase_coef * u
        ph = phase_coef * u
        with np.errstate(divide='ignore', invalid='ignore'):
            v = (np.sin(freq * hi + ph) - np.sin(freq * lo + ph)) / freq
        z = np.broadcast_to(np.abs(freq) < 1e-12, v.shape)
        return np.where(z, np.broadcast_to((hi - lo) * np.cos(ph), v.shape), v)
    return 0.5 * (part(om_m + om_l, om_l) + part(om_m - om_l, -om_l))

_T, _W = np.polynomial.legendre.leggauss(40)
def nodes(a, c, k):
    br = np.linspace(a, c, k + 1)
    xs = np.concatenate([(q - p) / 2 * _T + (q + p) / 2 for p, q in zip(br[:-1], br[1:])])
    ws = np.concatenate([(q - p) / 2 * _W for p, q in zip(br[:-1], br[1:])])
    return xs, ws

class Fam:
    def __init__(self, J, nu=24):
        self.J = J; self.Mx = 2 * J
        self.u, self.wu = nodes(0.0, 2.0, nu)
        self.C = Cmat_int(self.u, self.Mx)          # (N, Mx, Mx)
        self.C0 = Cmat_int(np.array([0.0]), self.Mx)[0]
    def g(self, a):
        c = h_coeffs(a)
        return np.einsum('ijk,j,k->i', self.C, c, c) / (c @ self.C0 @ c)
    def value(self, n, r1, w, a, b):
        g = self.g(a)
        x = b * self.u; wx = b * self.wu
        F = g * sech(w * x)
        X = 2 * b
        Ic = float(np.sum(wx * (1 - F) / (2 * np.sinh(x / 2)))) - math.log(math.tanh(X / 4))
        Ir = float(np.sum(wx * (1 - F) / (2 * np.cosh(x / 2)))) + math.pi / 2 - 2 * math.atan(math.tanh(X / 4))
        E = 4 * float(np.sum(wx * F * np.cosh(x / 2)))
        lB = EG + L8P - Ic; lA = lB + math.pi / 2 - Ir
        return (r1 * lA + (n - r1) * lB - E) / n

def nelder(f, x0, step, iters=400):
    d = len(x0); pts = [np.array(x0, float)] + [np.array(x0, float) + step * np.eye(d)[i] for i in range(d)]
    vals = [f(p) for p in pts]
    for _ in range(iters):
        o = np.argsort(vals)[::-1]; pts = [pts[i] for i in o]; vals = [vals[i] for i in o]
        c = sum(pts[:-1]) / d
        xr = c + (c - pts[-1]); fr = f(xr)
        if fr > vals[0]:
            xe = c + 2 * (c - pts[-1]); fe = f(xe)
            pts[-1], vals[-1] = (xe, fe) if fe > fr else (xr, fr)
        elif fr > vals[-2]:
            pts[-1], vals[-1] = xr, fr
        else:
            xc = c + 0.5 * (pts[-1] - c); fc = f(xc)
            if fc > vals[-1]: pts[-1], vals[-1] = xc, fc
            else:
                for i in range(1, d + 1):
                    pts[i] = pts[0] + 0.5 * (pts[i] - pts[0]); vals[i] = f(pts[i])
    i = int(np.argmax(vals)); return pts[i], vals[i]

def optimise(Fm, n, r1, w, b0):
    J = Fm.J
    def f(p):
        a = np.concatenate([[1.0], p[1:]]); b = math.exp(p[0])
        if b > 200 or b < 0.05: return -1e9
        return Fm.value(n, r1, w, a, b)
    p0 = np.concatenate([[math.log(b0)], np.zeros(J - 1)])
    best = None
    for rep in range(3):
        p, v = nelder(f, p0, 0.15 if rep == 0 else 0.05)
        if best is None or v > best[1]: best = (p, v)
        p0 = p
    p, v = best
    return math.exp(p[0]), np.concatenate([[1.0], p[1:]]), v

if __name__ == '__main__':
    J = int(sys.argv[1]) if len(sys.argv) > 1 else 3
    F1 = Fam(1)
    lA, lB, E = consts('H', 3.0, 0.5)
    print('sanity J=1 vs H (b=3, w=1/2, n=10, r1=0):', (10 * lB - E) / 10, F1.value(10, 0, 0.5, np.array([1.0]), 3.0))
    Fm = Fam(J)
    for w, lab in ((0.5, 'classical'), (15 / 32, 'ct')):
        for n, r1, b0 in ((8, 0, 2.1), (100, 0, 7.0), (100, 100, 8.6), (1000, 0, 15.0), (1000, 1000, 18.0), (2732, 0, 20.4), (10000, 0, 29.0), (10000, 10000, 34.0)):
            b, a, v = optimise(Fm, n, r1, w, b0)
            vH = max((r1 * consts('H', bb, w)[0] + (n - r1) * consts('H', bb, w)[1] - consts('H', bb, w)[2]) / n for bb in np.exp(np.linspace(math.log(b0/2), math.log(2*b0), 200)))
            print(f'{lab} n={n} r1={r1}: J={J} b={b:.4f} a={np.round(a,4)} rd>{math.exp(v):.4f}  (H: {math.exp(vH):.4f}, gain {100*(v-vH):+.3f}% in log)', flush=True)
