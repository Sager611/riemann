"""Falsification test: complex amplifiers. GapCeil keeps Re(h) of the Hermitian Taylor arrays h (real phi).
For complex P the main term is sesquilinear: P(w1) conj P(w2) with the same complex kernel, so the form is the
Hermitian operator built from h itself. Compare top eigenvalues sup_v lambda_max (real vs complex pencil)."""
import sys, time, resource
import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
from math import factorial
sys.path.insert(0, HERE)
import gcore as G
import gcore_r as GR
from scipy.optimize import minimize_scalar

def arrays(c, r, D, R=0.4, K=64):
    ang = np.exp(2j * np.pi * np.arange(K) / K)
    Z1, Z2 = np.meshgrid(R * ang, R * ang, indexing='ij')
    gs = []
    for g in GR.gbar_r(Z1, Z2, c, r):
        C = np.fft.fft2(g) / K ** 2
        sc = R ** -np.arange(K)
        gs.append((C * np.outer(sc, sc))[:D, :D])
    g0, g1, g2 = gs
    a = g0; b = G.shift_sum(g0) + g1
    cc = G.shift_sum(G.shift_sum(g0)) + 2 * G.shift_sum(g1) + g2
    return a, b, cc

def mats(c, th, r, N, D, cplx):
    a, b, cc = arrays(c, r, D)
    herm_err = max(np.abs(h - h.conj().T).max() / np.abs(h).max() for h in (a, b, cc))
    if cplx:
        a, b, cc = [0.5 * (h + h.conj().T) for h in (a, b, cc)]
    else:
        a, b, cc = [0.5 * (h + h.T).real for h in (a, b, cc)]
    N0 = N + D + 2
    B = G.Basis(th, N0, nq=N0 + 2 * r * r + 12)
    B.W = [(B.V * (B.w * B.X ** q)[:, None]).T @ B.V for q in range(r * r + 2)]
    Ipow = [np.eye(N0)]
    for k in range(D):
        Ipow.append(B.I @ Ipow[-1])
    Ip = [P[:, :N] for P in Ipow]
    s0 = 2 * r - 1
    fac = [float(factorial(k)) for k in range(D + s0 + 2)]
    q0 = r * r - 1
    def Gm(p, h):
        out = np.zeros((N, N), complex if cplx else float)
        WP = B.W[q0 + p]
        L = [WP @ Ip[k] for k in range(D)]
        for i in range(D):
            for j in range(D):
                hv = h[i, j]
                if hv == 0:
                    continue
                out = out + fac[i + s0] * fac[j + s0] * hv * (Ip[i].T @ L[j])
        return 0.5 * (out + out.conj().T)
    M0 = Gm(0, a); MS = Gm(1, a) + Gm(0, b); MSS = Gm(2, a) + 2 * Gm(1, b) + Gm(0, cc)
    return M0, MS, MSS, herm_err

def top(M0, MS, MSS, v):
    M1 = v * v * M0 - 2 * v * MS + MSS
    L = np.linalg.cholesky(M1)
    Li = np.linalg.inv(L)
    ev = np.linalg.eigvalsh(Li @ M0 @ Li.conj().T)
    return ev[-1]

def best(M0, MS, MSS):
    rr = minimize_scalar(lambda v: -top(M0, MS, MSS, v), bounds=(0.6, 2.2), method='bounded', options={'xatol': 1e-9})
    return -rr.fun, rr.x

def main():
  t0 = time.time()
  r = int(sys.argv[1]); ks = [float(x) for x in sys.argv[2].split(',')]; N = int(sys.argv[3]) if len(sys.argv) > 3 else 16
  for kap in ks:
      c = kap * np.pi
      out = []
      for cp in (False, True):
          M0, MS, MSS, he = mats(c, 0.25, r, N, 36, cp)
          lam, v = best(M0, MS, MSS)
          out.append((lam / kap ** 2, v))
      print('r=%d kappa=%.4f N=%d: real %.7f (v=%.4f)  complex %.7f (v=%.4f)  herm_err=%.1e [%.0fs]' % (r, kap, N, out[0][0], out[0][1], out[1][0], out[1][1], he, time.time() - t0), flush=True)
  print('peak %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))

if __name__ == "__main__":
  main()
