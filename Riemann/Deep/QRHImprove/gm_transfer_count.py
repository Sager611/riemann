#!/usr/bin/env python3
"""INDICATIVE ONLY: value of a Guth-Maynard-type large-value theorem for the sextic family (no such theorem exists).

Hypothesis H_GM: GM Theorem 1.1 (arXiv 2405.20552) holds verbatim for the family {psi_u : N(u) ~ U} with T -> U:
a polynomial of length U^s with bounded coefficients is >= U^w at <= U^{max(2s-2w, 18s/5-4w, 1+12s/5-4w)} rows.
It is applied to products M^i P^j Q of the inverse polynomial M (length r), the plain polynomial P (length m) and
selected prime sums Q (total length z <= L, each used once), next to Liu's moment counts e_I, e_S, e_L
(arXiv 2610.12234 (34.5)-(34.7); arXiv 2610.26835 (14.3)). Class count = min over cutoffs t of the max over
adversarial splits (r,m). Chen-Gupta-Li (arXiv 2507.08296) prove a GM analogue only for all characters mod one q.
Run: python3 gm_transfer_count.py   (numpy only; about 70 s, about 30 MB)
"""
import time
import numpy as np

T0 = time.time()
ALPHA = 5/6
RR = np.linspace(0.5, 1.5, 201)[:, None]
MM = np.linspace(0, 0.5, 101)[None, :]
TT = np.linspace(1, 1.5, 101)
PROD = [(i, j) for i in range(4) for j in range(5)]

def R143(d, x, k):
    c = 1/(3*k); D = 3-(1+2*c)*x; P = (2-2*c*x)*(1-x); J = (5/6-d)*D+d*P
    return 1-d+(5/6-d)*d*P/(2*J)

def count(d, x, k, L, gm=True):
    a = (1+d)/2; g = d*x; c = 1/(3*k)
    einv = np.where(RR <= 1, 1-d*RR-2*g*np.minimum((1-RR)/2, L), 1-ALPHA+(ALPHA-d)*RR)
    eS = 1-2*d*MM-2*g*np.minimum(np.maximum(1-2*MM, 0)*c/2, L)
    C = np.minimum(einv, eS)
    if gm:
        for i, j in PROD:
            for z in np.linspace(0, L, 8):
                s = i*RR+j*MM+z; w = a*(i*RR+j*MM)+z*(0.5+g)
                if i == 0 and j == 0 and z == 0:
                    continue
                C = np.minimum(C, np.maximum(np.maximum(2*s-2*w, 3.6*s-4*w), 1+2.4*s-4*w))
    best = 1.0
    for t in TT:
        W = (RR >= t-0.5-1e-12) & (RR <= t+1e-12) & (RR+MM >= t-1e-12)
        best = min(best, float(np.max(np.where(W, C, -np.inf))))
    return best

def E_full(d, x, b, l, k, R):          # Liu (14.8)
    c = 1/(3*k); h = (1+b+3*l)/2; ly = (1+b-l)/2; K = (3*R+1)/4
    return -1/4+b/6+5*l/4+d/2+d*x*l+d*x*c*(h/4-ly/2)-h+h*K

dg = np.linspace(1/8, 5/6, 72); xg = np.linspace(0, 0.5, 11)
DD, XX = np.meshgrid(dg, xg, indexing="ij")

def Rgrid(k, L, gm):
    return np.array([[count(d, x, k, L, gm) for x in xg] for d in dg])

def best_B(Rg, k):
    def maxE(b, l): return E_full(DD, XX, b, l, k, Rg).max()
    lo, hi = 0.15, 0.20
    for _ in range(30):
        l = (lo+hi)/2
        v = min(maxE(b, l) for b in np.linspace(0, 0.25, 501))
        lo, hi = (l, hi) if v <= 0 else (lo, l)
    return 11/12-lo/4, lo

if __name__ == "__main__":
    k0, L0 = 0.749773294, 0.21047
    dc = 0.469345
    print(f"check: minimax without GM at critical class = {count(dc, 0.5, k0, L0, False):.4f}"
          f"  (closed form (14.3) = {R143(dc, 0.5, k0):.4f})")
    print(f"critical class with H_GM: R = {count(dc, 0.5, k0, L0, True):.4f}")
    for gm in (False, True):
        k = k0
        B, l = 11/12-0.15/4, 0.15
        for _ in range(2):
            Rg = Rgrid(k, L0, gm)
            B, l = best_B(Rg, k)
            k = 2*B-1
        err = np.abs(Rgrid(k0, L0, False)-R143(DD, XX, k0)).max() if not gm else float("nan")
        print(f"{'H_GM ' if gm else 'no GM'}: coarse-grid B = {B:.5f}  l = {l:.5f}  (max |R_minimax-R_143| = {err:.4f})"
              f"  [{time.time()-T0:.0f}s]", flush=True)
