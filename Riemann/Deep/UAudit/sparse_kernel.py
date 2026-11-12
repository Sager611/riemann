#!/usr/bin/env python3
"""uaudit N2: sparse-kernel bound of DPCUncond PROOF 6c on actual zeta zeros.
K_a(u) = sum_{rho in S} k_{d_rho}(u - x_rho), k_d(u) = (2 pi)^-1 log(1 + d^2/u^2), x_rho = gamma - h/2,
S = subset of the zeros in [T,2T] (random / tightest local clusters / contiguous block), depths d in [a/L,(a+1)/L)
('mixed') or all = Y_a ('full', the kernel K'_a). Claim: ||K_a||_q^q << (a+2)^q sigma^(1-1/w) T, sigma = |S|/(T L).
Reports ||K||_q^q, the near-part bound C_q Y sum_S N_rho^(q-1) (N_rho = #S within 2Y), and the ratios
r1 = ||K||_q^q / ((a+2)^q sigma T)   (stronger than claimed),  r_near = ||K||_q^q / near-part bound.
Integration: scipy quad between consecutive singular points x_rho, kernels truncated at |u - x_rho| <= 40 Y
(the dropped tail is < 1e-3 relative).
"""
import sys, json, time, math
import numpy as np
from scipy.integrate import quad

ZF = sys.argv[2] if len(sys.argv) > 2 else "../ExploreHelio/results/zeta_zeros.npy"
Z = np.sort(np.load(ZF))
T = 2000.0
L = math.log(T)
PHI = 0.7
H = 2 * math.pi * PHI / L
Z1 = Z[(Z >= T) & (Z <= 2 * T)]
NZ = len(Z1)
rng = np.random.default_rng(12345)

def Cq(q):
    f = lambda v: (math.log1p(1.0 / (v * v)) / (2 * math.pi)) ** q
    return 2 * (quad(f, 0, 1, limit=200)[0] + quad(f, 1, np.inf, limit=200)[0])

def norm_q(x, d, q, Y):
    """int K^q, K = sum k_d(u - x), x sorted."""
    W = 40 * Y
    tot = 0.0
    n = len(x)
    # merged clusters of [x-W, x+W]
    i = 0
    while i < n:
        j = i
        while j + 1 < n and x[j + 1] - x[j] <= 2 * W:
            j += 1
        pts = [x[i] - W] + list(x[i:j + 1]) + [x[j] + W]
        for k in range(len(pts) - 1):
            lo, hi = pts[k], pts[k + 1]
            m0 = np.searchsorted(x, lo - W)
            m1 = np.searchsorted(x, hi + W)
            xs, ds = x[m0:m1], d[m0:m1]
            f = lambda u: (np.sum(np.log1p(ds * ds / ((u - xs) ** 2 + 1e-300))) / (2 * math.pi)) ** q
            tot += quad(f, lo, hi, limit=60, epsrel=1e-4, epsabs=1e-9)[0]
        i = j + 1
    return tot

def subset(kind, n, Y):
    if kind == "random":
        idx = np.sort(rng.choice(NZ, n, replace=False))
    elif kind == "tight":
        cnt = np.searchsorted(Z1, Z1 + 2 * Y) - np.searchsorted(Z1, Z1 - 2 * Y)
        idx = np.sort(np.argsort(-cnt + 1e-6 * rng.random(NZ), kind="stable")[:n])
    else:  # block
        idx = np.arange(n)
    return Z1[idx]

if __name__ == "__main__":
    t0 = time.time()
    out = dict(T=T, L=L, h=H, nzeros=NZ, rows=[])
    CQ = {q: Cq(q) for q in (2, 4)}
    out["Cq"] = CQ
    for a in (0, 2, 5):
        Y = (a + 1) / L
        for frac in (0.3, 0.05):
            n = int(round(frac * NZ))
            sigma = n / (T * L)
            for kind in ("random", "tight", "block"):
                g = subset(kind, n, Y)
                x = np.sort(g - H / 2)
                Nr = np.searchsorted(x, x + 2 * Y, side="right") - np.searchsorted(x, x - 2 * Y, side="left")
                for depth in ("mixed", "full"):
                    if depth == "mixed":
                        d = (a + rng.random(n)) / L
                    else:
                        d = np.full(n, Y)
                    for q in (2, 4):
                        Kq = norm_q(x, d, q, Y)
                        near = CQ[q] * Y * float(np.sum(Nr.astype(float) ** (q - 1)))
                        row = dict(a=a, frac=frac, n=n, sigma=sigma, kind=kind, depth=depth, q=q, Kq=Kq,
                                   near_bound=near, r1=Kq / ((a + 2) ** q * sigma * T), r_near=Kq / near,
                                   maxN=int(Nr.max()), meanN=float(Nr.mean()))
                        out["rows"].append(row)
                        print("a=%d frac=%.2f %-6s %-5s q=%d  ||K||^q=%.4g  r1=%.3g  r_near=%.3g  maxN=%d  %.0fs" % (
                            a, frac, kind, depth, q, Kq, row["r1"], row["r_near"], row["maxN"], time.time() - t0), flush=True)
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "sparse_kernel.json", "w"), indent=1)
    print("time %.1f s" % (time.time() - t0))
