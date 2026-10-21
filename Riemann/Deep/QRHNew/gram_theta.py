#!/usr/bin/env python3
"""Value of a smaller sixth-power loss in the additive Gram bound (low side of the quasi-RH probe).

OpenAI 7/8 Prop. 15.2: sum_m |A_m|^2 << (Q/Y')(1 + P_a^{1/6} + P_a^2/Y'); the P_a^{1/6} term comes from the
exceptional frequencies k = a^6 e. Suppose it were P_a^theta (theta in [0,1/6]). Then the direct bound is
Z^{lx/2 + theta*b/2}, so matching C_b(B) = B - (4+b)/6 gives
    B_low = 11/12 - l/4 - b*(1/12 - theta/2)          (balanced lengths lx=(1-b-l)/2, ly=(1+b-l)/2)
and the Gram side condition P_a^2/Y' <= P_a^theta needs ly - l >= (2-theta) b.
High side (unchanged): Liu arXiv 2610.26835 (14.8) relative to C_b(B):
    B_high = 11/12 - l/4 + max_{delta,x} E(delta,x,b,l,kappa),  kappa = 2B-1 (licensing fixed point).
Output: min over (b,l) of max(B_low, B_high).   Run: python3 gram_theta.py  (numpy; ~30 s, <30 MB)
"""
import numpy as np

DG = np.linspace(1/8, 5/6, 285)
XG = np.linspace(0, 0.5, 51)
DD, XX = np.meshgrid(DG, XG, indexing="ij")


def r143(d, x, k):
    c = 1/(3*k)
    dk = 3-(1+2*c)*x
    pk = (2-2*c*x)*(1-x)
    jk = (5/6-d)*dk+d*pk
    return 1-d+(5/6-d)*d*pk/(2*jk)


def max_e(b, ell, k):
    c = 1/(3*k)
    h = (1+b+3*ell)/2
    ly = (1+b-ell)/2
    rr = r143(DD, XX, k)
    e = -1/4+b/6+5*ell/4+DD/2+DD*XX*ell+DD*XX*c*(h/4-ly/2)-(3*h/4)*(1-rr)
    return float(e.max())


def solve(theta, k0=0.7498, iters=4):
    k = k0
    best_b = 9.0
    best_pt = (0.0, 0.0, 9.0, 9.0)
    for _ in range(iters):
        best_b = 9.0
        for ell in np.linspace(0.12, 0.20, 81):
            for b in np.linspace(0.0, 0.30, 121):
                lx, ly = (1-b-ell)/2, (1+b-ell)/2
                if lx <= ell or ly-ell < (2-theta)*b:
                    continue
                b_low = 11/12-ell/4-b*(1/12-theta/2)
                b_high = 11/12-ell/4+max_e(b, ell, k)
                bb = max(b_low, b_high)
                if bb < best_b:
                    best_b = bb
                    best_pt = (float(b), float(ell), b_low, b_high)
        k = max(13/18, 2*best_b-1)
    return best_b, best_pt, k


if __name__ == "__main__":
    for th in (1/6, 1/8, 1/12, 1/24, 0.0):
        bval, (b, ell, bl, bh), kap = solve(th)
        print(f"theta={th:.4f}: B={bval:.6f}  b={b:.4f} l={ell:.4f}  B_low={bl:.6f} B_high={bh:.6f} "
              f"kappa={kap:.4f}", flush=True)
