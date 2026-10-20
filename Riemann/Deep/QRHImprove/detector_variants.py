#!/usr/bin/env python3
"""Detector variants of the refined quasi-RH system (Liu 2610.26835 (14.7), averaged numerator).

Local factors from OpenAI 7/8 Lemma 7.1 (table (7.16)); exponents relative to C(B):
  standard slot (marks p | c, Z -> Zq_p, minus rescaled term): row factor q^(z-1)chi_p(u), h = 1-lx+l.
  v-marked slot (p | Gauss modulus s, Y -> Yq_p; t=0,k>=1 rows give W/(1-W)): row factor chi_p(u)/q_p,
     no z-coupling; direct side = standard with ly -> ly + l2; E_V = E_std(lx, ly, l1) + l2(g - 1/2).
  marked-only slot (no rescaled term): leftover scalar q^(z-1-w) gives E += l2(delta/2 - g) >= 0.
  mollified detector = Mellin/shape sum of standard products; shape of length ls < l costs
     E_std(ls) - (l - ls)/6 against the full principal.
Run: python3 detector_variants.py  (numpy only, ~170 s, 53 MB)
"""
import numpy as np
from refined_system import E_end, DD, XX, golden_min

KAPPA = 0.7497732880799279          # 2*B_flat - 1 (most favourable consistent moment parameter)
BFLAT = 0.874886647039964


def maxE_v(b, l1, l2, k=KAPPA):
    """two-piece: theta-marked l1 + v-marked l2; b = (ly + l2) - lx, lx + ly + l2 = 1 - l1."""
    B = 11/12 - l1/4
    lx, ly = (1 - l1 - b)/2, (1 - l1 + b)/2 - l2
    return float((E_end(DD, XX, lx, ly, l1, k, B) + l2*(DD*XX - 0.5)).max())


def best_l1(l2):
    lo, hi = 0.12, 0.17
    for _ in range(22):
        l1 = (lo + hi)/2
        _, v = golden_min(lambda b: maxE_v(b, l1, l2), 0.0, 0.25, 1e-6)
        lo, hi = (l1, hi) if v <= 0 else (lo, l1)
    return lo


if __name__ == "__main__":
    for l2 in (0.0, 0.01, 0.03):
        l1 = best_l1(l2)
        print(f"v-marked long piece l2={l2:.2f}: best l1={l1:.6f}  B={11/12 - l1/4:.9f}  (B_flat={BFLAT:.9f})")
    b_opt, _ = golden_min(lambda b: maxE_v(b, best_l1(0.0), 0.0), 0.0, 0.25, 1e-6)
    l = best_l1(0.0)
    lx, ly, B = (1 - l - b_opt)/2, (1 - l + b_opt)/2, 11/12 - l/4
    worst = max(float(E_end(DD, XX, lx, ly, s, KAPPA, B).max()) - (l - s)/6 for s in np.linspace(0, l, 9)[:-1])
    print(f"mollified shapes ls<l at the optimum: max E = {worst:.4f} < 0 (full shape: 0)")
