#!/usr/bin/env python3
"""Test B: cost bookkeeping of the numeric core at the explicit parameters of
Formalism/SizeWindowBrackets/SChain.lean SizesTextAt (all quantities in log-space; b = exp(x^eps) is huge).
Checks RestSize, CostWindows, SizeWindows (m <= x^(1/2+eps/2), x m <= k), and measures
errorBracketSizeWindows / (m nT b) against Main/(m nT b) ~ x log a, i.e. the relative error exponent.
Also: the m-balance (1/m vs m/x) and the hypothetical m^2-cost variant (the remark's 'stop at 2/3')."""
import math
def run(eps, Lx):          # Lx = log x
    x = math.exp(Lx) if Lx < 700 else None
    lg = lambda e: e * Lx  # log of x^e
    n = lg(2 + 8 * eps)                        # log n  (n ~ x^(2+8eps))
    m = lg(0.5 - eps / 16)                     # log m
    nT = lg(-eps / 8) + n - math.log(2)        # log nT
    nT2 = lg(-eps) + n
    Lb = math.exp(eps * Lx)                    # log b = x^eps
    La = Lx + Lb                               # log a, a = x b
    k = lg(4.5 + 10 * eps)
    f = math.log(2) + n                        # log f, f = |PhiSet| = |E0 x H| / 2 = 2n
    eta = math.exp(eps / 4 * Lx)               # log eta = x^(eps/4)
    theta = -math.exp(eps / 8 * Lx)            # log theta = -x^(eps/8)
    lx1 = Lx + math.log1p(math.exp(-Lx)) if Lx < 700 else Lx   # log(x+1)
    cost1 = math.log(2) + f + lx1 + theta <= -Lx               # 2 f (x+1) theta <= 1/x
    cost2 = math.log(6) + f + lx1 + max(eta, 0) + math.log1p(math.exp(-eta)) <= Lb - Lx   # 6 f (x+1)(eta+1) <= b/x
    rest = (m >= lg(0.5 - eps)) and (k <= 2 * n + lg(0.5 - eps))
    sizew = (m <= lg(0.5 + eps / 2)) and (Lx + m <= k) and (lg(3 * eps / 8) >= math.log(26))
    # error bracket / (m nT b), each term in log
    # every bracket term and Main carry exactly one factor b: cancel it analytically (b = exp(x^eps) is ~e^(1e239) here)
    t1 = nT + math.log(2) + Lx + math.log1p(math.exp(2 * m - k)) - (m + nT)
    t2 = m + math.log(2) + n + m + lg(eps) - (m + nT)
    t3 = m + nT2 + k + lg(eps) - (2 * n + lg(-eps)) - (m + nT)
    main = Lx + math.log(La)                   # Main/(m nT b) ~ x log a
    rel = max(t1, t2, t3) - main               # log relative error
    return cost1, cost2, rest, sizew, (t1 / Lx, t2 / Lx, t3 / Lx), rel / Lx
for eps in (0.04, 0.02, 0.01):
    thr = None
    for Lx in [10 ** (j / 8) for j in range(8, 80) if eps * 10 ** (j / 8) < 600]:
        c1, c2, r, s, ts, rel = run(eps, Lx)
        if c1 and c2 and r and s and thr is None: thr = Lx
    Lbig = 550 / eps
    c1, c2, r, s, ts, rel = run(eps, Lbig)
    print(f"eps={eps}: all constraints hold from log x ~ {thr:.3g}; at log x={Lbig:.3g} term exponents (t1,t2,t3)/log x = "
          f"{tuple(round(t, 4) for t in ts)}; relative error exponent = {rel:.4f} (expected -1/2 + O(eps))")
# m-balance: error/main ~ 1/m + m/x  -> optimum m = x^(1/2); with per-lift cost ~ m^2/x -> optimum x^(1/3), error x^(-1/3)
for name, cost in (("linear cost (Lean bracket)", lambda mu: mu - 1), ("quadratic cost (hypothetical)", lambda mu: 2 * mu - 1)):
    best = min((max(-mu, cost(mu)), mu) for mu in [i / 1000 for i in range(1001)])
    print(f"{name}: optimal chain exponent mu = {best[1]:.3f}, relative error x^{best[0]:.3f} -> Chebyshev error exponent {1 + best[0]:.3f}")
