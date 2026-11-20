# (L6) bookkeeping: from gapA <= x^(C eps) * errorBracketSizeWindows (ChTop.lean:37-44, ErrorFormsWeakWindowsFalse.lean:133-135)
# to G(x) <= x^(1/2+O(eps)), with the exact constants of RestSize (SChain.lean:53-58), ch_core (ChCore.lean:343-350)
# and the explicit sizes of SizesTextAt (SChain.lean:207-222). Normalisation of the thick part: N0 = m*nT*b*nT2/(2n)
# (|Phi| = 2n, |T2| thick places, w ~ b, a_t/w ~ x).
from mpmath import mp, mpf, log, exp, sqrt
from fractions import Fraction as Fr
mp.dps = 30
def EB(x, eps, n, nT, nT2, m, b, k, aE):
    return nT * aE * (1 + m**2 / k) + m * ((2 * n) * m * b * x**eps + nT2 * k * b * x**eps / (n * (n * x**(-eps))))
# 1. symbolic exponents (in units of eps above 1/2) of the three terms of EB/N0 under the worst case allowed by the hypotheses
terms = {"T1 = 2n nT aE(1+m^2/k)/(m nT b nT2)  [aE<=2xb, n/nT2<=x^eps, m>=x^(1/2-eps)]": (4, Fr(2)),
         "T2 = 4 n^2 m x^eps/(nT nT2)        [nT>=x^(-eps/8)n/2, nT2>=x^(-eps)n, m<=x^(1/2+eps/2)]": (8, Fr(1) + Fr(1, 8) + 1 + Fr(1, 2)),
         "T3 = 2 k x^(2eps)/(n nT)           [k<=n^2 x^(1/2-eps), nT>=x^(-eps/8)n/2]": (4, Fr(2) + Fr(1, 8) - 1)}
for t, (K, e) in terms.items(): print(f"{t}: <= {K} x^(1/2 + {e} eps)")
emax = max(e for _, e in terms.values()); print(f"worst exponent 1/2 + {emax} eps = 1/2 + {float(emax)} eps; total <= 20 x^(1/2+(C+{emax})eps) (T1 <= 4x^(1/2+2eps) + 4x^(5eps/2)) <= x^(1/2+(C+3)eps) once x^(3eps/8) >= 20 (Lean ch_core asks x^(3eps/8) >= 26)")
# 2. numbers at the explicit SizesTextAt parameters (log scale)
for eps, L in ((mpf('0.04'), mpf(3000)), (mpf('0.02'), mpf(5000)), (mpf('0.01'), mpf(12000))):
    x = exp(L); n = x**(2 + 8 * eps); m = x**(mpf(1) / 2 - eps / 16); b = exp(x**eps); a = x * b; k = x**(mpf(9) / 2 + 10 * eps)
    nT = x**(-eps / 8) * n / 2; nT2 = x**(-eps) * n; aE = a
    N0 = m * nT * b * nT2 / (2 * n)
    val = EB(x, eps, n, nT, nT2, m, b, k, aE) / N0
    parts = [nT * aE * (1 + m**2 / k) / N0, m * 2 * n * m * b * x**eps / N0, m * nT2 * k * b * x**eps / (n * n * x**(-eps)) / N0]
    ex = lambda v: float((log(v) / L - mpf(1) / 2) / eps)
    rs = [x**(mpf(1)/2 - eps) <= m, m <= x**(mpf(1)/2 + eps/2), x * m <= k, k <= n**2 * x**(mpf(1)/2 - eps), aE <= 2 * x * b]
    print(f"eps={float(eps)} log x={float(L)}: RestSize+ch_core size hyps hold: {all(rs)}; EB/N0 = x^(1/2 + {ex(val):.3f} eps); parts x^(1/2 + [{', '.join(f'{ex(p):.3f}' for p in parts)}] eps)")
# 3. optimisation in m with the Lean cost model (B0 = 2 n m b x^eps) and the alternative (per-lift B0 = 2 n m^2 b x^eps)
for model in ("Lean: B0 ~ m", "alt: B0 ~ m^2"):
    for L in (mpf(200), mpf(2000), mpf(20000)):
        x = exp(L)
        f = (lambda m: 2 * x / m + 4 * m) if model.startswith("Lean") else (lambda m: 2 * x / m + 4 * m**2)
        lo, hi = mpf(0), L   # golden search on log m
        for _ in range(200):
            m1, m2 = lo + (hi - lo) / 3, hi - (hi - lo) / 3
            if f(exp(m1)) < f(exp(m2)): hi = m2
            else: lo = m1
        mo = exp((lo + hi) / 2)
        print(f"{model}: log x={float(L)}: optimal m = x^{float(log(mo)/L):.4f}, min relative-error x * (2/m + cost) = x^{float(log(f(mo))/L):.4f}")
