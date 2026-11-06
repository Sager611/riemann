#!/usr/bin/env python3
"""Parameter bookkeeping for BV Thm 2 without RH (Riemann/Deep/BVIrr).
Illustrative: BV's unspecified absolute constants c_0, C_0 (Prop exp-root-2) and our C_3 (Prop 9') are set to 1;
tau = 1/2, kappa = tau/100 as in BV sec. 8.1. All logs natural. Pure python, no packages."""
import math
tau = 0.5; kappa = tau / 100; c0 = 1.0; C0 = 1.0
def lgamma_bell_upper(m):          # log B_m <= m log m  (crude, valid m >= 1)
    return m * math.log(m)
print("== Prop 19' margin and BV's window (condiX), part (1): m = 1, X1 = tau (c0 d/100)^(1/2) / log d")
print("%6s %12s %12s %8s %12s %10s" % ("log10d", "X1", "d^0.1", "condiX", "margin", "delta_min"))
for k in [12, 15, 20, 26, 30, 40, 50, 60, 80, 100]:
    d = 10.0 ** k; L = math.log(d)
    X = tau * math.sqrt(c0 * d / 100) / L
    lo = 100 * max(1 / kappa, 1 / tau, d ** 0.1); hi = kappa * tau / C0 * d / (L) ** 3
    margin = X / 160 - (2 * 1 * L + 2 * math.log(X))            # need > 0 (up to log(16 C3))
    dmin = (2.1 * L + 2 * math.log(X)) / X                        # strip width at which the error is O(1)
    print("%6d %12.4g %12.4g %8s %12.4g %10.3g" % (k, X, d ** 0.1, "yes" if lo < X < hi else "no", margin, dmin))
print()
print("== part (2): m = ceil(30 (log d)^2) as BV l.2463 (>= Wielandt 6 log d), X2^2 = c0 tau^2 d / (100 m (log md)^2)")
print("%6s %8s %12s %12s %8s %12s %12s %10s" % ("log10d", "m", "X2", "100m^2d^.1", "condiX", "margin19'", "margin_w", "delta_min"))
for k in [30, 40, 50, 60, 70, 80, 100, 150, 200]:
    d = 10.0 ** k; L = math.log(d)
    m = math.ceil(30 * L ** 2)
    X = math.sqrt(c0 * tau ** 2 * d / (100 * m * math.log(m * d) ** 2))
    lo = 100 * m ** 2 * max(1 / kappa, 1 / tau, d ** 0.1); hi = kappa * tau / C0 * d / (m * math.log(m * d)) ** 3
    margin = X / 160 - (2 * m * L + 2 * math.log(X) + math.log(4 * m))
    margin_w = X / 32 - lgamma_bell_upper(m) - 2 * math.log(X)     # B_m |1 - w| small
    dmin = (2 * m * L + 2 * math.log(X)) / X
    print("%6d %8d %12.4g %12.4g %8s %12.4g %12.4g %10.3g" % (k, m, X, lo, "yes" if lo < X < hi else "no", margin, margin_w, dmin))
print()
print("== smallest d with d^(1/10) >= a log d (a = 640: Prop 19' with only X >= m d^(1/10); a = 6.4: with condiX)")
for a in [640, 6.4]:
    k = 1.0
    while 10 ** (k / 10) < a * k * math.log(10): k += 0.01
    print("a = %6.1f : d >= 10^%.2f" % (a, k))
print()
print("== Varju-Xu with the strip: X = A log d, error d^2 X^2 log d e^{-X/32} = d^{2 - A/32} (log d)^3 * A^2")
for A in [5, 64, 80, 81, 96]:
    print("A = %3d : exponent 2 - A/32 = %+.5f  (RH: 2 - A/2 = %+.2f)  l > A/eps" % (A, 2 - A / 32, 2 - A / 2))
print()
print("== Varju-Xu Cor 1.3 range f(d) >= 2^{a pi(d)}: a = 1/(1-T) for a half-plane Re s > T for Dirichlet L-functions")
for name, T in [("RH (VX's range, a = 2)", 0.5), ("Cureton-Thorner 31/32 (via quadratic fields)", 31 / 32),
                ("7/8 (OpenAI zeta; Hu Hecke, unrefereed)", 7 / 8), ("Bebel Lean T* = B_2 + 2.2e-254", 0.87338196661320752)]:
    print("%-48s T = %.17f  a = %.6f" % (name, T, 1 / (1 - T)))
print()
print("== BGS 2502.17218 exponential step: x = e^{c'n}; need eps < c' delta / 2 (delta = strip width)")
for delta, name in [(0.5, "RH"), (1 / 32, "CT")]:
    print("%-3s: eps < c'/%g" % (name, 2 / delta))
print("== Eberhard m = 1: X = theta*c*n; errors e^{-delta X} and e^{X - cn}; best theta = 1/(1+delta)")
for delta, name in [(0.5, "RH"), (1 / 32, "CT")]:
    th = 1 / (1 + delta)
    print("%-3s: X = %.4f c n, both errors e^{-%.4f c n}" % (name, th, delta * th))
