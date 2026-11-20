# ceilaudit B2: (1) the 7/11 DPP's form factor vs Montgomery's (closed form, checked by direct overlap integral);
# (2) freedom 'non-uniform density slope' for Theorem U': hairpin mass S(L) = sum over depth bins a (depth a/L) of e^{q theta a} * n(a),
#     n(a) = min over known density bounds of the count of zeros at depth >= a/L, in units of T*L (N(T) ~ TL/2pi).
import math
def F_dpp(alpha, a):   # discrete sine DPP on a*Z, unit intensity, occupancy rho = a; F = 1 - overlap/(2 pi rho) of arcs on circle 2pi
    rho = a; l = 2 * math.pi * rho; om = 2 * math.pi * abs(alpha) * a; om = om % (2 * math.pi); om = min(om, 2 * math.pi - om)
    ov = max(l - om, 0) + max(l - (2 * math.pi - om), 0)
    return 1 - ov / l
a = 7 / 11; th = 4 / 7
print("7/11-DPP form factor F(alpha) vs Montgomery min(|alpha|,1):")
for al in [0.1, 0.3, 0.5, 4 / 7 - 1e-9, 0.6, 0.7, 0.8, 0.9, 1.0]:
    print(f"  alpha={al:.4f}: F_DPP={F_dpp(al, a):.6f}  Montgomery={min(al, 1):.6f}")
print("  => agrees exactly on |alpha| <= theta = 1/a - 1 = 4/7; equals 4/7 on [4/7, 1] (BGSTB's complex-gamma F ~ alpha there is the one input it violates)")
# (2) density bounds near the line, depth d = a/L, as fractions of T*L (constants set to 1, generous to the freedom):
#  Conrey (Baluyot thm): T^{1 - (8/7) d} log T  -> e^{-8a/7}
#  Ingham: T^{3(1-s)/(2-s)} (log T)^5, s = 1/2 + d -> exponent 1 - 2d/(3/2 - d) ~ 1 - (4/3) d ; factor L^4 relative to TL
#  Huxley: T^{12(1-s)/5} (log T)^{44}: useful only for s > 3/4 here (exponent < 1 needs d > 1/12)
#  Guth-Maynard: T^{30(1-s)/13 + eps}: exponent < 1 for d > 1/30; treat T^eps as L^{10}
def nbound(a, L):   # log-space minimum
    d = a / L
    if d >= 0.5: return 0.0
    lb = [-8 * a / 7, 4 * math.log(L) - L * 2 * d / (1.5 - d), 43 * math.log(L) + L * (12 * (0.5 - d) / 5 - 1), 9 * math.log(L) + L * (30 * (0.5 - d) / 13 - 1), 0.0]
    return math.exp(min(lb))
def lognb(a, L):
    d = a / L
    if d >= 0.5: return -1e300
    return min(-8 * a / 7, 4 * math.log(L) - L * 2 * d / (1.5 - d), 43 * math.log(L) + L * (12 * (0.5 - d) / 5 - 1), 9 * math.log(L) + L * (30 * (0.5 - d) / 13 - 1), 0.0)
for th_, q in [(0.56, 2.02), (0.58, 2.01), (0.60, 2.01), (0.66, 2.01)]:
    row = []
    for L in [1e2, 1e3, 1e4, 1e6, 1e8]:
        best, amax = -1e300, 0.0; a_ = 0.0
        while a_ < L / 2:
            lt = q * th_ * a_ + lognb(a_, L)
            if lt > best: best, amax = lt, a_
            a_ += 0.05 if a_ < 200 else a_ * 0.002
        row.append(f"L={L:.0e}: log10 max_a = {best/math.log(10):.2f} at a={amax:.0f}")
    print(f"theta={th_}, q={q} (q*theta={q*th_:.4f} vs 8/7={8/7:.4f}): " + "; ".join(row))
print("max_a e^{q theta a} n(a) must stay O(1) (sum over off-line zeros of X^{q(beta-1/2)} << T log T); growth with L = log T means the freedom fails.")
