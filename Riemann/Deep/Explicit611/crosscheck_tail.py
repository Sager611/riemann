# Float cross-check (not part of the certificate): direct per-block sums vs the closed-form tail bound of certify.py.
# For t in [1e10, 1e20], K = ceil(t^{1/3}), blocks r = 4..R (last one counted as full), each block = min(relaxed, exact Weyl).
# Run: uv run --no-project --with python-flint --with numpy python3 crosscheck_tail.py
import math, os, sys
sys.argv = ['certify.py', '0.568', '', '500', '39894', '1e10', '6000']        # NEX > 0 -> exact-Weyl tail (as in the certificate)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import certify as cf
from flint import arb
def direct(t):
    K = math.ceil(t ** (1 / 3)); n1 = math.floor(math.sqrt(t / (2 * math.pi))); R = n1 // K
    T1 = 2 * float((arb(0.5).zeta() - arb(0.5).zeta(arb(4 * K))).mid())     # 2 sum_{n <= 4K-1} n^{-1/2}
    s = sum(min(cf.blk_f(t, K, r, K)[0], cf.blk_exact_f(t, K, r, K)[0]) for r in range(4, R + 1))
    return T1 + s + 1.463 * t ** -0.25 + 0.127 * t ** -0.75, K, R
worst = -1
for Ts, tl in (("1e10", [1e10, 3e10, 1e11, 1e12]), ("1e13", [1e13, 1e14, 1e15]), ("8.97e17", [8.97e17, 2e18, 1e19, 1e20])):
    q, r1, ct, b1, b2, b3 = cf.best_tail(Ts)
    b1, b2, b3 = float(b1.upper()), float(b2.upper()), float(b3.upper())
    for t in tl:
        d, K, R = direct(t)
        f = b1 * t ** (1 / 6) * math.log(t) + b2 * t ** (1 / 6) + b3
        rr = t ** (1 / 6) * math.log(t)
        worst = max(worst, d / f)
        print("T=%-8s t=%.3g K=%d R=%d  direct ratio %.5f <= formula ratio %.5f  %s" % (Ts, t, K, R, d / rr, f / rr, "ok" if d <= f else "VIOLATION"))
print("max direct/formula = %.5f (must be <= 1)" % worst)
