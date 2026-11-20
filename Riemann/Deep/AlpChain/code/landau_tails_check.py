# Check of the tail lemma behind landau_tails (LandauTheorems/landauTails.lean:80-82) for delta = 1/2:
# exists sigma > 2, C with: for all U >= 1 there is k >= U such that
#   T(k,U) := int_{(0,inf) minus [U,1.5U]} u^k e^{(2-sigma)u}/k! du <= C / sigma^k.
# T(k,U) = (sigma-2)^{-k-1} * P(Gamma(k+1,1) not in [(sigma-2)U, 1.5(sigma-2)U]).
from mpmath import mp, mpf, gammainc, log, exp
mp.dps = 40
def tail_ratio(sigma, U, k):
    lam = sigma - 2
    lo, hi = lam * U, lam * mpf('1.5') * U
    pout = gammainc(k + 1, 0, lo, regularized=True) + gammainc(k + 1, hi, mp.inf, regularized=True)
    T = pout / lam**(k + 1)
    return T * mpf(sigma)**k   # must stay bounded in U (= C)
for sigma in (60, 110, 200):
    rows = []
    for U in (1, 3, 10, 30, 100, 300, 1000):
        lam = sigma - 2
        k = max(int(U) + 1, int(round(1.22 * lam * U)))   # centre of the window
        rows.append((U, k, float(log(tail_ratio(sigma, U, k)))))
    print(f"sigma={sigma}: log(T*sigma^k) at (U,k):", "; ".join(f"U={U},k={k}: {v:.2f}" for U, k, v in rows))
# Chernoff rate for the centred choice t=1.22: I(1/1.22) and I(1.5/1.22), I(t) = t - 1 - log t
I = lambda t: t - 1 - log(t)
print("rates:", float(I(1/mpf('1.22'))), float(I(mpf('1.5')/mpf('1.22'))), " threshold log(sigma/(sigma-2)) at sigma=110:", float(log(mpf(110)/108)))
