# t3_davenport.py — explore6 T3: the rotation family (contains the constant at alpha = 0).
# sup_alpha |S_x(alpha)|, S_x(alpha) = sum_{n<=x} mu(n) e(n alpha), on a 4x-oversampled grid (rfft), x = 2^k.
# Compares with |M(x)| = |S_x(0)| and with the proven exponents: Davenport (log-power), textbook split under a
# zero-free half-plane Re s > Theta for all Dirichlet L: (2+Theta)/3; Baker-Harman under GRH: 3/4.
import numpy as np, math, time
from fractions import Fraction
t0 = time.time()
LOGX = 22; X = 1 << LOGX
def mobius(N):
    mu = np.ones(N + 1, np.int8); acc = np.zeros(N + 1, np.float32)
    R = int(math.isqrt(N)) + 1
    sv = np.ones(R + 1, bool); sv[:2] = False
    for i in range(2, int(R**0.5) + 1):
        if sv[i]: sv[i*i::i] = False
    for p in np.nonzero(sv)[0].tolist():
        mu[::p] *= -1; acc[::p] += math.log(p); mu[::p*p] = 0
    n = np.arange(N + 1); big = acc < (np.log(np.maximum(n, 1)) - 1.0); mu[big] *= -1
    mu[0] = 0
    return mu
mu = mobius(X)
B2 = (777 - math.sqrt(1081)) / 852
print(f"# exponents: B_2 = {B2:.6f}; textbook split (2+B_2)/3 = {(2+B2)/3:.5f}; GRH textbook 5/6 = {5/6:.4f}; Baker-Harman GRH 3/4")
print("k   x        |M(x)|   sup|S|     sup/sqrt(x)  sup/sqrt(x log x)  alpha*        nearest a/q (q<=64)   |S| at 1/2, 1/3, 1/4, 1/6")
for k in (12, 14, 16, 18, 20, 22):
    x = 1 << k
    a = mu[1:x + 1].astype(np.float64)
    L = 4 * x
    F = np.fft.rfft(np.concatenate(([0.0], a)), L)      # F[j] = sum_{n=0}^{x} a_n e(-n j/L)
    absF = np.abs(F)
    j = int(np.argmax(absF)); s = absF[j]; al = j / L
    fr = Fraction(al).limit_denominator(64)
    def at(r):
        jj = int(round(r * L)); return absF[jj]
    print(f"{k:2d}  {x:<8d} {abs(a.sum()):7.0f}  {s:9.1f}  {s/math.sqrt(x):9.3f}    {s/math.sqrt(x*math.log(x)):9.3f}       {al:.6f}    {fr.numerator}/{fr.denominator} (off {abs(al-float(fr))*x:.2f}/x)   " +
          " ".join(f"{at(r)/math.sqrt(x):.3f}" for r in (0.5, 1/3, 0.25, 1/6)))
    del F, absF
print(f"seconds {time.time()-t0:.1f}")
