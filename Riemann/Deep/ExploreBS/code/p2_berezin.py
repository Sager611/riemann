# p2_berezin.py (explore57): phase-space (Weyl / Berezin-Lieb / CLR-type) bounds for Weil's form with its natural symbol.
# Q(f) = pole + (1/2pi) int sigma_c(t) |F(t)|^2 dt on L^2 of a window of length ell = ln c (wmat.py normalisation), where
# sigma_c(t) = Omega(t) - X_c(t), Omega(t) = Re psi(1/4 + it/2) - ln pi, X_c(t) = 2 sum_{n<c} Lambda(n) n^{-1/2} cos(t ln n).
# Berezin-Lieb (Tr phi(PAP) <= Tr P phi(A) P, phi convex, phi(0) = 0; phi = negative part):
#   sum_j lambda_j(T_c)_-  <=  B(c) := (ell/2pi) int_R (sigma_c)_-(t) dt     (T_c = pole-free part)
# Weyl count W(c) := (ell/2pi) |{t : sigma_c(t) < 0}|.  Bohr-Jessen (t log p equidistributed):  both ~ 2 ell E[exp X_c],
# E[exp X_c] = prod_p (1/2pi) int exp(2 ln p sum_{k: p^k < c} p^{-k/2} cos k th) dth.   Truth: N(Q_c) = 0 for c <= 47 (certified).
import sys, math, numpy as np
from scipy.special import digamma
def atoms(c):
    out = []
    for n in range(2, int(math.ceil(c))):
        m = n; p = None
        for q in range(2, n + 1):
            if m % q == 0: p = q; break
        while m % p == 0: m //= p
        if m == 1: out.append((n, math.log(p)/math.sqrt(n)))
    return out
def EexpX(c):
    tot = 1.0; th = np.linspace(0, 2*np.pi, 20001)[:-1]
    primes = sorted({n for n, w in atoms(c) if all(n % q for q in range(2, int(n**0.5) + 1))})
    for p in primes:
        s = np.zeros_like(th); k = 1
        while p**k < c:
            s += 2*math.log(p)*p**(-k/2)*np.cos(k*th); k += 1
        tot *= np.exp(s).mean()
    return tot
def grid(c, h=0.005):
    A = atoms(c); maxX = 2*sum(w for n, w in A)
    Tmax = 2*np.pi*math.exp(maxX + 0.5); out_m = 0.0; out_b = 0.0
    for a in np.arange(0, Tmax, 2000.0):
        t = np.arange(a, min(a + 2000.0, Tmax), h) + h/2
        om = digamma(0.25 + 0.5j*t).real - math.log(math.pi)
        X = np.zeros_like(t)
        for n, w in A: X += 2*w*np.cos(t*math.log(n))
        sg = om - X
        out_m += 2*h*np.count_nonzero(sg < 0); out_b += 2*h*np.clip(-sg, 0, None).sum()
    return out_m, out_b, maxX, Tmax
for c in [float(x) for x in sys.argv[1].split(',')]:
    ell = math.log(c); E = EexpX(c); A = atoms(c)
    line = f'c={c:g} ell={ell:.4f} atoms={len(A)} maxX={2*sum(w for n, w in A):.3f} E[expX]={E:.4g} BJ: W~B~2 ell E = {2*ell*E:.4g}'
    if c <= float(sys.argv[2]):
        m, b, mx, T = grid(c)
        line += f' | grid (|t| < {T:.3g}): |sigma<0| = {m:.4g}, W = {ell*m/(2*np.pi):.4g}, B = {ell*b/(2*np.pi):.4g}'
    print(line, flush=True)
