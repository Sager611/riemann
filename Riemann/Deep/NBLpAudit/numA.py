# nblpaudit numerics: C_Gamma, Lemma Gamma(a), Lemma H / Gamma(b), A4 residue identity + B(eps), Lemma A shape.
# usage: python3 numA.py [quick]
import sys, math, time
import numpy as np
from flint import acb, arb, ctx
import mpmath as mp
ctx.prec = 64
QUICK = len(sys.argv) > 1
LOGPI = math.log(math.pi)
def zeta(x, t): return acb(x, t).zeta()
def absf(z): return abs(complex(z.real.mid(), z.imag.mid()))
def gr_abs(x, t, e):  # |gamma(s+e)/gamma(s)|, gamma(s)=s(s-1)/2 pi^{-s/2} Gamma(s/2)
    s = acb(x, t); s2 = acb(x + e, t)
    lg = (s2 / 2).lgamma() - (s / 2).lgamma()
    return math.exp(float(lg.real.mid()) - e / 2 * LOGPI) * absf(s2 * (s2 - 1)) / absf(s * (s - 1))
t0 = time.time()
# 1. C_Gamma
CG = math.e**4 * 2 * mp.quad(lambda t: (8 + mp.log(2 + t)) * mp.sqrt(2 + t) / (0.25 + t * t), [0, 1, 10, 100, 1e3, 1e4, 1e6, mp.inf])
print("C_Gamma = %.4f" % CG)
# 2. Lemma Gamma(a): excess = Re gamma'/gamma(w) - 0.5 log(1+|w|)
def rgg(x, t):
    w = acb(x, t); v = 1 / w + 1 / (w - 1) + (w / 2).digamma() / 2
    return float(v.real.mid()) - LOGPI / 2 - 0.5 * math.log(1 + abs(complex(x, t)))
ts = np.concatenate([np.arange(0, 60 if QUICK else 200, 0.05), [500, 1e3, 1e4, 1e5]])
m1 = max(rgg(x, t) for x in np.arange(0.5, 0.999, 0.05) for t in ts if not (abs(x - 1) < 1e-9 and t == 0))
m2 = max(rgg(x, t) for x in np.arange(0.5, 3.0, 0.05) for t in ts if t >= 1)
print("LemmaGamma(a): max excess Re w<1: %.4f (claim <=4); |Im w|>=1, Re w<=3: %.4f (claim <=5)" % (m1, m2))
# 3. Lemma H / Gamma(b) on lines
cases = [(1.1, 0.1), (1.3, 0.2), (8 / 7, 0.05)]
tg = np.arange(0.0, 60 if QUICK else 300, 0.1)
for p, e in cases:
    sp = 1 / p - e / 2
    rH = max(absf(zeta(sp, t) / zeta(sp + e, t)) / gr_abs(sp, t, e) for t in tg)
    rG = max(gr_abs(sp, t, e) ** 2 / (math.exp(8 * e) * (2 + t) ** e) for t in tg)
    print("p=%.4f eps=%.2f sigma'=%.4f: max|R|/|Gr|=%.4f (<=1)  max|Gr|^2/(e^{8e}(2+t)^e)=%.4f (<=1)" % (p, e, sp, rH, rG))
for e in (0.2, 0.05):  # sharp configuration theta*=1/2 (RH zeros): sigma = 1/2 - e/2 (+-0.01)
    out = []
    for d in (0.0, 0.01, -0.01):
        x = 0.5 - e / 2 + d
        out.append(max(absf(zeta(x, t) / zeta(x + e, t)) / gr_abs(x, t, e) for t in np.arange(1.0, 60, 0.02)))
    print("theta*=1/2, eps=%.2f: max|xi(s)/xi(s+e)| at sigma=1/2-e/2: %.6f, +0.01: %.4f, -0.01 (outside lemma): %.4f" % (e, *out))
print("time %.1f" % (time.time() - t0)); sys.stdout.flush()
