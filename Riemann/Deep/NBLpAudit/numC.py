# nblpaudit: Lemma P (Riesz-Perron) rate check. P_N(s) = sum_{k<=N} mu(k)(1-k/N)k^{-s} vs 1/zeta(s).
import math
import numpy as np
from flint import acb, ctx
ctx.prec = 64
NMAX = 10**5
mu = np.ones(NMAX + 1, dtype=np.int8); mu[0] = 0
isp = np.ones(NMAX + 1, dtype=bool); isp[:2] = False
for q in range(2, NMAX + 1):
    if isp[q]:
        isp[2 * q::q] = False; mu[q::q] *= -1; mu[q * q::q * q] = 0
k = np.arange(1, NMAX + 1, dtype=float); lk = np.log(k); muk = mu[1:].astype(float)
def PN(N, s): return np.sum(muk[:N] * (1 - k[:N] / N) * np.exp(-s * lk[:N]))
def invz(s): z = acb(s.real, s.imag).zeta(); return 1 / complex(z.real.mid(), z.imag.mid())
Ns = [10**2, 10**3, 10**4, 10**5]
for sig in (0.8, 0.95):
    for t in (0.0, 5.0, 14.134725, 50.0, 300.0):
        s = complex(sig, t); iz = invz(s)
        e = [abs(PN(N, s) - iz) for N in Ns]
        sl = [math.log(e[i] / e[i + 1]) / math.log(10) for i in range(3)]
        print("sigma=%.2f t=%9.4f |P_N-1/zeta| N=1e2..1e5: %s  local exps %s (RH: -> sigma-1/2=%.2f)" % (sig, t, " ".join("%.2e" % x for x in e), " ".join("%.2f" % x for x in sl), sig - 0.5))
for eps in (0.1, 0.2):  # case sigma-eta>1 at s=1+eps, eta=eps/2: bound (1/(1-eta)+1/eta) N^{-eta}
    eta = eps / 2; iz = invz(complex(1 + eps, 0))
    print("s=1+%.1f: " % eps + " ".join("N=%d err=%.2e bound=%.2e" % (N, abs(PN(N, complex(1 + eps, 0)) - iz), (1 / (1 - eta) + 1 / eta) * N ** (-eta)) for N in Ns))
