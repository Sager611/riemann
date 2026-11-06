# sanity check of the main-term normalisations (q=1, d=1, delta=1): B = sum_n M(n) tau_{-1/2}(n)/n,  S(v) = sum_n M(n) tau(n)/n^{1+iv}
import numpy as np, math, cmath, sys
sys.path.insert(0, '.')
xi = int(sys.argv[1]) if len(sys.argv) > 1 else 10**6
th = 0.1
N = xi
tau = np.ones(N+1)
# multiplicative sieve via smallest prime factor
spf = np.zeros(N+1, dtype=np.int64)
for p in range(2, N+1):
    if spf[p] == 0:
        spf[p::p][spf[p::p] == 0] = p
c = [1.0]
for k in range(1, 40): c.append(c[-1]*(k - 1.5)/k)
tau[0] = 0
for n in range(2, N+1):
    p = spf[n]; m = n; k = 0
    while m % p == 0: m //= p; k += 1
    tau[n] = tau[m]*c[k]
n = np.arange(N+1, dtype=float); n[0] = 1
L = math.log(xi)
M = np.where(n <= xi**th, 1.0, np.log(xi/n)/((1-th)*L)); M[0] = 0
B = np.sum(M*tau/n)
pred = (2/math.sqrt(math.pi))*(1-math.sqrt(th))/((1-th)*math.sqrt(L))
print("xi=%d  B=%.6f  predicted main term %.6f  ratio %.4f" % (xi, B, pred, B/pred))
from mpmath import mp, quad, sqrt, exp, mpf
def D2(Y):
    Y = mpf(Y)
    return quad(lambda x: -2/sqrt(x) + quad(lambda t: (exp(-1j*t)-1)/t**1.5, [0, x]), [0, Y]) if Y > 0 else 0
for X in [0.25, 0.5, 1.0]:
    v = X/L
    S = np.sum(M*tau/n*np.exp(-1j*v*np.log(n)))
    Dl = complex(D2(X)) - complex(D2(X*th))
    predS = -Dl/((1-th)*2*math.sqrt(math.pi)*math.sqrt(v)*L)
    print("X=%.2f  S=%s  predicted %s  |ratio| %.4f" % (X, S, predS, abs(S)/abs(predS)))

