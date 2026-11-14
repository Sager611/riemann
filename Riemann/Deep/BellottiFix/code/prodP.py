# Rigorous upper bound for P = prod_p (1 + (sqrt(p) - 1)/p^2), used in |G(1/2+it)| <= 16 P W / log^2 W.
# Primes p <= N by sieve; float64 terms summed with math.fsum (error << 1e-9, added as margin);
# tail: sum_{p>N} (sqrt p - 1)/p^2 < sum_{p>N} p^(-3/2) <= (3/2)(1.25506)(2/(sqrt(N) log N))  [pi(u) <= 1.25506 u/log u, Rosser-Schoenfeld].
import math, numpy as np
N = 10**7
s = np.ones(N + 1, dtype=bool); s[:2] = False
for i in range(2, int(N**0.5) + 1):
    if s[i]: s[i*i::i] = False
p = np.nonzero(s)[0].astype(np.float64)
terms = np.log1p((np.sqrt(p) - 1.0) / (p * p))
logP = math.fsum(terms.tolist())
tail = 1.5 * 1.25506 * 2 / (math.sqrt(N) * math.log(N))
Pup = math.exp(logP + tail + 1e-9)
print('primes', len(p), 'logP(p<=N) %.12f' % logP, 'tail <= %.3e' % tail)
print('P in [%.8f, %.8f]' % (math.exp(logP), Pup))
print('16P <= %.6f' % (16 * Pup))

