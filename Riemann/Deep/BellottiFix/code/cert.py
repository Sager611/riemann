# Certified (ball arithmetic, python-flint arb) constants for the corrected |G(1/2+it)| bound.
# (1) P = prod_p (1 + (sqrt p - 1)/p^2): exact product over p <= 10^7 in arb, tail sum_{p>N} p^(-3/2) <= 3.765/(sqrt N log N)
#     (pi(u) <= 1.25506 u/log u, Rosser-Schoenfeld 1962, u > 1) and log(1+x) <= x.
# (2) c_tail = 24/pi^2 + 3 (L + L^2/4)/sqrt(Y0), L = log Y0, Y0 = 2*10^6, bounds S(y)/sqrt(y) for y >= Y0 (NOTE section 2).
# (3) c_small: sqfree.py's float64 sup over y <= Y0 (2.429643451) has relative error < 1e-9 (cumulative sums of < 1.3e6 positive
#     terms, n*u < 1.5e-10; log/sqrt/exp 1 ulp), so c_small < 2.4297 < c_tail: c = c_tail.
import numpy as np
from flint import arb, ctx
ctx.prec = 80
N = 10**7
s = np.ones(N + 1, dtype=bool); s[:2] = False
for i in range(2, int(N**0.5) + 1):
    if s[i]: s[i*i::i] = False
P = arb(1)
for p in np.nonzero(s)[0].tolist():
    P *= 1 + ((arb(p).sqrt() - 1) / (p * p))
tail = arb('3.765') / (arb(N).sqrt() * arb(N).log())
Pup = P * tail.exp()
Y0 = arb(2 * 10**6); L = Y0.log()
ctail = 24 / arb.pi()**2 + 3 * (L + L * L / 4) / Y0.sqrt()
G = ctail * ctail * Pup
print('P(p<=1e7) =', P.str(12), ' P <=', Pup.upper().str(10))
print('c_tail =', ctail.str(12), ' c_small(float) = 2.429643451 < c_tail')
print('c^2 P <=', G.upper().str(10), ' (used 9.7572):', bool(G < arb('9.7572')))
print('16 P <=', (16 * Pup).upper().str(10))
