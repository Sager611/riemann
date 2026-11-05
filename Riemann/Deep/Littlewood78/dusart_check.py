# Numerical sanity check (not a proof) of the two cited Dusart bounds on y <= 10^8:
#   S(y) = sum_{p<=y} 1/p <= loglog y + M + 1/(10 log^2 y) + 4/(15 log^3 y)   (y >= 10372, Thm 6.11)
#   R(y) = sum_{p<=y} log p/p <= log y + E + 0.2/log y + 0.2/log^2 y         (y >= 2974, Thm 6.12)
# The worst y in each gap between primes is the prime itself (S, R jump up there, the right sides increase).
# Run: uv run --no-project --with numpy python3 dusart_check.py
import numpy as np, math
N = 10**8
s = np.ones(N + 1, dtype=bool); s[:2] = False
for p in range(2, int(N**0.5) + 1):
    if s[p]: s[p*p::p] = False
P = np.nonzero(s)[0].astype(np.float64); del s
M = 0.2614972128476427837554268386; E = -1.3325822757332208817658287760
S = np.cumsum(1.0/P); R = np.cumsum(np.log(P)/P)
lg = np.log(P)
r1 = (S - np.log(lg) - M)/(1/(10*lg**2) + 4/(15*lg**3))
r2 = (R - lg - E)/(0.2/lg + 0.2/lg**2)
m1 = P >= 10372; m2 = P >= 2974
i1 = np.argmax(r1[m1]); i2 = np.argmax(r2[m2])
print(f"Thm 6.11 ratio max over primes in [10372, 1e8]: {r1[m1][i1]:.4f} at p = {int(P[m1][i1])}  (must be <= 1)")
print(f"Thm 6.12 ratio max over primes in [2974, 1e8]:  {r2[m2][i2]:.4f} at p = {int(P[m2][i2])}  (must be <= 1)")
