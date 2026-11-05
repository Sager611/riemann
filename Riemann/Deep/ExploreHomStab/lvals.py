# explore9: central values L(1/2, chi_{8d}) for odd squarefree d <= X (Soundararajan's family).
# Approximate functional equation (even real primitive chi mod q = 8d, root number 1):
#   L(1/2, chi) = 2 * sum_{n odd} chi(n) n^{-1/2} Q(1/4, pi n^2 / q),  Q = regularized upper incomplete gamma.
# chi_{8d}(p) = (2/p) (d/p) for odd primes p; chi(n) by complete multiplicativity; chi(even) = 0.
# Run: OMP_NUM_THREADS=1 uv run --no-project --with numpy --with scipy python3 lvals.py X out.npz
import sys, time
import numpy as np
from scipy.special import gammaincc

X = int(sys.argv[1]); out = sys.argv[2]
CUT = 60.0  # pi n^2/q > CUT  => Q(1/4, x) < 1e-27
t0 = time.time()

# odd squarefree d <= X
sq = np.ones(X + 1, dtype=bool); sq[0] = False
for p in range(3, int(X ** 0.5) + 2, 2):
    sq[p * p::p * p] = False
d = np.arange(X + 1)[sq & (np.arange(X + 1) % 2 == 1)].astype(np.int64)
ND = len(d)

NMAX = int(np.sqrt(CUT * 8 * X / np.pi)) + 2
# primes and smallest prime factor up to NMAX
spf = np.zeros(NMAX + 1, dtype=np.int64)
for i in range(2, NMAX + 1):
    if spf[i] == 0:
        spf[i::i][spf[i::i] == 0] = i
odd = np.arange(1, NMAX + 1, 2)            # odd n
idx = {int(n): k for k, n in enumerate(odd)}
primes = [p for p in range(3, NMAX + 1, 2) if spf[p] == p]

# Q(1/4, x) table in t = x^{1/4} (analytic in t)
tt = np.linspace(0.0, CUT ** 0.25, 2_000_001)
Qtab = gammaincc(0.25, tt ** 4)

L = np.zeros(ND)
B = 16000
for b0 in range(0, ND, B):
    db = d[b0:b0 + B]; nb = len(db)
    chip = {}
    for p in primes:
        T = np.zeros(p, dtype=np.int8)
        T[(np.arange(1, p) ** 2) % p] = 1
        T[1:][T[1:] == 0] = -1
        two = 1 if p % 8 in (1, 7) else -1
        chip[p] = (two * T[db % p]).astype(np.int8)
    chi = np.zeros((len(odd), nb), dtype=np.int8)
    chi[0] = 1
    S = np.zeros(nb)
    inv_q = np.pi / (8.0 * db)
    for k in range(len(odd)):
        n = int(odd[k])
        if n > 1:
            p = int(spf[n]); m = n // p
            chi[k] = chip[p] * chi[idx[m]]
        x = n * n * inv_q
        t = np.minimum(x, CUT) ** 0.25
        w = np.interp(t, tt, Qtab) / np.sqrt(n)
        S += chi[k] * w
    L[b0:b0 + nb] = 2.0 * S
    print(f"block {b0//B}: {b0+nb}/{ND} done, {time.time()-t0:.1f}s", flush=True)

np.savez(out, d=d, L=L)
print("ND", ND, "NMAX", NMAX, "time", round(time.time() - t0, 1))
