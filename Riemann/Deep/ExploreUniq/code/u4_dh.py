# explore37 u4 (adversarial): Davenport-Heilbronn analogue of the triangular Baez-Duarte picture.
# L_DH(s) = sum a_n n^-s (a periodic mod 5: 1, k, -k, -1, 0), entire; Mellin^{-1}[L(s)/s] = A(floor(1/x)), A(n) = sum_{i<=n} a_i,
# supported on (0,1]. BD-type problem: chi_(0,1] ~ sum_k c_k A(floor(1/(kx))); on I_j: A(floor(j/k)) (piecewise constant).
# Moebius analogue: chi = sum_k b_k A(floor(1/(k.))) on (1/(M+1),1], b = Dirichlet inverse of a. Floors eta_M(S)^2 as for zeta.
import math, sys, json, resource
import numpy as np
M = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
kap = (math.sqrt(10 - 2 * math.sqrt(5)) - 2) / (math.sqrt(5) - 1)
a = np.zeros(M + 1); per = [0.0, 1.0, kap, -kap, -1.0]  # index n mod 5 -> a_n (n%5==0 -> 0)
for n in range(1, M + 1): a[n] = per[n % 5]
A = np.cumsum(a)  # A[n] = sum_{i<=n} a_i
b = np.zeros(M + 1); b[1] = 1.0
for n in range(2, M + 1):
    s = 0.0
    for d in range(1, n):
        if n % d == 0: s += a[n // d] * b[d]
    b[n] = -s
w = np.sqrt(1.0 / (np.arange(1, M + 1) * (np.arange(1, M + 1) + 1.0)))  # sqrt|I_j|
js = np.arange(1, M + 1)
V = np.zeros((M, M))
for k in range(1, M + 1):
    jj = np.arange(k, M + 1); V[jj - 1, k - 1] = A[jj // k] * w[jj - 1]
chi = w.copy()
print(f'kappa={kap:.10f}; Moebius-analogue exactness: ||chi - V b|| = {np.linalg.norm(chi - V @ b[1:]):.2e}')
ess = [n for n in range(1, 61) if abs(b[n]) > 1e-12]
non = [n for n in range(1, 61) if abs(b[n]) <= 1e-12]
print('b_n, n=1..12:', ' '.join(f'{b[n]:+.4f}' for n in range(1, 13)))
print('essential n<=60 (b_n != 0):', ess)
print('non-essential n<=60 (b_n = 0):', non)
mu = np.ones(M + 1, dtype=np.int64); mu[0] = 0
for p in range(2, M + 1):
    if all(p % q for q in range(2, int(p ** 0.5) + 1)):
        mu[p::p] *= -1; mu[p * p::p * p] = 0
print('zeta essential (squarefree) n<=60 not DH-essential:', [n for n in range(1, 61) if mu[n] != 0 and abs(b[n]) <= 1e-12])
print('DH-essential n<=60 not squarefree:', [n for n in range(1, 61) if mu[n] == 0 and abs(b[n]) > 1e-12])
def eta2(S, Mc):
    cols = [k - 1 for k in S if k <= Mc]; Am = V[:Mc, cols]; y = chi[:Mc]
    c, *_ = np.linalg.lstsq(Am, y, rcond=None); r = y - Am @ c
    return float(r @ r)
primes = [p for p in range(2, M + 1) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
fam = {'N\\{1}': list(range(2, M + 1)), 'N\\{2}': [k for k in range(1, M + 1) if k != 2],
       'N\\{4}': [k for k in range(1, M + 1) if k != 4], 'N\\{5}': [k for k in range(1, M + 1) if k != 5],
       '{1}+primes': [1] + primes, 'squarefree': [k for k in range(1, M + 1) if mu[k] != 0],
       'supp b': [k for k in range(1, M + 1) if abs(b[k]) > 1e-12]}
print('DH family   eta_Mc^2 at Mc = 10, 100,', M)
out = {}
for nm, S in fam.items():
    v = [eta2(S, Mc) for Mc in (10, 100, M)]; out[nm] = v
    print(f'{nm:12s} ' + ' '.join(f'{x:.5e}' for x in v))
print('single removals m in supp b: m, |b_m|, eta^2(N\\{m}) at Mc=m, m^2 eta^2')
for m in (1, 2, 3, 4, 6, 7, 8, 9, 11, 16, 101, 199):
    if m <= M:
        v = eta2([k for k in range(1, m + 1) if k != m], m)
        print(f'{m:4d} {abs(b[m]):.4f} {v:.5e} {m*m*v:.4f}')
json.dump(out, open('/Users/adriansager/playground/riemann/.scratch/explore37/out/u4_dh_%d.json' % M, 'w'))
print(f'peak rss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20:.0f} MB')
