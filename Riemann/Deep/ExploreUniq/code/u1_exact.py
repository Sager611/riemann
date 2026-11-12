# explore37 u1: time-domain picture of Baez-Duarte (V1, L^2(0,inf)) after the unitary Mellin multiplier (s-1)/s:
#   chi_(0,1] -> e(x) = 1 + log x on (0,1];  rho_{1/k} -> -f(kx), f(x) = sum_{i<=1/x} (1 + log(i x)), supp f = (0,1].
#   d_N^2 = min_a || e - sum_{k<=N} a_k f(k.) ||^2_{L^2(0,1)}.
# On I_j = (1/(j+1), 1/j]: f(kx) = n(1 + log k + log x) + log n!, n = floor(j/k) (0 for j<k): 2-dim local space span{1, log x}.
# Lower bound (unconditional): for S not containing m, d_N(S)^2 >= eta_M(S)^2 := dist^2_{L^2(1/(M+1),1]}(e, span{f(k.): k in S, k<=M}), any M.
# usage: uv run --no-project --with mpmath --with numpy python3 u1_exact.py M
import sys, math, time, resource
import numpy as np, mpmath as mp
t0 = time.time()
M = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
mp.mp.dps = 40
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20  # MB on macOS (bytes)
# local orthonormalisation: variable l = log(x j) on I_j; Gram of {1, l}; upper Cholesky R_j
R00 = np.zeros(M + 1); R01 = np.zeros(M + 1); R11 = np.zeros(M + 1); LEN = np.zeros(M + 1)
for j in range(1, M + 1):
    a = mp.mpf(1) / (j + 1); b = mp.mpf(1) / j; Lj = mp.log(j)
    F1 = lambda x: x * (mp.log(x) + Lj) - x
    F2 = lambda x: (lambda lx: x * lx**2 - 2 * x * lx + 2 * x)(mp.log(x) + Lj)
    m0 = b - a; m1 = F1(b) - F1(a); m2 = F2(b) - F2(a)
    r00 = mp.sqrt(m0); r01 = m1 / r00; r11 = mp.sqrt(m2 - r01**2)
    R00[j], R01[j], R11[j], LEN[j] = float(r00), float(r01), float(r11), float(m0)
print(f'moments done {time.time()-t0:.1f}s', flush=True)
js = np.arange(1, M + 1)
lj = np.log(js.astype(float))
def block(alpha, beta, jj):  # coefficients in basis {1, l} -> orthonormal coords
    return np.stack([R00[jj] * alpha + R01[jj] * beta, R11[jj] * beta], axis=1).reshape(-1)
# e = 1 + log x = (1 - log j) + l on I_j
evec = block(1 - lj, np.ones(M), js)
V = np.zeros((2 * M, M))
lgam = np.array([math.lgamma(n + 1) for n in range(M + 2)])
for k in range(1, M + 1):
    jj = np.arange(k, M + 1); n = jj // k
    alpha = n * (1 + math.log(k) - np.log(jj.astype(float))) + lgam[n]; beta = n.astype(float)
    col = np.zeros(2 * M); col[2 * (k - 1):] = block(alpha, beta, jj); V[:, k - 1] = col
print(f'V built {V.shape} {time.time()-t0:.1f}s rss {rss():.0f} MB', flush=True)
np.save('/Users/adriansager/playground/riemann/.scratch/explore37/out/V_%d.npy' % M, V)
np.save('/Users/adriansager/playground/riemann/.scratch/explore37/out/e_%d.npy' % M, evec)
mu = np.ones(M + 1, dtype=np.int64); mu[0] = 0
for p in range(2, M + 1):
    if all(p % q for q in range(2, int(p ** 0.5) + 1)):
        mu[p::p] *= -1; mu[p * p::p * p] = 0
# check 1: Moebius exactness on (1/(M+1),1]: e = sum_{k<=M} mu(k) f(k.)
res = evec - V @ mu[1:].astype(float)
print(f'check Moebius exactness on (1/(M+1),1]: ||e - sum mu f(k.)|| = {np.linalg.norm(res):.3e} (||e||={np.linalg.norm(evec):.6f}, exact 1 - tail)', flush=True)
# check 2: Gram on (1/(M+1),1] vs Vasyunin values (tail O(log^2 M / M))
EG = 0.57721566490153286061
Gv = {(1, 1): math.log(2 * math.pi) - EG}
print(f'  <f,f>_(1/(M+1),1] = {V[:,0]@V[:,0]:.10f} ; Vasyunin G11 = {Gv[(1,1)]:.10f} ; diff*M/log^2M = {(Gv[(1,1)]-V[:,0]@V[:,0])*M/math.log(M)**2:.4f}')
print(f'  <e,f(k.)>, k=1..4: {[round(float(evec@V[:,k]),10) for k in range(4)]} ; b_k = {[round((math.log(k)+1-EG)/k,10) for k in range(1,5)]}')
print(f'  ||e||^2 on (1/(M+1),1] = {evec@evec:.10f} (1 - int_0^(1/(M+1)) (1+log x)^2 = {1 - (lambda x: x*(1+math.log(x))**2 - 2*x*math.log(x))(1/(M+1)):.10f})', flush=True)
def eta2(S, Mc):  # dist^2 on (1/(Mc+1),1] of e from span{f(k.): k in S, k <= Mc}
    rows = 2 * Mc; cols = [k - 1 for k in S if k <= Mc]
    if not cols: return float(evec[:rows] @ evec[:rows])
    A = V[:rows, cols]; y = evec[:rows]
    coef, *_ = np.linalg.lstsq(A, y, rcond=None)
    r = y - A @ coef
    return float(r @ r)
out = []
sqf = [m for m in range(1, M + 1) if mu[m] != 0]
# (a) S = N \ {m0}, squarefree m0: eta^2 = mu(m0)^2 delta'^2 at Mc = m0
lim = min(M, 300)
print('m0  eta2(N\\{m0}) [Mc=m0]   m0^2*eta2   m0^3*eta2', flush=True)
tab = []
for m0 in [m for m in sqf if m <= lim] + [m for m in (437, 501, 1001, 1499, 1999) if m <= M and mu[m] != 0]:
    S = [k for k in range(1, m0 + 1) if k != m0]
    v = eta2(S, m0); tab.append((m0, v))
    if m0 <= 30 or m0 in (47, 97, 101, 199, 211, 299, 293, 437, 501, 1001, 1499, 1999) or m0 == lim:
        print(f'{m0:5d} {v:.6e}  {m0**2*v:.5f}  {m0**3*v:.4f}', flush=True)
np.save('/Users/adriansager/playground/riemann/.scratch/explore37/out/eta_single_%d.npy' % M, np.array(tab))
# fit log eta2 vs log m0 for m0 in [20, lim]
arr = np.array([t for t in tab if 20 <= t[0] <= lim]); pfit = np.polyfit(np.log(arr[:, 0]), np.log(arr[:, 1]), 1)
print(f'fit eta2 ~ {math.exp(pfit[1]):.4f} m0^{pfit[0]:.4f} on [20,{lim}]', flush=True)
# (b) structured S: eta_M(S)^2 increasing in Mc
primes = [p for p in range(2, M + 1) if all(p % q for q in range(2, int(p ** 0.5) + 1))]
fam = {'N\\{1}': [k for k in range(2, M + 1)], 'N\\{2}': [k for k in range(1, M + 1) if k != 2],
       'odd': [k for k in range(1, M + 1, 2)], '{1}+primes': [1] + primes,
       'squarefree': sqf, 'N\\{4}': [k for k in range(1, M + 1) if k != 4]}
rng = np.random.default_rng(37)
for sd in range(3):
    fam[f'rand1/2 #{sd}'] = [1] + [k for k in range(2, M + 1) if rng.random() < 0.5]
print('family        eta_Mc^2 for Mc = 1,2,3,6,10,30,100,300' + (',1000' if M >= 1000 else ''), flush=True)
for nm, S in fam.items():
    vals = [eta2(S, Mc) for Mc in (1, 2, 3, 6, 10, 30, 100, 300, 1000) if Mc <= M]
    print(f'{nm:14s} ' + ' '.join(f'{v:.5e}' for v in vals), flush=True)
    out.append((nm, vals))
import json
json.dump({'fam': out, 'single': tab}, open('/Users/adriansager/playground/riemann/.scratch/explore37/out/u1_%d.json' % M, 'w'))
json.dump({k: v for k, v in fam.items()}, open('/Users/adriansager/playground/riemann/.scratch/explore37/out/fam_%d.json' % M, 'w'))
print(f'done {time.time()-t0:.1f}s peak rss {rss():.0f} MB', flush=True)
