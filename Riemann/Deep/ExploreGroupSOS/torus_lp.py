# torus_lp.py: optimal envelope on the S-unit torus T^4 (c = 10, primes 2,3,5,7) by LP with cutting planes (not a certificate):
#  minimize beta s.t. beta + Phi(theta) >= 0, Phi = -2 sum_{n<c} w_n cos(theta.v_n) + 2 sum_{v in F} e_v cos(theta.v),
#  F = S-units in an exponent box whose log lands in the gap [L, P-L] mod P ('sep': pure prime powers only). Then lattice scan.
import sys, math, itertools, numpy as np
from scipy.optimize import linprog
from weil_lat import atoms, m_sym, hplus
rng = np.random.default_rng(1)
c = float(sys.argv[5]) if len(sys.argv) > 5 else 10.0
L = math.log(c); fac = float(sys.argv[1]); mode = sys.argv[2]; Bmax = float(sys.argv[3])
P = 2 * L * fac; pr = np.array([p for p in range(2, int(c) + 1) if p < c and all(p % q for q in range(2, int(p**0.5) + 1))]); lg = np.log(pr); D = len(pr)
fixed = []
for i, p in enumerate(pr):
    k = 1
    while p**k < c:
        v = np.zeros(D, int); v[i] = k; fixed.append((v, math.log(p) / p**(k / 2))); k += 1
bx = [int(x) for x in sys.argv[4].split(',')]
box = [range(-m, m + 1) for m in bx]
F = []
for a in itertools.product(*box):
    a = np.array(a); nz = np.count_nonzero(a)
    if nz == 0: continue
    # one representative per +-pair
    first = a[np.nonzero(a)[0][0]]
    if first < 0: continue
    if nz == 1 and float(np.dot(a, lg)) < L - 1e-12: continue      # fixed prime power inside the window
    if mode == 'sep' and nz > 1: continue
    if mode == 'none': continue
    lv = abs(float(np.dot(a, lg)))
    if L - 1e-12 <= (lv % P) <= P - L + 1e-12: F.append(a)
F = np.array(F) if F else np.zeros((0, D), int)
nF = len(F)
Vfix = np.array([v for v, _ in fixed]); Wfix = np.array([w for _, w in fixed])
def phi_fixed(th): return -2 * np.cos(th @ Vfix.T) @ Wfix
def feats(th): return 2 * np.cos(th @ F.T) if nF else np.zeros((len(th), 0))
S = rng.uniform(-math.pi, math.pi, size=(20000, D)); S = np.vstack([np.zeros((1, D)), S])
for it in range(6):
    A = np.hstack([-np.ones((len(S), 1)), -feats(S)]); b = phi_fixed(S)
    res = linprog(np.concatenate([[1.0], np.zeros(nF)]), A_ub=A, b_ub=b,
                  bounds=[(None, None)] + [(-Bmax, Bmax)] * nF, method='highs')
    beta = res.x[0]; e = res.x[1:]
    worst = []; vmin = 0
    for _ in range(100):
        T = rng.uniform(-math.pi, math.pi, size=(10000, D))
        val = beta + phi_fixed(T) + (feats(T) @ e if nF else 0)
        vmin = min(vmin, val.min()); worst.append(T[np.argsort(val)[:40]])
    S = np.vstack([S] + worst)
    if vmin > -1e-3: break
np.savez(f'/Users/adriansager/playground/riemann/.scratch/explore8/tlp_c{c:g}_{fac}_{mode}.npz', F=F, e=e, beta=beta, Vfix=Vfix, Wfix=Wfix, lg=lg, P=P, L=L)
env = 2 * beta  # in the units of the table (2*beta vs 2 sigma): Phi >= -beta  <=> symbol prime part >= -beta
at = atoms(c)
sig = beta + max(0.0, -vmin)
r = 20.0
while True:
    envr = hplus(r) - 2 * sum(math.exp(-(2*n+0.5)*L) * ((2*n+0.5) + r) / ((2*n+0.5)**2 + r*r) for n in range(40)) - 4 * math.exp(L/2) / (r - 0.5) - sig
    if envr > 0: break
    r *= 1.02
K = int(r * P / (2 * math.pi)) + 2
k = np.arange(0, K + 1); rr = 2 * math.pi * k / P
w = m_sym(rr, L, at)
for a, ee in zip(F, e): w = w + 2 * ee * np.cos(rr * float(np.dot(a, lg)))
neg = int(np.sum(w < 0)); lastneg = float(rr[w < 0].max()) if neg else 0
print(f'c={c:g} P={fac}*2L mode={mode} box={bx} Bmax={Bmax}: atoms {nF}, torus envelope beta = {beta:.4f} (sampled min residual {vmin:.1e}); '
      f'2sigma = {2*Wfix.sum():.4f}, ||C_p|| = 2.510; sum|e| = {np.abs(e).sum():.2f}; lattice: nodes {K}, negative {neg}, last neg r {lastneg:.0f}, min w {w.min():.3f}', flush=True)
