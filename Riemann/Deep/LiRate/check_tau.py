# Checks of the Theorem 1 ingredients (numerical sanity, not proofs):
# (a) tau(theta) = sum_n phi(n/N) e^{-in theta} t(|n|) ~ log(2 pi theta)/theta^2 (phi = psi*psi, psi smooth bump);
# (b) Phi(x) = sum_n phi(n/N) e^{inx} >= 0 and Phi(0) ~ N int phi;
# (c) the two-sided weighted mean square ||R~||^2 vs the Theorem 1 bound Phi(0) N(T) (1 - tau_T/(Phi(0)N(T)))^2.
import math, numpy as np
g = 0.57721566490153286; kap = 1 + math.log(2 * math.pi) - g
def tr(m): return np.where(m > 0, m / 2 * (np.log(np.maximum(m, 1)) - kap), 0.0)
xs = np.linspace(-0.5, 0.5, 4001)[1:-1]; h = xs[1] - xs[0]
psi = np.exp(-1 / (1 - 4 * xs ** 2)); psi /= math.sqrt(np.sum(psi ** 2) * h)
phig = np.convolve(psi, psi) * h; xg = np.linspace(-1, 1, len(phig))
intphi = np.sum(phig) * (xg[1] - xg[0])
print(f'int phi = {intphi:.4f}, phi(0) = {phig.max():.4f}')
def phi(u): return np.interp(u, xg, phig, left=0, right=0)
lam = {}
for l in open('logs/lam_r999.txt'):
    if l.startswith('#'): continue
    a, b = l.split()[:2]; lam[int(a)] = float(b)
for N in (20000, 100000):
    n = np.arange(-N, N + 1); w = phi(n / N); tn = tr(np.abs(n).astype(float))
    print(f'N={N}: Phi(0)/N = {w.sum()/N:.4f}')
    for th in (0.003, 0.01, 0.03, 0.0707):
        tau = np.sum(w * tn * np.cos(n * th)); pred = math.log(2 * math.pi * th) / th ** 2
        Phi = np.sum(w * np.cos(n * th))
        print(f'  theta={th}: tau={tau:.2f}  log(2pi th)/th^2={pred:.2f}  diff={tau-pred:.2f}   Phi(theta)={Phi:.3e}')
    xx = np.linspace(0, 0.2, 2001); Pm = min(np.sum(w * np.cos(n * x)) for x in xx)
    print(f'  min Phi on [0,0.2] grid: {Pm:.3e}')
    Rt = np.array([lam[abs(k)] for k in n]) - tn
    lhs = np.sum(w * Rt ** 2)
    best = 0
    for T in np.linspace(20, 400, 381):
        L = math.log(T / (2 * math.pi)); NT = T / (2 * math.pi) * (L - 1) + 7 / 8
        tauT = T ** 3 * L ** 2 / (6 * math.pi); P0 = w.sum()
        if tauT < P0 * NT: best = max(best, P0 * NT * (1 - tauT / (P0 * NT)) ** 2)
    print(f'  sum phi R~^2 = {lhs:.1f};  Thm-1 bound (best T) = {best:.1f};  ratio {lhs/best:.2f}')
