# Smallest-case / random falsification of the explicit inequalities: prop:optimal (C(f_l)=C_l, coercivity constant 1-2l^2/pi^2,
# admissible eta_eps), lem:coeff sums, Lamzouri's Prop 2.1 (eq:finite) on random conjugation-invariant multisets.
import numpy as np
rng = np.random.default_rng(7)
def Cl(l): return l/2 + 1/(np.sqrt(2)*np.tan(l/np.sqrt(2)))
def c(th): return 2 - Cl(th)
print('c(1/2)=%.13f c(0.6)=%.7f c(0.75)=%.6f c(1)=%.10f' % (c(0.5), c(0.6), c(0.75), c(1.0)))
lo, hi = 0.5, 0.6
for _ in range(60):
    mid = (lo+hi)/2; lo, hi = (mid, hi) if c(mid) < 0 else (lo, mid)
Ct = Cl(0.25); q = int(np.floor(Ct)); print('theta0=%.10f; thm:distinct-all at 1/4: C=%.6f D=%.6f' % (lo, Ct, (2*q+1-Ct)/(q*(q+1))))
# --- C(f) by trapezoid on a fine grid (kernel |u-v| is Lipschitz: O(M^-2) error)
def Cfun(fv, u):
    du = u[1]-u[0]; wq = np.full(len(u), du); wq[0] = wq[-1] = du/2
    a = (wq*fv**2).sum()
    K = np.abs(u[:, None]-u[None, :])
    b = (wq*fv) @ K @ (wq*fv)
    return a + b, (wq*fv).sum()
M = 1601
for l in (0.3, 0.55, 0.8, 1.0):
    u = np.linspace(-l/2, l/2, M)
    fl = np.cos(np.sqrt(2)*u)/(np.sqrt(2)*np.sin(l/np.sqrt(2)))
    C0, m0 = Cfun(fl, u)
    worst = np.inf
    for k in range(300):  # random zero-mean perturbations h
        coef = rng.normal(size=8)*rng.uniform(0.01, 3)
        h = sum(cf*np.cos((j+1)*np.pi*(u/l+0.5)) for j, cf in enumerate(coef)) + rng.normal()*np.sin(np.pi*u/l)*0  # cos modes on J: zero mean
        h += rng.normal()*np.sin(2*np.pi*u/l)  # odd mode, zero mean
        C1, m1 = Cfun(fl+h, u)
        du = u[1]-u[0]; n2 = (h**2).sum()*du
        worst = min(worst, (C1-C0)/n2)
    print(f'lambda={l}: quad C(f_l)={C0:.10f} formula C_l={Cl(l):.10f} int f_l={m0:.10f}; min (C(f)-C_l)/||h||^2 over 300 random h = {worst:.4f} >= 1-2l^2/pi^2 = {1-2*l*l/np.pi**2:.4f}')
    # admissible eta_eps: chi_eps * sqrt(f_l), smooth cutoff; C(f_eps) -> C_l
    a = l/2
    for eps in (0.2*a, 0.05*a, 0.01*a):
        def chi(x):
            y = np.clip((a - np.abs(x))/eps, 0, 1)
            s = np.where((y > 0) & (y < 1), np.exp(-1/np.maximum(y,1e-300))/(np.exp(-1/np.maximum(y,1e-300)) + np.exp(-1/np.maximum(1-y,1e-300))), (y >= 1).astype(float))
            return s
        fe = chi(u)**2*fl; fe = fe/((fe*(u[1]-u[0])).sum())
        print(f'    eps={eps:.4f}: C(f_eps)={Cfun(fe,u)[0]:.8f}  (C_l={Cl(l):.8f})')
# --- lem:coeff sums
N = 2_000_000  # was 4e6: first run's max RSS 309 MB (footprint unmeasured) -> halved for the 300 MB cap
spf = np.zeros(N+1, dtype=np.int32)
for p in range(2, int(N**0.5)+1):
    if spf[p] == 0:
        blk = spf[p*p::p]; blk[blk == 0] = p
pr = np.nonzero(spf[2:] == 0)[0] + 2; spf[pr] = pr
Lam = np.zeros(N+1)
for p in pr:
    pk = p
    lp = np.log(p)
    while pk <= N:
        Lam[pk] = lp; pk *= p
n = np.arange(N+1, dtype=float)
print('Lambda check: sum_{n<=1e6} Lambda =', round(Lam[:1_000_001].sum()), '(psi(1e6)=999587)')
for x in (1.0, 2.0, 10.0, 100.0, 1e3, 1e4, 1e5):
    mn = np.minimum(n[2:]/x, x/n[2:])
    an2 = Lam[2:]**2/n[2:]*mn**2
    tail1 = x*x*(2*np.log(N)+1)/(4*N*N); tail2 = x*x*(np.log(N)**2+2*np.log(N)+2)/N
    s1 = an2.sum() + tail1; s2 = (n[2:]*an2).sum() + tail2
    print(f'x={x:g}: sum a_n^2 - log x = {s1-np.log(x):.4f};  sum n a_n^2/(x log^2 2x) = {s2/(x*np.log(2*x)**2):.4f}')
del spf, Lam, n
# --- Lamzouri Prop 2.1 (eq:finite) on random conjugation-invariant multisets, K = hat(f_l) at complex arguments
def Khat(z, l):
    a = l/2; b = np.sqrt(2); cc = 2*np.pi*z
    def part(wv):
        small = np.abs(wv) < 1e-9
        return np.where(small, a, np.sin(wv*a)/np.where(small, 1, wv))
    return (part(b-cc) + part(b+cc))/(np.sqrt(2)*np.sin(l/np.sqrt(2)))
worst = np.inf; worst_case = None
for trial in range(20000):
    l = rng.uniform(0.2, 1.0)
    reals = rng.uniform(0, rng.uniform(0.2, 6), size=rng.integers(0, 7))
    mult = rng.integers(1, 4, size=len(reals))
    Z = list(np.repeat(reals, mult).astype(complex))
    for k in range(rng.integers(0, 4)):
        zr = rng.uniform(0, 6); zi = rng.uniform(0.001, 1.5)
        Z += [complex(zr, zi), complex(zr, -zi)]
    if not Z: continue
    Z = np.array(Z)
    S = (Khat(Z[:, None]-Z[None, :], l)**2).sum()
    uniq, cnt = np.unique(reals, return_counts=True)
    simple_real = int(sum(1 for r_, m_ in zip(reals, mult) if m_ == 1 and np.sum(reals == r_) == 1))
    gap = simple_real - (2*len(Z) - S.real)
    if gap < worst: worst, worst_case = gap, (len(Z), simple_real, round(S.real,4), round(abs(S.imag),12))
print(f'Lamzouri (eq:finite): min over 20000 random sets of LHS-RHS = {worst:.6f} (case |Z|, N_s, S_re, |S_im| = {worst_case})')
