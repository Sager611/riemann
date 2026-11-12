# N2: distribution and L^2-tail of |R_T|^2 over [T,2T] (DPC coefficients), T = 1e4, 1e6; theta = 0.45, 0.55, 0.65.
import numpy as np, sys
rng = np.random.default_rng(1)
def bump01(x):
    x = np.clip(x, 0, 1); g = lambda y: np.where(y > 0, np.exp(-1/np.maximum(y, 1e-300)), 0.0)
    return g(x)/(g(x)+g(1-x))
def coeffs(T, th, eta=0.1, b=None):
    L = np.log(T); b = b or th+0.07; N = int(T**th)+1
    rho = lambda u: (1-bump01((u-th)/(b-th)))*bump01(u/eta-1)
    f = lambda x: np.exp(1-1/(1-(x/th)**2)) if x < th else 0.0
    spf = np.zeros(N+1, dtype=np.int64)            # smallest prime factor sieve
    for p in range(2, N+1):
        if spf[p] == 0: spf[p::p] = np.where(spf[p::p] == 0, p, spf[p::p])
    ns, cs = [], []
    for n in range(1, N):
        m, w, mu, ok = n, 1.0, 1, True
        while m > 1:
            p = spf[m]; m //= p
            if m % p == 0: ok = False; break
            w *= rho(np.log(p)/L); mu = -mu
            if w == 0: ok = False; break
        if not ok: continue
        c = mu*w*f(np.log(n)/L)/np.sqrt(n)
        if c != 0: ns.append(n); cs.append(c)
    return np.array(ns, float), np.array(cs)
def run(T, th, M=200000, chunk=1000):
    ns, cs = coeffs(T, th); S = np.sum(cs**2); lg = np.log(ns)
    t = rng.uniform(T, 2*T, M); R2 = np.empty(M)
    for i in range(0, M, chunk):
        tt = t[i:i+chunk]; R = np.exp(-1j*np.outer(tt, lg)) @ cs; R2[i:i+chunk] = np.abs(R)**2/S
    m2 = R2.mean(); m4 = (R2**2).mean()
    tails = [ (R2*(R2 > W)).mean()/m2 for W in (2,4,8,16,32,64) ]
    q = np.quantile(R2, [0.5, 0.9, 0.99, 0.999])
    print(f"T={T:.0e} th={th} terms={len(ns)} S={S:.3f} E|R|^2/S={m2:.3f} E|R|^4/(E|R|^2)^2={m4/m2**2:.2f} max={R2.max():.1f} q50,90,99,99.9={np.round(q,2)} tailmass(W=2..64)={np.round(tails,4)}", flush=True)
for T in (1e4, 1e6):
    for th in (0.45, 0.55, 0.65): run(T, th)

