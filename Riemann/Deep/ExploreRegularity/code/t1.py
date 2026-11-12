# T1: DPC's rough Mobius polynomial R_T and its prime sum S_1 (2610.28527 Lemma 5), tested on actual zeros.
# (a) explicit formula  S_L(t) := (1/L) sum_n Lambda(n) n^{-1/2-it} Fhat(log n/L) = -sum_gamma F(L(gamma-t)/2pi) + (negligible)
#     with Fhat(u) = rho(u)/u, F(x) = int Fhat(u) e^{2 pi i x u} du  (unconditional; off-line zeros would enter at complex gamma)
# (b) moments of R_T at theta = 0.45, 0.55, 0.60 vs T (sampling t), normalised by (sum |c_n|^2)^{q/2}
import sys, json, numpy as np
rng = np.random.default_rng(7)
def S(x):
    x = np.clip(x, 1e-12, 1-1e-12); a = np.exp(-1/x); b = np.exp(-1/(1-x)); return a/(a+b)
def smooth_step(x):
    out = np.where(x <= 0, 0.0, np.where(x >= 1, 1.0, 0.0))
    m = (x > 0) & (x < 1); out = out.astype(float); out[m] = S(x[m]); return out
def rho(u, eta, theta, b):
    u = np.asarray(u, float)
    return (1 - smooth_step((u - theta)/(b - theta))) * smooth_step(u/eta - 1)
def f(x, theta):
    x = np.asarray(x, float); return np.where(x < theta, (1 - x/theta)**2, 0.0)
def primes_upto(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
def mobius_rough(N, pmin):
    # squarefree n <= N all of whose prime factors are > pmin; returns n, mu(n), list of prime factors via smallest-factor table
    spf = np.zeros(N+1, np.int64)
    for p in primes_upto(N):
        m = spf[p::p] == 0; idx = np.arange(p, N+1, p)[m]; spf[idx] = p
    ns, mus, facs = [1], [1], [[]]
    for n in range(2, N+1):
        m, fl, ok = n, [], True
        while m > 1:
            p = spf[m]; m //= p
            if p <= pmin or (m % p == 0): ok = False; break
            fl.append(p)
        if ok: ns.append(n); mus.append((-1)**len(fl)); facs.append(fl)
    return np.array(ns), np.array(mus), facs
def explicit_test(T, zeros, eta=0.1, theta=0.5, b=0.6, npts=60, half=None, centers=None):
    L = np.log(T)
    uu = np.linspace(eta/2, b, 4001); w = np.gradient(uu)
    Fh = rho(uu, eta, theta, b)/uu
    # prime-power side
    nmax = int(np.exp(b*L)) + 1
    ps = primes_upto(nmax)
    nn, lam = [], []
    for p in ps:
        q = p
        while q <= nmax:
            nn.append(q); lam.append(np.log(p)); q *= p
    nn = np.array(nn, float); lam = np.array(lam)
    coef = lam*rho(np.log(nn)/L, eta, theta, b)/(np.log(nn)/L)/L/np.sqrt(nn)
    keep = coef != 0; nn, coef = nn[keep], coef[keep]
    ts = centers if centers is not None else np.sort(rng.uniform(zeros.min()+200, zeros.max()-200, npts))
    res = []
    for t in ts:
        SL = np.sum(coef*np.exp(-1j*t*np.log(nn)))
        g = zeros[np.abs(zeros - t) < (half if half else 300/L*2*np.pi)]
        x = L*(g - t)/(2*np.pi)
        Fx = (np.exp(2j*np.pi*np.outer(x, uu)) * (Fh*w)).sum(axis=1)
        Z = -Fx.sum()
        res.append((t, SL, Z))
    return res, len(nn)
def rt_moments(T, theta, eta=0.05, nsamp=20000, qs=(2, 3, 4, 6), diag=False):
    L = np.log(T); N = int(np.exp(theta*L)); b = theta + 0.1
    n, mu, facs = mobius_rough(N, np.exp(eta*L/2))
    w = np.array([np.prod(rho(np.log(np.array(fl, float))/L, eta, theta, b)) if fl else 1.0 for fl in facs])
    c = mu/np.sqrt(n)*f(np.log(n)/L, theta)*w
    keep = c != 0; n, c = n[keep].astype(float), c[keep]
    G = np.sum(c**2)
    D4 = np.nan
    if diag and len(n) <= 1500:
        prod = np.outer(n.astype(np.int64), n.astype(np.int64)).ravel(); wv = np.outer(c, c).ravel()
        u, inv = np.unique(prod, return_inverse=True); bsum = np.bincount(inv, weights=wv)
        D4 = float(np.sum(bsum**2)/G**2); del prod, wv, u, inv, bsum
    ts = rng.uniform(T, 2*T, nsamp)
    vals = np.empty(nsamp, complex); ln = np.log(n)
    ch = max(1, int(2e6/len(n)))
    for i in range(0, nsamp, ch):
        vals[i:i+ch] = (np.exp(-1j*np.outer(ts[i:i+ch], ln))*c).sum(axis=1)
    a = np.abs(vals)
    return dict(T=T, theta=theta, N=N, terms=len(n), G=G, kmax=int(theta/(eta/2)),
                M={q: float(np.mean(a**q)/G**(q/2)) for q in qs}, mv4=1 + T**(2*theta - 1), maxR=float(a.max()/np.sqrt(G)), D4=D4)
if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'a':
        import os
        z = np.load(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'ExploreHelio', 'results', 'zeta_zeros.npy'))
        for th, b in ((0.5, 0.6), (0.8, 1.0)):
            res, nterms = explicit_test(1.0e4, z, eta=0.1, theta=th, b=b, npts=40)
            d = np.array([abs(r[1] - r[2]) for r in res]); s = np.array([abs(r[1]) for r in res])
            print('T=1e4 supp rho=(%.2f,%.2f)xL: %d prime powers; |S_L| mean %.3f max %.3f; |S_L - zero side| max %.2e mean %.2e'
                  % (0.05, b, nterms, s.mean(), s.max(), d.max(), d.mean()))
    if mode == 'b2':
        for T, th in ((1e4,0.45),(1e4,0.6),(1e6,0.45),(1e6,0.55),(1e6,0.6),(1e8,0.45),(1e8,0.55),(1e8,0.6)):
            r = rt_moments(T, th, diag=True)
            print('T=%.0e th=%.2f N=%d terms=%d G=%.3f M2=%.3f M3=%.3f M4=%.3f (diag %.3f) M6=%.3f MVfactor=%.2f max|R|/sqrtG=%.2f' % (T, th, r['N'], r['terms'], r['G'], r['M'][2], r['M'][3], r['M'][4], r['D4'], r['M'][6], r['mv4'], r['maxR']), flush=True)
    if mode == 'b':
        for T in (1e4, 1e6, 1e8, 1e10):
            for th in (0.45, 0.55, 0.6):
                if T == 1e10 and th > 0.55: continue
                r = rt_moments(T, th)
                print('T=%.0e th=%.2f N=%d terms=%d G=%.3f kmax=%d M2=%.3f M3=%.3f M4=%.3f M6=%.3f MVfactor4=%.2f max|R|/sqrtG=%.2f'
                      % (T, th, r['N'], r['terms'], r['G'], r['kmax'], r['M'][2], r['M'][3], r['M'][4], r['M'][6], r['mv4'], r['maxR']), flush=True)
