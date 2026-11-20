# hdx_ceiling.py (explore64): rigorous lower bounds (python-flint arb) for the Kronecker ceiling of local-to-global certificates.
#  Lam_univ(c, c') = sup_{u>=0} sum_{n<c} 2 w_n max(0, <u, T_{ln n} u>/|u|^2 - cos(pi/(ceil(ln c'/ln n)+1)))   (any links of diameter <= ln c')
#  Lam_long(c, c') = sup_u sum_{c'<=n<c} 2 w_n <u, T_{ln n} u>/|u|^2                                            (atoms no link can see)
#  T_s = truncated shift on L^2[-a,a]; u piecewise constant on N cells (exact band formula, two diagonals per atom).
#  Any u >= 0 gives a valid lower bound; u = Perron vector of the active-atom matrix (float), then re-evaluated in arb with rational u.
import sys, math, numpy as np
from flint import arb, fmpq, ctx
ctx.prec = 128

def ppow(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n; e = 0
            while m % p == 0: m //= p; e += 1
            return (p, e) if m == 1 else None
    return None

def atoms(cq):  # prime powers n < c (c rational): (n, p)
    return [(n, ppow(n)[0]) for n in range(2, int(math.ceil(float(cq))) + 1) if fmpq(n) < cq and ppow(n)]

def band_float(u, s, h):
    k0 = int(math.floor(s/h)); N = len(u); tot = 0.0
    for k in (k0, k0 + 1):
        wgt = h - abs(k*h - s)
        if wgt > 0 and k < N: tot += wgt*float(np.dot(u[:N-k], u[k:]))
    return tot

def ceil_terms(cq, clq):
    lc, lcl = math.log(float(cq)), math.log(float(clq)); out = []
    for n, p in atoms(cq):
        s = math.log(n)
        if fmpq(n) >= clq: m = 1
        else: m = math.ceil(lcl/s - 1e-12)
        mw = math.ceil(lc/s - 1e-12)          # block size on the window: <u,T_s u> <= cos(pi/(mw+1))|u|^2
        wl = 0.0 if m == 1 else math.cos(math.pi/(m + 1))
        if mw > m: out.append((n, p, m, wl))
    return out

def optimize(cq, clq, N, longonly=False):
    a = math.log(float(cq))/2; h = 2*a/N; terms = ceil_terms(cq, clq)
    if longonly: terms = [t for t in terms if t[2] == 1]
    if not terms: return None, 0.0, [], None
    W = {n: math.log(p)/math.sqrt(n) for n, p, _, _ in terms}
    def mat(S):
        from scipy.sparse import diags
        dg = {0: np.zeros(N)}
        for n, p, m, wl in terms:
            if n not in S: continue
            s = math.log(n); k0 = int(math.floor(s/h))
            for k in (k0, k0 + 1):
                wgt = h - abs(k*h - s)
                if wgt > 0 and k < N: dg[k] = dg.get(k, np.zeros(N - k)) + W[n]*wgt/h
            dg[0] = dg[0] - 2*W[n]*wl
        ks = sorted(dg); vals = [dg[k] for k in ks] + [dg[k] for k in ks if k > 0]
        return diags(vals, ks + [-k for k in ks if k > 0], shape=(N, N), format='csr')
    def top(A):
        from scipy.sparse.linalg import eigsh
        w_, V = eigsh(A, k=1, which='LA', tol=1e-10, maxiter=20000); return V[:, 0]
    def F(u):
        nu = h*float(np.dot(u, u)); tot = 0.0
        for n, p, m, wl in terms: tot += 2*W[n]*max(0.0, band_float(u, math.log(n), h)/nu - wl)
        return tot
    S = set(n for n, _, _, _ in terms); best = (None, -1.0, None)
    for it in range(8):
        u = np.abs(top(mat(S))); val = F(u)
        if val > best[1]: best = (u, val, sorted(S))
        nu = h*float(np.dot(u, u))
        S2 = set(n for n, p, m, wl in terms if band_float(u, math.log(n), h)/nu > wl)
        if S2 == S or not S2: break
        S = S2
    return best[0], best[1], terms, F

def certify(cq, clq, u, terms, N, den=10**6):
    # rigorous lower bound of sum_n 2 w_n max(0, <u,T u>/|u|^2 - wl) with u rationalized (u >= 0)
    a = arb(cq).log()/2; h = 2*a/N
    q = [fmpq(int(round(x*den)), den) for x in u/np.max(u)]
    nu = h*sum(arb(x)*arb(x) for x in q); tot = arb(0)
    for n, p, m, wl in terms:
        s = arb(n).log(); k0 = int(math.floor(float(s.mid())/float(h.mid()))); acc = arb(0)
        for k in (k0 - 1, k0, k0 + 1, k0 + 2):
            if k < 0 or k >= N: continue
            wgt = h - abs(arb(k)*h - s)
            lo = wgt.lower()
            if lo <= 0: continue   # dropping a nonnegative term keeps a lower bound
            acc += arb(lo)*sum(arb(q[i])*arb(q[i + k]) for i in range(N - k))
        wlA = arb(0) if m == 1 else (arb.pi()/(m + 1)).cos()
        r = acc/nu - wlA
        if r.lower() > 0: tot += 2*arb(p).log()/arb(n).sqrt()*arb(r.lower())
    return tot

if __name__ == '__main__':
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 240
    pairs = [(10, '10^0.5'), (10, '10^0.75'), (10, '10^0.9'), (20, '20^0.5'), (30, '30^0.5'), (47, '47^0.5'), (47, '47^0.75'), (47, '47^0.9'),
             (49, 47), (50, 47), (53, 47), (60, 47), (71, 47), (100, 47), (72, 71), (80, 71), (100, 71),
             (fmpq(83, 2), fmpq(227, 6)), (45, fmpq(227, 6)), (50, fmpq(227, 6))]
    print(f'# N={N} cells. c, c_link, Lam_long (float, arb lower bound), Lam_univ (float, arb lower bound), active atoms (n:m)')
    for cq, cl in pairs:
        if isinstance(cl, str):
            b, e = cl.split('^'); clf = float(b)**float(e); clq = fmpq(int(round(clf*10**6)), 10**6)
        else: clq = fmpq(cl); clf = float(clq)
        cq = fmpq(cq)
        u1, v1, t1, F1 = optimize(cq, clq, N, longonly=True)
        c1 = certify(cq, clq, u1, [t for t in t1 if t[2] == 1], N) if u1 is not None else arb(0)
        u2, v2, t2, F2 = optimize(cq, clq, N)
        if u1 is not None and u2 is not None and F2(u1) > v2: u2, v2 = u1, F2(u1)
        c2 = certify(cq, clq, u2, t2, N) if u2 is not None else arb(0)
        act = ' '.join(f'{n}:{m}' for n, p, m, wl in t2)
        print(f'{float(cq):7.3f} {clf:8.3f}  long {v1:.4f} >= {float(c1.lower()):.4f}   univ {v2:.4f} >= {float(c2.lower()):.4f}   [{act}]', flush=True)
