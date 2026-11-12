
# explore40: explore34's lattice dual LP for Lamzouri's Thm Main2 + realizability (hypermetric) rows.
# Dual: mu = delta_0 + sum_{k!=0} h(1+eta_|k|) delta_{kh}, eta_k=-1 (k<k0), P(al) >= |al| 1_theta + t on [0,1/(2h)].
# Added (valid for any simple point set with hard core lam): triangle  mu_p+mu_q-mu_{p+q}<=1, mu_p+mu_q-mu_{|q-p|}<=1,
# packing sum_{k=p}^{p+k0-1} mu_k <= 1.   t*>0: augmented method cannot prove D(lam)>0 (on this lattice).
import numpy as np, sys, time, threading, os, resource, json
def _wd(cap=285e6):
    while True:
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > cap:
            print("MEMCAP exceeded, abort", flush=True); os._exit(3)
        time.sleep(0.1)
threading.Thread(target=_wd, daemon=True).start()
from scipy.optimize import linprog
import scipy.sparse as sp
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6

def Pvals(al, eta, h, chunk=1500):
    K = len(eta); kh = np.arange(1, K+1)*h
    out = np.empty(len(al))
    for i in range(0, len(al), chunk):
        out[i:i+chunk] = (1-h) + 2*h*(np.cos(2*np.pi*np.outer(al[i:i+chunk], kh)) @ eta)
    return out

def real_rows(K, k0, h, tri=True, pack=True):
    # returns sparse A (rows x nf) and b for constraints on eta_free (k=k0..K) ; column index k-k0
    nf = K-k0+1; rows=[]; cols=[]; vals=[]; b=[]; r=0
    def col(k): return k-k0
    if tri:
        for p in range(k0, K+1):
            for q in range(p, K+1):
                # T1: eta_p+eta_q-eta_{p+q} <= 1/h-1  (eta_{p+q}=0 if p+q>K)
                cs={col(p):1.0}; cs[col(q)]=cs.get(col(q),0)+1.0
                if p+q<=K: cs[col(p+q)]=cs.get(col(p+q),0)-1.0
                for c,v in cs.items(): rows.append(r); cols.append(c); vals.append(v)
                b.append(1/h-1); r+=1
                if q>p:
                    d=q-p; cs={col(p):1.0, col(q):1.0}; rhs=1/h-1
                    if d<k0: rhs=1/h-2
                    else: cs[col(d)]=cs.get(col(d),0)-1.0
                    for c,v in cs.items(): rows.append(r); cols.append(c); vals.append(v)
                    b.append(rhs); r+=1
    if pack:
        for p in range(1, K+1):
            ks=[k for k in range(p,p+k0) if k0<=k<=K]; nout=sum(1 for k in range(p,p+k0) if k>K)
            if not ks: continue
            for k in ks: rows.append(r); cols.append(col(k)); vals.append(h)
            b.append(1-h*len(ks)-h*nout); r+=1
    A=sp.csr_matrix((vals,(rows,cols)),shape=(r,nf))
    return A, np.array(b)

def solve(al, theta, h, K, k0, B, RA, Rb):
    ks = np.arange(1, K+1); free = ks >= k0; nf = int(free.sum())
    tgt = np.where(al <= theta, al, 0.0)
    C = 2*h*np.cos(2*np.pi*np.outer(al, ks*h))
    base = (1-h) - C[:, ~free].sum(axis=1)
    A = np.empty((len(al), nf+1)); A[:, :nf] = -C[:, free]; A[:, -1] = 1.0
    del C
    c = np.zeros(nf+1); c[-1] = -1.0
    if RA is not None:
        Afull = sp.vstack([sp.csr_matrix(A), sp.hstack([RA, sp.csr_matrix((RA.shape[0],1))])]).tocsr()
        bfull = np.concatenate([base-tgt, Rb])
    else:
        Afull = A; bfull = base-tgt
    ub = min(B, 1/h-1)
    for meth in ('highs-ds', 'highs-ipm'):
        r = linprog(c, A_ub=Afull, b_ub=bfull, bounds=[(-1.0, ub)]*nf+[(None, 1.0)], method=meth)
        if r.status == 0: break
    if r.status != 0: return None, None
    eta = np.full(K, -1.0); eta[free] = r.x[:nf]
    return r.x[-1], eta

def dual(lam, theta, h, X, B=20.0, fine=0.0002, iters=40, tol=1e-7, band=0.003, newmax=250, rowmax=1100, real=True, verbose=False):
    k0 = int(round(lam/h)); assert abs(k0*h-lam) < 1e-9
    K = int(round(X/h)); amax = 0.5/h
    RA, Rb = real_rows(K, k0, h) if real else (None, None)
    alf = np.arange(0, amax+1e-12, fine); bf = np.where(alf <= theta, alf, 0.0)
    keep = np.array([0.0, theta, min(theta+fine, amax)])
    al = np.unique(np.concatenate([np.linspace(0, amax, 600), keep, np.linspace(0, 0.02, 11), np.linspace(theta-0.02, theta, 11)]))
    t = None
    for it in range(iters):
        t, eta = solve(al, theta, h, K, k0, B, RA, Rb)
        if t is None: return None, None, None, it
        v = Pvals(alf, eta, h) - bf; mn = v.min()
        if verbose: print(f"   it{it} rows={len(al)} t={t:.5e} finemin={mn:.5e} rss={rss():.0f}", flush=True)
        if mn >= t - tol: return t, mn, eta, it
        lmask = (v <= np.r_[np.inf, v[:-1]]) & (v <= np.r_[v[1:], np.inf])
        lm = np.where(lmask & (v < t + band))[0]
        if len(lm) > newmax: lm = lm[np.argsort(v[lm])[:newmax]]
        vr = Pvals(al, eta, h) - np.where(al <= theta, al, 0.0)
        sl = vr - t
        old = al[sl < band]
        if len(old) > rowmax - len(lm): old = al[np.argsort(sl)[:max(0, rowmax - len(lm))]]
        al = np.unique(np.concatenate([old, keep, alf[lm]]))
    return t, mn, eta, iters

if __name__ == '__main__':
    theta = eval(sys.argv[1]); k0 = int(sys.argv[2])
    Xs = [float(x) for x in sys.argv[3].split(',')]
    lams = [float(x) for x in sys.argv[4].split(',')]
    tag = sys.argv[5]; real = '-noreal' not in sys.argv; verbose = '-v' in sys.argv
    for X in Xs:
        for lam in lams:
            t0 = time.time(); h = lam/k0
            t, mn, eta, it = dual(lam, theta, h, X, real=real, verbose=verbose)
            print(f"real={real} theta={theta:.7f} k0={k0} h={h:.5f} X={X} K={int(round(X/h))} lam={lam:.5f} t*={t} finemin={mn} iters={it} ({time.time()-t0:.1f}s) maxrss={rss():.0f}MB", flush=True)
            if t is not None:
                fn = f"logs/real_{tag}_X{int(X)}_lam{lam}.json"
                json.dump({'theta': theta, 'lam': lam, 'h': h, 'X': X, 'k0': k0, 't': t, 'finemin': mn, 'eta': eta.tolist(), 'real': real}, open(fn, 'w'))
