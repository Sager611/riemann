
# explore40: small-gap dual LP (explore34 / Lamzouri Thm Main2) + EXACT local realizability:
# the Palm weights mu_k (k < R) must be the pair statistics of a stationary hard-core 0/1 lattice process,
# encoded by a shift-consistent probability vector p on hard-core words of length R (1D: consistent R-marginals
# extend to a stationary (R-1)-Markov process, so this is necessary and sufficient for the range-(R-1) pair data).
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

def words(R, k0):
    out=[]
    def rec(pos, last, w):
        if pos==R: out.append(w); return
        rec(pos+1, last, w)                      # 0 at pos
        if last is None or pos-last>=k0: rec(pos+1, pos, w|(1<<pos))
    rec(0,None,0); return out

def pattern_block(R, k0, h, K):
    W = words(R, k0); nw=len(W)
    m = (1<<(R-1))-1
    pre = {}; idx=0
    W1 = words(R-1, k0); id1 = {u:i for i,u in enumerate(W1)}
    rows=[]; cols=[]; vals=[]; r=0
    # stationarity
    for j,w in enumerate(W):
        rows.append(id1[w & m]); cols.append(j); vals.append(1.0)
        rows.append(id1[w >> 1]); cols.append(j); vals.append(-1.0)
    r = len(W1); beq=[0.0]*r
    # linking k0..R-1 : sum_w p_w [b0 & bk] - h^2 eta_k = h^2 ; eta_k column = nw + (k-k0)
    for k in range(k0, R):
        for j,w in enumerate(W):
            if (w & 1) and ((w>>k)&1): rows.append(r); cols.append(j); vals.append(1.0)
        rows.append(r); cols.append(nw+(k-k0)); vals.append(-h*h); beq.append(h*h); r+=1
    # density, normalisation
    for j,w in enumerate(W):
        if w & 1: rows.append(r); cols.append(j); vals.append(1.0)
    beq.append(h); r+=1
    for j,w in enumerate(W): rows.append(r); cols.append(j); vals.append(1.0)
    beq.append(1.0); r+=1
    nf = K-k0+1
    Aeq = sp.csr_matrix((vals,(rows,cols)), shape=(r, nw+nf+1))
    return nw, Aeq, np.array(beq)

def Pvals(al, eta, h, chunk=1500):
    K = len(eta); kh = np.arange(1, K+1)*h
    out = np.empty(len(al))
    for i in range(0, len(al), chunk):
        out[i:i+chunk] = (1-h) + 2*h*(np.cos(2*np.pi*np.outer(al[i:i+chunk], kh)) @ eta)
    return out

def solve(al, theta, h, K, k0, B, PB):
    ks = np.arange(1, K+1); free = ks >= k0; nf = int(free.sum())
    nw, Aeq, beq = PB if PB is not None else (0, None, None)
    tgt = np.where(al <= theta, al, 0.0)
    C = 2*h*np.cos(2*np.pi*np.outer(al, ks*h))
    base = (1-h) - C[:, ~free].sum(axis=1)
    A = np.zeros((len(al), nw+nf+1)); A[:, nw:nw+nf] = -C[:, free]; A[:, -1] = 1.0
    del C
    c = np.zeros(nw+nf+1); c[-1] = -1.0
    bounds = [(0.0,1.0)]*nw + [(-1.0, min(B, 1/h-1))]*nf + [(None, 1.0)]
    Aub = sp.csr_matrix(A)
    for meth in ('highs-ds', 'highs-ipm'):
        r = linprog(c, A_ub=Aub, b_ub=base-tgt, A_eq=Aeq, b_eq=beq, bounds=bounds, method=meth)
        if r.status == 0: break
    if r.status != 0: print('   LP status', r.status, r.message, flush=True); return None, None, None
    eta = np.full(K, -1.0); eta[free] = r.x[nw:nw+nf]
    return r.x[-1], eta, r.x[:nw]

def dual(lam, theta, h, X, R, B=20.0, fine=0.0002, iters=40, tol=1e-7, band=0.003, newmax=250, rowmax=1100, verbose=False):
    k0 = int(round(lam/h)); assert abs(k0*h-lam) < 1e-9
    K = int(round(X/h)); amax = 0.5/h
    PB = pattern_block(R, k0, h, K) if R>0 else None
    if PB is not None and verbose: print(f"   words={PB[0]} eqrows={PB[1].shape[0]} rss={rss():.0f}", flush=True)
    alf = np.arange(0, amax+1e-12, fine); bf = np.where(alf <= theta, alf, 0.0)
    keep = np.array([0.0, theta, min(theta+fine, amax)])
    al = np.unique(np.concatenate([np.linspace(0, amax, 600), keep, np.linspace(0, 0.02, 11), np.linspace(theta-0.02, theta, 11)]))
    t = None; p=None
    for it in range(iters):
        t, eta, p = solve(al, theta, h, K, k0, B, PB)
        if t is None: return None, None, None, it, None
        v = Pvals(alf, eta, h) - bf; mn = v.min()
        if verbose: print(f"   it{it} rows={len(al)} t={t:.5e} finemin={mn:.5e} rss={rss():.0f}", flush=True)
        if mn >= t - tol: return t, mn, eta, it, p
        lmask = (v <= np.r_[np.inf, v[:-1]]) & (v <= np.r_[v[1:], np.inf])
        lm = np.where(lmask & (v < t + band))[0]
        if len(lm) > newmax: lm = lm[np.argsort(v[lm])[:newmax]]
        vr = Pvals(al, eta, h) - np.where(al <= theta, al, 0.0)
        sl = vr - t
        old = al[sl < band]
        if len(old) > rowmax - len(lm): old = al[np.argsort(sl)[:max(0, rowmax - len(lm))]]
        al = np.unique(np.concatenate([old, keep, alf[lm]]))
    return t, mn, eta, iters, p

if __name__ == '__main__':
    theta = eval(sys.argv[1]); k0 = int(sys.argv[2]); X=float(sys.argv[3]); R=int(sys.argv[4])
    lams = [float(x) for x in sys.argv[5].split(',')]; tag=sys.argv[6]; verbose='-v' in sys.argv
    for lam in lams:
        t0=time.time(); h=lam/k0
        t, mn, eta, it, p = dual(lam, theta, h, X, R, verbose=verbose)
        print(f"R={R} theta={theta:.7f} k0={k0} h={h:.5f} X={X} K={int(round(X/h))} lam={lam:.5f} t*={t} finemin={mn} iters={it} ({time.time()-t0:.1f}s) maxrss={rss():.0f}MB", flush=True)
        if t is not None:
            json.dump({'theta':theta,'lam':lam,'h':h,'X':X,'k0':k0,'R':R,'t':t,'finemin':mn,'eta':eta.tolist(),'p':(p.tolist() if p is not None else None)}, open(f"logs/pat_{tag}_R{R}_X{int(X)}_lam{lam}.json",'w'))
