
import numpy as np, sys, time, threading, os, resource, json
def _wd(cap=285e6):
    while True:
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > cap:
            print("MEMCAP exceeded, abort", flush=True); os._exit(3)
        time.sleep(0.1)
threading.Thread(target=_wd, daemon=True).start()
from scipy.optimize import linprog
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6

TAIL = '-tail' in sys.argv
def Th(al, h):
    # sum_m Delta(al - m/h) for al in [0, 1/(2h)] (h < 1): m = 0, 1 only
    p = 1.0/h
    return np.clip(1-np.abs(al), 0, None) + np.clip(1-np.abs(al-p), 0, None)
def tailcol(al, h, K):
    kh = np.arange(1, K+1)*h; s2 = np.sinc(kh)**2
    out = np.empty(len(al))
    for i in range(0, len(al), 1500):
        a = al[i:i+1500]
        out[i:i+1500] = -Th(a, h) + h + 2*h*(np.cos(2*np.pi*np.outer(a, kh)) @ s2)
    return out
def Pvals(al, eta, h, chunk=1500, ctail=0.0):
    K = len(eta); kh = np.arange(1, K+1)*h
    out = np.empty(len(al))
    for i in range(0, len(al), chunk):
        out[i:i+chunk] = (1-h) + 2*h*(np.cos(2*np.pi*np.outer(al[i:i+chunk], kh)) @ eta)
    if ctail != 0.0: out += ctail*tailcol(al, h, K)
    return out

def solve(al, theta, h, K, k0, B):
    ks = np.arange(1, K+1); free = ks >= k0; nf = int(free.sum())
    b = np.where(al <= theta, al, 0.0)
    C = 2*h*np.cos(2*np.pi*np.outer(al, ks*h))
    base = (1-h) - C[:, ~free].sum(axis=1)
    nt = 1 if TAIL else 0
    A = np.empty((len(al), nf+1+nt)); A[:, :nf] = -C[:, free]; A[:, -1] = 1.0
    if TAIL: A[:, nf] = -tailcol(al, h, K)
    del C
    c = np.zeros(nf+1+nt); c[-1] = -1.0
    cmax = 0.99*(np.pi*K*h)**2
    for meth in ('highs-ds', 'highs-ipm'):
        r = linprog(c, A_ub=A, b_ub=base-b, bounds=[(-1.0, B)]*nf+[(0.0, cmax)]*nt+[(None, 1.0)], method=meth)
        if r.status == 0: break
    if r.status != 0: return None, None
    eta = np.full(K, -1.0); eta[free] = r.x[:nf]
    global CT
    CT = r.x[nf] if TAIL else 0.0
    return r.x[-1], eta

KEEPALL = '-keep' in sys.argv
CT = 0.0
BDEF = float(os.environ.get('DUAL_B', '20'))
FINEDEF = float(os.environ.get('DUAL_FINE', '0.0002'))
def dual(lam, theta, h, X, B=None, fine=FINEDEF, iters=40, tol=1e-7, band=0.003, newmax=250, rowmax=1100, verbose=False):
    if B is None: B = BDEF
    k0 = int(round(lam/h)); assert abs(k0*h-lam) < 1e-9
    K = int(round(X/h)); amax = 0.5/h
    alf = np.arange(0, amax+1e-12, fine); bf = np.where(alf <= theta, alf, 0.0)
    keep = np.array([0.0, theta, min(theta+fine, amax)])
    al = np.unique(np.concatenate([np.linspace(0, amax, 600), keep, np.linspace(0, 0.02, 11), np.linspace(theta-0.02, theta, 11)]))
    t = None
    for it in range(iters):
        t, eta = solve(al, theta, h, K, k0, B)
        if t is None: return None, None, None, it
        v = Pvals(alf, eta, h, ctail=CT) - bf; mn = v.min()
        if verbose: print(f"   it{it} rows={len(al)} t={t:.5e} finemin={mn:.5e} rss={rss():.0f}", flush=True)
        if mn >= t - tol: return t, mn, eta, it
        lmask = (v <= np.r_[np.inf, v[:-1]]) & (v <= np.r_[v[1:], np.inf])
        lm = np.where(lmask & (v < t + band))[0]
        if len(lm) > newmax: lm = lm[np.argsort(v[lm])[:newmax]]
        vr = Pvals(al, eta, h, ctail=CT) - np.where(al <= theta, al, 0.0)
        sl = vr - t
        old = al if KEEPALL else al[sl < band]
        if len(old) > rowmax - len(lm): old = al[np.argsort(sl)[:max(0, rowmax - len(lm))]]
        al = np.unique(np.concatenate([old, keep, alf[lm]]))
    return t, mn, eta, iters

if __name__ == '__main__':
    theta = eval(sys.argv[1]); k0 = int(sys.argv[2])
    Xs = [float(x) for x in sys.argv[3].split(',')]
    lams = [float(x) for x in sys.argv[4].split(',')]
    tag = sys.argv[5]
    verbose = '-v' in sys.argv
    for X in Xs:
        for lam in lams:
            t0 = time.time(); h = lam/k0
            t, mn, eta, it = dual(lam, theta, h, X, verbose=verbose)
            print(f"theta={theta:.7f} k0={k0} h={h:.5f} X={X} K={int(round(X/h))} lam={lam:.5f} ctail={CT:.5f} t*={t} finemin={mn} iters={it} ({time.time()-t0:.1f}s) maxrss={rss():.0f}MB", flush=True)
            if t is not None and mn is not None and mn > 0:
                fn = f"logs/cert_{tag}_X{int(X)}_lam{lam}.json"
                json.dump({'theta': theta, 'lam': lam, 'h': h, 'X': X, 'k0': k0, 't': t, 'finemin': mn, 'eta': eta.tolist(), 'ctail': CT}, open(fn, 'w'))
                print("saved", fn, flush=True)
