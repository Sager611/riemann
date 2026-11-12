
import numpy as np, sys, time, threading, os, resource, json
def _wd(cap=280e6):
    while True:
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > cap:
            print("MEMCAP exceeded, abort", flush=True); os._exit(3)
        time.sleep(0.2)
threading.Thread(target=_wd, daemon=True).start()
from scipy.optimize import linprog

def Pvals(al, eta, h):
    K = len(eta); kh = np.arange(1, K+1)*h
    out = np.empty(len(al))
    for i in range(0, len(al), 2000):
        out[i:i+2000] = (1-h) + 2*h*(np.cos(2*np.pi*np.outer(al[i:i+2000], kh)) @ eta)
    return out

def solve(al, lam, theta, h, K, k0, B):
    ks = np.arange(1, K+1); free = ks >= k0; nf = int(free.sum())
    b = np.where(al <= theta, al, 0.0)
    C = 2*h*np.cos(2*np.pi*np.outer(al, ks*h))
    base = (1-h) - C[:, ~free].sum(axis=1)
    A = np.hstack([-C[:, free], np.ones((len(al),1))])
    c = np.zeros(nf+1); c[-1] = -1.0
    for meth in ('highs-ds', 'highs-ipm', 'highs'):
        r = linprog(c, A_ub=A, b_ub=base-b, bounds=[(-1.0, B)]*nf+[(None, 1.0)], method=meth)
        if r.status == 0: break
    if r.status != 0: return None, None
    eta = np.full(K, -1.0); eta[free] = r.x[:nf]
    return r.x[-1], eta

def dual(lam, theta, h, X, B=20.0, fine=0.0002, iters=40, tol=1e-7, band=0.015, maxrows=400, verbose=False):
    k0 = int(round(lam/h)); assert abs(k0*h-lam) < 1e-9, (k0, h, lam)
    K = int(round(X/h)); amax = 0.5/h
    al = np.unique(np.concatenate([np.linspace(0, amax, 500), [theta], np.linspace(0, 0.02, 11), np.linspace(theta-0.02, theta, 11)]))
    alf = np.arange(0, amax+1e-12, fine); bf = np.where(alf <= theta, alf, 0.0)
    keep_always = np.array([0.0, theta])
    for it in range(iters):
        t, eta = solve(al, lam, theta, h, K, k0, B)
        if t is None: return None, None, None
        v = Pvals(alf, eta, h) - bf; mn = v.min()
        if verbose: print(f"   it{it} rows={len(al)} t={t:.4e} finemin={mn:.4e}", flush=True)
        if mn >= t - tol: break
        # active set: local minima of v within band of t, plus current rows with small slack
        lmask = (v <= np.r_[np.inf, v[:-1]]) & (v <= np.r_[v[1:], np.inf])
        lm = np.where(lmask & (v < t + band))[0]
        vr = Pvals(al, eta, h) - np.where(al <= theta, al, 0.0)
        old = al[vr < t + band]
        new = alf[lm]
        if len(new) > maxrows: new = alf[lm[np.argsort(v[lm])[:maxrows]]]
        al = np.unique(np.concatenate([old, new, keep_always, alf[lm] - fine, alf[lm] + fine]))
        al = al[(al >= 0) & (al <= amax)]
    return t, mn, eta

if __name__ == '__main__':
    theta = eval(sys.argv[1]); k0 = int(sys.argv[2]); X = float(sys.argv[3])
    lams = [float(x) for x in sys.argv[4].split(',')]
    save = sys.argv[5] if len(sys.argv) > 5 else None
    best = None
    for lam in lams:
        t0 = time.time(); h = lam/k0
        t, mn, eta = dual(lam, theta, h, X)
        print(f"theta={theta:.7f} k0={k0} h={h:.5f} X={X} lam={lam:.5f} t*={t} finemin={mn} ({time.time()-t0:.1f}s) maxrss={resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f}MB", flush=True)
        if t is not None and mn is not None and mn > 0 and (best is None or lam > best[0]): best = (lam, h, t, mn, eta.tolist())
    if save and best:
        json.dump({'theta': theta, 'lam': best[0], 'h': best[1], 'X': X, 'k0': k0, 't': best[2], 'finemin': best[3], 'eta': best[4]}, open(save, 'w'))
        print("saved", save, best[0])
