
# Kink-exact dual: sigma = sum_j c_j psi_{a_j} + sum_{k>=k0} h eta_k delta_{kh}, psi_a = pair measure of the a-discrete sine (a = m h >= lam),
# sum c_j = 1, c_j >= 0. P(al) = 1 - sum_j c_j T_j(al) + 2h sum_k eta_k cos(2 pi al k h), T_a(al) = sum_m Delta(al - m/a).
import numpy as np, sys, time, threading, os, resource, json
def _wd(cap=285e6):
    while True:
        if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > cap:
            print("MEMCAP exceeded, abort", flush=True); os._exit(3)
        time.sleep(0.1)
threading.Thread(target=_wd, daemon=True).start()
from scipy.optimize import linprog
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6
def tri_comb(al, a):
    # sum_m Delta(al - m/a), al >= 0 array
    out = np.zeros_like(al); p = 1.0/a
    mmax = int(np.ceil((al.max()+1)/p)) + 1
    for m in range(-1, mmax+1):
        out += np.clip(1 - np.abs(al - m*p), 0, None)
    return out
def cols(al, h, ms, k0, K):
    Tc = np.stack([-tri_comb(al, m*h) for m in ms], axis=1)
    ks = np.arange(k0, K+1)
    Ec = 2*h*np.cos(2*np.pi*np.outer(al, ks*h))
    return Tc, Ec
def Pvals(al, c, eta, h, ms, k0, K, chunk=2000):
    out = np.empty(len(al))
    for i in range(0, len(al), chunk):
        Tc, Ec = cols(al[i:i+chunk], h, ms, k0, K)
        out[i:i+chunk] = 1 + Tc @ c + Ec @ eta
    return out
def solve(al, theta, h, ms, k0, K, B):
    J = len(ms); ks = np.arange(k0, K+1); nk = len(ks)
    Tc, Ec = cols(al, h, ms, k0, K)
    b = np.where(al <= theta, al, 0.0)
    # rows: -(Tc c + Ec eta) + t <= 1 - b
    A1 = np.hstack([-Tc, -Ec, np.ones((len(al),1))]); b1 = 1 - b
    # weights W_k >= 0: -(sum_j c_j [m_j|k] a_j (1-sinc^2(kh)) + h eta_k) <= 0
    Wc = np.zeros((nk, J))
    for j, m in enumerate(ms):
        sel = (ks % m == 0); Wc[sel, j] = m*h*(1 - np.sinc(ks[sel]*h)**2)
    A2 = np.hstack([-Wc, -h*np.eye(nk), np.zeros((nk,1))]); b2 = np.zeros(nk)
    Aeq = np.hstack([np.ones((1,J)), np.zeros((1,nk)), np.zeros((1,1))]); beq = [1.0]
    c = np.zeros(J+nk+1); c[-1] = -1
    bounds = [(0, None)]*J + [(-B, B)]*nk + [(None, 1.0)]
    for meth in ('highs-ds', 'highs-ipm'):
        r = linprog(c, A_ub=np.vstack([A1, A2]), b_ub=np.r_[b1, b2], A_eq=Aeq, b_eq=beq, bounds=bounds, method=meth)
        if r.status == 0: break
    if r.status != 0: return None, None, None
    return r.x[-1], r.x[:J], r.x[J:J+nk]
def dual(lam, theta, k0, X, mmax, B=20.0, fine=0.0002, iters=60, tol=1e-8, verbose=False):
    h = lam/k0; K = int(round(X/h)); amax = 0.5/h; ms = list(range(k0, mmax+1))
    alf = np.arange(0, amax+1e-12, fine); bf = np.where(alf <= theta, alf, 0.0)
    al = np.unique(np.concatenate([np.linspace(0, amax, 800), [theta], np.linspace(0, 0.02, 11)]))
    for it in range(iters):
        t, c, eta = solve(al, theta, h, ms, k0, K, B)
        if t is None: return None
        v = Pvals(alf, c, eta, h, ms, k0, K) - bf; mn = v.min()
        if verbose: print(f"   it{it} rows={len(al)} t={t:.5e} finemin={mn:.5e} rss={rss():.0f}", flush=True)
        if mn >= t - tol: break
        lmask = (v <= np.r_[np.inf, v[:-1]]) & (v <= np.r_[v[1:], np.inf])
        lm = np.where(lmask & (v < t + 0.003))[0]
        al = np.unique(np.concatenate([al, alf[lm], alf[np.argsort(v)[:50]]]))
    return dict(t=t, finemin=mn, c=c.tolist(), eta=eta.tolist(), h=h, k0=k0, K=K, ms=ms, lam=lam, theta=theta, X=X, iters=it)
if __name__ == '__main__':
    theta = eval(sys.argv[1]); k0 = int(sys.argv[2]); X = float(sys.argv[3]); mmax = int(sys.argv[4])
    lams = [float(x) for x in sys.argv[5].split(',')]; tag = sys.argv[6]
    for lam in lams:
        t0 = time.time(); d = dual(lam, theta, k0, X, mmax, verbose='-v' in sys.argv)
        if d is None: print(f"lam={lam} LP failed"); continue
        print(f"theta={theta:.7f} k0={k0} X={X} mmax={mmax} lam={lam} t*={d['t']:.6e} finemin={d['finemin']:.6e} iters={d['iters']} c>1e-6 at m={[m for m,cc in zip(d['ms'],d['c']) if cc>1e-6]} ({time.time()-t0:.1f}s) maxrss={rss():.0f}MB", flush=True)
        if d['finemin'] > 0:
            fn = f"logs/dpp_{tag}_lam{lam}.json"; json.dump(d, open(fn, 'w')); print("saved", fn, flush=True)
