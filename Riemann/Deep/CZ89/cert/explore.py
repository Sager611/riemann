
import numpy as np, time, sys
b = 1/np.sqrt(2)
a = np.array([2/25, 1/25, -1/200])
c0 = b/(2*np.sin(b))
W = np.zeros((7,7))
rows = {0:[0.306,0.645,1.000,1.000,1.000,2.000],1:[0.347,0.047,0.000,0.000,1.000],2:[0.347,0.616,0.000,1.000],3:[0.347,0.047,1.000],4:[0.347,0.645],5:[0.306]}
for i,v in rows.items():
    for t,val in enumerate(v):
        W[i,i+1+t] = val
pairs = [(i,j) for i in range(7) for j in range(i+1,7)]
w = np.array([W[i,j] for i,j in pairs])
PI = np.array([i for i,j in pairs]); PJ = np.array([j for i,j in pairs])
A = np.zeros((21,6))
for p,(i,j) in enumerate(pairs):
    A[p, i:j] = 1.0
def s(z):
    return np.sinc(z/np.pi)
def sp(z):
    out = np.empty_like(z)
    small = np.abs(z) < 1e-4
    zs = z[small]; out[small] = -zs/3 + zs**3/30
    zb = z[~small]; out[~small] = (zb*np.cos(zb) - np.sin(zb))/zb**2
    return out
def K(x):
    pi = np.pi
    r = c0*(s(pi*x-b)+s(pi*x+b))
    for k in (1,2,3):
        r = r + a[k-1]/2*(s(pi*(x-k))+s(pi*(x+k)))
    return r
def Kp(x):
    pi = np.pi
    r = pi*c0*(sp(pi*x-b)+sp(pi*x+b))
    for k in (1,2,3):
        r = r + pi*a[k-1]/2*(sp(pi*(x-k))+sp(pi*(x+k)))
    return r
H = 1e-4; XMAX = 40.0
xs = np.arange(0, XMAX+2*H, H)
kk = K(xs); kpp = Kp(xs)
QT = kk*kk; QPT = 2*kk*kpp
def qinterp(d):
    t = d/H
    i = np.minimum(t.astype(np.int64), len(xs)-2)
    f = t - i
    q = QT[i]*(1-f) + QT[i+1]*f
    qp = QPT[i]*(1-f) + QPT[i+1]*f
    return q, qp
def S_and_grad(g):
    y = np.concatenate([np.zeros((g.shape[0],1)), np.cumsum(g,axis=1)], axis=1)
    d = y[:, PJ] - y[:, PI]
    q, qp = qinterp(d)
    S = q @ w + g.sum(axis=1)/200
    G = (qp*w) @ A + 1/200
    return S, G
def S_exact(g):
    g = np.atleast_2d(g)
    y = np.concatenate([np.zeros((g.shape[0],1)), np.cumsum(g,axis=1)], axis=1)
    d = y[:, PJ] - y[:, PI]
    k = K(d)
    return (k*k) @ w + g.sum(axis=1)/200
def optimize(g, iters=3000, lo=0.0, hi=6.0):
    M = g.shape[0]
    step = np.full(M, 0.05)
    S, G = S_and_grad(g)
    for it in range(iters):
        gn = np.clip(g - step[:,None]*G, lo, hi)
        Sn, Gn = S_and_grad(gn)
        acc = Sn < S
        g = np.where(acc[:,None], gn, g)
        S = np.where(acc, Sn, S); G = np.where(acc[:,None], Gn, G)
        step = np.where(acc, step*1.3, step*0.5)
        step = np.clip(step, 1e-12, 1.0)
    return g, S
if __name__ == "__main__":
    seed = int(sys.argv[1]); nb = int(sys.argv[2]); M = int(sys.argv[3]); tmax = float(sys.argv[4])
    rng = np.random.default_rng(seed)
    t0 = time.time()
    best = []
    for batch in range(nb):
        g0 = rng.uniform(0.0, 3.0, size=(M,6))
        g, S = optimize(g0, iters=600)
        idx = np.argsort(S)[:30]
        for i in idx:
            best.append((S[i], tuple(np.round(g[i],6))))
        print("batch", batch, "min", S.min(), "time", round(time.time()-t0,1), flush=True)
        if time.time()-t0 > tmax: break
    best.sort()
    seen = []
    for Sv, gv in best:
        if all(np.abs(np.array(gv)-np.array(x)).max() > 2e-2 for _,x in seen):
            seen.append((Sv,gv))
        if len(seen) >= 15: break
    for Sv, gv in seen:
        print("%.10f exact %.10f" % (Sv, S_exact(np.array(gv))[0]), gv, "sum=%.5f" % sum(gv))

