# uc_main.py - explore2: dBN flow near Lehmer's pair via Cauchy jets; tests T1 (doubling), T2 (Poon frequency), T3 (windowed Laguerre)
import mpmath as mp, time, sys
mp.mp.dps = 80
T0 = time.time()
I = mp.mpc(0, 1)
X0 = mp.mpf('14010.163')            # z = 2*gamma; Lehmer pair gamma ~ 7005.063, 7005.101
RHO = mp.mpf(7); M = 384; NJ = 150; NREL = 120; NMAX = 70; KMAX = 25
def H0(z):
    s = (1 + I*z)/2
    return s*(s-1)/2 * mp.power(mp.pi, -s/2) * mp.gamma(s/2) * mp.zeta(s) / 8
def log(*a):
    print(*a, flush=True)
zs = [X0 + RHO*mp.expj(2*mp.pi*m/M) for m in range(M)]
raw = [H0(z) for z in zs]
S = max(abs(v) for v in raw)
vals = [v/S for v in raw]
log("circle done %.1fs, log10 S = %s" % (time.time()-T0, mp.nstr(mp.log10(S), 8)))
W = [mp.expj(-2*mp.pi*m/M) for m in range(M)]
def jets(vs):
    c = []
    for j in range(NJ+1):
        acc = mp.mpc(0)
        for m in range(M):
            acc += vs[m]*W[(j*m) % M]
        c.append(acc/M/RHO**j)
    return c
cH = jets(vals)
imr = max(abs(mp.im(x)) for x in cH[:NREL]) / max(abs(mp.re(x)) for x in cH[:NREL])
cH = [mp.re(x) for x in cH]
log("jets done %.1fs; max|Im c|/max|Re c| = %s" % (time.time()-T0, mp.nstr(imr, 3)))
A = [[mp.factorial(n+2*k)/(mp.factorial(n)*mp.factorial(k)) for k in range(KMAX+1)] for n in range(NMAX+1)]
def heat(c, t):
    out = []
    for n in range(NMAX+1):
        acc = mp.mpf(0); tk = mp.mpf(1)
        for k in range(KMAX+1):
            if n+2*k > NREL: break
            acc += A[n][k]*tk*c[n+2*k]; tk *= -t
        out.append(acc)
    return out
def pv(c, h):
    acc = mp.mpf(0) if not isinstance(h, mp.mpc) else mp.mpc(0)
    for x in reversed(c): acc = acc*h + x
    return acc
def der(c): return [(n+1)*c[n+1] for n in range(len(c)-1)]
def shift(c, h):
    d = len(c); q = []
    for n in range(d):
        acc = mp.mpf(0); hp = mp.mpf(1)
        for m in range(n, d):
            acc += c[m]*mp.binomial(m, n)*hp; hp *= h
        q.append(acc)
    return q
# --- validation
c0 = heat(cH, 0)
for h in [mp.mpf(-4), mp.mpf('0.3'), mp.mpc(2.5, 1)]:
    d = H0(X0+h)/S; p = pv(c0, h)
    log("valid P0 h=%s rel.err %s" % (mp.nstr(h,3), mp.nstr(abs(p-d)/abs(d), 3)))
# vertical Gaussian formula H_t(x) = E[H_0(x + i sqrt(2t) Y)], Y ~ N(0,1); Gauss-Hermite via numpy nodes
import numpy as np
xs, ws = np.polynomial.hermite.hermgauss(60)
t = mp.mpf('0.1'); h = mp.mpf('0.5')
acc = mp.mpc(0)
for xi, wi in zip(xs, ws):
    acc += mp.mpf(wi) * H0(X0 + h + I*mp.sqrt(2*t)*mp.sqrt(2)*mp.mpf(xi))
acc = acc/mp.sqrt(mp.pi)/S
ct = heat(cH, t)
log("valid heat t=0.1 h=0.5: jets %s  vertical-GH %s  rel.err %s" % (mp.nstr(pv(ct,h),12), mp.nstr(mp.re(acc),12), mp.nstr(abs(pv(ct,h)-acc)/abs(acc),3)))
# --- real zeros at t=0
def real_zeros(c, lo=-6.5, hi=6.5, step=0.002):
    zs_ = []; hs = [mp.mpf(lo) + k*mp.mpf(step) for k in range(int((hi-lo)/step)+1)]
    prev = pv(c, hs[0])
    for k in range(1, len(hs)):
        cur = pv(c, hs[k])
        if prev == 0 or prev*cur < 0:
            r = mp.findroot(lambda x: pv(c, x), (hs[k-1], hs[k]), solver='anderson')
            zs_.append(r)
        prev = cur
    return zs_
Z0 = real_zeros(c0)
log("zeros of H_0 in [-6.5,6.5]: " + ", ".join("%s(g=%s)" % (mp.nstr(z,8), mp.nstr((X0+z)/2,10)) for z in Z0))
log("time %.1fs" % (time.time()-T0))
import pickle
pickle.dump({'cH': [str(x) for x in cH], 'S': str(S), 'X0': str(X0), 'raw': [str(v) for v in raw], 'Z0': [str(z) for z in Z0]}, open('jets.pkl', 'wb'))
log("saved")
