
# B2: scan t in [T0,T1]: real zeros of the AFE main term Z_main (Spira's natural range N=floor(sqrt(t/2pi)))
#     vs Riemann-Siegel with C0 correction (proxy for true Z). Deficits => off-line zeros of F_N; verify with mpmath.
import numpy as np, sys, time, mpmath as mp
T0, T1, H = float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])
t0c = time.time()
def theta(t):
    return t/2*np.log(t/(2*np.pi)) - t/2 - np.pi/8 + 1/(48*t) + 7/(5760*t**3)
def C0(p):
    c = np.cos(2*np.pi*p); num = np.cos(2*np.pi*(p*p - p - 1/16))
    bad = np.abs(c) < 1e-7
    out = np.where(bad, 0.0, num/np.where(bad, 1.0, c))
    if bad.any():  # removable singularity: use nearby value
        pp = p[bad] + 1e-6; out[bad] = np.cos(2*np.pi*(pp*pp-pp-1/16))/np.cos(2*np.pi*pp)
    return out
cands = []; nsc_main = 0; nsc_rs = 0
CH = 400000
ts_all = np.arange(T0, T1, H)
for i0 in range(0, len(ts_all), CH):
    t = ts_all[i0:i0+CH+1]
    N = np.floor(np.sqrt(t/(2*np.pi))).astype(int)
    th = theta(t)
    Zm = np.zeros_like(t)
    for n in range(1, N.max()+1):
        msk = N >= n
        Zm[msk] += 2*np.cos(th[msk] - t[msk]*np.log(n))/np.sqrt(n)
    a = np.sqrt(t/(2*np.pi)); p = a - N
    Zr = Zm + (-1.0)**(N-1)*a**(-0.5)*C0(p)
    sm = np.signbit(Zm); sr = np.signbit(Zr)
    chm = np.nonzero(sm[1:] != sm[:-1])[0]; chr_ = np.nonzero(sr[1:] != sr[:-1])[0]
    nsc_main += len(chm); nsc_rs += len(chr_)
    # candidate: RS has two sign changes within 0.3 where main has none (pair lost) -- compare locally
    setm = t[chm]; 
    for j in range(len(chr_)-1):
        ta, tb = t[chr_[j]], t[chr_[j+1]]
        if tb - ta < 0.3:
            # main-term sign changes in [ta-0.05, tb+0.05]
            k = np.sum((setm > ta-0.05) & (setm < tb+0.05))
            if k == 0: cands.append((float(ta), float(tb), int(N[chr_[j]])))
print('grid [%g,%g] step %g: sign changes main=%d, RS(C0)=%d, lost-pair candidates=%d  (%.0fs)' % (T0,T1,H,nsc_main,nsc_rs,len(cands),time.time()-t0c), flush=True)
mp.mp.dps = 25
def F(s, N):
    chi = mp.power(2, s)*mp.power(mp.pi, s-1)*mp.sin(mp.pi*s/2)*mp.gamma(1-s)
    return mp.fsum(mp.power(n, -s) for n in range(1, N+1)) + chi*mp.fsum(mp.power(n, s-1) for n in range(1, N+1))
for (ta, tb, N) in cands[:40]:
    tm = (ta+tb)/2
    zt = [mp.siegelz(ta+0.001), mp.siegelz(tm), mp.siegelz(tb+0.001)]
    best = None
    for d0 in (0.02, 0.05, 0.1):
        try:
            r = mp.findroot(lambda s: F(s, N), mp.mpc(0.5+d0, tm))
            if abs(F(r, N)) < 1e-12 and abs(r.imag - tm) < 0.3 and abs(r.real-0.5) > 1e-6: best = r; break
        except Exception: pass
    print('pair %.4f-%.4f N=%d  Z(mid)=%s  offline zero of F_N: %s' % (ta, tb, N, mp.nstr(zt[1],4), ('sigma-1/2=%s t=%s' % (mp.nstr(best.real-0.5,6), mp.nstr(best.imag,10))) if best is not None else 'not found'), flush=True)
