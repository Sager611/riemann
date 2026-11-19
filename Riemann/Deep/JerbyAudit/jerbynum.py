# JerbyAudit: which numerical Z reproduces Jerby's own quoted viscosities, and are our counterexamples robust under it?
# V0 exact (mpmath siegelz, derivative=1); V1 Riemann-Siegel main sum 2*sum_{k<=a} only (no remainder);
# V2 Jerby's Z_N(t;1) = cos theta + sum_{k=1}^{N} (k+1)^(-1/2) cos(theta - ln(k+1) t), N = [t/2] (2406.06548 abstract, l. 235, 395).
import sys, time, numpy as np
import mpmath as mp
mp.mp.dps = 30
PI = mp.pi
def gram(n):
    x0 = 2*PI*mp.exp(1 + mp.lambertw((8*n + 1)/(8*mp.e)).real)
    return mp.findroot(lambda t: mp.siegeltheta(t) - n*PI, x0)
def partial(t, K, start):   # sum_{m=start}^{K} m^(-1/2) cos(th - t ln m), and its t-derivative; phases reduced in high precision per chunk
    th = mp.siegeltheta(t); thp = mp.diff(mp.siegeltheta, t)
    tf = float(t); S = 0.0; dS = 0.0; CH = 2_000_000
    for lo in range(start, K + 1, CH):
        hi = min(K, lo + CH - 1)
        # phase(m) = th - t ln m ; write t = t0 + e (t0 = float), reduce th - t0*ln(lo) exactly, then add -t0*ln(m/lo) - e*ln m
        t0 = mp.mpf(tf); e = float(t - t0)
        base = float(mp.fmod(th - t0*mp.log(lo), 2*PI))
        m = np.arange(lo, hi + 1, dtype=np.float64)
        lr = np.log1p((m - lo)/lo)                       # ln(m/lo), accurate
        ph = base - (tf*lr - np.floor(tf*lr/(2*np.pi))*2*np.pi) - e*np.log(m)
        w = 1/np.sqrt(m); S += float(np.sum(w*np.cos(ph))); dS += float(np.sum(-w*np.sin(ph)*(float(thp) - np.log(m))))
    return S, dS, th, thp
def variants(t):
    Z0 = mp.siegelz(t); dZ0 = mp.siegelz(t, derivative=1)
    a = int(mp.floor(mp.sqrt(t/(2*PI))))
    S1, dS1, th, thp = partial(t, a, 1)
    N = int(mp.floor(t/2))
    S2, dS2, _, _ = partial(t, N + 1, 2)
    Z2 = float(mp.cos(th)) + S2; dZ2 = float(-mp.sin(th)*thp) + dS2
    return (float(Z0), float(dZ0)), (2*S1, 2*dS1), (Z2, dZ2)
which = sys.argv[1]
pts = {'jerby': [(730119, 4.46023), (9807962, 0.0750883)], 'ours': [(1565734083, None)]}[which]
for n, quoted in pts:
    t0 = time.time(); g = gram(n)
    out = variants(g)
    s = '  '.join(f'{nm}: Z={Z:.10g} dZ={dZ:.10g} mu={abs(dZ/Z):.7f}' for nm, (Z, dZ) in zip(['V0 exact', 'V1 RS main sum', 'V2 Jerby Z_[t/2]'], out))
    print(f'n={n} g={mp.nstr(g, 15)} quoted mu={quoted}: {s}  [{time.time() - t0:.1f} s]', flush=True)
