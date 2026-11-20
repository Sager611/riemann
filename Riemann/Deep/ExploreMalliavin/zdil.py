# explore62 (c): zeta's own dilution at finite height. On a window of ~W zeros at height t0 (unfolded grid, 16 points per unit),
# cut Z to its large values: E = {|F| >= tau}; on each run of E^c keep the sign if both ends agree (f = +-tau, no zero),
# else join linearly (one zero). Report: zeros kept, relative L2 error, worst second-order-mean error, fourth-order errors, M1.
# F(v) = Z(t)/sqrt(2 theta'(t)), v = theta(t)/pi (explore60's unfolding). zeta via Riemann-Siegel (main sum + C0), DH via oscz copy.
import sys, numpy as np
from scipy.special import loggamma
PI = np.pi
def th(t): return loggamma(0.25 + 0.5j*t).imag - 0.5*t*np.log(PI)
def thp(t): return 0.5*np.log(t/(2*PI)) + 1/(48*t*t)   # theta'(t) to O(t^-2)
def zrs(t, CH=None):
    out = np.empty(len(t)); CH = CH or max(8, int(4e5/np.sqrt(t.max()/(2*PI))))
    for i in range(0, len(t), CH):
        tt = t[i:i+CH]; a = np.sqrt(tt/(2*PI)); N = int(a.max()); n = np.arange(1, N+1)
        ph = th(tt)[:, None] - tt[:, None]*np.log(n)[None, :]
        mask = n[None, :] <= np.floor(a)[:, None]
        s = 2*(np.cos(ph)*mask/np.sqrt(n)[None, :]).sum(1)
        Nn = np.floor(a); p = a - Nn
        C0 = np.cos(2*PI*(p*p - p - 1/16))/np.cos(2*PI*p)
        out[i:i+CH] = s + (-1)**(Nn-1)*a**-0.5*C0
    return out
def vgrid(t0, W, step=1/16):
    v0 = th(np.array([t0]))[0]/PI; v = v0 + np.arange(0, W, step); t = np.full(len(v), float(t0))
    for _ in range(6): t = t - (th(t) - PI*v)/thp(t)
    return v, t
def stats(F, S=16):
    r2 = (F*F).mean(); r4 = (F**4).mean()
    C = np.array([(F[:len(F)-u]*F[u:]).mean() for u in range(0, 8*S+1)])
    K = {x: (F[:len(F)-int(x*S)]**2*F[int(x*S):]**2).mean() for x in (0.5, 1, 2, 4)}
    return r2, r4, C, K, np.abs(F).mean()/np.sqrt(r2)
def nz(F): return int((np.sign(F[:-1])*np.sign(F[1:]) < 0).sum())
def cut(F, tau):
    f = F.copy(); small = np.abs(F) < tau; n = len(F); i = 0
    while i < n:
        if not small[i]: i += 1; continue
        j = i
        while j+1 < n and small[j+1]: j += 1
        l = F[i-1] if i > 0 else None; r = F[j+1] if j+1 < n else None
        if l is None: l = r
        if r is None: r = l
        if np.sign(l) == np.sign(r): f[i:j+1] = np.sign(l)*tau
        else: f[i:j+1] = l + (r - l)*(np.arange(i, j+1) - (i-1))/((j+1) - (i-1))
        i = j+1
    return f
def report(name, F, S=16):
    r2, r4, C, K, m1 = stats(F, S); z = nz(F); rms = np.sqrt(r2)
    print(f"{name}: points={len(F)} zeros={z} (per unit {z*S/len(F):.4f}) <F^2>={r2:.4f} M4={r4/r2**2:.3f} M1={m1:.4f} "
          f"G(.5,1,2,4)=" + ",".join(f"{K[x]/r4:.4f}" for x in K) + f" K(4)={K[4]/r2**2:.3f}", flush=True)
    for a in (0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5):
        tau = a*rms; f = cut(F, tau); E = np.abs(F) >= tau
        r2f, r4f, Cf, Kf, m1f = stats(f, S)
        e2 = ((f-F)**2).mean()/r2; eC = np.abs(Cf - C).max()/r2
        eK = max(abs(Kf[x]/r4f - K[x]/r4) for x in K)
        print(f"  tau={a:.1f}rms: zeros kept {nz(f)/z:.4f}; |E| {E.mean():.4f}; E carries {(F[E]**2).sum()/(F*F).sum():.4f} of <F^2>, "
              f"{(F[E]**4).sum()/(F**4).sum():.4f} of <F^4>; L2 err {e2:.4f}; max|C_f-C_F|/<F^2> {eC:.4f}; M4 {r4f/r2f**2/(r4/r2**2)-1:+.4f} rel; "
              f"max|G_f-G_F| {eK:.4f}; M1 {m1f:.4f}", flush=True)
if __name__ == '__main__':
    W = float(sys.argv[2]); which = sys.argv[1]
    if which == 'zeta':
        for t0 in [float(x) for x in sys.argv[3].split(',')]:
            v, t = vgrid(t0, W); z = zrs(t); report(f"zeta t0={t0:.0e} L={np.log(t0/(2*PI)):.2f}", z/np.sqrt(2*thp(t)))
    else:
        sys.path.insert(0, '.'); import oscz_copy as O
        which = {'zetaem': 'zeta'}.get(which, which)  # Euler-Maclaurin zeta (accurate at low t)
        t0 = float(sys.argv[3]); v = None
        tt = np.array([t0]); v0 = O.theta(which, tt)[0]/PI; vv = v0 + np.arange(0, W, 1/16); t = np.full(len(vv), t0)
        thw = lambda t: (O.theta(which, t + 1e-5) - O.theta(which, t - 1e-5))/2e-5  # this function's theta'
        for _ in range(12): t = t - (O.theta(which, t) - PI*vv)/thw(t)
        print(f"{which}: unfolding residual max|theta(t)/pi - v| = {np.abs(O.theta(which, t)/PI - vv).max():.1e}")
        z = np.empty(len(t)); im = np.empty(len(t))
        for k in range(0, len(t), 512): z[k:k+512], im[k:k+512] = O.Z(which, t[k:k+512], chunk=128)
        print(f"{which}: max|Im|/max|Z| = {np.abs(im).max()/np.abs(z).max():.1e}")
        report(f"{which} t0={t0:.0e}", z/np.sqrt(2*thw(t)))
