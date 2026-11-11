# T2b: odd sector = 'roll' problem with NEGATIVE lift: Q_odd = -2 s s^T + R, s_n = <p_n, sinh(x/2)>, R = A - P.
# Munk: sup_f 2<f,sinh>^2 / R(f) = 2 s^T R^{-1} s (maximal lift^2/drag, optimiser R^{-1}s); Q_odd >= 0 <=> R > 0 and 2 s^T R^{-1} s <= 1.
# Also: tip value of the even ground state vs basis size K (Kutta-type edge behaviour).
import numpy as np, sys
sys.path.insert(0, '.')
from weilleg import WeilLeg
from scipy.special import spherical_in
class WeilLegOdd(WeilLeg):
    def __init__(s, a, K, **kw):
        super().__init__(a, K, **kw); s.deg = 2*np.arange(K) + 1; s.nrm = np.sqrt((2*s.deg + 1)/(2*a))
        dmax = s.deg[-1]; s.nin = dmax + 40; s.nout = 2*dmax + 60
        from numpy.polynomial import legendre as Lg; s.ti, s.wi = Lg.leggauss(s.nin)
    def pole(s): return s.nrm*2*s.a*spherical_in(s.deg, s.a/2)     # <p_n, sinh(x/2)> for odd n
    def build(s):
        super().build(); s.Q = -2*np.outer(s.l, s.l) + s.A - s.P; return s
print(" a      c     lam0_odd    lam0(R)    Munk eff 2s^T R^-1 s   1-eff      cos(opt, x*sqrt(1-x^2/a^2))")
for a in [0.1, 0.2, 0.3, np.log(2)/2, 0.4, 0.45, 0.5, np.log(3)/2, 0.6, 0.65]:
    W = WeilLegOdd(a, 20).build()
    R = W.A - W.P; eR = np.linalg.eigvalsh(R); eff = 2*W.l @ np.linalg.solve(R, W.l)
    opt = np.linalg.solve(R, W.l); opt /= np.linalg.norm(opt)
    xg, wg = np.polynomial.legendre.leggauss(300); x = a*xg; wx = a*wg
    roll = (W.vals(x)*wx) @ (x*np.sqrt(1 - (x/a)**2)); roll /= np.linalg.norm(roll)
    print(f"{a:5.3f} {np.exp(2*a):6.3f} {np.linalg.eigvalsh(W.Q)[0]:11.4e} {eR[0]:10.3e}   {eff:.12f}   {1-eff:10.3e}   {abs(opt @ roll):.6f}")
print("Even ground state tip value f0(a)/f0(0) vs number of even Legendre modes K:")
for a in [0.1, 0.3, 0.45]:
    out = []
    for K in [6, 12, 24, 48, 96]:
        W = WeilLeg(a, K).build(); e, V = np.linalg.eigh(W.Q); c0 = V[:, 0]
        f_tip = c0 @ W.nrm; f_mid = c0 @ W.vals(np.array([0.0]))[:, 0]
        out.append(f"K={K}: {f_tip/f_mid:.5f} (lam0 {e[0]:.8e})")
    print(f"  a={a}: " + '; '.join(out))
