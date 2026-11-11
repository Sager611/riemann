# Certify (or refute) a screened counterexample to Jerby's Gram-point repulsion conjecture
# (2406.06548 conj:rep: g_n bad, g_{n-1}, g_{n+1} good  =>  |Z'(g_n)| > 4|Z(g_n)|).
# Arb (python-flint) ball arithmetic at 256 bits: interval-Newton enclosure of each Gram point g_k (theta(g_k) = k pi, theta from
# the continuous log-gamma branch), enclosures of Z(g_k) for k = n-2..n+2 and of Z'(g_n), mu = |Z'/Z|; independent mpmath
# Riemann-Siegel recheck (siegelz, derivative=1) at 30 digits.  usage: certify_cx.py n [n ...]
import sys
from flint import arb, acb, ctx, acb_series
import mpmath
ctx.prec = 256
PI = arb.pi()
def theta(t): return acb(arb(0.25), t/2).lgamma().imag - t/2*PI.log()
def dtheta(t): return acb(arb(0.25), t/2).digamma().real/2 - PI.log()/2
def gram_ball(k, t0):
    t = arb(t0)
    for _ in range(8):
        m = arb(t.mid()); t = arb((m - (theta(m) - k*PI)/dtheta(m)).mid())
    e = arb('1e-45'); X = arb(t.mid(), e); m = arb(t.mid())
    N = m - (theta(m) - k*PI)/dtheta(X)
    assert N.rad() < e and X.contains(N), 'interval Newton failed'
    return N
def ZdZ(t):
    s = acb_series([acb(arb(0.5), t), acb(0, 1)], prec=2); z = s.zeta().coeffs()
    th = theta(t); thp = dtheta(t); e = acb(th.cos(), th.sin())
    return (e*z[0]).real, (e*(acb(0, thp)*z[0] + z[1])).real, (e*z[0]).imag
mpmath.mp.dps = 30
for a in sys.argv[1:]:
    n = int(a); R = {}
    for k in range(n-2, n+3):
        g = gram_ball(k, float(mpmath.grampoint(k)))
        Z, dZ, imZ = ZdZ(g)
        good = bool(Z > 0) if k % 2 == 0 else bool(Z < 0); bad = bool(Z < 0) if k % 2 == 0 else bool(Z > 0)
        R[k] = (g, Z, dZ, good, bad)
        print(f'k={k} g={g.str(30)} Z={Z.str(16)} good={good} bad={bad} |Im e^(i theta) zeta|<={float(abs(imZ).upper()):.1e}')
    g, Z, dZ, good, bad = R[n]
    mu = abs(dZ)/abs(Z)
    iso = R[n-1][3] and R[n][4] and R[n+1][3]
    print(f"  Z'(g_{n}) in {dZ.str(16)};  4|Z(g_n)| in {(4*abs(Z)).str(16)}")
    print(f'  mu(g_{n}) in {mu.str(16)};  certified mu<4: {bool(mu < 4)};  certified mu>4: {bool(mu > 4)};  isolated bad (certified): {iso}')
    gm = mpmath.grampoint(n); Zm = mpmath.siegelz(gm); dZm = mpmath.siegelz(gm, derivative=1)
    Zl = mpmath.siegelz(mpmath.grampoint(n-1)); Zr = mpmath.siegelz(mpmath.grampoint(n+1))
    print(f'  mpmath: g_n={mpmath.nstr(gm, 25)} Z(g_n-1)={mpmath.nstr(Zl, 12)} Z(g_n)={mpmath.nstr(Zm, 15)} Z(g_n+1)={mpmath.nstr(Zr, 12)} Zp(g_n)={mpmath.nstr(dZm, 15)} mu={mpmath.nstr(abs(dZm/Zm), 12)}', flush=True)
