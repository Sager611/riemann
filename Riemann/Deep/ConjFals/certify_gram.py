# Certify a Gram-point viscosity value in ball arithmetic (python-flint/Arb) and cross-check with mpmath.
# usage: certify.py n [n ...]   prints, for g_{n-1}, g_n, g_{n+1}: certified Z, and for g_n: Z', mu = |Z'/Z|
import sys
from flint import arb, acb, ctx, acb_series
import mpmath
ctx.prec = 200
PI = arb.pi()
def theta(t):            # t: arb
    return acb(arb(0.25), t/2).lgamma().imag - t/2*PI.log()
def dtheta(t):
    return acb(arb(0.25), t/2).digamma().real/2 - PI.log()/2
def gram_ball(n, t0):
    t = arb(t0)
    for _ in range(6):   # Newton on midpoints
        t = arb(t.mid()) - (theta(arb(t.mid())) - n*PI)/dtheta(arb(t.mid()))
        t = arb(t.mid())
    # interval Newton step on [t - e, t + e]
    e = arb('1e-40')
    X = arb(t.mid(), e)
    N = arb(t.mid()) - (theta(arb(t.mid())) - n*PI)/dtheta(X)
    assert N.rad() < e and X.contains(N), 'interval Newton failed'
    return N
def Z_and_dZ(t):         # t: arb ball; returns enclosures of Z(t), Z'(t)
    s = acb_series([acb(arb(0.5), t), acb(0, 1)], prec=2)   # s(x) = 1/2 + i(t + x)
    z = s.zeta().coeffs()          # zeta(s(x)) = z0 + z1 x
    th = theta(t); thp = dtheta(t)
    e = acb(th.cos(), th.sin())
    Z = (e*z[0]).real
    dZ = (e*(acb(0, thp)*z[0] + z[1])).real
    imZ = (e*z[0]).imag
    return Z, dZ, imZ
mpmath.mp.dps = 40
for arg in sys.argv[1:]:
    n = int(arg)
    out = []
    for k in (n-1, n, n+1):
        t0 = float(mpmath.grampoint(k))
        g = gram_ball(k, t0)
        Z, dZ, imZ = Z_and_dZ(g)
        good = (Z > 0) if k % 2 == 0 else (Z < 0)
        bad = (Z < 0) if k % 2 == 0 else (Z > 0)
        out.append((k, g, Z, dZ, good, bad))
        print(f'n={k} g={g.str(25)} Z={Z.str(15)} dZ={dZ.str(15)} good={good} bad={bad} |Im|<={float(abs(imZ).upper()):.1e}')
    k, g, Z, dZ, good, bad = out[1]
    mu = abs(dZ)/abs(Z)
    print(f'  mu(g_{n}) in {mu.str(15)}; certified mu<4: {bool(mu < 4)}; isolated bad: {bool(out[0][4] and out[1][5] and out[2][4])}')
    # independent cross-check (mpmath Riemann-Siegel)
    gm = mpmath.grampoint(n)
    Zm = mpmath.siegelz(gm); dZm = mpmath.siegelz(gm, derivative=1)
    print(f'  mpmath: g={mpmath.nstr(gm,25)} Z={mpmath.nstr(Zm,15)} dZ={mpmath.nstr(dZm,15)} mu={mpmath.nstr(abs(dZm/Zm),12)}')
