# certify_conf.py -- jerbyconf: Arb (python-flint, 200 bits) certificate for a confined isolated bad Gram point:
# enclosures of g_{n-1}, g_n, g_{n+1} (interval Newton on theta(g) = k pi), Z there (isolated bad), Z'(g_n), mu = |Z'/Z| and
# s*mu = 2 pi mu / log(g_n / 2 pi) (Jerby 2610.21952 l. 657-663); plus Arb signs of Z at given interior points (midpoints between the
# two screened zeros), certifying >= 2 zeros in that Gram interval.  mpmath siegelz recheck.  usage: certify_conf.py n[:t1[:t2...]] ...
import sys, os
from flint import arb, acb, ctx, acb_series
import mpmath
ctx.prec = 200
PI = arb.pi()
def theta(t): return acb(arb(0.25), t/2).lgamma().imag - t/2*PI.log()
def dtheta(t): return acb(arb(0.25), t/2).digamma().real/2 - PI.log()/2
def gram_ball(k, t0):
    t = arb(t0)
    for _ in range(8):
        m = arb(t.mid()); t = arb((m - (theta(m) - k*PI)/dtheta(m)).mid())
    e = arb('1e-40'); X = arb(t.mid(), e); m = arb(t.mid())
    N = m - (theta(m) - k*PI)/dtheta(X)
    assert N.rad() < e and X.contains(N), 'interval Newton failed'
    return N
def ZdZ(t):
    s = acb_series([acb(arb(0.5), t), acb(0, 1)], prec=2); z = s.zeta().coeffs()
    th = theta(t); thp = dtheta(t); e = acb(th.cos(), th.sin())
    return (e*z[0]).real, (e*(acb(0, thp)*z[0] + z[1])).real
def Zpt(t):
    th = theta(t); return (acb(th.cos(), th.sin())*acb(arb(0.5), t).zeta()).real
mpmath.mp.dps = 30
for a in sys.argv[1:]:
    parts = a.split(':'); n = int(parts[0]); R = {}
    for k in (n-1, n, n+1):
        g = gram_ball(k, float(mpmath.grampoint(k))); Z, dZ = ZdZ(g)
        good = bool(Z > 0) if k % 2 == 0 else bool(Z < 0); bad = bool(Z < 0) if k % 2 == 0 else bool(Z > 0)
        R[k] = (g, Z, dZ, good, bad)
    g, Z, dZ, good, bad = R[n]
    mu = abs(dZ)/abs(Z); smu = 2*PI*mu/(g/(2*PI)).log()
    iso = R[n-1][3] and R[n][4] and R[n+1][3]
    print(f'n={n} g_n={g.str(25)} Z(g_n-1)={R[n-1][1].str(12)} Z(g_n)={Z.str(15)} Z(g_n+1)={R[n+1][1].str(12)} isolated_bad(certified)={iso}')
    print(f"  Z'(g_n)={dZ.str(15)} mu={mu.str(15)} s*mu={smu.str(15)}  certified s*mu<2: {bool(smu < 2)}  certified mu>4: {bool(mu > 4)}")
    for p in parts[1:]:
        tp = arb(p); zv = Zpt(tp)
        side = 'left' if bool(tp < g) else 'right'
        opp = bool(zv*Z < 0)
        print(f'  Z({p})={zv.str(12)} ({side} interval) sign opposite to Z(g_n): {opp}' + ('  => >= 2 zeros in that Gram interval (certified)' if opp else ''))
    if not os.environ.get('SKIP_MPMATH'):
        gm = mpmath.grampoint(n); Zm = mpmath.siegelz(gm); dZm = mpmath.siegelz(gm, derivative=1)
        print(f'  mpmath: Z(g_n)={mpmath.nstr(Zm, 15)} Zp={mpmath.nstr(dZm, 15)} mu={mpmath.nstr(abs(dZm/Zm), 12)}')
    sys.stdout.flush()
