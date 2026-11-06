# explore18 common: completed functions with a functional equation F(s) = F(1-s), real on the critical line.
#   xi   : Riemann xi.
#   mv   : xi with the 6th/7th zero pair moved off the line to 1/2 +- 0.15 + i*39.2524 (FE and reality kept, no Euler product).
#   dh   : Davenport-Heilbronn Lambda_f(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) f(s) (FE, real on the line, zeros off the line).
# Each name X has X(s) and Xp(s) = d/ds X(s).
import mpmath as mp
mp.mp.dps = 25

def h(s):
    return s * (s - 1) / 2 * mp.power(mp.pi, -s / 2) * mp.gamma(s / 2)

def hlogd(s):
    return 1 / s + 1 / (s - 1) - mp.log(mp.pi) / 2 + mp.digamma(s / 2) / 2

def xi(s):
    return h(s) * mp.zeta(s)

def xip(s):
    return h(s) * (hlogd(s) * mp.zeta(s) + mp.zeta(s, 1, 1))

GA = mp.mpf('37.586178158825671257217763481')   # 6th zero of zeta
GB = mp.mpf('40.918719012147495187398126914')   # 7th zero
GBAR = (GA + GB) / 2
DELTA = mp.mpf('0.15')
NUM = [mp.mpc(0.5 + DELTA, GBAR), mp.mpc(0.5 + DELTA, -GBAR), mp.mpc(0.5 - DELTA, GBAR), mp.mpc(0.5 - DELTA, -GBAR)]
DEN = [mp.mpc(0.5, GA), mp.mpc(0.5, -GA), mp.mpc(0.5, GB), mp.mpc(0.5, -GB)]

def Q(s):
    return mp.fprod([s - r for r in NUM]) / mp.fprod([s - r for r in DEN])

def Qlogd(s):
    return sum(1 / (s - r) for r in NUM) - sum(1 / (s - r) for r in DEN)

def mv(s):
    return xi(s) * Q(s)

def mvp(s):
    q = Q(s)
    return xip(s) * q + xi(s) * q * Qlogd(s)

KAP = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
ACOEF = {1: mp.mpf(1), 2: KAP, 3: -KAP, 4: mp.mpf(-1)}

def dh_f(s):
    return mp.power(5, -s) * sum(ACOEF[r] * mp.zeta(s, mp.mpf(r) / 5) for r in ACOEF)

def dh_fp(s):
    return mp.power(5, -s) * sum(ACOEF[r] * (mp.zeta(s, mp.mpf(r) / 5, 1) - mp.log(5) * mp.zeta(s, mp.mpf(r) / 5)) for r in ACOEF)

def dh_g(s):
    return mp.power(5 / mp.pi, (s + 1) / 2) * mp.gamma((s + 1) / 2)

def dh(s):
    return dh_g(s) * dh_f(s)

def dhp(s):
    return dh_g(s) * ((mp.log(5 / mp.pi) / 2 + mp.digamma((s + 1) / 2) / 2) * dh_f(s) + dh_fp(s))

FUNCS = {'xi': (xi, xip), 'mv': (mv, mvp), 'dh': (dh, dhp)}

# ---- argument principle: total change of arg g along a polyline, adaptive bisection ----
def _seg(g, z0, z1, g0, g1, tol, depth):
    d = mp.im(mp.log(g1 / g0))
    if depth >= 40:
        return d
    zm = (z0 + z1) / 2
    gm = g(zm)
    d1 = mp.im(mp.log(gm / g0)); d2 = mp.im(mp.log(g1 / gm))
    if abs(d) < tol and abs(d1 + d2 - d) < 1e-8:
        return d
    return _seg(g, z0, zm, g0, gm, tol, depth + 1) + _seg(g, zm, z1, gm, g1, tol, depth + 1)

def darg(g, z0, z1, step=0.25, tol=0.3):
    n = max(1, int(mp.ceil(abs(z1 - z0) / step)))
    pts = [z0 + (z1 - z0) * k / n for k in range(n + 1)]
    vals = [g(z) for z in pts]
    return sum(_seg(g, pts[k], pts[k + 1], vals[k], vals[k + 1], tol, 0) for k in range(n))

def zeros_in_rect(g, s0, s1, t0, t1, step=0.25):
    # counterclockwise boundary of [s0, s1] x [t0, t1]
    c = [mp.mpc(s0, t0), mp.mpc(s1, t0), mp.mpc(s1, t1), mp.mpc(s0, t1), mp.mpc(s0, t0)]
    w = sum(darg(g, c[k], c[k + 1], step) for k in range(4))
    return w / (2 * mp.pi)

def line_sign_changes(F, t0, t1, step=0.02):
    n = int(mp.ceil((t1 - t0) / step)); cnt = 0; prev = None; locs = []
    for k in range(n + 1):
        t = t0 + (t1 - t0) * mp.mpf(k) / n
        v = mp.re(F(mp.mpc(0.5, t)))
        if prev is not None and v * prev < 0:
            cnt += 1; locs.append(float(t))
        prev = v
    return cnt, locs
