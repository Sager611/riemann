# explore18 test T2: class AIII/BDI winding (Toeplitz-type index) of a companion E with Re E = F on the critical line.
#   Exact companion  E_c = F + c F'       (on the line: A - i c A', A = F(1/2+it) real)
#   Levinson's       E_L = (1 - lam/L) F + F'/L,  lam = gamma-factor log-derivative (= h * (zeta + zeta'/L) for xi)
# Index bound: N_0(window) >= |W| - 1, W = (1/pi) * change of arg E along the line.
# Bulk side: W = N_total - 2*kappa + b, kappa = zeros of E in Re s > 1/2 (window), b = horizontal-edge (eta-like) term.
import sys, time
from topo_common import *
from topo_common import _seg
mp.mp.dps = 20

def lam_xi(s):
    return hlogd(s)

def lam_dh(s):
    return mp.log(5 / mp.pi) / 2 + mp.digamma((s + 1) / 2) / 2

def sample(F, Fp, lam, pts):
    out = []
    for z in pts:
        out.append((F(z), Fp(z), lam(z)))
    return out

def edge_pts(z0, z1, step):
    n = max(1, int(mp.ceil(abs(z1 - z0) / step)))
    return [z0 + (z1 - z0) * k / n for k in range(n + 1)]

def winding(E, F, Fp, lam, pts, vals, tol=0.8):
    tot = mp.mpf(0); nref = 0
    ev = [E(v[0], v[1], v[2]) for v in vals]
    for k in range(len(pts) - 1):
        d = mp.im(mp.log(ev[k + 1] / ev[k]))
        if abs(d) > tol:
            g = lambda z: E(F(z), Fp(z), lam(z))
            d = _seg(g, pts[k], pts[k + 1], ev[k], ev[k + 1], 0.3, 0); nref += 1
        tot += d
    return tot, nref

def run(name, t0, t1, cs, Lq):
    F, Fp = FUNCS[name]
    lam = lam_dh if name == 'dh' else lam_xi
    t0 = mp.mpf(t0); t1 = mp.mpf(t1)
    L = mp.log(Lq * (t0 + t1) / 2 / (2 * mp.pi))
    # line (upward), then the rest of the right rectangle [1/2, 3] x [t0, t1], counterclockwise = down the line
    line = edge_pts(mp.mpc(0.5, t0), mp.mpc(0.5, t1), 0.01)
    top = edge_pts(mp.mpc(0.5, t1), mp.mpc(3, t1), 0.05)
    right = edge_pts(mp.mpc(3, t1), mp.mpc(3, t0), 0.05)
    bot = edge_pts(mp.mpc(3, t0), mp.mpc(0.5, t0), 0.05)
    T = time.time()
    vl = sample(F, Fp, lam, line)
    vt = sample(F, Fp, lam, top); vr = sample(F, Fp, lam, right); vb = sample(F, Fp, lam, bot)
    # on-line zeros of F: sign changes of Re F
    re = [mp.re(v[0]) for v in vl]
    n0 = sum(1 for k in range(len(re) - 1) if re[k] * re[k + 1] < 0)
    # zeros of F in Re s > 1/2: rectangle [0.52, 3] x [t0, t1]
    nr = zeros_in_rect(F, 0.52, 3, t0, t1, 0.05)
    ntot = n0 + 2 * int(mp.nint(nr))
    print(f"== {name} window [{float(t0)}, {float(t1)}]  N_0 = {n0}  zeros of F right of line = {mp.nstr(nr, 6)}  N_total = {ntot}  L = {mp.nstr(L, 5)}  ({time.time() - T:.0f} s)", flush=True)
    comps = [(f"F+{c}F'", (lambda c: (lambda a, b, l: a + c * b))(mp.mpf(c))) for c in cs]
    comps.append(("Levinson E_L", lambda a, b, l: (1 - l / L) * a + b / L))
    for nm, E in comps:
        wl, r1 = winding(E, F, Fp, lam, line, vl)
        W = wl / mp.pi
        other = mp.mpf(0); r2 = 0
        for pts, vals in ((top, vt), (right, vr), (bot, vb)):
            w, rr = winding(E, F, Fp, lam, pts, vals); other += w; r2 += rr
        kap = -(other + wl) / (2 * mp.pi)  # path = line up, top rightward, right edge down, bottom leftward = clockwise
        print(f"   {nm:14s} W = {mp.nstr(W, 6):>9s}  kappa(right) = {mp.nstr(kap, 6):>9s}  N_total - 2 kappa - W = {mp.nstr(ntot - 2 * mp.nint(kap) - W, 4):>7s}  index bound N_0 >= W - 1: {n0 >= W - 1}  refinements {r1 + r2}", flush=True)

if __name__ == '__main__':
    which = sys.argv[1]
    cs = ['0.05', '0.2', '1', '5']
    if which == 'xi':
        run('xi', 10, 60, cs, 1)
    elif which == 'mv':
        run('mv', 10, 60, cs, 1)
    elif which == 'dh1':
        run('dh', 83, 88, cs, 5)
    elif which == 'dh2':
        run('dh', 112, 116.5, cs, 5)
