#!/usr/bin/env python3
"""q2_decay.py -- explore16, follow-up of q1_margins.py: decay of the Jensen-Laguerre margin s*(d,0) in d for Xi and
Davenport-Heilbronn (d <= 64), a fine grid for the Rogers-Szego threshold ratio eps_RS / (2 s*), and the near-circle zeros of
the Rogers-Szego-Jensen polynomials for DH.
Run: uv run --no-project --with python-flint python3 -u q2_decay.py
"""
import time, math
from flint import arb, acb, acb_poly, ctx
from q1_margins import gam_xi, gam_dh, hyp, margin, rs, real_count

def main():
    t0 = time.time()
    N = 140
    ctx.prec = 6000; ctx.cap = 2 * N + 2
    gx = gam_xi(N); gd = gam_dh(N)
    assert all(x > 0 for x in gx) and all(x > 0 for x in gd)
    ctx.prec = 1600
    print(f'gammas ready ({time.time() - t0:.0f} s)')
    print(' d    Xi s*(d,0)   2d s*     DH s*(d,0)   2d s*')
    for d in (24, 32, 40, 48, 56, 64):
        sx = margin(gx, d, 0); sd = margin(gd, d, 0)
        print(f'{d:3d}  {sx:11.4e}  {2 * d * sx:7.4f}   {sd:11.4e}  {2 * d * sd:7.4f}   ({time.time() - t0:.0f} s)')
    print('\nfine ratio eps_RS/(2 s*): bisection on the first loss of hyperbolicity of RS_{e^-eps}')
    for name, g in (('Xi', gx), ('DH', gd)):
        for (d, n) in ((3, 0), (8, 0), (24, 0), (8, 20), (24, 40)):
            s = margin(g, d, n)
            hi = 2.2 * s; lo = 1.9 * s
            assert hyp(rs(g, d, n, lo)) and not hyp(rs(g, d, n, hi))
            for _ in range(30):
                mid = (lo + hi) / 2
                if hyp(rs(g, d, n, mid)): lo = mid
                else: hi = mid
            print(f'  {name} d={d:2d} n={n:2d}: s* = {s:.5e}, eps_RS = {lo:.5e}, ratio = {lo / (2 * s):.4f}')
    print('\nRS_q zeros for DH, d = 12, n = 0:')
    for eps in (0.5, 2.0):
        c = rs(gd, 12, 0, eps)
        rr = acb_poly([acb(a.mid()) for a in c]).roots(tol=arb(2) ** (-300), maxprec=6400)
        mods = [math.log(abs(complex(float(r.real.mid()), float(r.imag.mid())))) for r in rr]
        angs = [abs(math.atan2(float(r.imag.mid()), -float(r.real.mid()))) for r in rr]
        print(f'  eps = {eps}: log|X| in [{min(mods):.3f}, {max(mods):.3f}], |arg(-X)| in [{min(angs):.3f}, {max(angs):.3f}], real: {real_count(c)}/12')
    print(f'done ({time.time() - t0:.0f} s)')

if __name__ == '__main__':
    main()
