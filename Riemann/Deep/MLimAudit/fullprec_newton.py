#!/usr/bin/env python3
"""mlimaudit: fresh full-precision check of a few Weisz zeros (own code). Box B = [1/2 +- r] + i[gamma_W +- r], r = 2(rad + ulp);
one complex interval Newton step N(B) = s0 - zeta(s0)/zeta'(B) at P bits. N(B) inside B => exactly one zero in B (convex box,
integral mean value form + Brouwer), on the line by the symmetry of B about Re s = 1/2. Then Im N(B) inside Weisz's gamma ball,
and 1/(N(B) zeta'(N(B))) inside his c balls.
Usage: python3 fullprec_newton.py DATADIR n1 n2 ..."""
import sys, os, time
sys.set_int_max_str_digits(0)
from flint import arb, acb, acb_series, ctx
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from npzstream import NpzRows
NPZ = os.path.join(sys.argv[1], 'zeta_zeros_14400.npz')
want = sorted(int(x) for x in sys.argv[2:])
rows = {}
st = [NpzRows(NPZ, nm) for nm in ('gamma', 'gamma_rad', 'c_re', 'c_re_rad', 'c_im', 'c_im_rad')]
for i, r in enumerate(zip(*st)):
    if i + 1 in want: rows[i + 1] = r
def ulp(s):
    m, _, e = s.strip().lower().partition('e'); return (int(e) if e else 0) - len(m.partition('.')[2])
def wball(mid, rad):
    return arb(mid) + arb(0, (arb(rad) + arb('1e%d' % ulp(mid))).upper())
ctx.prec = 33100
for n in want:
    t0 = time.time()
    gs, gr, crs, crr, cis, cir = rows[n]
    gm = arb(gs)
    r = (2 * (arb(gr) + arb('1e%d' % ulp(gs))) + arb(gm.rad())).upper()
    B = acb(arb('0.5') + arb(0, r), gm.mid() + arb(0, r))
    s0 = acb(arb('0.5'), gm.mid())
    z0 = s0.zeta()
    zp = acb_series([B, 1], prec=2).zeta().coeffs()[1]
    NB = s0 - z0 / zp
    inside = bool(B.contains(NB))
    g_ok = bool(wball(gs, gr).contains(NB.imag))
    zq = acb_series([NB, 1], prec=2).zeta().coeffs()[1]
    c = 1 / (NB * zq)
    c_ok = bool(wball(crs, crr).contains(c.real)) and bool(wball(cis, cir).contains(c.imag))
    print('zero %d (gamma %s): N(B) in B: %s; Im N(B) in Weisz gamma ball: %s; c in Weisz c balls: %s; |zeta(s0)| ~ %s; %.0f s' %
          (n, gs[:9], inside, g_ok, c_ok, abs(z0).str(3), time.time() - t0), flush=True)
