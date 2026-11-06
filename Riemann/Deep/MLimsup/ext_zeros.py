#!/usr/bin/env python3
"""Zeros 14401..N(T1max) of zeta with Arb at low precision: gamma balls and rigorous upper bounds for |c| = 1/|rho zeta'(rho)|.
Their phases at Weisz's points are unknown (gamma known only to ~25 digits), so certify.py bounds each term by 2|k(gamma)||c|.
Usage: python3 ext_zeros.py T1max out.npz"""
import sys, time
import numpy as np
from flint import arb, acb, acb_series, ctx
T1max = int(sys.argv[1]); out = sys.argv[2]
t0 = time.time()
ctx.prec = 96
NT = int(arb(T1max).zeta_nzeros().unique_fmpz())
n0 = 14401
cnt = NT - n0 + 1
gam_mid, gam_rad, cup = [], [], []
B = 500
for s in range(n0, NT + 1, B):
    zs = acb.zeta_zeros(s, min(B, NT + 1 - s))
    for z in zs:
        zp = acb_series([z, 1], prec=2).zeta().coeffs()[1]
        cu = (1 / (abs(z) * abs(zp))).upper()
        gam_mid.append(float(z.imag.mid())); gam_rad.append(float(z.imag.rad())); cup.append(float(cu) * (1 + 1e-12))
    print('  up to zero %d: %.0f s' % (s + len(zs) - 1, time.time() - t0), flush=True)
np.savez(out, gam=np.array(gam_mid), gam_rad=np.array(gam_rad), cup=np.array(cup), NT=NT, T1max=T1max)
print('N(%d) = %d; %d zeros; sum 2|c| = %.5f; %.0f s' % (T1max, NT, cnt, 2 * sum(cup), time.time() - t0))
