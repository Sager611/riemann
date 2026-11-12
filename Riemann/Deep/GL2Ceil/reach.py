# gl2ceil reach.py (float, untrusted): fixed-dual gap vs offsets in R and lam at the optimum
import sys, time
sys.argv = ['x', '25', '256', '0.07006']
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ceil as cc
from flint import arb, fmpq
import numpy as np
Rc = fmpq(40420, 4096); Rf = float(Rc.p)/float(Rc.q); lam_i = 50500; lamf = lam_i*cc.HL
sc, mc, beta, J = cc.build(Rf, lamf)
C, B = cc.CB_arb(lam_i)
def Dval(q, C, B): return q['K0'] - q['K1']*cc.XI/C - q['K2']/(cc.XI*B)
print("Rc", Rf, "lam", lamf, "J", J)
for dk in [1, 4, 16, 64, 256]:
    for sgn in (-1, 1):
        R2 = Rc + sgn*fmpq(dk, 4096); R2f = float(R2.p)/float(R2.q)
        q = cc.rig(sc, mc, R2, R2)
        D = float(Dval(q, C, B).mid()); _, _, _, J2 = cc.build(R2f, lamf)
        print(f"dR={sgn*dk/4096:+.5f}: m={J2:.10e} D_fixed={D:.10e} relgap={(J2-D)/J2:.3e}")
q = cc.rig(sc, mc, Rc, Rc)
for dl in [64, 256, 1024, 4096, 16384]:
    for sgn in (-1, 1):
        i2 = lam_i + sgn*dl; C2, B2 = cc.CB_arb(i2)
        D = float(Dval(q, C2, B2).mid()); _, _, _, J2 = cc.build(Rf, i2*cc.HL)
        print(f"dlam={sgn*dl*cc.HL:+.5f}: m={J2:.10e} D_fixed={D:.10e} relgap={(J2-D)/J2:.3e}")
