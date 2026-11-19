#!/usr/bin/env python3
"""zcmp.py (ppraudit): compare own Lemma Z coefficients (mpmath quad) with ExplProg closed forms at the tail point,
and own sup_{L>=L1} Z/L with EP's Z_over_L_majorant; quad accuracy check against an exact antiderivative."""
import sys, importlib.util
sys.argv = [sys.argv[0]]; sys.dont_write_bytecode = True
import indep_kappa as IK
M = IK.M
spec = importlib.util.spec_from_file_location("ep", "../ExplProg/cert.py"); EP = importlib.util.module_from_spec(spec); spec.loader.exec_module(EP)
from flint import arb
from fractions import Fraction as Fr
# exact test: int_5^inf rho0(t)/t^2 dt
f = lambda t: IK.rho0(t) / t ** 2
ex = lambda t: -(M.log(t / (2 * M.pi)) + 1) / (M.pi * t) + IK.CB * (-1 / (2 * t) + M.log((t + 2) / t) / 4)
print("quad test: quad = %s, exact = %s, diff = %s" % (M.nstr(M.quad(f, [5, M.inf]), 20), M.nstr(-ex(M.mpf(5)), 20), M.nstr(M.quad(f, [5, M.inf]) + ex(M.mpf(5)), 5)))
P = IK.PR[:10001]; L1 = M.fsum(M.log(q) for q in P)
beta = M.mpf(34) / 10
X = beta * M.log(M.mpf(391) / 200 * L1) * IK.IB; Y = X - IK.ell; eY = M.e * Y
for name, (N0, K1, K2) in (("right", (IK.n0, IK.n1 + IK.n0 / eY, IK.n2 + 2 * IK.n1 / eY + (2 / eY) ** 2 * IK.n0)),
                           ("left", (IK.m0, IK.m1 + IK.m0 / eY, IK.m2 + 2 * IK.m1 / eY + (2 / eY) ** 2 * IK.m0))):
    mine = IK.zco(N0, K1, K2)
    tA = lambda v: arb(M.nstr(v, 60))
    ep = EP.zsum_coeffs(tA(N0), tA(K1), tA(K2))
    print("%s: own (za, zb, ze) = %s; EP = %s" % (name, [M.nstr(v, 12) for v in mine], [e.str(12) for e in ep]))
    L1a = arb(M.nstr(L1, 40))
    print("   sup Z/L own = %s; EP majorant = %s" % (M.nstr(IK.supZL(mine, N0, L1), 14), EP.Z_over_L_majorant(ep, tA(N0), L1a).str(14)))
