#!/usr/bin/env python3
"""gl2ceil cross-check against independent code: the exact record values of Degree2Mol/cert.py (read-only import) must be >= our Arb
dual lower bound at the same (xi, R, lam), and close to it (the record is near-optimal). Not part of the proof; a test of the dual code.
usage: uv run --no-project --with numpy --with scipy --with python-flint --with mpmath python3 xcheck.py"""
import sys, json, os
HERE = os.path.dirname(os.path.abspath(__file__)); D2 = os.path.join(HERE, '..', 'Degree2Mol')
prm = json.load(open(os.path.join(D2, 'params.json')))
key = sys.argv[1] if len(sys.argv) > 1 else 'AT21_uncond'
xn, xd = prm[key]['xi'].split('/')
sys.argv = ['x', xn, xd, '0.07006']; sys.path.insert(0, HERE); sys.path.insert(0, D2)
import ceil as cc
import cert
from flint import fmpq, arb
from fractions import Fraction
F = lambda s: (lambda fr: fmpq(fr.numerator, fr.denominator))(Fraction(s))
R = F(prm[key]['R']); lam = F(prm[key]['lam']); h = [F(s) for s in prm[key]['h']]; xi = F(prm[key]['xi'])
Q = cert.Q_cheb(h); P = cert.P_sinh(lam)
c_rec, k_rec = cert.c_kappa(P, Q, R, xi)          # exact c of the record polynomials (cert.py, 4000 bits)
cc.ctx.prec = 256
Rf = float(R.p)/float(R.q); lf = float(lam.p)/float(lam.q)
sc, mc, beta, J = cc.build(Rf, lf)
q = cc.rig(sc, mc, R, R)
l = arb(lam); s = l.sinh(); s2 = (2*l).sinh()/(4*l); C = l*l*(arb(1)/2 + s2)/(s*s); B = (s2 - arb(1)/2)/(s*s)
D = q['K0'] - q['K1']*arb(xi)/C - q['K2']/(arb(xi)*B)
print(f"{key}: xi={xi} R={float(Rf):.8f} lam={lf:.8f}")
print(f"  record c (cert.py exact)  = {c_rec.str(18)}  kappa = {k_rec.str(14)}")
print(f"  dual lower bound D (Arb)  = {D.str(18)}  kappa_D = {(1 - D.log()/arb(R)).str(14)}  float primal J = {J:.12e}")
print(f"  c_rec - D >= 0: {bool(c_rec - D >= 0) or bool(c_rec > D)};  (c_rec - D)/c_rec = {float(((c_rec - D)/c_rec).mid()):.3e}")
