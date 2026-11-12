#!/usr/bin/env python3
"""gl2ceil implementation test: for sample cells, the Taylor lower bound lower(+-d) used by check_lam must lie below the dual value D
evaluated directly at R = Rc +- d (same sigma, mu; Arb), and the gap must shrink like d^3 (tests K0', K1', K1'' and the remainders).
Also checks D(sigma, mu) against an mpmath quadrature of the dual integrals at the cell centre (tests the moment formulas).
usage: uv run --no-project --with numpy --with scipy --with python-flint --with mpmath python3 taylorcheck.py"""
import sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = ['x', '25', '256', '0.0700564']; sys.path.insert(0, HERE)
import ceil as cc
from flint import fmpq, arb
import mpmath as mp
mp.mp.dps = 40
def Dval(q, x, z): return q['K0'] - x*q['K1'] - z*q['K2']
for Rc, i in ((fmpq(40421, 4096), 50500), (fmpq(20*4096, 4096), 70000), (fmpq(6*4096 + 77, 4096), 30000)):
    Rf = float(Rc.p)/float(Rc.q); sc, mc, beta, J = cc.build(Rf, i*cc.HL)
    C, B = cc.CB_arb(i); x = cc.XI/C; z = 1/(cc.XI*B)
    for dk in (1, 4, 16):
        d = fmpq(dk, 4096); q = cc.rig(sc, mc, Rc, Rc - d); da = arb(d)
        a2 = -x*q['K1pp']/2; a2l = a2.lower(); a2c = a2l if a2l < 0 else arb(0)
        out = []
        for sgn in (-1, 1):
            low = (q['K0'] + sgn*q['K0p']*da - x*(q['K1'] + sgn*q['K1p']*da) - z*q['K2'] + a2c*da*da - da*da/2*q['F0'] - x*q['G1']*da*da*da/6)
            R2 = Rc + sgn*d; q2 = cc.rig(sc, mc, R2, R2); Dd = Dval(q2, x, z)
            out.append((bool(Dd > low), float(((Dd - low)/Dd).mid())))
        print(f"Rc={Rf:.5f} lam={i*cc.HL:.4f} d={float(d):.2e}: direct D >= Taylor lower at -d, +d: {out[0][0]}, {out[1][0]}; rel gaps {out[0][1]:.2e}, {out[1][1]:.2e}")
    # mpmath quadrature of the dual at Rc (independent of the moment code)
    q = cc.rig(sc, mc, Rc, Rc); Dm = Dval(q, x, z)
    K = cc.KP
    def sig(y, der=0):
        k = min(int(y*K), K-1); xx = mp.mpf(y)*K - k; c = [mp.mpf(float(t)) for t in sc[k]]
        if k > 0:  # continuity fix as in rig: constant term = previous end value
            c[0] = sum(sig_end[k-1])
        if der == 0: return mp.polyval(c[::-1], xx)
        return K*mp.polyval([j*c[j] for j in range(len(c)-1, 0, -1)], xx)
    sig_end = []
    for k in range(K):
        c = [mp.mpf(float(t)) for t in sc[k]]
        if k > 0: c[0] = sum(sig_end[k-1])
        sig_end.append(c)
    def mu(y):
        z_ = y if y >= 0.5 else 1 - y; j = min(int((z_ - 0.5)*K), K//2 - 1); xx = (mp.mpf(z_) - mp.mpf(1)/2)*K - j
        return mp.polyval([mp.mpf(float(t)) for t in mc[j]][::-1], xx)
    Rm = mp.mpf(Rf); xm = mp.mpf(float(x.mid())); zm = mp.mpf(float(z.mid()))
    pts = [mp.mpf(k)/K for k in range(K+1)]
    I1 = mp.quad(lambda y: (sig(y, 1) + mu(y)*mp.e**(-Rm*y))**2, pts)
    I2 = mp.quad(lambda y: sig(y)**2, pts); M = mp.quad(mu, pts)
    s0 = sig(mp.mpf(0)); s1 = sig(mp.mpf(1) - mp.mpf(10)**-30)
    Dq = mp.mpf(1)/2 + M - 2*s0 - (2*s1 + M*mp.e**(-Rm))**2/2 - xm*I1 - zm*I2
    print(f"   quadrature D = {mp.nstr(Dq, 15)}  vs Arb D = {Dm.str(15)}  rel diff {float((Dq - mp.mpf(float(Dm.mid())))/Dq):.1e}")
