#!/usr/bin/env python3
"""gl2ceil float ceilings (untrusted, for the report): kappa*(xi) = max_{R,lam} 1 - log M(R,lam)/R with M from the exact Euler-Lagrange
ODE (ceil.build), and the Ritz table: max over (R, lam) with Q = 1 + sum_{n<=N} h_n (T_{2n-1}(1-2y) - 1) (degree 2N-1), P = P_lam.
usage: uv run --no-project --with numpy --with scipy --with python-flint python3 floatceil.py  > logs/floatceil.log"""
import sys, os, numpy as np
from scipy.optimize import minimize
from numpy.polynomial import chebyshev as Ch
from numpy.polynomial import legendre as Lg
HERE = os.path.dirname(os.path.abspath(__file__))
sys.argv = ['x', '25', '256', '0.07006']; sys.path.insert(0, HERE)
import ceil as cc
xg, wg = Lg.leggauss(600); yg = (xg + 1)/2; wg = wg/2
def ritz(R, lam, xi, N):      # Degree2Mol/opt2.py basis (w = e^{Ry}Q, beta free)
    z = 1 - 2*yg; E = np.exp(R*yg); Phi = [E]; dPhi = [R*E]
    for n in range(1, N+1):
        c = np.zeros(2*n); c[2*n-1] = 1.0; T = Ch.chebval(z, c); dT = Ch.chebval(z, Ch.chebder(c))*(-2.0)
        Phi.append(E*(T - 1)); dPhi.append(R*E*(T - 1) + E*dT)
    Phi = np.array(Phi); dPhi = np.array(dPhi); end = np.array([np.exp(R)] + [-2*np.exp(R)]*N)
    C, B = cc.PAB_f(lam); H = (C/xi)*((Phi*wg) @ Phi.T) + xi*B*((dPhi*wg) @ dPhi.T) + 0.5*np.outer(end, end)
    h = np.linalg.solve(H[1:, 1:], -H[1:, 0]); v = np.concatenate([[1.0], h])
    return 1 - np.log(0.5 + v @ H @ v)/R
def kap_ode(R, lam, xi):
    cc.xi = xi; sc, mc, beta, J = cc.build(R, lam); return 1 - np.log(J)/R, beta
def argmax(f, x0, tol):
    r = minimize(lambda v: -f(v[0], abs(v[1])), x0, method='Nelder-Mead', options=dict(xatol=tol, fatol=1e-13, maxiter=500))
    return -r.fun, r.x[0], abs(r.x[1])
for name, xi, x0, table in (('GL2, nu < 25/128: xi = 25/256', 25/256, [9.87, 0.77], True), ('GL2, Ramanujan: xi = 1/8', 1/8, [7.69, 0.77], True),
                            ('GL2 at the record length xi = 25/256 - 1e-7', 25/256 - 1e-7, [9.87, 0.77], False),
                            ('zeta, Conrey theta = 4/7', 4/7, [1.29, 1.0], False)):
    k, R, lam = argmax(lambda R, l: kap_ode(R, l, xi)[0], x0, 1e-7); beta = kap_ode(R, lam, xi)[1]
    print(f"{name}: kappa* = {k:.10f} at R = {R:.6f}, lam = {lam:.6f}, Q(1) = beta - 1 = {beta-1:.3e}", flush=True)
    if table:
        for N in (1, 2, 3, 4, 6, 8, 10, 14, 20, 30):
            kN, RN, lN = argmax(lambda R, l: ritz(R, l, xi, N), [R, lam], 1e-6)
            print(f"   Ritz deg Q = {2*N-1:2d}: max kappa = {kN:.10f} at R = {RN:.5f}, lam = {lN:.5f}", flush=True)
