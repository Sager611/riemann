#!/usr/bin/env python3
"""Exact limit of DPC's two-line method on [0, theta]: Im C, C = e^{pi A_phi} on L^2[0,theta].
Kernel (audit D7): Im C(x,y) = sin(a) delta(x-y) + a e^{ia(x-y)} e^{ia sgn(x-y)}, a = pi*phi.
Claim: Im C is PSD iff phi*(1+theta) <= 1 (for 1/2 < phi < 1, phi*theta < 1/2); i.e. the method gives exactly lambda > 1/(1+theta).
Also checks <1,(I+aS)^{-1}1> = tan(a theta)/a, S F(x) = i int_0^theta sgn(x-y) F(y) dy."""
import numpy as np, json, sys
def lam_min(phi, theta, n=1200):
    a = np.pi*phi; h = theta/n; x = (np.arange(n)+0.5)*h
    D = x[:,None]-x[None,:]
    K = a*np.exp(1j*a*D)*np.exp(1j*a*np.sign(D))
    K[np.arange(n),np.arange(n)] = a*np.exp(0j)  # sgn(0)=0 on the diagonal (measure zero)
    M = np.sin(a)*np.eye(n) + K*h
    M = (M+M.conj().T)/2
    return float(np.linalg.eigvalsh(M)[0])/h  # eigenvalue of the operator, normalised by grid weight
def quad_form(phi, theta, n=1200):
    a = np.pi*phi; h = theta/n; x = (np.arange(n)+0.5)*h
    S = 1j*np.sign(x[:,None]-x[None,:])*h
    one = np.ones(n)
    v = np.linalg.solve(np.eye(n)+a*S, one)
    return float(np.real(np.vdot(one, v))*h), float(np.tan(a*theta)/a)
out = {}
for theta in (0.5, 0.45, 4/7, 1.0):
    phic = 1/(1+theta)
    rows = []
    for dphi in (-0.02, -0.005, 0.005, 0.02):
        phi = phic + dphi
        rows.append((phi, lam_min(phi, theta, 900)))
    qf = quad_form(phic+0.01, theta, 900) if (phic+0.01)*theta < 0.5 else None
    out[str(theta)] = dict(phi_c=phic, lam_min=rows, quad_form_vs_tan=qf)
    print("theta=%.4f  1/(1+theta)=%.5f  " % (theta, phic) + "  ".join("phi=%.4f:min eig(ImC)=%+.4e" % r for r in rows) + ("  <1,(I+aS)^-1 1>=%.6f vs tan(a th)/a=%.6f" % qf if qf else ""))
json.dump(out, open(sys.argv[1] if len(sys.argv)>1 else "threshold.json","w"), indent=1)
