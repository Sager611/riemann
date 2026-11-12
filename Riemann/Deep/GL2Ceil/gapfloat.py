# gl2ceil gapfloat.py (float, untrusted): ODE-based optimum S (2508.11108 BVP, DERIVATION s4) -> piecewise-polynomial Fenchel dual, float64 check
import numpy as np
from scipy.integrate import solve_ivp
from numpy.polynomial import chebyshev as Ch
from numpy.polynomial import legendre as Lg
def PAB(lam):
    s = np.sinh(lam)
    return lam**2*(0.5 + np.sinh(2*lam)/(4*lam))/s**2, (np.sinh(2*lam)/(4*lam) - 0.5)/s**2
def solve_S(R, lam, xi, tq):
    Cl, Bl = PAB(lam)
    c0 = Cl/xi - xi*Bl*R*R; c1 = 4*xi*Bl*R*R; c = -c0/c1
    def f(t, Y):  # Y = [up, up', phi, phi']
        th = np.tanh(t); g = 1/(1+np.exp(2*t))
        return [Y[1], c*g - th*Y[1] - c*Y[0], Y[3], -th*Y[3] - c*Y[2]]
    ts = np.unique(np.concatenate([[0.0, R], tq]))
    sol = solve_ivp(f, (0, R), [0.5, 0.0, 0.0, 1.0], method='DOP853', rtol=1e-13, atol=1e-14, t_eval=ts)
    Y = sol.y; up, dup, ph, dph = Y
    a = (1 - up[-1])/ph[-1]; u = up + a*ph; du = dup + a*dph
    v = -ph/ph[-1]; dv = -dph/ph[-1]
    return ts, u, du, v, dv, c, Cl, Bl
def J_of_beta(R, lam, xi, NG=400):
    # J(beta) = (1+e^{2R}(beta-1)^2)/2 + (C/xi) W0 + xi B W1, S = beta u + v, quadratic in beta
    xg, wg = Lg.leggauss(NG); tq = (xg+1)/2*R; wq = wg/2*R
    ts, u, du, v, dv, c, Cl, Bl = solve_S(R, lam, xi, tq)
    idx = np.searchsorted(ts, tq); u, du, v, dv = u[idx], du[idx], v[idx], dv[idx]
    et = np.exp(tq); emt = np.exp(-tq); pre = np.exp(R)/(2*R)
    def J(beta):
        S = beta*u + v; dS = beta*du + dv
        W0 = pre*np.sum(wq*(et*S**2 + emt*(beta - S)**2))
        W1 = pre*R*R*np.sum(wq*(et*(S + 2*dS)**2 + emt*(beta - S + 2*dS)**2))
        return 0.5*(1 + np.exp(2*R)*(beta-1)**2) + (Cl/xi)*W0 + xi*Bl*W1
    j0, j1, j2 = J(0.0), J(1.0), J(2.0)
    A2 = (j2 - 2*j1 + j0)/2; A1 = (j1 - j0) - A2
    beta = -A1/(2*A2)
    return beta, J(beta), (Cl, Bl, c)
def dual_pieces(R, lam, xi, beta, K=24, deg=14):
    Cl, Bl = PAB(lam); al = Cl/xi; ga = xi*Bl
    # nodes in y: pieces on [0,1] for sigma, on [1/2,1] for m
    yk = np.linspace(0, 1, K+1); ym = np.linspace(0.5, 1, K//2+1)
    xc = (np.cos(np.pi*(np.arange(deg+1)+0.5)/(deg+1)) + 1)/2   # Chebyshev pts on [0,1]
    ys_s = np.concatenate([yk[i] + (yk[i+1]-yk[i])*xc for i in range(K)])
    ys_m = np.concatenate([ym[i] + (ym[i+1]-ym[i])*xc for i in range(len(ym)-1)])
    tt = np.abs(R*(2*np.concatenate([ys_s, ys_m]) - 1))
    ts, u, du, v, dv, c, _, _ = solve_S(R, lam, xi, tt)
    S = beta*u + v; dS = beta*du + dv
    d2S = c*beta/(1+np.exp(2*ts)) - np.tanh(ts)*dS - c*S
    def Qd(y):
        t = R*(2*y-1); i = np.searchsorted(ts, np.abs(t))
        up = t >= 0
        Q = np.where(up, S[i], beta - S[i]); Q1 = 2*R*dS[i]; Q2 = np.where(up, 4*R*R*d2S[i], -4*R*R*d2S[i])
        return Q, Q1, Q2
    Q, Q1, Q2 = Qd(ys_s); sig = ga*np.exp(R*ys_s)*(R*Q + Q1)
    Q, Q1, Q2 = Qd(ys_m); L = al*Q - ga*(R*R*Q + 2*R*Q1 + Q2); mm = np.exp(2*R*ys_m)*L
    # fit local polys in x in [0,1]
    def fit(vals, nodes):
        out = []
        for i in range(len(nodes)-1):
            cc = Ch.chebfit(2*xc-1, vals[i*(deg+1):(i+1)*(deg+1)], deg)
            # convert to poly in x on [0,1]: z = 2x-1
            p = Ch.cheb2poly(cc); # poly in z
            # substitute z = 2x - 1
            P = np.polynomial.Polynomial(p)(np.polynomial.Polynomial([-1, 2]))
            out.append(P.coef)
        return out
    sp = fit(sig, yk); mp_ = fit(mm, ym)
    # enforce continuity of sigma (float version: adjust constant terms)
    for i in range(1, K): sp[i][0] += np.polyval(sp[i-1][::-1], 1.0) - sp[i][0]
    return yk, sp, ym, mp_, al, ga
def D_float(R, yk, sp, ym, mp_, al, ga, NG=40):
    xg, wg = Lg.leggauss(NG); x = (xg+1)/2; w = wg/2
    P = np.polynomial.Polynomial
    def sig(y):  # evaluate sigma and sigma' piecewise
        i = np.minimum(np.searchsorted(yk, y, side='right')-1, len(sp)-1); h = yk[1]-yk[0]
        xx = (y - yk[i])/h
        return np.array([P(sp[j])(xx[k]) for k, j in enumerate(i)]), np.array([P(sp[j]).deriv()(xx[k])/h for k, j in enumerate(i)])
    def mu(y):
        z = np.where(y >= 0.5, y, 1-y); i = np.minimum(np.searchsorted(ym, z, side='right')-1, len(mp_)-1); h = ym[1]-ym[0]
        xx = (z - ym[i])/h
        return np.array([P(mp_[j])(xx[k]) for k, j in enumerate(i)])
    # quadrature over [0,1] piecewise (pieces of sigma grid; mu grid is a refinement-compatible? use fine uniform pieces)
    nodes = np.linspace(0, 1, 2*len(sp)+1)
    yy = np.concatenate([nodes[i] + (nodes[i+1]-nodes[i])*x for i in range(len(nodes)-1)])
    ww = np.concatenate([(nodes[i+1]-nodes[i])*w for i in range(len(nodes)-1)])
    s, ds = sig(yy); m = mu(yy)
    p = ds + m*np.exp(-R*yy)
    M = np.sum(ww*m)
    s0 = P(sp[0])(0.0); s1 = P(sp[-1])(1.0)
    return 0.5 + M - 2*s0 - (2*s1 + M*np.exp(-R))**2/2 - np.sum(ww*p*p)/al - np.sum(ww*s*s)/ga
xi = 25/256
for (R, lam) in [(9.8678, 0.77059), (6.0, 0.77), (14.0, 1.2), (20.0, 2.0), (2.0, 0.3), (9.8678, 0.6)]:
    beta, J, (Cl, Bl, c) = J_of_beta(R, lam, xi)
    for K, deg in ((16, 10), (24, 14), (40, 16)):
        yk, sp, ym, mp_, al, ga = dual_pieces(R, lam, xi, beta, K, deg)
        D = D_float(R, yk, sp, ym, mp_, al, ga)
        print(f"R={R} lam={lam} K={K} deg={deg}: beta={beta:.8f} J={J:.12e} D={D:.12e} gap/J={(J-D)/J:.2e} kappa_J={1-np.log(J)/R:.10f} kappa_D={1-np.log(D)/R:.10f}", flush=True)
