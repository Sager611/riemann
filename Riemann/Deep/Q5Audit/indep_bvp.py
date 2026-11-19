#!/usr/bin/env python3
"""q5audit: independent kappa(theta) - 2theta/3 for 2508.11108 Question 5 from the paper's formulas (no Degree2Mol code).
BVP (mainODEbdry, TeX l.433): S'' + tanh(t) S' - S + 1/(1+e^{2t}) = 0, S(0)=1/2, S(R)=0, R = sqrt(3/5)/theta (simplipara l.372),
solved by scipy solve_bvp (linear: S, S', J with J' = e^{-t}S, J(0)=0), max 5000 nodes, RSS guard 230 MB.
c from (crRdef l.728): c = 1/2 - 1/sqrt15 + (2/sqrt15) e^R I, I = 1 - e^{-R} - S'(0) - J(R); check: I from the functional by Simpson on sol.sol.
kappa = 1 - log(c)/R evaluated as -log((2/sqrt15) I + (1/2 - 1/sqrt15) e^{-R})/R.
R > 40 (theta < 0.0194): I_R = Istar - e^{-R} + (exp. small), Istar := lim (I_R + e^{-R}) = 1 - S_inf'(0) - J_inf, taken from R = 40
after checking that I_R + e^{-R} is constant in R = 20, 30, 40.
Run: uv run --no-project --with numpy --with scipy python3 indep_bvp.py
"""
import numpy as np, resource
from scipy.integrate import solve_bvp, simpson
s15 = np.sqrt(15.0); a0 = 2 / s15; c00 = 0.5 - 1 / s15; s35 = np.sqrt(0.6); PHI = (1 + np.sqrt(5)) / 2
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e6   # macOS: bytes -> MB
def r(t): e = np.exp(-2 * t); return e / (1 + e)          # 1/(1+e^{2t}) for t >= 0
def solve(R, tol=1e-9):
    if rss() > 230: raise SystemExit("RSS guard: %.0f MB" % rss())
    t = np.linspace(0, R, int(max(200, 10 * R)))
    S0 = 0.5 * np.exp(-t / PHI) * (1 - t / R); dS0 = np.gradient(S0, t)
    y0 = np.vstack([S0, dS0, np.zeros_like(t)])
    f = lambda t, y: np.vstack([y[1], -np.tanh(t) * y[1] + y[0] - r(t), np.exp(-t) * y[0]])
    bc = lambda ya, yb: np.array([ya[0] - 0.5, yb[0], ya[2]])
    sol = solve_bvp(f, bc, t, y0, tol=tol, max_nodes=5000, bc_tol=1e-12)
    Sp0 = sol.y[1, 0]; J = sol.y[2, -1]
    I = 1 - np.exp(-R) - Sp0 - J
    tt = np.linspace(0, R, 20001); Y = sol.sol(tt); S, P = Y[0], Y[1]
    Gd = simpson(np.exp(tt) * (S ** 2 + P ** 2) + np.exp(-tt) * ((1 - S) ** 2 + P ** 2), x=tt)
    return sol, Sp0, J, I, Gd
def kappa_from_I(R, I): return -np.log(a0 * I + c00 * np.exp(-R)) / R
# large-R constant
ist = {}
for R in [20.0, 30.0, 40.0]:
    sol, Sp0, J, I, Gd = solve(R); ist[R] = I + np.exp(-R)
    print(f"R={R}: status {sol.status} nodes {sol.x.size}: S'(0)={Sp0:.13f} J={J:.13f} I+e^-R={I+np.exp(-R):.13f} rss={rss():.0f}MB", flush=True)
Istar = ist[40.0]
print(f"Istar = {Istar:.13f};  spread over R=20,30,40: {max(ist.values())-min(ist.values()):.2e};  paper l.~768: lim(1-S'(0)-int) <= 1.1421")
thetas = [2/3, 0.5, 0.49999, 0.45, 0.4, 0.38, 0.37997, 0.3799, 0.35, 0.3, 0.25, 0.2285, 0.22845, 0.2284, 0.2, 1/6, 0.15, 1/8, 0.1,
          5/54, 0.05, 0.03, 0.02, 0.01, 0.005, 0.002, 0.001, 1e-4]
print("theta        R          S'(0)          I(formula)      I(direct)       kappa           kappa-2th/3     rel.margin   nodes st")
rows = []
for th in thetas:
    R = s35 / th
    if R <= 40:
        sol, Sp0, J, I, Gd = solve(R); info = f"{sol.x.size:<6d} {sol.status}"
    else:
        Sp0 = float('nan'); I = Istar - np.exp(-R); Gd = float('nan'); info = "asympt"
    k = kappa_from_I(R, I); m = k - 2 * th / 3
    rows.append((th, m, m / (2 * th / 3)))
    print(f"{th:<12.6g} {R:<10.5f} {Sp0:<14.10f} {I:<15.12f} {Gd:<15.12f} {k:<15.11f} {m:<15.4e} {m/(2*th/3):<12.4e} {info} rss={rss():.0f}MB", flush=True)
for th in [0.5, 0.38, 0.2285]:
    R = s35 / th; _, a, J1, I1, _ = solve(R, 1e-7); _, b2, J2, I2, _ = solve(R, 1e-9)
    print(f"tol check theta={th}: S'(0) tol1e-7 {a:.12f} tol1e-9 {b2:.12f}; kappa diff {abs(kappa_from_I(R,I1)-kappa_from_I(R,I2)):.2e}")
inr = [x for x in rows if x[0] <= 0.5]
mn = min(inr, key=lambda x: x[2]); print(f"grid (0,1/2]: all margins > 0: {all(x[1] > 0 for x in inr)}; smallest relative margin theta={mn[0]:.6g} rel={mn[2]:.4e}")
