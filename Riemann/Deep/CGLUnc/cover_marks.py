"""One-window ceiling over the whole window family (escape-route cap for Q-side refinements).
For on-line S-periodic patterns X_i with masses {1,2} (density 1, simple fraction nu_i) put
Phi_i = (2/S^2) sum_{1<=m<S} |sum_j m_j e^{2 pi i m x_j/S}|^2 Tsym(m/S)   and   Psi = 2 int_0^1 Tsym(a) a da,
so that HS_W(X) - Khat(0) = tr(Phi G)/(int psi)^2 and R_W - Khat(0) = tr(Psi G)/(int psi)^2 for the coefficient Gram G.
If lambda (simplex) has Psi - sum lambda_i Phi_i PSD, the block mixture Xbar has HS_W <= R_W for EVERY window of the
family (single, mixed, complex); after a window-wise move of points to HS_W = R_W it is an on-line configuration with
N0s/N = sum lambda_i nu_i matching that window's trace and HS data, so no certificate using one window's data proves more."""
import sys, json, itertools
import numpy as np
from scipy.optimize import minimize
import cvxpy as cp
from winopt import Ttensor
from tenwin import TEN_C, MT_C
from marks import necklaces_marks

Smax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
xg, wg = np.polynomial.legendre.leggauss(80)
A01 = 0.5 * xg + 0.5; W01 = 0.5 * wg
Psi = sum(w * a * Ttensor(a) for a, w in zip(A01, W01)) * 2
Psi = 0.5 * (Psi + Psi.T)
TS = {}
def Tset(S):
    if S not in TS:
        T = np.array([Ttensor(m / S) for m in range(1, S)]) if S > 1 else np.zeros((0, 17, 17))
        TS[S] = 0.5 * (T + np.transpose(T, (0, 2, 1)))
    return TS[S]

def relax(c, S, xs, mw):
    T = Tset(S); ms = np.arange(1, S)
    kh = np.einsum('i,mij,j->m', c, T, c) if S > 1 else np.zeros(0)
    def f(x):
        ph = np.exp(2j * np.pi * np.outer(ms, x) / S); mu = ph @ mw
        val = 2 * np.sum(kh * np.abs(mu) ** 2) / S ** 2
        g = (4 * np.real(np.conj(mu)[:, None] * (2j * np.pi * ms[:, None] / S) * ph * mw[None, :]) * kh[:, None]).sum(axis=0) / S ** 2
        return val, g
    if S == 1:
        return xs
    r = minimize(f, xs, jac=True, method="L-BFGS-B", options={"maxiter": 2000, "gtol": 1e-14})
    return r.x if r.fun <= f(xs)[0] else xs

def Phimat(S, x, mw):
    if S == 1:
        return np.zeros((17, 17))
    T = Tset(S); ms = np.arange(1, S)
    mu = np.exp(2j * np.pi * np.outer(ms, x) / S) @ mw
    return 2 * np.einsum('m,mij->ij', np.abs(mu) ** 2, T) / S ** 2

def margin_check(lib, Ph, nus, nu_cap):
    lam = cp.Variable(len(lib), nonneg=True); t = cp.Variable()
    M = Psi - sum(lam[i] * Ph[i] for i in range(len(lib)))
    cons = [cp.sum(lam) == 1, 0.5 * (M + M.T) - t * np.eye(17) >> 0]
    if nu_cap is not None:
        cons.append(nus @ lam <= nu_cap)
    prob = cp.Problem(cp.Maximize(t), cons)
    prob.solve(solver=cp.SCS, eps=1e-9, max_iters=200000)
    return t.value, lam.value

if __name__ == "__main__":
    cs = [TEN_C, MT_C, np.array(json.load(open("winopt_best.json"))["c"])]
    lib = []
    for S in range(1, Smax + 1):
        for nk in necklaces_marks(S):
            xs0 = np.array([i for i in range(S) if nk[i] > 0], float)
            mw = np.array([nk[i] for i in range(S) if nk[i] > 0], float)
            nu = float(np.sum(mw == 1)) / S
            if nu > 0.85:
                continue
            for c in cs:
                x = relax(c, S, xs0, mw)
                lib.append((S, nu, x, mw, ''.join('.12'[v] for v in nk)))
    print("library", len(lib), flush=True)
    Ph = [Phimat(S, x, mw) for (S, nu, x, mw, _) in lib]
    nus = np.array([q[1] for q in lib])
    t0, lv = margin_check(lib, Ph, nus, None)
    print("max PSD margin without nu cap: t* =", t0, " nu_bar at that lambda =", float(nus @ lv), flush=True)
    for cap in [0.70, 0.69, 0.685, 0.6818, 0.68, 0.679]:
        tc, lc = margin_check(lib, Ph, nus, cap)
        print(f"  cap nu_bar <= {cap}: t* = {tc}", flush=True)
    class P: pass
    prob = P(); prob.value = float(nus @ lv)
    used = sorted([(lv[i], lib[i][4], lib[i][1]) for i in range(len(lib)) if lv[i] > 1e-6], reverse=True)[:8]
    print("top patterns:", [(round(a, 4), b, round(c_, 4)) for a, b, c_ in used])
    json.dump({"nu_bar": prob.value, "lam": lv.tolist(), "lib": [(q[0], q[1], q[2].tolist(), q[3].tolist(), q[4]) for q in lib]},
              open(f"cover_marks_{Smax}.json", "w"))
