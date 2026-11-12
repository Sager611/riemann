"""Single-window adversaries with doubles (escape-route test for Q-side refinements).
On-line S-periodic configuration: sites x_j with masses m_j in {1,2}, total mass S (density 1). Window data:
HS/length = Khat(0) + 2 S^-2 sum_{1<=m<S} Khat(m/S)|sum_j m_j e^{2 pi i m x_j/S}|^2 must not exceed R_W = Khat(0) + 2 int_0^1 Khat a da
(then move points to reach equality). Any certificate using only window W's trace/HS data (and any configuration-wise
lemmas) has kappa <= simple fraction of such a configuration (or of block mixtures: LP over the found patterns)."""
import sys, json, itertools
import numpy as np
from scipy.optimize import minimize, linprog
from tenwin import TEN_C, MT_C, Khat, H_of, autoc, int_psi
from periodic import setup

def necklaces_marks(S):
    seen = set()
    for seq in itertools.product((0, 1, 2), repeat=S):
        if sum(seq) != S:
            continue
        r = min(seq[i:] + seq[:i] for i in range(S))
        if r in seen:
            continue
        seen.add(r)
        yield r

def phi_target(c):
    return H_target(c)

def H_target(c):
    # 2 int_0^1 Khat(a) a da = R_W - Khat(0)
    xg, wg = np.polynomial.legendre.leggauss(200)
    a = 0.5 * xg + 0.5; w = 0.5 * wg
    return 2 * np.sum(w * a * Khat(a, c))

def relaxPhi(c, S, xs, ms_, cache={}):
    key = (id(c), S)
    if key not in cache:
        cache[key] = setup(c, S)
    ms, kh, k0 = cache[key]
    mw = np.array(ms_, float)
    def f(x):
        ph = np.exp(2j * np.pi * np.outer(ms, x) / S)
        mu = ph @ mw
        val = 2 * np.sum(kh * np.abs(mu) ** 2) / S ** 2
        g = (4 * np.real(np.conj(mu)[:, None] * (2j * np.pi * ms[:, None] / S) * ph * mw[None, :]) * kh[:, None]).sum(axis=0) / S ** 2
        return val, g
    r = minimize(f, xs, jac=True, method="L-BFGS-B", options={"maxiter": 3000, "gtol": 1e-14})
    return min(r.fun, f(xs)[0])

if __name__ == "__main__":
    which = sys.argv[1]; Smax = int(sys.argv[2])
    c = TEN_C if which == "TEN" else (MT_C if which == "MT" else np.array(json.load(open(which))["c"]))
    tgt = H_target(c); H = H_of(c)
    pts = []
    for S in range(1, Smax + 1):
        for nk in necklaces_marks(S):
            xs = np.array([i for i in range(S) if nk[i] > 0], float)
            mw = [nk[i] for i in range(S) if nk[i] > 0]
            nu = sum(1 for v in mw if v == 1) / S
            if nu > 0.8:
                continue
            ph = relaxPhi(c, S, xs, mw)
            pts.append((nu, ph, S, ''.join('.12'[v] for v in nk)))
    # LP: minimise sum l_i nu_i s.t. sum l_i phi_i <= tgt, sum l_i = 1, l >= 0
    nu = np.array([p[0] for p in pts]); ph = np.array([p[1] for p in pts])
    res = linprog(nu, A_ub=[ph], b_ub=[tgt], A_eq=[np.ones(len(pts))], b_eq=[1], bounds=(0, None), method="highs")
    used = [(pts[i][3], round(res.x[i], 4), round(pts[i][0], 4), pts[i][1]) for i in range(len(pts)) if res.x[i] > 1e-9]
    print(which, "H_W =", H, "target 2int Khat a =", tgt, "patterns", len(pts))
    print("single-window double-adversary ceiling nu_double =", res.fun)
    print("mixture:", used)
