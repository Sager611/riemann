# Where do the kernel constraints bind, and what would the exact K_h (instead of Lemma 10's lower bound) give? (model only)
import numpy as np
from scipy.optimize import linprog
from dbn import scalars_f, CG, lp_best, Bf, pf
B = Bf
def cons(A2, exact):
    U = np.concatenate([np.linspace(0, 0.999, 500), 1 - np.logspace(-3, -9, 40)])
    V = np.linspace(0, 1, 201)
    U, V = np.meshgrid(U, V, indexing='ij'); U = U.ravel(); V = V.ravel()
    D = B * U / (1 - U); r = A2 * V
    w = D + B - r; Dl = w * w + 4 * D * r
    s = (1 - U) ** 2  # scaling (N/Dl^2 ~ 1/D at large D)
    a = s * 2 * w * Dl / Dl ** 2
    c = s * (8 * B * w * w - (2 * w + 4 * B) * Dl) / Dl ** 2
    if exact:
        b = 1 / 25
        f = lambda v: np.where(v * v + 4 * D * r > 0, 2 * v / np.maximum(v * v + 4 * D * r, 1e-300), np.inf)
        lhs = np.minimum(f(D + np.maximum(r, b) - r), f(D + A2 - r))
    else:
        lhs = 2 / (D + A2 + 3 * r)
    rr = s * lhs; keep = np.isfinite(rr)
    return a[keep], c[keep], rr[keep], U[keep], V[keep]
for (tm, tp, q) in [(0.0, 0.00025, 0.5576563720), (0.06, 0.06025, 0.235), (0.118, 0.11825, 0.0405)]:
    sc = scalars_f(tm, tp); C, G, _ = CG(sc)
    al = C / pf - 2 / (B - q); ga = G + 2 * (B + q) / (B - q) ** 2
    out = []
    for ex in (False, True):
        a, c, rhs, U, V = cons(q, ex)
        res = linprog([-al, ga], A_ub=np.stack([a, c], 1), b_ub=rhs, bounds=[(0, None), (0, None)], method='highs')
        lam, mu = res.x; k = 4 * (al * lam - ga * mu)
        slack = rhs - a * lam - c * mu
        idx = np.argsort(slack)[:3]
        out.append("exact=%s lam=%.5f mu=%.5f k=%.5f  active at (U,V)=%s" % (ex, lam, mu, k, [(round(U[j], 3), round(V[j], 3)) for j in idx]))
    print("q=%.4f t=%.3f" % (q, tm)); [print("  " + o) for o in out]
lam, mu, k = lp_best(0.5576563720, *CG(scalars_f(0.0, 0.00025))[:2], 16, 1); print("Bernstein-LP 16x1 at q0:", lam, mu, k)
