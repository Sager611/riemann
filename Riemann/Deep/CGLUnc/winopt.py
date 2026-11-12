"""Float search for the TEN-class ceiling over the 16-harmonic window family:
maximise nu*(c) = least root of nu = H(c) + ebar_c(nu), ebar_c = lower convex hull of relaxed periodic patterns."""
import sys, json, itertools
import numpy as np
from scipy.optimize import minimize
from tenwin import TEN_C, MT_C, OM, psi, quad, H_of
from necklace import necklaces, envelope

TH0 = 1 / np.sqrt(2)
I1SQ = 2 * np.sin(TH0) ** 2
WJ = np.array([0.0] + [(0.5 - 1 / (4 * np.pi ** 2 * j ** 2)) / I1SQ for j in range(1, 17)])
H_MT = 1.5 - 1 / np.tan(TH0) / np.sqrt(2)
_xg, _wg = np.polynomial.legendre.leggauss(200)


def Hc(c):
    return H_MT - np.sum(WJ * c ** 2)


def Ttensor(alpha):
    a = abs(alpha)
    if a >= 1:
        return np.zeros((17, 17))
    lo, hi = -0.5, 0.5 - a
    s = 0.5 * (hi - lo) * _xg + 0.5 * (hi + lo)
    w = 0.5 * (hi - lo) * _wg
    A = np.cos(np.outer(s, OM))
    B = np.cos(np.outer(s + a, OM))
    return (A * w[:, None]).T @ B


def lib_patterns(Smax=16, numin=0.6, numax=0.8):
    pats = []
    for S in range(1, Smax + 1):
        for nk in necklaces(S):
            p = sum(nk)
            if numin <= p / S <= numax:
                pats.append("".join("x" if b else "." for b in nk))
    return pats


class Evaluator:
    def __init__(self, pats):
        self.pats = pats
        Ss = sorted(set(len(p) for p in pats))
        self.T = {}
        for S in Ss:
            ms = np.arange(1, S)
            self.T[S] = (ms, np.array([Ttensor(m / S) for m in ms]))
        self.T0 = Ttensor(0.0)
        self.x = {p: np.array([i for i, ch in enumerate(p) if ch == 'x'], float) for p in pats}

    def energies(self, c, relax=True):
        k0 = c @ self.T0 @ c / I1SQ
        out = []
        for pat in self.pats:
            S = len(pat)
            ms, TT = self.T[S]
            kh = np.einsum('i,mij,j->m', c, TT, c) / I1SQ
            x0 = self.x[pat]
            p = len(x0)

            def f(x):
                ph = np.exp(2j * np.pi * np.outer(ms, x) / S)
                mu = ph.sum(axis=1)
                e = (k0 * p * p + 2 * np.sum(kh * np.abs(mu) ** 2)) / S ** 2 - p / S
                g = (4 * np.real(np.conj(mu)[:, None] * (2j * np.pi * ms[:, None] / S) * ph) * kh[:, None]).sum(axis=0) / S ** 2
                return e, g
            if relax:
                r = minimize(f, x0, jac=True, method="L-BFGS-B", options={"maxiter": 500, "gtol": 1e-13})
                e, xx = r.fun, r.x
                self.x[pat] = xx
            else:
                e = f(x0)[0]
            out.append((p / S, e))
        return out

    def root(self, c):
        H = Hc(c)
        pts = self.energies(c)
        hull = envelope(pts)
        for (a, b), (a2, b2) in zip(hull, hull[1:]):
            fa, fb = a - H - b, a2 - H - b2
            if fa < 0 <= fb:
                return a + (-fa / (fb - fa)) * (a2 - a)
        return np.nan


if __name__ == "__main__":
    pats = lib_patterns(int(sys.argv[1]) if len(sys.argv) > 1 else 12)
    print("library size", len(pats), flush=True)
    ev = Evaluator(pats)
    for name, c in [("MT", MT_C.copy()), ("TEN", TEN_C.copy())]:
        print(name, "H =", Hc(c), "nu* =", ev.root(c), flush=True)
    best = [ev.root(TEN_C), TEN_C.copy()]
    it = [0]

    def obj(z):
        c = np.concatenate([[1.0], z])
        v = ev.root(c)
        it[0] += 1
        if v > best[0]:
            best[0], best[1] = v, c.copy()
        if it[0] % 100 == 0:
            print("it", it[0], "best", best[0], flush=True)
        return -v if np.isfinite(v) else 1.0
    for start in [TEN_C[1:], np.zeros(16)]:
        r = minimize(obj, start, method="Powell", options={"maxfev": int(sys.argv[2]) if len(sys.argv) > 2 else 1500, "xtol": 1e-6, "ftol": 1e-10})
        print("start done:", -r.fun, flush=True)
    print("BEST nu* =", best[0], "H =", Hc(best[1]))
    print("c =", json.dumps([round(v, 9) for v in best[1].tolist()]))
    json.dump({"nu": best[0], "c": best[1].tolist()}, open("winopt_best.json", "w"))
