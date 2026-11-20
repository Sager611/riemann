"""gapjinf: dictionary bound by cutting planes (LP). Each cut is the moment vector of ONE genuine function
F(x) = sum_n d_n e^{i(th + 2 pi n/P) x} vanishing on a periodic zero pattern Z (gaps <= pi), d = N u (u = bottom
eigenvector of N^* diag(P(xi^2)) N), or a sinusoid sin(w x), w >= 1 (Lemma C of REPORT.md).
h_dict(J) = max { -a0 m0 : sum_{j>=1} a_j m_j = 1, sum_j a_j mu_t(j) >= 0 for all dictionary functions t }  >= h_J(true).
"""
import numpy as np
from scipy.optimize import linprog

PI = np.pi


def pattern_fibres(gaps, nth, Xi, tol=1e-9):
    gaps = np.asarray(gaps, float)
    z = np.concatenate([[0.0], np.cumsum(gaps)[:-1]]); P = gaps.sum(); base = 2 * PI / P
    out = []
    for t in range(nth):
        th = (t + 0.5) / nth * base
        n = np.arange(int(np.floor((-Xi - th) / base)) - 1, int(np.ceil((Xi - th) / base)) + 2)
        xi = th + base * n
        k = np.abs(xi) <= Xi
        n, xi = n[k], xi[k]
        E = np.exp(2j * PI * np.outer(z, n) / P)
        u, s, vh = np.linalg.svd(E)
        rk = int((s > tol * max(1, s.max())).sum())
        N = vh[rk:].conj().T
        if N.shape[1] == 0:
            continue
        out.append((xi, N))
    return out


class DictLP:
    def __init__(self, m, J, fibres, X0=6.25, xmax=1e4):
        self.J = J; self.X0 = X0
        self.ms = np.asarray(m[:J + 1], float) / X0 ** np.arange(J + 1)
        self.fib = [(xi ** 2 / X0, N) for xi, N in fibres]
        self.cuts = []
        for w in np.concatenate([np.linspace(1, 3, 41), np.geomspace(3, np.sqrt(xmax), 40)]):
            self.add_point(w * w / X0)
        self.c1 = 1.0 / X0; self.xmax = xmax / X0

    def add_point(self, x):
        v = x ** np.arange(self.J + 1.0)
        self.cuts.append(v / np.abs(v).max())

    def add_vec(self, x, w):
        v = np.array([np.dot(w, x ** j) for j in range(self.J + 1)])
        self.cuts.append(v / np.abs(v).max())

    def solve(self):
        J = self.J
        A = -np.array(self.cuts)
        cobj = np.zeros(J + 1); cobj[0] = self.ms[0]                     # minimise a0 m0  (= maximise -a0 m0)
        res = linprog(cobj, A_ub=A, b_ub=np.zeros(len(A)), A_eq=self.ms[None, :] * (np.arange(J + 1) >= 1),
                      b_eq=[1.0], bounds=[(-1e4, 1e4)] * (J + 1), method='highs')
        if res.status != 0:
            raise RuntimeError('LP status %d: %s' % (res.status, res.message))
        return -res.fun, res.x

    def separate(self, a, tol=1e-10):
        nnew = 0
        # positivity on [1, xmax]
        pc = np.polynomial.Polynomial(a)
        r = pc.deriv().roots(); r = r[np.abs(r.imag) < 1e-9].real
        xs = np.concatenate([[self.c1], r[(r > self.c1) & (r < self.xmax)]])
        vals = pc(xs)
        for x, val in zip(xs, vals):
            if val < -tol:
                self.add_point(x); nnew += 1
        worst = 0.0
        for x, N in self.fib:
            Pv = pc(x)
            M = N.conj().T @ (Pv[:, None] * N)
            ev, U = np.linalg.eigh(M)
            if ev[0] < -tol * np.abs(Pv).max():
                d = N @ U[:, 0]
                self.add_vec(x, np.abs(d) ** 2); nnew += 1
                worst = min(worst, ev[0] / np.abs(Pv).max())
        return nnew, worst

    def run(self, rounds=60, verbose=False):
        for it in range(rounds):
            val, a = self.solve()
            nnew, worst = self.separate(a)
            if verbose:
                print('   round %d: h=%.6f cuts=%d new=%d worst=%.2e' % (it, val, len(self.cuts), nnew, worst), flush=True)
            if nnew == 0:
                break
        return val, a * self.X0 ** -np.arange(self.J + 1.0), it
