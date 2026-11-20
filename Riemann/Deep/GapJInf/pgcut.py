"""gapjinf: TRUE per-gap class of order J by cutting planes; exact separation via pgschur (all modes k <= KM eliminated).
LP: max -a0 m0 s.t. sum_{j>=1} a_j m_j = 1, <x, w(piece)> >= 0 for all cut pieces (balanced-mixture dual)."""
import numpy as np
from scipy.optimize import linprog
import pgschur as PS
from validate53 import basis as polybasis

PI = np.pi
GG = PI * np.concatenate([np.linspace(1, 0.96, 9), np.linspace(0.95, 0.5, 19), np.linspace(0.47, 0.05, 15)])


class PerGap:
    def __init__(self, m, J, X0=6.25, ggrid=GG, KM=3000):
        self.J = J; self.X0 = X0; self.KM = KM
        self.gp = [(g, par) for g in ggrid for par in (1, 0)]
        D0 = PS.mode_data(PI, J, 1, 10); self.nv = D0['nv']; self.tkeys = D0['tkeys']
        self.scale = np.concatenate([X0 ** -np.arange(J + 1.0), [X0 ** -(mm + l + 0.5) for (mm, l) in self.tkeys]])
        self.ms = np.asarray(m[:J + 1], float) / X0 ** np.arange(J + 1)
        self.cuts = []
        for w in np.concatenate([np.linspace(1, 3, 21), np.geomspace(3, 100, 30)]):
            self.add_point(w * w)
        # second oracle (explore53's Legendre bubble basis, degree NP): small gaps, where the mode sums lose accuracy
        NP = 30
        self.PG, self.PB = polybasis(NP, J)
        self.pkeys = [(2 * mm + 1, 2 * l) for (mm, l) in self.tkeys]
        self.gsmall = PI * np.geomspace(1e-4, 0.3, 40)

    def add_point(self, x):
        v = np.zeros(self.nv); v[:self.J + 1] = x ** np.arange(self.J + 1.0)
        self.add_raw(v)

    def add_raw(self, v):
        v = v * (1 / self.scale) ** -1 if False else v / self.scale ** -1
        self.cuts.append(v / np.abs(v).max())

    def solve(self):
        A = -np.array(self.cuts)
        c = np.zeros(self.nv); c[0] = self.ms[0]
        Aeq = np.zeros((1, self.nv)); Aeq[0, 1:self.J + 1] = self.ms[1:]
        res = linprog(c, A_ub=A, b_ub=np.zeros(len(A)), A_eq=Aeq, b_eq=[1.0], bounds=[(-1e3, 1e3)] * self.nv, method='highs')
        if res.status != 0:
            raise RuntimeError('LP: %s' % res.message)
        self.res = res
        return -res.fun, res.x

    def form(self, x):
        return x * self.scale

    def separate(self, x, tol=1e-12):
        xt = self.form(x); nnew = 0; worst = 0.0
        a = xt[:self.J + 1]
        pc = np.polynomial.Polynomial(a)
        r = pc.deriv().roots(); r = r[np.abs(r.imag) < 1e-9].real
        for xx in np.concatenate([[1.0], r[r > 1]]):
            if pc(xx) < -tol * np.abs(a).max() * max(1, xx) ** self.J:
                self.add_point(xx); nnew += 1
        for g, par in self.gp:
            D = PS.mode_data(g, self.J, par, self.KM)
            Amin, R, Bt = PS.reduce(D, xt)
            Av = (xt @ D['A']) / (np.abs(xt) @ np.abs(D['A']))
            if Av.min() < -tol:
                k = int(np.argmin(Av)); xx = D['w'][k] ** 2
                self.add_point(xx); nnew += 1; worst = min(worst, Av.min()); continue
            if R is None:
                continue
            ev, U = np.linalg.eigh(R)
            scl = np.abs(np.diag(R)).max() + 1e-300
            if ev[0] < -tol * scl:
                self.add_raw(PS.piece_vector(D, U[:, 0], Bt)); nnew += 1
                worst = min(worst, ev[0] / scl)
        for g in self.gsmall:
            # basis function n has exact degree n+2: scale by g^(min(n+2, J) - 1/2) so that all entries stay O(1) as g -> 0
            sc = g ** (np.minimum(np.arange(len(self.PG[0])) + 2, self.J) - 0.5)
            S2 = np.outer(sc, sc)
            M = sum(xt[j] * g ** (1 - 2 * j) * self.PG[j] * S2 for j in range(self.J + 1))
            for i, pq in enumerate(self.pkeys):
                M = M + xt[self.J + 1 + i] * g ** (-pq[0] - pq[1]) * self.PB[pq] * S2
            Nrm = sum(g ** (1 - 2 * j) * self.PG[j] * S2 for j in range(self.J + 1))
            Lc = np.linalg.cholesky(Nrm); Li = np.linalg.inv(Lc)
            ev, U = np.linalg.eigh(Li @ M @ Li.T)
            if ev[0] < -1e-10:
                c = sc * (Li.T @ U[:, 0])
                v = np.array([g ** (1 - 2 * j) * c @ self.PG[j] @ c for j in range(self.J + 1)] +
                             [g ** (-pq[0] - pq[1]) * c @ self.PB[pq] @ c for pq in self.pkeys])
                self.add_raw(v); nnew += 1; worst = min(worst, ev[0])
        return nnew, worst

    def prune(self, x, keep=1e-6):
        """drop cuts with large slack (memory/LP size); the first 51 sinusoid cuts are kept"""
        A = np.array(self.cuts); sl = A @ x
        idx = [i for i in range(len(self.cuts)) if i < 51 or sl[i] < keep]
        self.cuts = [self.cuts[i] for i in idx]

    def run(self, rounds=200, verbose=False, prune_every=25):
        hist = []
        for it in range(rounds):
            val, x = self.solve()
            if prune_every and it and it % prune_every == 0:
                self.prune(x)
            nnew, worst = self.separate(x)
            hist.append((val, nnew, worst))
            if verbose and it % 10 == 0:
                print('   round %d: h=%.7f cuts=%d new=%d worst=%.2e' % (it, val, len(self.cuts), nnew, worst), flush=True)
            if nnew == 0:
                break
        return val, self.form(x), it, hist
