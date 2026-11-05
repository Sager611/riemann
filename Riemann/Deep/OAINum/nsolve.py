
# nsolve.py -- faithful implementation of Theorem norm:solve (Section 4 of the factorization paper)
# over k = F_p with zeta in F_p.  Promised norm equation N_{L/k(X)}(h) = G, L = k(X)[Y]/(Y^q - F).
from ffpoly import *
from unitroot import unitroot
import numpy as np

# ---------------- rational functions ----------------
class Frac:
    __slots__ = ('n', 'd', 'p')
    def __init__(self, n, d, p, reduce=True):
        self.p = p
        n = norm(list(n)); d = norm(list(d))
        if not d:
            raise ZeroDivisionError
        if not n:
            self.n, self.d = [], [1]
            return
        if reduce:
            g = pgcd(n, d, p)
            if g != [1]:
                n = pdiv_exact(n, g, p); d = pdiv_exact(d, g, p)
        inv = pow(d[-1], p - 2, p)
        self.n = pscal(n, inv, p); self.d = pscal(d, inv, p)
    def __add__(s, o):
        return Frac(padd(pmul(s.n, o.d, s.p), pmul(o.n, s.d, s.p), s.p), pmul(s.d, o.d, s.p), s.p)
    def __sub__(s, o):
        return Frac(psub(pmul(s.n, o.d, s.p), pmul(o.n, s.d, s.p), s.p), pmul(s.d, o.d, s.p), s.p)
    def __mul__(s, o):
        return Frac(pmul(s.n, o.n, s.p), pmul(s.d, o.d, s.p), s.p)
    def inv(s):
        return Frac(s.d, s.n, s.p)
    def __truediv__(s, o):
        return s * o.inv()
    def iszero(s):
        return not s.n
    def __eq__(s, o):
        return s.n == o.n and s.d == o.d

def fdet(Mx, p):
    n = len(Mx)
    A = [row[:] for row in Mx]
    det = Frac([1], [1], p)
    for c in range(n):
        piv = None
        for r in range(c, n):
            if not A[r][c].iszero():
                piv = r; break
        if piv is None:
            return Frac([], [1], p)
        if piv != c:
            A[c], A[piv] = A[piv], A[c]
            det = det * Frac([p - 1], [1], p)
        det = det * A[c][c]
        inv = A[c][c].inv()
        for r in range(c + 1, n):
            if A[r][c].iszero():
                continue
            f = A[r][c] * inv
            A[r] = [A[r][j] - f * A[c][j] for j in range(n)]
    return det

def fsolve(Aug, nunk, p):
    # Aug: list of rows (lists of Frac) length nunk+1 ; return solution list or None
    A = [row[:] for row in Aug]
    rows = len(A)
    piv_cols = []
    r = 0
    for c in range(nunk):
        piv = None
        for i in range(r, rows):
            if not A[i][c].iszero():
                piv = i; break
        if piv is None:
            continue
        A[r], A[piv] = A[piv], A[r]
        inv = A[r][c].inv()
        A[r] = [x * inv for x in A[r]]
        for i in range(rows):
            if i != r and not A[i][c].iszero():
                f = A[i][c]
                A[i] = [A[i][j] - f * A[r][j] for j in range(nunk + 1)]
        piv_cols.append(c)
        r += 1
    for i in range(r, rows):
        if not A[i][nunk].iszero():
            return None
    sol = [Frac([], [1], p) for _ in range(nunk)]
    for i, c in enumerate(piv_cols):
        sol[c] = A[i][nunk]
    return sol

# ---------------- power series ----------------
def ps_trunc(a, n):
    return norm(list(a[:n]))
def ps_mul(a, b, n, p):
    return ps_trunc(pmul(a, b, p), n)
def ps_inv(a, n, p):
    assert a and a[0] % p
    b = [pow(a[0], p - 2, p)]
    k = 1
    while k < n:
        k = min(2 * k, n)
        # b <- b (2 - a b)
        ab = ps_mul(a, b, k, p)
        two_minus = psub([2], ab, p)
        b = ps_mul(b, two_minus, k, p)
    return ps_trunc(b, n)
def ps_root(f, q, n, p):
    # S^q = f, S(0) = 1, f(0) = 1
    assert f[0] % p == 1
    S = [1]
    k = 1
    qinv = pow(q, p - 2, p)
    while k < n:
        k = min(2 * k, n)
        Sq1 = [1]
        for _ in range(q - 1):
            Sq1 = ps_mul(Sq1, S, k, p)
        Sq = ps_mul(Sq1, S, k, p)
        num = psub(Sq, ps_trunc(f, k), p)
        den = pscal(Sq1, q, p)
        S = psub(S, ps_mul(num, ps_inv(den, k, p), k, p), p)
    return ps_trunc(S, n)

def hensel_root(r0, u, q, m, p):
    # R^q = u mod m (m = H^prec), R = r0 mod H ; Newton iteration
    R = pmod(r0, m, p)
    for _ in range(64):
        Rq1 = ppowmod(R, q - 1, m, p)
        Rq = pmod(pmul(Rq1, R, p), m, p)
        diff = psub(Rq, pmod(u, m, p), p)
        if not diff:
            return R
        R = pmod(psub(R, pmul(diff, pinvmod(pscal(Rq1, q, p), m, p), p), p), m, p)
    raise RuntimeError("hensel failed")

# ---------------- matrices of polys mod m (or series) ----------------
def mat_zero(q):
    return [[[] for _ in range(q)] for _ in range(q)]
def mat_id(q):
    Mx = mat_zero(q)
    for i in range(q):
        Mx[i][i] = [1]
    return Mx
def mat_mul_mod(A, B, red):
    q = len(A)
    C = mat_zero(q)
    for i in range(q):
        for k in range(q):
            if not A[i][k]:
                continue
            for j in range(q):
                if B[k][j]:
                    C[i][j] = padd(C[i][j], pmul(A[i][k], B[k][j], P_), P_)
    for i in range(q):
        for j in range(q):
            C[i][j] = red(C[i][j])
    return C
P_ = None

def basis_mats(Ys, Vs, q, red):
    # Ys = (scale, Mat) for Y, Vs = (scale, Mat) for V. returns dict (a,b) -> (scale, Mat) of Y^a V^b
    out = {}
    Ypow = [(0, mat_id(q))]
    for a in range(1, q):
        s, Mx = Ypow[-1]
        Ypow.append((s + Ys[0], mat_mul_mod(Mx, Ys[1], red)))
    Vpow = [(0, mat_id(q))]
    for b in range(1, q):
        s, Mx = Vpow[-1]
        Vpow.append((s + Vs[0], mat_mul_mod(Mx, Vs[1], red)))
    for a in range(q):
        for b in range(q):
            out[(a, b)] = (Ypow[a][0] + Vpow[b][0], mat_mul_mod(Ypow[a][1], Vpow[b][1], red))
    return out

def coeffs_mod(poly, n):
    v = list(poly) + [0] * (n - len(poly))
    return v[:n]

class NormSolver:
    def __init__(self, p, q, zeta, omega, F, A, B, M=None, verbose=False):
        global P_
        P_ = p
        self.p, self.q, self.zeta, self.omega = p, q, zeta, omega
        self.F = monic(F, p)
        self.N = deg(self.F)
        assert self.N % q == 0
        g = pgcd(A, B, p)
        A = pdiv_exact(A, g, p); B = pdiv_exact(B, g, p)
        lcB = B[-1]
        self.A = pscal(A, pow(lcB, p - 2, p), p)
        self.B = monic(B, p)
        self.D = max(deg(self.A), deg(self.B))
        self.M = M if M is not None else self.N + self.D + 2
        self.verbose = verbose
        self.info = {}

    def groups(self):
        p = self.p
        F, A, B = self.F, self.A, self.B
        gr = []
        used = [1]
        for (Aj, j) in squarefree_decomp(A, p) if deg(A) > 0 else []:
            g1 = pgcd(Aj, F, p)
            if deg(g1) > 0: gr.append((g1, 1, j))
            g2 = pdiv_exact(Aj, g1, p)
            if deg(g2) > 0: gr.append((g2, 0, j))
            used = pmul(used, Aj, p)
        for (Bj, j) in squarefree_decomp(B, p) if deg(B) > 0 else []:
            g1 = pgcd(Bj, F, p)
            if deg(g1) > 0: gr.append((g1, 1, -j))
            g2 = pdiv_exact(Bj, g1, p)
            if deg(g2) > 0: gr.append((g2, 0, -j))
            used = pmul(used, Bj, p)
        Fr = pdiv_exact(F, pgcd(F, used, p), p)
        if deg(Fr) > 0: gr.append((Fr, 1, 0))
        return gr

    def local_model(self, Hi, eps, w, prec):
        p, q, zeta = self.p, self.q, self.zeta
        m = ppow(Hi, prec, p)
        red = lambda f: pmod(f, m, p)
        F, A, B = self.F, self.A, self.B
        if eps == 1:
            Fp_ = pdiv_exact(F, Hi, p)
            if w >= 0:
                Ap = pdiv_exact(A, ppow(Hi, w, p), p)
                u = red(pmul(Ap, pinvmod(red(pmul(B, ppow(Fp_, w, p), p)), m, p), p))
            else:
                Bp = pdiv_exact(B, ppow(Hi, -w, p), p)
                u = red(pmul(pmul(A, ppow(Fp_, -w, p), p), pinvmod(red(Bp), m, p), p))
            r0 = unitroot(pmod(u, Hi, p), Hi, q, zeta, self.omega, p)
            R = hensel_root(r0, u, q, m, p)
            Ym = mat_zero(q)
            for j in range(q - 1):
                Ym[j + 1][j] = [1]
            Ym[0][q - 1] = red(F)
            Yi = mat_zero(q)
            for j in range(q - 1):
                Yi[j][j + 1] = red(Hi)
            Yi[q - 1][0] = pinvmod(red(Fp_), m, p)
            Dm = mat_zero(q)
            for j in range(q):
                Dm[j][j] = red(pscal(R, pow(zeta, j, p), p))
            if w >= 0:
                V = Dm
                for _ in range(w):
                    V = mat_mul_mod(Ym, V, red)
                Vs = (0, V)
            else:
                V = Dm
                for _ in range(-w):
                    V = mat_mul_mod(Yi, V, red)
                Vs = (-w, V)
            return basis_mats((0, Ym), Vs, q, red)
        else:
            assert w % q != 0
            S0 = unitroot(pmod(F, Hi, p), Hi, q, zeta, self.omega, p)
            S = hensel_root(S0, F, q, m, p)
            zinv = pow(zeta, p - 2, p)
            Ym = mat_zero(q)
            for j in range(q):
                Ym[j][j] = red(pscal(S, pow(zinv, j, p), p))
            V = mat_zero(q)
            if w > 0:
                for j in range(q - 1):
                    V[j + 1][j] = [1]
                V[0][q - 1] = red(pmul(A, pinvmod(red(B), m, p), p))
                Vs = (0, V)
            else:
                Bp = pdiv_exact(B, ppow(Hi, -w, p), p)
                Hw = red(ppow(Hi, -w, p))
                for j in range(q - 1):
                    V[j + 1][j] = Hw
                V[0][q - 1] = red(pmul(A, pinvmod(red(Bp), m, p), p))
                Vs = (-w, V)
            return basis_mats((0, Ym), Vs, q, red)

    def infinity_model(self, prec):
        p, q, zeta = self.p, self.q, self.zeta
        N = self.N
        red = lambda f: ps_trunc(f, prec)
        Finf = [self.F[N - k] for k in range(N + 1)]
        Sinf = ps_root(Finf, q, prec, p)
        zinv = pow(zeta, p - 2, p)
        Ym = mat_zero(q)
        for j in range(q):
            Ym[j][j] = red(pscal(Sinf, pow(zinv, j, p), p))
        Ys = (N // q, Ym)
        A, B = self.A, self.B
        dA, dB = deg(A), deg(B)
        At = [A[dA - k] for k in range(dA + 1)]
        Bt = [B[dB - k] for k in range(dB + 1)]
        winf = dB - dA
        Gt = ps_mul(At, ps_inv(Bt, prec, p), prec, p)
        V = mat_zero(q)
        if winf >= 0:
            for j in range(q - 1):
                V[j + 1][j] = [1]
            V[0][q - 1] = red([0] * winf + Gt)
            Vs = (0, V)
        else:
            for j in range(q - 1):
                V[j + 1][j] = [0] * (-winf) + [1]
            V[0][q - 1] = Gt
            Vs = (-winf, V)
        return basis_mats(Ys, Vs, q, red), winf

    def sections(self):
        p, q = self.p, self.q
        M = self.M
        gr = self.groups()
        H = [1]
        for (Hi, e, w) in gr:
            H = pmul(H, Hi, p)
        self.H = H
        dH = deg(H)
        Dmax = M * (dH + 1)
        nper = Dmax + 1
        nunk = q * q * nper
        idx = lambda a, b, n: ((a * q) + b) * nper + n
        rows = []
        self.info['groups'] = [(deg(Hi), e, w) for (Hi, e, w) in gr]
        for (Hi, eps, w) in gr:
            di = deg(Hi)
            if eps == 0 and w % q == 0:
                for a in range(q):
                    for b in range(q):
                        need = M - (b * w) // q
                        if need <= 0:
                            continue
                        mod = ppow(Hi, need, p)
                        L = need * di
                        block = np.zeros((L, nunk), dtype=np.int64)
                        cur = [1]
                        for n in range(nper):
                            r = pmod(cur, mod, p)
                            block[:, idx(a, b, n)] = coeffs_mod(r, L)
                            cur = pmul(cur, [0, 1], p) if len(cur) < 2 * L else pmod(pmul(cur, [0, 1], p), mod, p)
                        rows.append(block)
                continue
            # matrix model: need max scale first
            if eps == 1:
                Smax = (q - 1) * max(0, -w)
            else:
                Smax = (q - 1) * max(0, -w)
            prec = M + Smax
            mats = self.local_model(Hi, eps, w, prec)
            mod = ppow(Hi, prec, p)
            L = prec * di
            for r_ in range(q):
                for c_ in range(q):
                    block = np.zeros((L, nunk), dtype=np.int64)
                    for a in range(q):
                        for b in range(q):
                            s_ab, Mx = mats[(a, b)]
                            ent = Mx[r_][c_]
                            if not ent:
                                continue
                            mult = pmod(pmul(ppow(Hi, Smax - s_ab, p), ent, p), mod, p)
                            cur = mult
                            for n in range(nper):
                                block[:, idx(a, b, n)] = coeffs_mod(cur, L)
                                cur = pmod(pmul(cur, [0, 1], p), mod, p)
                    rows.append(block)
        # infinity
        # scales: s_ab = a*N/q + b*max(0,-winf); compute winf first
        winf = deg(self.B) - deg(self.A)
        Sinf = (q - 1) * (self.N // q) + (q - 1) * max(0, -winf)
        prec = M + Sinf + 1
        mats, _ = self.infinity_model(prec)
        Ht = [H[dH - k] for k in range(dH + 1)]
        HtM = ps_inv(ppow_trunc(Ht, M, prec, p), prec, p)
        evrows = {}
        for r_ in range(q):
            for c_ in range(q):
                block = np.zeros((M + Sinf, nunk), dtype=np.int64)
                evrow = np.zeros(nunk, dtype=np.int64)
                for a in range(q):
                    for b in range(q):
                        s_ab, Mx = mats[(a, b)]
                        ent = Mx[r_][c_]
                        if not ent:
                            continue
                        T = ps_mul(HtM, [0] * (Sinf - s_ab) + ent, prec, p)
                        T = coeffs_mod(T, prec)
                        for n in range(nper):
                            sh = Dmax - n
                            if sh >= prec:
                                continue
                            col = idx(a, b, n)
                            # coefficient of z^m is T[m - sh]
                            for mm in range(sh, M + Sinf):
                                block[mm, col] = (block[mm, col] + T[mm - sh]) % p
                            mm = M + Sinf
                            if mm - sh >= 0:
                                evrow[col] = (evrow[col] + T[mm - sh]) % p
                rows.append(block)
                evrows[(r_, c_)] = evrow
        Amat = np.concatenate(rows, axis=0) if rows else np.zeros((0, nunk), dtype=np.int64)
        del rows
        self.info['unknowns'] = nunk
        self.info['equations'] = Amat.shape[0]
        ker = nullspace(Amat, nunk, p)
        self.info['dimS'] = len(ker)
        self.ker = ker
        self.nper, self.Dmax = nper, Dmax
        evs = []
        for v in ker:
            Em = np.zeros((q, q), dtype=np.int64)
            for (r_, c_), er in evrows.items():
                Em[r_, c_] = int((er @ v) % p) if True else 0
            evs.append(Em % p)
        self.evs = evs
        return ker, evs

    def element_from_vec(self, v):
        q, p, nper = self.q, self.p, self.nper
        P = [[norm([int(x) for x in v[((a * q) + b) * nper:((a * q) + b + 1) * nper]]) for b in range(q)] for a in range(q)]
        return (ppow(self.H, self.M, p), P)

    def run(self):
        p, q = self.p, self.q
        ker, evs = self.sections()
        if not ker:
            raise RuntimeError('no sections')
        # B = span of evs
        flat = np.array([E.flatten() for E in evs], dtype=np.int64) % p
        R, piv = rref(flat, p)
        Bbasis = [R[i].reshape(q, q) for i in range(R.shape[0])]
        self.info['dimB'] = len(Bbasis)
        nb = len(Bbasis)
        Gram = np.zeros((nb, nb), dtype=np.int64)
        for i in range(nb):
            for j in range(nb):
                Gram[i, j] = int(np.trace((Bbasis[i] @ Bbasis[j]) % p)) % p
        nk = nullspace(Gram, nb, p)
        Nmats = []
        for c in nk:
            Mx = np.zeros((q, q), dtype=np.int64)
            for i in range(nb):
                Mx = (Mx + int(c[i]) * Bbasis[i]) % p
            Nmats.append(Mx)
        self.info['dimN'] = len(Nmats)
        if Nmats:
            stack = np.concatenate(Nmats, axis=0)
            U = nullspace(stack, q, p)
        else:
            U = [np.eye(q, dtype=np.int64)[i] for i in range(q)]
        self.info['dimU'] = len(U)
        vv = U[0] % p
        k0 = int(np.nonzero(vv)[0][0])
        wrow = np.zeros(q, dtype=np.int64)
        wrow[k0] = pow(int(vv[k0]), p - 2, p)
        a0 = np.outer(vv, wrow) % p
        # check a0 in B
        test = np.concatenate([flat, a0.flatten()[None, :]], axis=0)
        _, piv2 = rref(test, p)
        self.info['a0_in_B'] = (len(piv2) == len(piv))
        # lift a0 to S
        Lsys = np.array([E.flatten() for E in evs], dtype=np.int64).T % p
        x = solve_lin(Lsys, a0.flatten(), p)
        if x is None:
            raise RuntimeError('lift failed')
        vec = np.zeros(len(ker[0]), dtype=np.int64)
        for xi, kv in zip(x, ker):
            if int(xi):
                vec = (vec + int(xi) * kv) % p
        b = self.element_from_vec(vec)
        e = b
        for _ in range(q - 1):
            e = self.emul(e, b)
        self.e = e
        # checks: idempotent and reduced trace 1
        ee = self.emul(e, e)
        self.info['idempotent'] = self.eeq(ee, e)
        c00 = Frac(e[1][0][0], e[0], p)
        self.info['trd_e'] = (c00 * Frac([q], [1], p))
        h = self.extract_h(e)
        return h

    # ----- cyclic algebra arithmetic, elements (den, num[a][b]) -----
    def emul(self, x, y):
        p, q, zeta = self.p, self.q, self.zeta
        F, A, B = self.F, self.A, self.B
        dx, X = x; dy, Yy = y
        out = [[[] for _ in range(q)] for _ in range(q)]
        for a in range(q):
            for b in range(q):
                if not X[a][b]:
                    continue
                for c in range(q):
                    for d in range(q):
                        if not Yy[c][d]:
                            continue
                        t = pmul(X[a][b], Yy[c][d], p)
                        t = pscal(t, pow(zeta, b * c, p), p)
                        aa, bb = a + c, b + d
                        if aa >= q:
                            t = pmul(t, F, p); aa -= q
                        if bb >= q:
                            t = pmul(t, A, p); bb -= q
                        else:
                            t = pmul(t, B, p)
                        out[aa][bb] = padd(out[aa][bb], t, p)
        den = pmul(pmul(dx, dy, p), B, p)
        return self.ereduce((den, out))

    def ereduce(self, x):
        p = self.p
        den, X = x
        g = den
        for row in X:
            for c in row:
                if c:
                    g = pgcd(g, c, p)
                    if g == [1]:
                        break
        if deg(g) > 0:
            den = pdiv_exact(den, g, p)
            X = [[pdiv_exact(c, g, p) if c else [] for c in row] for row in X]
        inv = pow(den[-1], p - 2, p)
        return (pscal(den, inv, p), [[pscal(c, inv, p) for c in row] for row in X])

    def eeq(self, x, y):
        p = self.p
        for a in range(self.q):
            for b in range(self.q):
                if psub(pmul(x[1][a][b], y[0], p), pmul(y[1][a][b], x[0], p), p):
                    return False
        return True

    def extract_h(self, e):
        p, q, zeta = self.p, self.q, self.zeta
        F, A, B = self.F, self.A, self.B
        den, Ecoef = e
        Ef = [[Frac(Ecoef[a][b], den, p) for b in range(q)] for a in range(q)]
        G = Frac(A, B, p)
        Fq = Frac(F, [1], p)
        # Y^a e
        Ya = []
        for a in range(q):
            Mx = [[Frac([], [1], p) for _ in range(q)] for _ in range(q)]
            for c in range(q):
                for d in range(q):
                    if Ef[c][d].iszero():
                        continue
                    t = Ef[c][d]
                    cc = a + c
                    if cc >= q:
                        t = t * Fq; cc -= q
                    Mx[cc][d] = Mx[cc][d] + t
            Ya.append(Mx)
        Ve = [[Frac([], [1], p) for _ in range(q)] for _ in range(q)]
        for c in range(q):
            for d in range(q):
                if Ef[c][d].iszero():
                    continue
                t = Ef[c][d] * Frac([pow(zeta, c, p)], [1], p)
                dd = d + 1
                if dd >= q:
                    t = t * G; dd -= q
                Ve[c][dd] = Ve[c][dd] + t
        Aug = []
        for c in range(q):
            for d in range(q):
                Aug.append([Ya[a][c][d] for a in range(q)] + [Ve[c][d]])
        sol = fsolve(Aug, q, p)
        if sol is None:
            raise RuntimeError('Ve = he inconsistent')
        return sol

    def norm_of(self, h):
        p, q = self.p, self.q
        Fq = Frac(self.F, [1], p)
        Mx = [[Frac([], [1], p) for _ in range(q)] for _ in range(q)]
        for j in range(q):
            for a in range(q):
                if h[a].iszero():
                    continue
                c = a + j
                t = h[a]
                if c >= q:
                    t = t * Fq; c -= q
                Mx[c][j] = Mx[c][j] + t
        return fdet(Mx, p)

def ppow_trunc(a, e, n, p):
    r = [1]
    for _ in range(e):
        r = ps_mul(r, a, n, p)
    return r

