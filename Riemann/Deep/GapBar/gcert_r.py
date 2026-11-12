"""gapceil: Arb certificate that BM's criterion with general r (here r = 2) fails at a given kappa for EVERY amplifier
(polynomial P of any degree, any v, any theta < 1/4).  Extends gcert.py (r = 1); see DERIVATION.md §7.

For r >= 2 the forms carry the weight X^b, b = r^2-1:  Q_j = sum_ij Gam^(r)_ij <I^i phi, X^(b+p) I^j phi>, phi = I^(2r) f,
so the natural space is L2(X^b dX) on [0,1/4]. Basis: orthonormal Jacobi polynomials q_n of L2(X^b dX) (Cholesky of the
weighted Gram matrix of the Legendre P~_n). Key facts (DERIVATION.md §7):
  * I maps W_m = span{q_n : n >= m} into W_{m-1} and is tridiagonal there:
    I q_n = th (a_n q_{n-1} + b_n q_n + c_n q_{n+1}),  |a_n| <= 1/(2n+b-1), |b_n| <= b/(2n+b)^2, |c_n| <= 1/(2n+b+1),
    so |I w| <= tau_m |w| on W_m with tau_m = th (1/(2m+b-1) + b/(2m+b)^2 + 1/(2m+b+1)).
  * B p = X^b g_p with g_p a polynomial, so <w, B p> = <w, g_p>_b and |P_W g_p|^2 = sum_{k>=N} <q_k, B p>^2 (exact window).
  * High modes: <w, B w> >= (delta - C)|w|_b^2, delta = min d(X)/X^b, C from
      (0,1)+(1,0), p = 0: 2|Gam_01| t_1;  p >= 1: (b+p) th^(p-1) |Gam_01| t_1^2 (integration by parts; X^(b+p) (Iw)^2 -> 0 at 0)
      (0,j)+(j,0), j >= 2: 2|Gam_0j| th^p t_j;   (i,j >= 1): |Gam_ij| th^p t_i t_j;   t_k = prod_{l<k} tau_{N-l}.
    Series tail: |Gam^(r)_ij| <= 4(i+j+8)^2((i+2r)(j+2r))^(2r) gam^(i+j) (certr.py), |I^k| <= th^k/k! in L2(X^b).
Run: uv run --no-project --with python-flint --with numpy --with scipy python3 runcert_r.py <kappa> [N] [r]"""
import time, math
from flint import arb, arb_mat, fmpq
import gcert as GC
import certr as CR


class OpsR:
    def __init__(self, kap, N, D, gam, r, Nw=None):
        t0 = time.time()
        self.kap, self.N, self.D, self.gam, self.r = kap, N, D, gam, r
        self.b = r * r - 1
        th = GC.th
        self.c = GC.A_(kap) * arb.pi()
        self.Gam = CR.gamma_arrays_r(self.c, D, gam, r)
        self.t_gam = time.time() - t0
        self.Nw = Nw or (N + 2 * D + 4)
        n2 = self.Nw + D + 6
        big = n2 + self.b + 6
        I, T, G0 = GC.legendre_ops(big)
        G0d = arb_mat([[G0[i] if i == j else arb(0) for j in range(big)] for i in range(big)])
        Tp = arb_mat([[arb(1) if i == j else arb(0) for j in range(big)] for i in range(big)])
        self.Wq = {}
        for q in range(self.b + 3):
            if q >= self.b:
                self.Wq[q - self.b] = GC.sub(G0d * Tp, n2, n2)
            Tp = Tp * T
        Is = GC.sub(I, n2, n2)
        self.IsT = Is.transpose()
        Y = arb_mat([[arb(1) if i == j else arb(0) for j in range(N)] for i in range(n2)])
        self.Y = [Y]
        for k in range(1, D):
            self.Y.append(Is * self.Y[-1])
        self.n2 = n2
        # orthonormal basis of L2(X^b dX) on [0, th]: q_k = sqrt((2k+b+1)/th^(b+1)) P_k^(0,b)(2X/th - 1), built in exact rationals
        # from the Jacobi recurrence (DLMF 18.9.2, alpha = 0) in the Legendre coefficient space (x = 2X/th - 1 acts as (2/th) T - Id)
        n = self.Nw; b = self.b
        thq = GC.TH
        Tq = [[fmpq(0)] * (n + 2) for _ in range(n + 2)]
        for k in range(n + 2):
            Tq[k][k] = thq / 2
            if k + 1 < n + 2: Tq[k + 1][k] = thq / 2 * fmpq(k + 1, 2 * k + 1)
            if k >= 1: Tq[k - 1][k] = thq / 2 * fmpq(k, 2 * k + 1)
        def xmul(vec):
            out = [fmpq(0)] * (n + 2)
            for k, cf in enumerate(vec):
                if cf == 0: continue
                for kk in (k - 1, k, k + 1):
                    if 0 <= kk < n + 2:
                        out[kk] += Tq[kk][k] * cf * 2 / thq
            return [o - c for o, c in zip(out, vec)]
        R = [[fmpq(1)] + [fmpq(0)] * (n + 1)]
        xv = xmul(R[0])
        R.append([((b + 2) * xv[k] - b * R[0][k]) / 2 for k in range(n + 2)])
        for m in range(1, n - 1):
            xv = xmul(R[m])
            A1 = 2 * (m + 1) * (m + b + 1) * (2 * m + b)
            R.append([(((2 * m + b + 1) * ((2 * m + b + 2) * (2 * m + b) * xv[k] - b * b * R[m][k])) - 2 * m * (m + b) * (2 * m + b + 2) * R[m - 1][k]) / A1 for k in range(n + 2)])
        Li = []
        for k in range(n):
            nk = (arb(2 * k + b + 1) / GC.A_(thq) ** (b + 1)).sqrt()
            Li.append([GC.A_(R[k][l]) * nk for l in range(n)])
        self.QT = arb_mat(Li)              # rows = coefficients of q_k in the P~ basis (q_k = sum_l Li[k][l] P~_l)
        self.QTN = GC.sub(self.QT, N, N)
        self.cache = {}
        self.t_leg = time.time() - t0

    def A(self, h, p):
        """Nw x N block of <q_k, G_p^(h) q_l> (unweighted inner product; G_p^(h) = sum Gam_ij (I^i)^T X^(b+p) I^j)."""
        key = (h, p)
        if key in self.cache: return self.cache[key]
        G = self.Gam[h]; D = self.D; W = self.Wq[p]
        H = None
        for i in range(D - 1, -1, -1):
            Z = arb_mat(self.n2, self.N)
            for j in range(D):
                Z = Z + G[i][j] * self.Y[j]
            WZ = W * Z
            H = WZ if H is None else WZ + self.IsT * H
        AP = GC.sub(H, self.Nw, self.N)                 # <P~_k, B P~_l>
        out = self.QT * AP * self.QTN.transpose()       # <q_k, B q_l>
        self.cache[key] = out
        return out


def tvec_b(N, b, kmax):
    th = GC.th
    K0 = N // 2
    tau = lambda m: th * (arb(1) / (2 * m + b - 1) + arb(b) / (2 * m + b) ** 2 + arb(1) / (2 * m + b + 1))
    t = [arb(1)]
    for k in range(1, kmax + 1):
        if k <= K0:
            t.append(t[-1] * tau(N - k + 1))
        else:
            t.append(t[K0] * th ** (k - K0) / GC.FAC[k - K0])
    return t


_TC = {}
def tail_gam(i, j, gam, r):
    return 4 * (i + j + 8) ** 2 * (arb(i + 2 * r) * (j + 2 * r)) ** (2 * r) * arb(gam) ** (i + j)

def tail_norm_r(D, gam, r, p):
    key = ('eps', D, gam, r, p)
    if key not in _TC:
        th = GC.th; q = arb(gam) * th; tot = arb(0)
        for i in range(340):
            for j in range(340):
                if i >= D or j >= D:
                    tot += 4 * (i + j + 8) ** 2 * (arb(i + 2 * r) * (j + 2 * r)) ** (2 * r) * q ** (i + j) / (GC.FAC[i] * GC.FAC[j])
        _TC[key] = (tot + arb(10) ** -60) * th ** p
    return _TC[key]

def C_eff_b(G, p, D, gam, r, cfs, N, b):
    th = GC.th
    key = ('t', N, b)
    if key not in _TC:
        _TC[key] = tvec_b(N, b, 345)
    t = _TC[key]
    a = lambda i, j: abs(G[i][j])
    s = arb(0)
    if p == 0:
        s += 2 * a(0, 1) * t[1]
    else:
        s += a(0, 1) * (b + p) * th ** (p - 1) * t[1] ** 2
    for j in range(2, D):
        s += 2 * a(0, j) * th ** p * t[j]
    for i in range(1, D):
        for j in range(1, D):
            s += a(i, j) * th ** p * t[i] * t[j]
    tk = ('tail', N, D, gam, r, b)
    if tk not in _TC:
        tt = arb(0)
        for i in range(340):
            for j in range(340):
                if i >= D or j >= D:
                    if i == 0 or j == 0:
                        w = 2 * t[max(i, j)]          # (0,j) with j >= D >= 2: 2|Gam_0j| th^p t_j, th^p <= 1
                    else:
                        w = t[i] * t[j]
                    tt += tail_gam(i, j, gam, r) * w
        _TC[tk] = tt + arb(10) ** -60
    return s + _TC[tk] * cfs


def certify_r(ops, terms, label=''):
    N, Nw, D, gam, r, b = ops.N, ops.Nw, ops.D, ops.gam, ops.r, ops.b
    M = arb_mat(Nw, N)
    for (h, p, cf) in terms:
        M = M + cf * ops.A(h, p)
    B = M.tolist()
    eps = sum((abs(cf) * tail_norm_r(D, gam, r, p) for (h, p, cf) in terms), arb(0))
    dD = None; nsub = 256
    for tt in range(nsub):
        lo = GC.A_(GC.TH) * tt / nsub; hi = GC.A_(GC.TH) * (tt + 1) / nsub
        Xb = (lo + hi) / 2 + arb(0, ((hi - lo) / 2).upper())
        dX = sum((cf * ops.Gam[h][0][0] * Xb ** p for (h, p, cf) in terms), arb(0))
        l = dX.lower()
        dD = l if dD is None or l < dD else dD
    dD = arb(dD)
    byp = {}
    for (h, p, cf) in terms:
        G = ops.Gam[h]
        if p not in byp:
            byp[p] = [[cf * G[m][n] for n in range(D)] for m in range(D)]
        else:
            byp[p] = [[byp[p][m][n] + cf * G[m][n] for n in range(D)] for m in range(D)]
    C = arb(0)
    for p, G in byp.items():
        cfs = sum((abs(cf) for (h, pp, cf) in terms if pp == p), arb(0))
        C += C_eff_b(G, p, D, gam, r, cfs, N, b)
    dW = dD - C
    info = {'label': label, 'deltaD': dD.str(5), 'C': C.str(5), 'deltaW': dW.str(5), 'eps_tail': eps.str(3)}
    if not dW > 0:
        info['fail'] = 'deltaW <= 0'
        return False, info
    S = [[B[k][l] + arb(0, eps.upper()) for l in range(N)] for k in range(N)]
    for k in range(N):
        for l in range(N):
            s = arb(0)
            for w in range(N, Nw):
                s += B[w][k] * B[w][l]
            S[k][l] -= (1 + GC.ETA) * s / dW
    for k in range(N):
        S[k][k] -= (1 + 1 / GC.ETA) * eps ** 2 / dW
    ok, kk, d = GC.cholesky_ok(arb_mat(S))
    info['cholesky'] = 'OK' if ok else ('FAIL at %d (%s)' % (kk, d.str(5)))
    return ok, info
