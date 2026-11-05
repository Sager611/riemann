#!/usr/bin/env python3
"""Independent checker for leaves of the reduced Krzyz critical-stratum search (N = n-1).

Written for the audit of Zenodo 23964169 (Wikstrom).  It does NOT import the author's code.
For a parameter cell P = theta_1..theta_N x q and a mass box L it tries to prove that no
(theta, q, lam) in P x L satisfies the necessary conditions for an n-1 atom extremal:
   (M)  sum_j (lam_j - 1/k) e^{i k theta_j} = T_k(q)/k,  k = 1..n   (Thm kkt, rotated)
   (S)  lam_j >= 0,  T_lo <= t = sum lam_j <= T_hi(n)
   (V)  e^{-t} C(theta, q) >= 2/e   (|a_n| on the moment locus, Prop. central-coeff)
   (D)  theta_1 <= ... <= theta_N,  theta_1 + theta_N <= 2 pi,  q >= 1.
Method (different from the production kernel): tight per-entry interval ranges of the
trigonometric matrix, LP relaxation with exact (Arb-verified) Farkas/duality certificates
(scipy-HiGHS only proposes multipliers), LP-certified mass and total-mass bounds, a value cut
t <= ln C_max + 1 - ln 2, a Krawczyk operator on the full (2n-1)-dimensional system
(theta, q, lam), and bisection.  All proof arithmetic is python-flint Arb ball arithmetic.
"""
import os
for _k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    os.environ[_k] = '1'
import math
import numpy as np
from scipy.optimize import linprog
from flint import arb, acb, arb_mat, ctx

ctx.prec = 128
PI = arb.pi()
TWO_PI = 2 * PI
TWO_OVER_E = 2 / arb(1).exp()
LN_TWO_OVER_E = TWO_OVER_E.log()
TLO_BALL = -(1 - TWO_OVER_E).log() / 2
TLO = TLO_BALL.lower()                      # rigorous lower bound of T_lo
THI = {3: 3.54, 4: 4.69, 5: 5.84, 6: 6.99, 7: 8.15}   # binary64 cutoffs used by production


def iv(lo, hi):
    """Ball containing [lo, hi] (lo, hi exact arb points or floats)."""
    return arb(lo).union(arb(hi))


def lower(x):
    return x.lower()


def upper(x):
    return x.upper()


def certainly_lt(a, b):
    return bool(arb(a).upper() < arb(b).lower())


# ---------------------------------------------------------------- trig ranges
def trig_range(a, b, k):
    """Exact-endpoint enclosure of cos(k th), sin(k th) for th in [a, b] (a <= b exact)."""
    A = k * a
    Bv = k * b
    ca, sa = A.cos(), A.sin()
    cb, sb = Bv.cos(), Bv.sin()
    clo = ca.lower().min(cb.lower())
    chi = ca.upper().max(cb.upper())
    slo = sa.lower().min(sb.lower())
    shi = sa.upper().max(sb.upper())
    m0 = int(math.floor(float((2 * A / PI).lower()))) - 1
    m1 = int(math.ceil(float((2 * Bv / PI).upper()))) + 1
    for m in range(m0, m1 + 1):
        t = m * PI / 2
        if t.upper() < A.lower() or t.lower() > Bv.upper():
            continue
        r = m % 4
        if r == 0:
            chi = arb(1)
        elif r == 1:
            shi = arb(1)
        elif r == 2:
            clo = arb(-1)
        else:
            slo = arb(-1)
    return iv(clo.max(arb(-1)), chi.min(arb(1))), iv(slo.max(arb(-1)), shi.min(arb(1)))


def cheb_T(q, k):
    a, b = arb(1), q
    if k == 0:
        return a
    for _ in range(k - 1):
        a, b = b, 2 * q * b - a
    return b


def cheb_U(q, k):  # U_k
    a, b = arb(1), 2 * q
    if k == 0:
        return a
    for _ in range(k - 1):
        a, b = b, 2 * q * b - a
    return b


def T_range(qlo, qhi, k):
    """T_k increasing on [1, inf): range over [qlo, qhi], qlo >= 1."""
    return iv(cheb_T(qlo, k).lower(), cheb_T(qhi, k).upper())


def dT_range(qlo, qhi, k):
    """T_k' = k U_{k-1}, increasing on [1, inf)."""
    if k == 1:
        return arb(1)
    return iv((k * cheb_U(qlo, k - 1)).lower(), (k * cheb_U(qhi, k - 1)).upper())


# ---------------------------------------------------------------- value C(theta, q)
def poly_coeffs(alphas):
    p = [acb(1)]
    for a in alphas:
        nxt = [acb(0)] * (len(p) + 1)
        for k, c in enumerate(p):
            nxt[k] -= a * c
            nxt[k + 1] += c
        p = nxt
    return p


def C_of(alphas, q):
    p = poly_coeffs(alphas)
    s0 = sum((c.real * c.real + c.imag * c.imag for c in p), arb(0))
    s1 = sum(((p[k + 1] * p[k].conjugate()).real for k in range(len(p) - 1)), arb(0))
    return 2 * q * s0 - 2 * s1, p, s0


def C_range(P, N):
    """Natural extension intersected with a mean-value form over the box P."""
    al = []
    for j in range(N):
        c, s = trig_range(P[j][0], P[j][1], 1)
        al.append(acb(c, s))
    qb = iv(P[N][0], P[N][1])
    Cnat, p, s0 = C_of(al, qb)
    # gradient enclosure
    grads = []
    for l in range(N):
        rest = poly_coeffs([al[j] for j in range(N) if j != l])
        dp = [(-acb(0, 1) * al[l]) * c for c in rest] + [acb(0)]
        d0 = 2 * sum(((dp[k] * p[k].conjugate()).real for k in range(len(p))), arb(0))
        d1 = sum(((dp[k + 1] * p[k].conjugate() + p[k + 1] * dp[k].conjugate()).real for k in range(len(p) - 1)), arb(0))
        grads.append(2 * qb * d0 - 2 * d1)
    grads.append(2 * s0)
    ctr = [(P[j][0] + P[j][1]) / 2 for j in range(N + 1)]
    rad = [(P[j][1] - P[j][0]) / 2 for j in range(N + 1)]
    alc = [acb(ctr[j].cos(), ctr[j].sin()) for j in range(N)]
    C0 = C_of(alc, ctr[N])[0]
    spread = sum((abs(g).upper() * rad[k] for k, g in enumerate(grads)), arb(0))
    mv = C0 + arb(0, spread.upper())
    lo = Cnat.lower().max(mv.lower())
    hi = Cnat.upper().min(mv.upper())
    return lo, hi


# ---------------------------------------------------------------- moment system
def system_intervals(P, n):
    """Interval A (2n x N) and b (2n) over P (tight per entry)."""
    N = n - 1
    A = [[None] * N for _ in range(2 * n)]
    b = [None] * (2 * n)
    for k in range(1, n + 1):
        sc = sr = arb(0)
        for j in range(N):
            c, s = trig_range(P[j][0], P[j][1], k)
            A[2 * k - 2][j] = c
            A[2 * k - 1][j] = s
            sc += c
            sr += s
        Tk = T_range(P[N][0], P[N][1], k)
        b[2 * k - 2] = (sc + Tk) / k
        b[2 * k - 1] = sr / k
    return A, b


def lp_rows(A, b, n, tlo, tcap):
    """Rows G lam <= g, as exact arb points, valid for lam >= 0."""
    N = n - 1
    G, g = [], []
    for i in range(2 * n):
        G.append([A[i][j].lower() for j in range(N)])
        g.append(b[i].upper())
        G.append([-A[i][j].upper() for j in range(N)])
        g.append(-b[i].lower())
    G.append([arb(-1)] * N)
    g.append(-tlo)
    G.append([arb(1)] * N)
    g.append(tcap)
    return G, g


def _floats(M):
    return np.array([[float(x) for x in row] for row in M])


def certify_infeasible(G, g, l, u):
    """Phase-1 LP proposes y >= 0; Arb verifies min_{l<=lam<=u} y^T G lam > y^T g."""
    Gf = _floats(G)
    gf = np.array([float(x) for x in g])
    m, N = Gf.shape
    c = np.concatenate([np.zeros(N), np.ones(m)])
    Aub = np.hstack([Gf, -np.eye(m)])
    bounds = [(float(l[j]), float(u[j])) for j in range(N)] + [(0, None)] * m
    try:
        res = linprog(c, A_ub=Aub, b_ub=gf, bounds=bounds, method='highs')
    except Exception:
        return False
    if res.status != 0 or res.fun <= 1e-13:
        return False
    y = [max(0.0, -float(v)) for v in res.ineqlin.marginals]
    return verify_farkas(G, g, l, u, y)


def verify_farkas(G, g, l, u, y):
    N = len(G[0])
    yb = [arb(v) for v in y]  # exact doubles
    cvec = [sum((yb[i] * G[i][j] for i in range(len(G))), arb(0)) for j in range(N)]
    lhs = arb(0)
    for j in range(N):
        lhs += (cvec[j] * l[j]).min(cvec[j] * u[j])
    rhs = sum((yb[i] * g[i] for i in range(len(G))), arb(0))
    return bool(lhs.lower() > rhs.upper())


def certified_max(G, g, l, u, cobj):
    """Rigorous upper bound of cobj^T lam over {G lam <= g, l <= lam <= u} via LP duality."""
    Gf = _floats(G)
    gf = np.array([float(x) for x in g])
    N = Gf.shape[1]
    bounds = [(float(l[j]), float(u[j])) for j in range(N)]
    try:
        res = linprog(-np.array([float(x) for x in cobj]), A_ub=Gf, b_ub=gf, bounds=bounds, method='highs')
    except Exception:
        return None
    if res.status != 0:
        return None
    y = [max(0.0, -float(v)) for v in res.ineqlin.marginals]
    yb = [arb(v) for v in y]
    val = sum((yb[i] * g[i] for i in range(len(G))), arb(0))
    for j in range(N):
        cj = arb(cobj[j]) - sum((yb[i] * G[i][j] for i in range(len(G))), arb(0))
        val += (cj * l[j]).max(cj * u[j])
    return val.upper()


# ---------------------------------------------------------------- Krawczyk
def F_point(u, n):
    N = n - 1
    out = []
    for k in range(1, n + 1):
        re = -cheb_T(u[N], k) / k
        im = arb(0)
        for j in range(N):
            a = u[N + 1 + j] - arb(1) / k
            re += a * (k * u[j]).cos()
            im += a * (k * u[j]).sin()
        out += [re, im]
    return out


def J_box(U, n):
    """Interval Jacobian of F over the box U (list of balls)."""
    N = n - 1
    D = 2 * N + 1
    J = [[arb(0)] * D for _ in range(2 * n)]
    qlo, qhi = U[N].lower(), U[N].upper()
    for k in range(1, n + 1):
        for j in range(N):
            c, s = trig_range(U[j].lower(), U[j].upper(), k)
            a = U[N + 1 + j] - arb(1) / k
            J[2 * k - 2][j] = -a * k * s
            J[2 * k - 1][j] = a * k * c
            J[2 * k - 2][N + 1 + j] = c
            J[2 * k - 1][N + 1 + j] = s
        J[2 * k - 2][N] = -dT_range(qlo, qhi, k) / k
        J[2 * k - 1][N] = arb(0)
    return J


def J_point_float(u, n):
    N = n - 1
    D = 2 * N + 1
    J = np.zeros((2 * n, D))
    for k in range(1, n + 1):
        for j in range(N):
            th, lam = u[j], u[N + 1 + j]
            a = lam - 1.0 / k
            J[2 * k - 2, j] = -a * k * math.sin(k * th)
            J[2 * k - 1, j] = a * k * math.cos(k * th)
            J[2 * k - 2, N + 1 + j] = math.cos(k * th)
            J[2 * k - 1, N + 1 + j] = math.sin(k * th)
        q = u[N]
        # T_k'(q) = k U_{k-1}(q)
        Ua, Ub = 1.0, 2 * q
        if k == 1:
            dT = 1.0
        else:
            for _ in range(k - 2):
                Ua, Ub = Ub, 2 * q * Ub - Ua
            dT = k * (Ub if k >= 2 else 1.0)
        J[2 * k - 2, N] = -dT / k
    return J


def krawczyk(U, n, iters=6):
    """Returns ('empty', None) | ('unique', K) | ('open', K) with K = contracted box."""
    N = n - 1
    D = 2 * N + 1
    status = 'open'
    for _ in range(iters):
        mid = [x.mid() for x in U]
        midf = [float(x) for x in mid]
        Jm = J_point_float(midf, n)
        try:
            Y = np.linalg.pinv(Jm)
        except Exception:
            return 'open', U
        if not np.all(np.isfinite(Y)):
            return 'open', U
        Yb = arb_mat(D, 2 * n, [arb(float(v)) for v in Y.ravel()])
        Fm = arb_mat(2 * n, 1, F_point(mid, n))
        JU = arb_mat(2 * n, D, [x for row in J_box(U, n) for x in row])
        Mm = arb_mat(D, D, [arb(1) if i == j else arb(0) for i in range(D) for j in range(D)]) - Yb * JU
        dU = arb_mat(D, 1, [U[i] - mid[i] for i in range(D)])
        K = arb_mat(D, 1, mid) - Yb * Fm + Mm * dU
        newU = []
        inside = True
        for i in range(D):
            ki = K[i, 0]
            lo = ki.lower().max(U[i].lower())
            hi = ki.upper().min(U[i].upper())
            if lo > hi:
                return 'empty', None
            if not (ki.lower() > U[i].lower() and ki.upper() < U[i].upper()):
                inside = False
            newU.append(iv(lo, hi))
        shrink = max(float(newU[i].rad() / U[i].rad()) if U[i].rad() > 0 else 0.0 for i in range(D))
        U = newU
        if inside:
            status = 'unique'
        if shrink > 0.95 and not inside:
            break
    return status, U


def residual_excludes(U, n):
    """Mean-value enclosure of z^T F over U for z spanning the left kernel of J(mid)."""
    N = n - 1
    D = 2 * N + 1
    mid = [x.mid() for x in U]
    Jm = J_point_float([float(x) for x in mid], n)
    Usvd, S, Vt = np.linalg.svd(Jm)
    z = Usvd[:, D]
    zb = [arb(float(v)) for v in z]
    Fm = F_point(mid, n)
    JU = J_box(U, n)
    val = sum((zb[i] * Fm[i] for i in range(2 * n)), arb(0))
    for c in range(D):
        coef = sum((zb[i] * JU[i][c] for i in range(2 * n)), arb(0))
        val += coef * (U[c] - mid[c])
    return not val.contains(0)


# ---------------------------------------------------------------- one box
def domain_contract(P, N):
    P = [list(x) for x in P]
    for j in range(1, N):
        P[j][0] = P[j][0].max(P[j - 1][0])
    for j in range(N - 2, -1, -1):
        P[j][1] = P[j][1].min(P[j + 1][1])
    if (P[0][0] + P[N - 1][0]).lower() > TWO_PI.upper():
        return None
    P[0][1] = P[0][1].min((TWO_PI - P[N - 1][0]).upper())
    P[N - 1][1] = P[N - 1][1].min((TWO_PI - P[0][0]).upper())
    P[N][0] = P[N][0].max(arb(1))
    for j in range(N + 1):
        if P[j][0] > P[j][1]:
            return None
    return P


def test_box(P, L, n, stats):
    """Returns (reason, None, None) if excluded, else (None, P', L')."""
    N = n - 1
    T = arb(THI[n])
    P = domain_contract(P, N)
    if P is None:
        return 'domain', None, None
    Clo, Chi = C_range(P, N)
    tcap = (Chi.log() - TWO_OVER_E.log()).upper().min(T)
    if tcap < TLO:
        return 'value_cut', None, None
    l = [L[j][0].max(arb(0)) for j in range(N)]
    u = [L[j][1].min(T) for j in range(N)]
    if any(l[j] > u[j] for j in range(N)):
        return 'mass_box', None, None
    A, b = system_intervals(P, n)
    G, g = lp_rows(A, b, n, TLO, tcap)
    if certify_infeasible(G, g, l, u):
        return 'lp', None, None
    # LP-certified bounds on masses and on t
    for j in range(N):
        e = [arb(1) if i == j else arb(0) for i in range(N)]
        hi = certified_max(G, g, l, u, e)
        if hi is not None and hi < u[j]:
            u[j] = hi
        lo = certified_max(G, g, l, u, [-x for x in e])
        if lo is not None and -lo > l[j]:
            l[j] = -lo
        if l[j] > u[j]:
            return 'lp_bounds', None, None
    tmin = certified_max(G, g, l, u, [arb(-1)] * N)
    tmin = -tmin if tmin is not None else TLO
    tmin = tmin.max(sum(l, arb(0)))
    if (-tmin).exp().upper() * Chi < TWO_OVER_E.lower():
        return 'lp_value', None, None
    wmax = max(float(P[j][1] - P[j][0]) for j in range(N))
    if wmax < 0.35:
        U = [iv(P[j][0], P[j][1]) for j in range(N + 1)] + [iv(l[j], u[j]) for j in range(N)]
        st, K = krawczyk(U, n)
        if st == 'empty':
            return 'kraw_empty', None, None
        lamK = K[N + 1:]
        if any(x.upper() < 0 for x in lamK):
            return 'kraw_negmass', None, None
        tK = sum(lamK, arb(0))
        if tK.upper() < TLO or tK.lower() > T:
            return 'kraw_mass', None, None
        PK = [[K[j].lower(), K[j].upper()] for j in range(N + 1)]
        PKc = domain_contract(PK, N)
        if PKc is None:
            return 'kraw_domain', None, None
        Clo2, Chi2 = C_range(PKc, N)
        if (-tK.lower().max(TLO)).exp().upper() * Chi2 < TWO_OVER_E.lower():
            return 'kraw_value', None, None
        if residual_excludes(K, n):
            return 'kraw_residual', None, None
        P = PKc
        l = [K[N + 1 + j].lower().max(arb(0)) for j in range(N)]
        u = [K[N + 1 + j].upper().min(T) for j in range(N)]
    return None, P, [[l[j], u[j]] for j in range(N)]


def split(P, n):
    N = n - 1
    best, bi = -1.0, -1
    for j in range(N + 1):
        w = float(P[j][1] - P[j][0]) * (n * n if j == N else n)
        if w > best:
            best, bi = w, j
    mid = (P[bi][0] + P[bi][1]) / 2
    a = [list(x) for x in P]
    b = [list(x) for x in P]
    a[bi][1] = mid
    b[bi][0] = mid
    return a, b


def check_cell(n, lo, hi, ml=None, mh=None, limit=3000, stats=None):
    """Prove exclusion of a logged cell. If ml/mh is None, masses range over [0, T_hi(n)]."""
    N = n - 1
    P = [[arb(lo[j]), arb(hi[j])] for j in range(N + 1)]
    if ml is None:
        L = [[arb(0), arb(THI[n])] for _ in range(N)]
    else:
        L = [[arb(ml[j]), arb(mh[j])] for j in range(N)]
    stack = [(P, L, 0)]
    nodes = 0
    reasons = {}
    maxdepth = 0
    while stack:
        P, L, d = stack.pop()
        nodes += 1
        maxdepth = max(maxdepth, d)
        why, P2, L2 = test_box(P, L, n, stats)
        if why:
            reasons[why] = reasons.get(why, 0) + 1
            continue
        if nodes >= limit:
            return False, nodes, reasons, maxdepth
        c1, c2 = split(P2, n)
        stack.append((c1, L2, d + 1))
        stack.append((c2, L2, d + 1))
    return True, nodes, reasons, maxdepth
