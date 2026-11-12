"""gapceil: Arb certificate that Bui-Milinovich's large-gap criterion (arXiv 1410.3635v3, r = 1) fails at a given kappa
for EVERY admissible amplifier (polynomial P of any degree, any v in R, any theta < 1/4).  See DERIVATION.md.

Reduction (DERIVATION.md 2-3): with X_i = theta x_i, tau1 = a t1, tau2 = b t2 (ab = Jacobian), f(U) = P(1 - U/theta) on [0, theta],
  c_j ∝ Q_j[f] = ∫ f(X+X1+X2) f(X+X3+X4) rho_j,  rho_j = ∫_0^{1-X1-X3} ∫_0^{1-X2-X4} (v-S)^{2j} U conj(U') dτ1 dτ2,
  U conj U' = e^{ic(X1-X3)} + e^{ic(X2-X4)} - e^{ic(X1-X4)} e^{ic(τ1-τ2)} - e^{ic(X2-X3)} e^{-ic(τ1-τ2)},  c = kappa pi,
  S = X + X1+X2+X3+X4 + τ1+τ2; theta enters only through supp f, so f in L2[0,1/4] covers every theta <= 1/4.
  phi = I^2 f, I g(X) = ∫_X^{1/4} g:  Q_0 = <phi, M0 phi>,  Q_1(v) = <phi, (v^2 M0 - 2v MS + MSS) phi>, each operator
  sum_{i,j>=0} Gam_ij (I^i)^T X^p I^j  (Gam = fiber-averaged Taylor data of rho_j, arrays 'a','b','c').
Certified: kappa^2 Q_1(v) - Q_0 >= 0 for all phi in L2[0,1/4] and all v  ==>  h1 = c0/(kappa^2 c1) <= 1 for every amplifier.
Each operator B = D(X) (multiplication, the (0,0) term) + K; PSD via V = span{P_n : n < N}, W = V^perp:
  <w,Bw> >= (delta_D - C)|w|^2, C from integration by parts ((Iw)(0) = <w,1> = 0, (Iw)(1/4) = 0) and |I^k w| <= t_k |w|;
  <phi,B phi> >= <p,Bp> - |P_W B p|^2/delta_W, |P_W B p|^2 = p^T E^T Lam E p (+ tail), tested by an Arb Cholesky.
"""
import sys, time, json, resource, math
from math import comb, factorial
from flint import arb, acb, arb_mat, acb_poly, fmpq, ctx

TH = fmpq(1, 4)
def A_(q):
    q = fmpq(q)
    return arb(q.p) / arb(q.q)
def dec(s):
    neg = s.startswith('-'); s = s.lstrip('+-'); a, b = (s.split('.') + [''])[:2]
    q = fmpq(int(a + b), 10 ** len(b)); return -q if neg else q
def set_prec(bits):
    """set working precision and rebuild the precision-dependent constants."""
    global th, FAC, ETA
    ctx.prec = bits
    th = A_(TH)
    FAC = [arb(factorial(i)) for i in range(700)]
    ETA = arb(10) ** -6
set_prec(128)

# ---------------- Part 1: Gamma arrays ----------------
def h_series(k, z, zero, Q):
    """ordinary Taylor coefficients (q < Q) of h(y) = m_k(1-y, z) = ∫_0^{1-y} t^k e^{zt} dt."""
    if zero:
        # (1-y)^{k+1}/(k+1)
        out = [acb(0)] * Q
        for j in range(k + 2):
            if j < Q:
                out[j] = acb(comb(k + 1, j) * (-1) ** j) / (k + 1)
        return out
    ez = z.exp()
    e = [acb(1)]
    for n in range(1, Q):
        e.append(e[-1] * (-z) / n)
    g = [acb(0)] * Q
    for i in range(k + 1):
        a = acb(comb(k, i) * (-1) ** i)
        for n in range(Q - i):
            g[i + n] += a * e[n]
    g = [x * ez for x in g]
    m0 = (ez - 1) / z; m1 = (ez - m0) / z; m2 = (ez - 2 * m1) / z
    m = [m0, m1, m2][k]
    return [m] + [-g[q - 1] / q for q in range(1, Q)]

def uhat(al1, al3, k, z, zero, D, gam):
    """scaled derivative array  uhat[a][b] / gam^(a+b),  uhat = d^a_X1 d^b_X3 [e^{i al1 X1 + i al3 X3} h(X1+X3)] at 0."""
    h = h_series(k, z, zero, 2 * D)
    I = acb(0, 1)
    # EGF form: uhat_ab = sum C(a,n1) C(b,n2) (i al1)^{a-n1} (i al3)^{b-n2} hhat_{n1+n2}, hhat_q = q! h_q ; scale by gam^-(a+b)
    gi = 1 / arb(gam)
    hh = [h[q] * FAC[q] * gi ** q for q in range(2 * D)]
    p1 = [(I * al1 * gi) ** n for n in range(D)]
    p3 = [(I * al3 * gi) ** n for n in range(D)]
    # first along b: R[n1][b] = sum_{n2<=b} C(b,n2) p3[b-n2] hh[n1+n2]
    R = [[sum((comb(b, n2) * p3[b - n2] * hh[n1 + n2] for n2 in range(b + 1)), acb(0)) for b in range(D)] for n1 in range(D)]
    U = [[sum((comb(a, n1) * p1[a - n1] * R[n1][b] for n1 in range(a + 1)), acb(0)) for b in range(D)] for a in range(D)]
    return U

def conv2(U, W, D):
    """2D convolution (truncated to D x D) via Kronecker substitution."""
    L = 2 * D - 1
    pu = [acb(0)] * (D * L); pw = [acb(0)] * (D * L)
    for a in range(D):
        for b in range(D):
            pu[a * L + b] = U[a][b]; pw[a * L + b] = W[a][b]
    P = acb_poly(pu) * acb_poly(pw)
    cf = P.coeffs()
    cf = cf + [acb(0)] * (D * L * 2 - len(cf))
    return [[cf[M * L + N] for N in range(D)] for M in range(D)]

def gamma_arrays(c, D, gam):
    """real symmetric arrays Gam^(a), Gam^(b), Gam^(c) (true scale), index [i][j] for (I^i, I^j)."""
    I = acb(0, 1)
    ic = I * c
    terms = [(1, (c, -c), (0, 0), (acb(0), acb(0)), True),     # TA: phase c(X1-X3); u-side (X1,X3), w-side (X2,X4)
             (1, (0, 0), (c, -c), (acb(0), acb(0)), True),     # TB: phase c(X2-X4)
             (-1, (c, 0), (0, -c), (ic, -ic), False),           # TC: phase c(X1-X4), eps=+1
             (-1, (0, -c), (c, 0), (-ic, ic), False)]           # TD: phase c(X2-X3), eps=-1
    g = [[[acb(0)] * D for _ in range(D)] for _ in range(3)]
    for s, (a1, a3), (a2, a4), (z1, z2), zero in terms:
        us = [uhat(arb(a1) if not isinstance(a1, arb) else a1, arb(a3) if not isinstance(a3, arb) else a3, k, z1, zero, D, gam) for k in range(3)]
        ws = [uhat(arb(a2) if not isinstance(a2, arb) else a2, arb(a4) if not isinstance(a4, arb) else a4, k, z2, zero, D, gam) for k in range(3)]
        combos = [[(0, 0, 1)], [(1, 0, 1), (0, 1, 1)], [(2, 0, 1), (1, 1, 2), (0, 2, 1)]]
        for gi, cl in enumerate(combos):
            for (k, l, mult) in cl:
                C = conv2(us[k], ws[l], D)
                for M in range(D):
                    for N in range(D):
                        g[gi][M][N] += s * mult * C[M][N]
    # unscale: true Gam_MN = gam^(M+N) * scaled
    gp = [arb(gam) ** n for n in range(2 * D)]
    g = [[[x[M][N] * gp[M + N] for N in range(D)] for M in range(D)] for x in g]
    def sh(G):
        return [[(G[M - 1][N] * (M + 1) if M > 0 else 0) + (G[M][N - 1] * (N + 1) if N > 0 else 0) for N in range(D)] for M in range(D)]
    a = g[0]
    b = [[x + y for x, y in zip(r1, r2)] for r1, r2 in zip(sh(a), g[1])]
    s1 = sh(g[1]); ss = sh(sh(a))
    cc = [[ss[M][N] + 2 * s1[M][N] + g[2][M][N] for N in range(D)] for M in range(D)]
    def sym(G):
        return [[(G[M][N].real + G[N][M].real) / 2 for N in range(D)] for M in range(D)]
    return {'a': sym(a), 'b': sym(b), 'c': sym(cc)}

def gam_bound(i, j, gam):
    """|Gam^(h)_ij| <= 4 (i+j+4)^2 (i+1)(j+1) gam^(i+j) for h in a,b,c."""
    return 4 * (i + j + 4) ** 2 * (i + 1) * (j + 1) * arb(gam) ** (i + j)

# ---------------- Part 2: Legendre (unnormalized P_n(2X/th-1)) ----------------
def legendre_ops(n):
    I = [[arb(0)] * n for _ in range(n)]
    T = [[arb(0)] * n for _ in range(n)]
    h2 = th / 2
    I[0][0] = h2
    if n > 1: I[1][0] = -h2
    for k in range(1, n):
        I[k - 1][k] = h2 / (2 * k + 1)
        if k + 1 < n: I[k + 1][k] = -h2 / (2 * k + 1)
    for k in range(n):
        T[k][k] = h2
        if k + 1 < n: T[k + 1][k] = h2 * (k + 1) / (2 * k + 1)
        if k >= 1: T[k - 1][k] = h2 * k / (2 * k + 1)
    G0 = [th / (2 * k + 1) for k in range(n)]
    return arb_mat(I), arb_mat(T), G0

def sub(M, r, c):
    L = M.tolist()
    return arb_mat([[L[i][j] for j in range(c)] for i in range(r)])

class Ops:
    """Gram-type blocks A_(h,p)[k,l] = <P_k, G_p^(h) P_l>, k < Nw, l < N, for the truncated (i,j < D) series.
    Horner: A = sum_i (I^T)^i W_p Z_i, Z_i = sum_j Gam_ij I^j E_N."""
    def __init__(self, kap, N, D, gam, Nw=None):
        t0 = time.time()
        self.kap, self.N, self.D, self.gam = kap, N, D, gam
        self.c = A_(kap) * arb.pi()
        self.Gam = gamma_arrays(self.c, D, gam)
        self.t_gam = time.time() - t0
        self.Nw = Nw or (N + 2 * D + 4)
        n2 = self.Nw + D + 6
        I, T, G0 = legendre_ops(n2 + 3)
        G0d = arb_mat([[G0[i] if i == j else arb(0) for j in range(n2 + 3)] for i in range(n2 + 3)])
        Wp = [G0d, G0d * T, G0d * T * T]
        self.Wp = [sub(W, n2, n2) for W in Wp]
        Is = sub(I, n2, n2)
        self.IsT = Is.transpose()
        Y = arb_mat([[arb(1) if i == j else arb(0) for j in range(N)] for i in range(n2)])
        self.Y = [Y]
        for k in range(1, D):
            self.Y.append(Is * self.Y[-1])
        self.n2 = n2
        self.cache = {}
        self.t_leg = time.time() - t0
    def A(self, h, p):
        key = (h, p)
        if key in self.cache: return self.cache[key]
        G = self.Gam[h]; D = self.D; W = self.Wp[p]
        H = None
        for i in range(D - 1, -1, -1):
            Z = arb_mat(self.n2, self.N)
            for j in range(D):
                Z = Z + G[i][j] * self.Y[j]
            WZ = W * Z
            H = WZ if H is None else WZ + self.IsT * H
        out = sub(H, self.Nw, self.N)
        self.cache[key] = out
        return out

# ---------------- Part 3: certificate ----------------
_TC = {}
def tail_norm(D, gam, p, maxi=320):
    if (D, gam, p) in _TC: return _TC[(D, gam, p)]
    """bound on || sum_{i or j >= D} Gam_ij (I^i)^T X^p I^j ||  using |Gam| bound and |I^k| <= th^k/k!."""
    q = arb(gam) * th
    tot = arb(0)
    for i in range(maxi):
        for j in range(maxi):
            if i >= D or j >= D:
                tot += 4 * (i + j + 4) ** 2 * (i + 1) * (j + 1) * q ** (i + j) / (FAC[i] * FAC[j])
    _TC[(D, gam, p)] = (tot + arb(10) ** -80) * th ** p
    return _TC[(D, gam, p)]

def cholesky_ok(S):
    L = S.tolist(); n = len(L)
    for k in range(n):
        d = L[k][k]
        for j in range(k):
            d -= L[k][j] * L[k][j]
        if not d > 0:
            return False, k, d
        d = d.sqrt(); L[k][k] = d
        for i in range(k + 1, n):
            s = L[i][k]
            for j in range(k):
                s -= L[i][j] * L[k][j]
            L[i][k] = s / d
    return True, n, None

def certify(ops, terms, label=''):
    """terms: list of (h, p, coef) defining the symmetric operator B = sum coef * G_p^(h). Returns (ok, info)."""
    N, Nw, D, gam = ops.N, ops.Nw, ops.D, ops.gam
    M = arb_mat(Nw, N)
    for (h, p, cf) in terms:
        M = M + cf * ops.A(h, p)
    B = M.tolist()
    # operator-norm bound of the series tail (i or j >= D)
    eps = sum((abs(cf) * tail_norm(D, gam, p) for (h, p, cf) in terms), arb(0))
    # multiplication part d(X) = sum coef Gam^(h)_00 X^p, lower bound on [0, 1/4] by 256 subintervals
    dD = None
    nsub = 256
    for t in range(nsub):
        lo = A_(TH) * t / nsub; hi = A_(TH) * (t + 1) / nsub
        Xb = (lo + hi) / 2 + arb(0, ((hi - lo) / 2).upper())
        dX = sum((cf * ops.Gam[h][0][0] * Xb ** p for (h, p, cf) in terms), arb(0))
        l = dX.lower()
        dD = l if dD is None or l < dD else dD
    dD = arb(dD)
    # C: |<w, K w>| <= C |w|^2 on W
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
        C += C_eff(G, p, D, gam, cfs, N)
    dW = dD - C
    info = {'label': label, 'deltaD': dD.str(5), 'C': C.str(5), 'deltaW': dW.str(5), 'eps_tail': eps.str(3)}
    if not dW > 0:
        info['fail'] = 'deltaW <= 0'
        return False, info
    # S = A_VV - ((1+eta) E^T Lam E + (1+1/eta) eps^2 G_V)/delta_W ; tail radius on A_VV
    Gl = [th / (2 * l + 1) for l in range(Nw)]
    lam = [(2 * k + 1) / th for k in range(Nw)]
    S = [[B[k][l] + arb(0, (eps * (Gl[k] * Gl[l]).sqrt()).upper()) for l in range(N)] for k in range(N)]
    for k in range(N):
        for l in range(N):
            s = arb(0)
            for w in range(N, Nw):
                s += B[w][k] * lam[w] * B[w][l]
            S[k][l] -= (1 + ETA) * s / dW
    for k in range(N):
        S[k][k] -= (1 + 1 / ETA) * eps ** 2 * Gl[k] / dW
    ok, kk, d = cholesky_ok(arb_mat(S))
    info['cholesky'] = 'OK' if ok else ('FAIL at %d (%s)' % (kk, d.str(5)))
    return ok, info

def tvec(N, kmax):
    """t_k >= sup_{w in W} |I^k w|/|w|, W = span{P_n : n >= N}: I maps span{n>=m} into span{n>=m-1} with norm <= 2 beta_m,
    beta_m = (th/2)/sqrt((2m-1)(2m+1)); for k > K0 = N//2 use |I^k w| <= |I^(k-K0)| |I^K0 w|, |I^j| <= th^j/j!."""
    K0 = N // 2
    beta = lambda m: (th / 2) / ((arb(2 * m - 1) * arb(2 * m + 1)).sqrt())
    t = [arb(1)]
    for k in range(1, kmax + 1):
        if k <= K0:
            t.append(t[-1] * 2 * beta(N - k + 1))
        else:
            t.append(t[K0] * th ** (k - K0) / FAC[k - K0])
    return t

def C_eff(G, p, D, gam, cfs, N, tcache={}):
    """bound C with |<w, K_p w>| <= C |w|^2 for w in W (all (i,j) != (0,0) terms of sum Gam_ij (I^i)^T X^p I^j):
       (0,1)+(1,0): Gam_01 p <Iw, X^{p-1} Iw>               -> |Gam_01| p th^{p-1} t1^2
       (0,j)+(j,0), j>=2: 2 Gam_0j <Iw, p X^{p-1} I^j w - X^p I^{j-1} w>  -> 2|Gam_0j| (p th^{p-1} t1 tj + th^p t1 t_{j-1})
       (i,j), i,j>=1: |Gam_ij| th^p t_i t_j.   (integration by parts uses (Iw)(0) = <w,1> = 0, (Iw)(th) = 0.)
       tail (i or j >= D) with |Gam_ij| <= 4(i+j+4)^2(i+1)(j+1) gam^(i+j), weights 6 t1 t_{j-1} (i=0) / t_i t_j (i,j>=1)."""
    key = (N, 330)
    if key not in tcache:
        tcache[key] = tvec(N, 330)
    t = tcache[key]
    s = arb(0)
    a = lambda i, j: abs(G[i][j])
    if p >= 1:
        s += a(0, 1) * p * th ** (p - 1) * t[1] ** 2
    for j in range(2, D):
        s += 2 * a(0, j) * (p * th ** max(p - 1, 0) * t[1] * t[j] + th ** p * t[1] * t[j - 1])
    for i in range(1, D):
        for j in range(1, D):
            s += a(i, j) * th ** p * t[i] * t[j]
    tk = ('tail', N, D, gam)
    if tk not in tcache:
        tt = arb(0)
        for i in range(320):
            for j in range(320):
                if i >= D or j >= D:
                    if i == 0 or j == 0:
                        k = max(i, j)
                        w = 2 * 6 * t[1] * t[k - 1] if k >= 2 else arb(0)
                    else:
                        w = t[i] * t[j]
                    tt += 4 * (i + j + 4) ** 2 * (i + 1) * (j + 1) * arb(gam) ** (i + j) * w
        tcache[tk] = tt + arb(10) ** -80
    return s + tcache[tk] * cfs

