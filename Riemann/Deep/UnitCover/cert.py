#!/usr/bin/env python3
"""
cert.py -- certificate for: every modulus q <= Q is covered by the subset products of O(log Q) small primes
(all moduli, explicit constants), from the Lean-checked zero-free half-plane Re s > B for every Dirichlet L-function.
Run (one command, one thread, < 100 MB):
    uv run --no-project --with python-flint --with mpmath python3 cert.py all
Sections: core | delta (Lemma D) | cover (Theorem 1, table) | small (exact closure, q <= 2000) | all.
Every inequality used as a proof step is a certified sign of an arb ball (pos()); floats only choose parameters.
Proofs and notation: DERIVATION.md (section numbers in the comments).
Reused (imported, not modified): Riemann/Deep/ExplProg/cert.py -- the half-plane constant B, the kernel KER
(explprog DERIVATION sec 3.1), the zero sums Z(L) (sec 2, BMOR zero count) and the zeta sum (Lemma Zzeta).
"""
import sys, os, math, time, importlib.util
from fractions import Fraction as Fr
sys.dont_write_bytecode = True          # no __pycache__ in this folder or in ../ExplProg

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("explprog_cert", os.path.join(HERE, "..", "ExplProg", "cert.py"))
ep = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(ep)
from flint import arb, ctx
ctx.prec = 200

A, pos, amax = ep.A, ep.pos, ep.amax
B, b, INV_b, LOG2 = ep.B, ep.b, ep.INV_b, ep.LOG2
K = ep.KER                                  # triangle kernel, ell = 7, rise 49/10 (explprog sec 3.1)
L9 = arb(9).log()
RS_THETA = A(Fr(101624, 100000))            # theta(x) < 1.01624 x (Rosser-Schoenfeld 1962, Thm 9)

def up(x):
    return float(x.upper())

def ceil_str(x, d):
    """decimal string >= every point of the ball x, with d decimals (upward rounding for statements)"""
    v = math.ceil(float(x.upper()) * 10 ** d + 1e-9) / 10 ** d
    return ("%%.%df" % d) % v

# ---------------------------------------------------------------- sec 0: inputs
def core():
    ep.core()
    print("[core] b = 1 - B = %s; zeta coefficient z1/k2 = (sum 1/gamma^2)/2 = %s; BMOR constant CB = %s (explprog)"
          % (b.str(16), (ep.SUM_INV_G2 / 2).str(6), ep.CB.str(5)))
    cinf = ep.limit_constant(K)
    print("[core] kernel: W(1) = %s, W(0) = %s, ||w||_inf <= %s; c_inf = za/(2W(1)) = %s; c_inf/delta0 = %s (delta0 = 1/1000)"
          % (K.W1.str(10), K.W0.str(10), K.winf.str(8), cinf.str(6), (cinf * 1000).str(5)))
    print("[core] limits of the method: k = 4 slope as delta -> 0: 5/(4 log 2) = %s; 1/log 2 = %s; log2(phi(1999)) = %s"
          % ((5 / (4 * LOG2)).str(5), (1 / LOG2).str(7), (arb(1998).log() / LOG2).str(5)))
    # remark (REPORT caveat): y = (820 log Q)^E exceeds Q exactly while E log(820 L) > L; find the crossover on a 0.1 grid
    Lx = 50.0
    while not (INV_b * (820 * A(Fr(Lx).limit_denominator(10))).log() < A(Fr(Lx).limit_denominator(10))):
        Lx += 0.1
    print("[core] y = (820 log Q)^E < Q from log Q = %.1f on (Q > 10^%.1f)" % (Lx, Lx / math.log(10)))

# ---------------------------------------------------------------- sec 2: Lemma D (character means of the prime law)
def pp_bound(X):
    """sum over prime powers p^k <= y, k >= 2, of log p  <=  psi(y) - theta(y) <= 1.01624 (sqrt y + y^(1/3) log y / log 2)"""
    return RS_THETA * ((X / 2).exp() + (X / 3).exp() * X / LOG2)

def delta_direct(cp, L, N, L0=None):
    """sec 2: certified upper bound for delta(y, Q, N) with y = (cp L0)^(1/b), L = max(log Q, log 9) >= log 9, N draws.
    delta >= Re E[chi(g_j) | F_{j-1}] for every q <= Q, every nonprincipal chi mod q, every j <= N.
    L0 <= L (default L0 = L): the bound is non-increasing in y and non-decreasing in L, so on a cell L in [L0, L1]
    delta_direct(cp, L1, N, L0) bounds delta at every L of the cell with y = (cp L)^(1/b) >= (cp L0)^(1/b)."""
    L = A(L); pos(L - L9 + arb(10) ** -30, "L >= log 9")
    L0 = L if L0 is None else A(L0)
    X = (A(cp) * L0).log() / b
    Y = X - K.ell
    pos(Y - 10, "Y = log(y e^-7) >= 10")
    (n0, k1, k2), (m0, k1L, k2L) = K.ks(Y)
    Z = ep.zsum_dir(L, n0, k1, k2); ZL = ep.zsum_dir(L, m0, k1L, k2L)
    z1 = ep.SUM_INV_G2 * k2 / 2; z1L = ep.SUM_INV_G2 * k2L / 2
    y = X.exp(); yB = (B * X).exp(); yh = (X / 2).exp()
    triv = K.W0 / (1 - (-2 * Y).exp())
    PP = K.winf * pp_bound(X)
    num = (yB * Z + yh * ZL) / 2 + triv + PP + 2 * K.winf * ((A(N) - 1) * X + L)
    den = y * K.W1 - (yB * z1 + yh * z1L) - triv - PP
    pos(den, "M > 0")
    return num / den, X

def delta_tail(cp, d0, Lm, C1, C2):
    """sec 2.3: certified margin of  delta(y(L), L, N) <= d0  for ALL L >= Lm, with y = (cp L)^(1/b) and N - 1 <= C1 L + C2."""
    Lm = A(Lm); pos(Lm - L9 + arb(10) ** -30, "Lm >= log 9")
    X = (A(cp) * Lm).log() / b
    Y = X - K.ell
    pos(Y - 10, "Y >= 10"); pos(X - 2, "X > 1/(B - 1/3)"); pos(X - 1 / B, "X > 1/B")
    (n0, k1, k2), (m0, k1L, k2L) = K.ks(Y)
    co = ep.zsum_coeffs(n0, k1, k2); coL = ep.zsum_coeffs(m0, k1L, k2L)
    yB = (B * X).exp(); xhB = ((A(Fr(1, 2)) - B) * X).exp(); x3B = ((A(Fr(1, 3)) - B) * X).exp()
    z1 = ep.SUM_INV_G2 * k2 / 2; z1L = ep.SUM_INV_G2 * k2L / 2
    d0 = A(d0)
    lhs = ep.Z_over_L_majorant(co, n0, Lm) / 2 + xhB * ep.Z_over_L_majorant(coL, m0, Lm) / 2
    lhs += (1 + d0) * (K.W0 / ((1 - (-2 * Y).exp()) * yB) + K.winf * RS_THETA * (xhB + x3B * X / LOG2)) / Lm
    lhs += 2 * K.winf * ((A(C1) + A(C2) / Lm) * X + 1) / yB
    lhs += d0 * (z1 + xhB * z1L) / Lm
    rhs = d0 * A(cp) * K.W1
    return (rhs - lhs) / rhs

def delta_cells(cp, d0, C1, C2, Lt):
    """sec 2.3: cells [L0, L1] covering [log 9, Lt]; on each, N - 1 <= C1 L1 + C2; returns worst certified delta"""
    L0 = L9; worst = None; n = 0
    while L0 < A(Lt):
        L1 = L0 + amax(A(Fr(1, 100)), L0 / 200)
        dd, _ = delta_direct(cp, L1, A(C1) * L1 + A(C2) + 1, L0)
        pos(A(d0) - dd, "Lemma D cell L0=%s" % L0.str(6))
        worst = dd if worst is None or dd > worst else worst
        L0 = L1; n += 1
    return worst, n, L0

def delta_section(cp=Fr(820), d0=Fr(1, 1000), C1=Fr(4), C2=Fr(200), Lt=Fr(60)):
    w, n, Le = delta_cells(cp, d0, C1, C2, Lt)
    m = delta_tail(cp, d0, Le, C1, C2)
    pos(m, "Lemma D tail")
    print("[delta] Lemma D: for every Q >= 3, y = (%s max(log Q, log 9))^(1/(1-B)), N <= %s log Q + %s draws: delta <= %s"
          % (cp, C1, C2, d0))
    print("[delta]   %d cells on [log 9, %s] (worst delta %s), tail L >= %s margin %s" % (n, Le.str(5), w.str(5), Le.str(5), m.str(4)))
    for (Lf, N) in [(L9, 30), (A(Fr(6908, 1000)), 40), (A(Fr(27631, 1000)), 100), (A(Fr(230258, 1000)), 600), (A(Fr(2302585, 1000)), 6000)]:
        dd, X = delta_direct(cp, Lf, N)
        print("[delta]   direct at log Q = %.3f, N = %d: delta <= %s, log10 y <= %s" % (float(Lf.mid()), N, dd.str(5), ceil_str(X / arb(10).log(), 3)))
    return m

# ---------------------------------------------------------------- sec 3: the covering lemma (moments of the collision)
# sec 3.3: closed relation sets.  For a k-tuple a of nonprincipal characters, R(a) = {e in {-1,0,1}^k : prod chi_i^e_i = 1}
# is "closed": R = span_Z(R) cap {-1,0,1}^k and span_Z(R) contains no unit vector.  pi(R) = sum_{e in R} 2^(-k-|e|).
# #{a : R(a) = R} <= |ker(Ghat^k -> Ghat^R)| = H^(k - rho) prod_i |Ghat[d_i]|  (rho = rank, d_i = elementary divisors).
import itertools

def hnf(rows, k):
    """Hermite normal form (row style, positive pivots, strictly increasing pivot columns) of the Z-span of rows"""
    M = [list(r) for r in rows if any(r)]
    basis = []; col = 0
    while M and col < k:
        if not any(r[col] for r in M):
            col += 1; continue
        while sum(1 for r in M if r[col] != 0) > 1:
            nz = sorted([r for r in M if r[col] != 0], key=lambda r: abs(r[col]))
            piv = nz[0]
            for r in nz[1:]:
                qq = r[col] // piv[col]
                for j in range(k): r[j] -= qq * piv[j]
            M = [r for r in M if any(r)]
        piv = [r for r in M if r[col] != 0][0]
        if piv[col] < 0: piv[:] = [-x for x in piv]
        basis.append(piv); M = [r for r in M if r is not piv]; col += 1
    return [tuple(r) for r in basis]

def in_lattice(v, basis, k):
    v = list(v)
    for r in basis:
        pc = next(j for j in range(k) if r[j] != 0)
        if v[pc] % r[pc] != 0: return False
        qq = v[pc] // r[pc]
        v = [x - qq * c for x, c in zip(v, r)]
    return not any(v)

def det(m):
    if len(m) == 1: return m[0][0]
    return sum(((-1) ** j) * m[0][j] * det([row[:j] + row[j + 1:] for row in m[1:]]) for j in range(len(m)) if m[0][j])

def elem_divisors(basis, k):
    """elementary divisors > 1 of the lattice with this basis: d_j = D_j / D_(j-1), D_j = gcd of j x j minors"""
    rho = len(basis); D = [1]
    for j in range(1, rho + 1):
        g = 0
        for rs in itertools.combinations(range(rho), j):
            for cs in itertools.combinations(range(k), j):
                g = math.gcd(g, det([[basis[r][c] for c in cs] for r in rs]))
        D.append(abs(g))
    return [D[j] // D[j - 1] for j in range(1, rho + 1) if D[j] // D[j - 1] > 1]

def closed_sets(k):
    """all closed relation sets (breadth-first closure from {0}); returns list of (pi, rho, elementary divisors, R)"""
    cube = list(itertools.product((-1, 0, 1), repeat=k))
    units = [tuple(1 if i == j else 0 for i in range(k)) for j in range(k)]
    def close(gens):
        Bs = hnf(gens, k)
        if any(in_lattice(u, Bs, k) for u in units): return None, None
        return frozenset(e for e in cube if in_lattice(e, Bs, k)), Bs
    R0, B0 = close([])
    seen = {R0: B0}; frontier = [R0]
    while frontier:
        nxt = []
        for R in frontier:
            Bs = seen[R]
            for v in cube:
                if v in R: continue
                R2, B2 = close(list(Bs) + [v])
                if R2 is None or R2 in seen: continue
                seen[R2] = B2; nxt.append(R2)
        frontier = nxt
    out = []
    for R, Bs in seen.items():
        pi = sum(Fr(1, 2 ** (k + sum(1 for x in e if x))) for e in R)
        out.append((pi, len(Bs), tuple(elem_divisors(Bs, k)), R))
    return out

_CLS = {}
def classes(k):
    """aggregated classes: list of (pi, m = k - rho, n = #nontrivial d_i, D = prod d_i, multiplicity)"""
    if k in _CLS: return _CLS[k]
    agg = {}
    for (pi, rho, ed, R) in closed_sets(k):
        D = 1
        for d in ed: D *= d
        key = (pi, k - rho, len(ed), D)
        agg[key] = agg.get(key, 0) + 1
    out = sorted([(pi, m, n, D, c) for (pi, m, n, D), c in agg.items()], key=lambda x: (-x[1], x[0]))
    assert out[0][:4] == (Fr(1, 2 ** k), k, 0, 1) and out[0][4] == 1      # the generic class R = {0}
    _CLS[k] = out
    return out

def qsum_cls(m, n, D, Q):
    """sec 3.4: upper bound for sum over 3 <= q <= Q of phi(q)^m prod_i |G_q[d_i]|  (arb)"""
    Q = A(Q)
    if n == 0:
        return Q if m == 0 else Q ** (m + 1) / (m + 1)       # sum (q-1)^m <= Q^(m+1)/(m+1)
    return arb(2) ** n * Q ** (m + 1) * (1 + Q.log()) ** (D - 1)   # |G[d]| <= 2 d^omega(q), D^omega <= tau_D, sum tau_D <= Q(1+log Q)^(D-1)

def union_bound(Q, k, d, al, r, s):
    """sec 3.5: certified upper bound for sum over 3 <= q <= Q of Pr(S_{r+s}(q) != (Z/qZ)^x)"""
    d = A(d); al = A(al); t = al / (1 - al); beta = d + (1 - d) * al
    tot = arb(0)
    for (pi, m, n, D, c) in classes(k):
        pa = A(pi); rho = pa + d * (1 - pa)
        tot += c * qsum_cls(m, n, D, Q) * rho ** r
    tot = tot / t ** k
    tot += al * A(Q) ** 2 / 2 * beta ** s
    return tot

def best_rs(lQ, k, d, al):
    """float search for the least r + s (parameter choice only; certified afterwards by union_bound)"""
    lQ = float(lQ); d = float(d); al = float(al); t = al / (1 - al); beta = d + (1 - d) * al
    def lq(m, n, D):
        if n == 0: return lQ if m == 0 else (m + 1) * lQ - math.log(m + 1)
        return n * math.log(2) + (m + 1) * lQ + (D - 1) * math.log(1 + lQ)
    terms = [(math.log(c) + lq(m, n, D), float(pi) + d * (1 - float(pi))) for (pi, m, n, D, c) in classes(k)]
    def F(r):
        xs = [a + r * math.log(p) for (a, p) in terms]; mx = max(xs)
        return mx + math.log(sum(math.exp(x - mx) for x in xs)) - k * math.log(t)
    best = None
    for f1 in [0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.99]:
        lo, hi = 0, 10 ** 7
        while lo < hi:
            mid = (lo + hi) // 2
            if F(mid) <= math.log(f1 * 0.999): hi = mid
            else: lo = mid + 1
        s = max(0, math.ceil((math.log(al) + 2 * lQ - math.log(2) - math.log((1 - f1) * 0.999)) / math.log(1 / beta)))
        if best is None or lo + s < best[0] + best[1]: best = (lo, s)
    return best

def table_row(lQ_lo, Qexact, k, d, alphas):
    """least certified N = r + s for one Q (Qexact = arb upper bound for Q) over a grid of alpha"""
    best = None
    for al in alphas:
        r, s = best_rs(lQ_lo, k, d, al)
        while r > 0 and union_bound(Qexact, k, d, al, r - 1, s) < 1: r -= 1
        while s > 0 and union_bound(Qexact, k, d, al, r, s - 1) < 1: s -= 1
        U = union_bound(Qexact, k, d, al, r, s)
        if U < 1 and (best is None or r + s < best[0] + best[1]):
            best = (r, s, al, U)
    r, s, al, U = best
    pos(1 - U, "union bound < 1")
    return r, s, al, U

# ---------------------------------------------------------------- sec 4: Theorem 1 for every Q >= Q0 (linear bound)
def sup_poly_exp(f, eps, Lmin):
    """sup over L >= Lmin of (1+L)^f e^(-eps L), f >= 0, eps > 0 (arb); maximum at 1 + L = f/eps"""
    f = A(f)
    if f == 0: return (-eps * Lmin).exp()
    Ls = f / eps - 1
    if Ls <= Lmin: return (1 + Lmin) ** f * (-eps * Lmin).exp()
    if Ls >= Lmin: return (f / eps) ** f * (eps - f).exp()
    return amax((1 + Lmin) ** f * (-eps * Lmin).exp(), (f / eps) ** f * (eps - f).exp())

def linear_theorem(k, d, al, a0, b0, Lmin=arb(3).log()):
    """sec 4: r = ceil(a1 log Q + a0), s = ceil(b1 log Q + b0), a1 = (k+1)/log(1/rho_gen), b1 = 2/log(1/beta).
    Returns (a1, b1, certified sup over log Q >= Lmin of the union bound)."""
    d = A(d); al = A(al); t = al / (1 - al); beta = d + (1 - d) * al
    cl = classes(k)
    pg = A(cl[0][0]); rho_g = pg + d * (1 - pg)
    a1 = (k + 1) / (1 / rho_g).log()
    b1 = 2 / (1 / beta).log()
    tot = arb(0)
    for idx, (pi, m, n, D, c) in enumerate(cl):
        pa = A(pi); rho = pa + d * (1 - pa)
        if n == 0:
            e, f, coef = (1, 0, arb(1)) if m == 0 else (m + 1, 0, arb(1) / (m + 1))
        else:
            e, f, coef = m + 1, D - 1, arb(2) ** n
        if idx == 0:
            mfac = arb(1)                         # generic class: Q^(k+1) rho_g^(a1 log Q) = 1 exactly
        else:
            eps = a1 * (1 / rho).log() - e        # Q^e rho^(a1 log Q) = Q^(-eps)
            pos(eps, "class exponent eps > 0 (pi=%s, m=%d)" % (pi, m))
            mfac = sup_poly_exp(f, eps, Lmin)
        tot += c * coef * mfac * rho ** A(a0)       # rho^r <= rho^(a1 log Q + a0)
    tot = tot / t ** k
    tot += al / 2 * beta ** A(b0)                  # alpha Q^2/2 beta^s <= alpha/2 beta^b0
    return a1, b1, tot

def cover_section(d0=Fr(1, 1000), plan=None):
    for k in (1, 2, 3, 4):
        cl = classes(k)
        print("[cover] k=%d: %d closed relation sets in %d classes (pi, k-rank, #d_i>1, prod d_i): max pi = %s"
              % (k, sum(c for *_, c in cl), len(cl), max(x[0] for x in cl)), flush=True)
        if k >= 2:   # DERIVATION sec 9: least slack log2(1/pi_R) - (k - rho + 1) k/(k + 1) over non-generic classes (delta = 0)
            sl = min(((1 / A(pi)).log() / LOG2 - A(Fr(k * (m + 1), k + 1)), pi, m) for (pi, m, n, D, c) in cl[1:])
            print("[cover]   generic-dominance slack (bits, delta = 0): least %s at pi = %s, k - rho = %d" % (sl[0].str(4), sl[1], sl[2]))
    res = {}
    if plan is None:
        plan = [(1, Fr(1, 10 ** 4), arb(3).log()), (2, Fr(1, 10 ** 4), arb(3).log()), (3, Fr(1, 10 ** 4), arb(3).log()),
                (4, Fr(1, 10 ** 4), arb(3).log()), (4, Fr(1, 10 ** 4), arb(10) ** 2 * arb(10).log())]
    print("[cover] delta = %s:" % d0)
    for (k, al, Lmin) in plan:
        a0 = None
        for a0f in range(0, 2000):
            a1, b1, tot = linear_theorem(k, d0, al, Fr(a0f), Fr(0), Lmin)
            if tot - A(al) / 2 < A(Fr(9, 10)):
                a0 = Fr(a0f); break
        beta = A(d0) + (1 - A(d0)) * A(al)
        b0 = Fr(0)
        while not (A(al) / 2 * beta ** A(b0) < A(Fr(9, 100))): b0 += 1
        a1, b1, tot = linear_theorem(k, d0, al, a0, b0, Lmin)
        pos(1 - tot, "linear theorem k=%d" % k)
        C1 = a1 + b1; C2 = A(a0) + A(b0) + 2
        Q0 = "3" if Lmin < 2 else "10^%d" % round(float((Lmin / arb(10).log()).mid()))
        print("[cover] k=%d, alpha=%s, all Q >= %s: r = ceil(%s log Q + %s), s = ceil(%s log Q + %s); sup union bound <= %s < 1"
              % (k, al, Q0, a1.str(6), a0, b1.str(6), b0, tot.str(4)))
        print("[cover]   => |P| <= %s log Q + %s = (1/log 2 + %s) log Q + %s  (rounded up; slope %s log2 Q; hole share b1 <= %s)"
              % (ceil_str(C1, 4), C2.str(4), ceil_str(C1 - 1 / LOG2, 4), C2.str(4), ceil_str(C1 * LOG2, 4), ceil_str(b1, 6)), flush=True)
        print("[cover]   upper ends: a1 <= %s, b1 <= %s, C = a1 + b1 <= %s" % (ceil_str(a1, 6), ceil_str(b1, 6), ceil_str(C1, 6)))
        res[(k, Q0)] = (C1, C2)
    return res

def min_cp(L, N, d0):
    """least cp (two significant figures, rounded up) with certified delta_direct(cp, L, N) <= d0"""
    lo, hi = 1.0, 1e9
    for _ in range(80):
        mid = (lo + hi) / 2
        dd, _ = delta_direct(Fr(mid).limit_denominator(10 ** 6), L, N)
        if dd < A(d0): hi = mid
        else: lo = mid
    e = 10 ** (len(str(int(hi))) - 2)
    cp = Fr((int(hi) // e + 1) * e)
    dd, _ = delta_direct(cp, L, N)
    pos(A(d0) - dd, "min_cp")
    return cp

def table_section(d0=Fr(1, 1000), rows=(3, 6, 12, 30, 100, 1000), ks=(1, 2, 3, 4)):
    alphas = [Fr(1, 10 ** j) for j in range(1, 9)] + [Fr(3, 10 ** j) for j in range(2, 9)]
    for e10 in rows:
        Qx = arb(10) ** e10
        lQ = float((Qx.log()).lower())
        L = Qx.log()
        line = []; bestN = None
        for k in ks:
            r, s, al, U = table_row(lQ, Qx, k, d0, alphas)
            if bestN is None or r + s < bestN: bestN = r + s
            line.append("k=%d: %d (r=%d, s=%d, alpha=%s)" % (k, r + s, r, s, al))
        Nmax = max(int(x.split(": ")[1].split(" ")[0]) for x in line)
        cpr = min_cp(L, Nmax, d0)
        dd, X = delta_direct(cpr, L, Nmax)
        pos(A(d0) - dd, "delta at Q = 10^%d" % e10)
        lgy = (L * A(cpr)).log() / b / arb(10).log()
        print("[table] Q = 10^%d (log2 Q = %.1f, delta = %s): N = %d;  %s;  c' = %s, max P <= (c' log Q)^(1/(1-B)) <= 10^%s"
              % (e10, e10 * math.log2(10), d0, bestN, "; ".join(line), cpr, ceil_str(lgy, 2)), flush=True)

# ---------------------------------------------------------------- sec 5: exact closure for small Q
def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def need(q, seq):
    """least m such that the subset products of {p in seq[:m] : p does not divide q} are all of (Z/qZ)^x (exact);
    returns (m, phi(q))"""
    units = [x for x in range(1, q) if math.gcd(x, q) == 1]
    ph = len(units)
    if ph == 1: return 0, 1
    S = bytearray(q); S[1] = 1; mem = [1]
    for i, p in enumerate(seq):
        if q % p == 0: continue
        pm = p % q
        new = []
        for x in mem:
            z = (x * pm) % q
            if not S[z]:
                S[z] = 1; new.append(z)
        mem += new
        if len(mem) == ph: return i + 1, ph
    return None, ph

def small_section(Qs=2000):
    t0 = time.time()
    pr = primes_upto(20000)
    first = pr[:60]
    above = [p for p in pr if p > Qs][:60]
    res = {}
    for name, seq in (("first primes", first), ("primes > %d" % Qs, above)):
        worst = (0, None); hard = []
        for q in range(2, Qs + 1):
            m, ph = need(q, seq)
            if m is None: raise SystemExit("CERT FAIL: %s do not cover q = %d" % (name, q))
            if m > worst[0]: worst = (m, q)
            hard.append((m, q))
        hard.sort(reverse=True)
        m = worst[0]
        P = seq[:m]
        res[name] = P
        print("[small] %s: the %d primes %s ... %s cover every 2 <= q <= %d (p | q excluded); hardest q: %s"
              % (name, m, P[:6], P[-3:], Qs, ", ".join("%d(m=%d)" % (q, mm) for (mm, q) in hard[:6])), flush=True)
        print("[small]   P = %s" % P)
    # adversarial moduli, recomputed with the printed sets
    adv = [3, 4, 6, 8, 30, 210, 1024, 1155, 1365, 1001, 1980, 1999, 1997, 1993, 2000]
    for name, P in res.items():
        oks = []
        for q in adv:
            m, ph = need(q, P)
            oks.append("%d:%s" % (q, "ok" if m is not None else "FAIL"))
            if m is None: raise SystemExit("CERT FAIL adversarial")
        print("[small]   %s adversarial check: %s" % (name, " ".join(oks)))
    # least m for smaller Q (same two sequences)
    for Qs2 in (100, 300, 1000):
        ab2 = [p for p in pr if p > Qs2][:60]
        m1 = max(need(q, first)[0] for q in range(2, Qs2 + 1))
        m2 = max(need(q, ab2)[0] for q in range(2, Qs2 + 1))
        print("[small] Q = %d: the first %d primes (up to %d) cover every 2 <= q <= Q; so do the %d least primes > Q (up to %d)"
              % (Qs2, m1, first[m1 - 1], m2, ab2[m2 - 1]))
    print("[small] done in %.1f s" % (time.time() - t0))
    return res

if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if mode in ("core", "all"): core()
    if mode in ("delta", "all"):
        delta_section()
        delta_section(cp=Fr(820000), d0=Fr(1, 10 ** 6))
    if mode in ("cover", "all"):
        cover_section(); table_section()
        cover_section(d0=Fr(1, 10 ** 6), plan=[(3, Fr(1, 10 ** 7), arb(3).log()), (4, Fr(1, 10 ** 7), arb(3).log())])
        table_section(d0=Fr(1, 10 ** 6), rows=(12, 100, 1000), ks=(3, 4))
    if mode in ("small", "all"): small_section()
    print("[cert] ALL CERTIFIED (%s) in %.1f s" % (mode, time.time() - t0))
