"""Exact evaluation of c(0,0), c_1(u), c_2(u) of Bui-Hall-Subira Jorge (arXiv 2511.14415, Thm 1.1 l.75-100,
c_1 l.~617-640, c_2 l.~648-663) for r in N, polynomial P, rational theta.

Reduction (see DERIVATION.md):
 * slack w_k = argument of the k-th P factor; constraint k holds one (x,z) pair:
   C1: (x3,z3,w1)  C2: (x4,z4,w2)  C3: (x1,z1,w3)  C4: (x2,z2,w4)
   R1 = 1-t1-t2, R2 = 1-t3-t4, R3 = 1-t1-t3, R4 = 1-t2-t4.
 * the integrand (without P and weights) depends on x,z,t only through
   X12=x1+x2, Z12=z1+z2, X34=x3+x4, Z34=z3+z4, w1..w4.
 * v1..v4 integrated exactly (fmpq_mpoly), then a moment functional:
   Dirichlet integrals on each 2-simplex, then the 4-dim t-integral G (closed form, two-case split).
"""
from flint import fmpq, fmpq_mpoly_ctx, fmpq_mat
from math import comb, factorial
from functools import lru_cache

NAMES = ('v1', 'v2', 'v3', 'v4', 'X12', 'Z12', 'X34', 'Z34', 'w1', 'w2', 'w3', 'w4', 'u')


def poly_parts(theta):
    """theta: fmpq. Returns dict name -> fmpq_mpoly in (X12..w4,u) after integrating v1..v4 over [0,1]^4."""
    ctx = fmpq_mpoly_ctx.get(NAMES, 'lex')
    v1, v2, v3, v4, X12, Z12, X34, Z34, w1, w2, w3, w4, u = ctx.gens()
    th = fmpq(theta)
    Sx = X12 + X34
    Sz = Z12 + Z34
    a = 1 - th * Sx                      # 1 - theta*sum x_j
    b = 1 - th * Sz                      # 1 - theta*sum z_j
    D12 = v1 * a - v2 * b + th * (X12 - Z12)   # v1-v2+theta(x1+x2-z1-z2-v1 Sx+v2 Sz)
    D34 = v1 * a - v2 * b + th * (X34 - Z34)   # v1-v2+theta(x3+x4-z3-z4-v1 Sx+v2 Sz)
    # t1+t2+x1+x2+x3+z3 = (1-w1)+X12 etc.
    F1 = u - th * (1 - w1 + X12) - v1 * a + v3 * D12
    F2 = u - th * (1 - w2 + Z12) - v2 * b - v3 * D12
    F3 = u - th * (1 - w3 + X34) - v1 * a + v4 * D34
    F4 = u - th * (1 - w4 + Z34) - v2 * b - v4 * D34
    base = a * b * D12 * D34

    def vint(p):
        for i in range(4):
            q = p.integral(i)
            p = q.subs({NAMES[i]: 1}) - q.subs({NAMES[i]: 0})
        return p
    return {'c': vint(base), 'c1': vint(base * F1 * F2 * F3 * F4), 'c2': vint(base * F1 * F3)}


def split_terms(p):
    """fmpq_mpoly -> dict {(gA, gB, k): coeff}, gA=(a3,a4,b1,b2), gB=(a1,a2,b3,b4), k = power of u."""
    out = {}
    for mon, cf in zip(p.monoms(), p.coeffs()):
        assert mon[0] == mon[1] == mon[2] == mon[3] == 0
        a1, a2, a3, a4, b1, b2, b3, b4, k = mon[4:]
        out[((a3, a4, b1, b2), (a1, a2, b3, b4), k)] = fmpq(cf)
    return out


# ---------- exact special integrals ----------
@lru_cache(maxsize=None)
def fact(n):
    return factorial(n)


def Beta(p, q):
    """B(p,q) for positive integers, exact."""
    return fmpq(fact(p - 1) * fact(q - 1), fact(p + q - 1))


def dirichlet(p, q, m, r):
    """int over {x,z>=0, x+z<=R} of x^(p+r-1) z^(q+r-1) (R-x-z)^m = coefficient * R^(p+q+m+2r)."""
    return fmpq(fact(p + r - 1) * fact(q + r - 1) * fact(m), fact(p + q + m + 2 * r))


def Hcase(n1, n2, n3, n4, e1, e2, e3, e4):
    """Case R1 <= R3 of int t^e R^n over the t-polytope (coordinates S,Delta,R1,R4)."""
    s = fmpq(0)
    for j in range(n3 + 1):
        for k in range(n2 + 1):
            P_ = n1 + j + e1 + 1 + n4 + k + e4 + 1
            Q_ = n3 - j + n2 - k
            s += (comb(n3, j) * comb(n2, k) * Beta(n1 + j + 1, e1 + 1) * Beta(n4 + k + 1, e4 + 1)
                  * Beta(Q_ + 1, e3 + 1) * Beta(P_ + 1, e2 + e3 + Q_ + 2))
    return s


@lru_cache(maxsize=None)
def G(n1, n2, n3, n4, e):
    """int over {t>=0, R_k>=0} of (t1t2t3t4)^e R1^n1 R2^n2 R3^n3 R4^n4 dt (exact)."""
    return Hcase(n1, n2, n3, n4, e, e, e, e) + Hcase(n3, n4, n1, n2, e, e, e, e)


def G_exact_table(r, n0, nl, return_H=False):
    """Exact dict {(n1,n2,n3,n4): G} for n_i in [n0, n0+nl-1], via the Hankel form of Hcase:
    H(n1,n2,n3,n4) = sum_{j'<=n3, k'<=n2} u_{n1,n3}[j'] gamma_N(j'+k') u_{n4,n2}[k'],
    u_{a,c}[j'] = C(c,j') B(a+c-j'+1, e+1), gamma_N(s) = B(s+1,e+1) B(N-s+2e+3, 2e+s+2), N = n1+n2+n3+n4;
    G(n1,n2,n3,n4) = H(n1,n2,n3,n4) + H(n3,n4,n1,n2)  (same terms as Hcase/G above, regrouped)."""
    e = r * r - 1
    nmax = n0 + nl - 1
    J = nmax + 1
    groups = {}
    for a in range(n0, nmax + 1):
        for c in range(n0, nmax + 1):
            groups.setdefault(a + c, []).append((a, c))
    U = {A: fmpq_mat([[fmpq(comb(c, jp)) * Beta(a + c - jp + 1, e + 1) if jp <= c else fmpq(0) for jp in range(J)]
                      for (a, c) in prs]) for A, prs in groups.items()}
    Ut = {A: m.transpose() for A, m in U.items()}
    gam = {}
    H = {}
    for A, prsA in groups.items():
        for B_, prsB in groups.items():
            N = A + B_
            if N not in gam:
                # entries with s > N - 2*n0 never meet nonzero u-entries (j' <= n3, k' <= n2); set them to 0
                gv = [Beta(s + 1, e + 1) * Beta(N - s + 2 * e + 3, 2 * e + s + 2) if s <= N - 2 * n0 else fmpq(0)
                      for s in range(2 * J - 1)]
                gam[N] = fmpq_mat(J, J, [gv[jp + kp] for jp in range(J) for kp in range(J)])
            M = U[A] * gam[N] * Ut[B_]
            for i, (n1, n3) in enumerate(prsA):
                for j, (n4, n2) in enumerate(prsB):
                    H[(n1, n2, n3, n4)] = M[i, j]
    if return_H:
        return H
    return {k: H[k] + H[(k[2], k[3], k[0], k[1])] for k in H}


# ---------- assembling ----------
def phi(p, q, b, r, P):
    """dict power->coef: int over 2-simplex of size R of x^(p+r-1) z^(q+r-1) w^b P(w)."""
    out = {}
    for c, pc in enumerate(P):
        if pc == 0:
            continue
        n = p + q + b + c + 2 * r
        out[n] = out.get(n, fmpq(0)) + fmpq(pc) * dirichlet(p, q, b + c, r)
    return out


def V(g, r, P):
    """group monomial (a,a',b,b') -> dict {(n_first, n_second): coef}."""
    aa, ab, b1, b2 = g
    out = {}
    for i in range(aa + 1):
        for k in range(ab + 1):
            cf = comb(aa, i) * comb(ab, k)
            f1 = phi(i, k, b1, r, P)
            f2 = phi(aa - i, ab - k, b2, r, P)
            for n1, c1 in f1.items():
                for n2, c2 in f2.items():
                    key = (n1, n2)
                    out[key] = out.get(key, fmpq(0)) + cf * c1 * c2
    return out


def evaluate(parts_terms, r, P, uval):
    """Exact values {name: fmpq} for P (list of fmpq coefficients, P(x)=sum P[c] x^c) and u=uval (fmpq)."""
    up = evaluate_upoly(parts_terms, r, P)
    uv = fmpq(uval)
    return {name: sum((cf * uv ** k for k, cf in enumerate(cs)), fmpq(0)) for name, cs in up.items()}


def evaluate_upoly(parts_terms, r, P, Gtab=None):
    """Exact {name: [coef of u^0, u^1, ...]} for P (list of fmpq coefficients, P(x)=sum P[c] x^c).
    Gtab: optional exact dict from G_exact_table (same values as G, faster)."""
    e = r * r - 1
    gAs = sorted({k[0] for t in parts_terms.values() for k in t})
    gBs = sorted({k[1] for t in parts_terms.values() for k in t})
    VA = {g: V(g, r, P) for g in gAs}
    VB = {g: V(g, r, P) for g in gBs}
    colsA = sorted({n for d in VA.values() for n in d})
    colsB = sorted({n for d in VB.values() for n in d})
    iA = {n: i for i, n in enumerate(colsA)}
    iB = {n: i for i, n in enumerate(colsB)}
    MA = fmpq_mat(len(gAs), len(colsA))
    for i, g in enumerate(gAs):
        for n, cf in VA[g].items():
            MA[i, iA[n]] = cf
    MB = fmpq_mat(len(colsB), len(gBs))
    for j, g in enumerate(gBs):
        for n, cf in VB[g].items():
            MB[iB[n], j] = cf
    GM = fmpq_mat(len(colsA), len(colsB))
    for i, (n1, n2) in enumerate(colsA):
        for j, (n3, n4) in enumerate(colsB):
            GM[i, j] = Gtab[(n1, n2, n3, n4)] if Gtab is not None else G(n1, n2, n3, n4, e)
    M = MA * GM * MB
    rA = {g: i for i, g in enumerate(gAs)}
    rB = {g: j for j, g in enumerate(gBs)}
    out = {}
    for name, terms in parts_terms.items():
        deg = max(k for (_, _, k) in terms)
        cs = [fmpq(0)] * (deg + 1)
        for (gA, gB, k), cf in terms.items():
            cs[k] += cf * M[rA[gA], rB[gB]]
        out[name] = cs
    return out
