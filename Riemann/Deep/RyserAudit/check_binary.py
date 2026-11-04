#!/usr/bin/env python3
"""Check of Section 5 ('binary coefficient identity') of the preprint.

For h in {+-1}^{C_4 x P}: c = h(Z=1)/2, d = h(Z=-1)/2, g = h(Z=i)/2, b = (H1+H3)/2, H2 + i H3 = (1+i)J + 2U.
Exact group-ring identities: d - c = -2b,  g - c = (1+i)(i b - J - (1-i) U).
Consequence used by the paper: at the prime (1+i) over 2 (unique here),
   lambda_2(R) + lambda_{1+i}(W) = u^2 / c_empty  (mod (1+i))   whenever all c_S, d_S, g_S are 2-units,
R = prod_S (d_S/c_S)^{eps_S},  W = prod_S (g_S/c_S)^{eps_S}.  So R = W = 1 is impossible for ANY sign array;
the proof then needs R, W of odd order, which comes from the character-comparison lemma (P cyclic).
Tests: (1) order 4 example; (2) random sign arrays on C_4 x C_{u^2}, u = 3, 5, 9 (|I| = 1, prime (1+i) unique),
    (2') same for n = 36, 100, 196, 324, 900 with the full alternating product over all S subset I, computed in the
    residue field of one prime over 2 (also when 2 splits: u = 7, 15);
(3) near miss: a Menon-Hadamard difference set in C_4 x C_3 x C_3 (perfect binary array 12 x 3):
    all norms hold, but P is not cyclic; we show exactly which step of the argument fails there.
"""
import random, sys, warnings
warnings.filterwarnings('ignore', category=DeprecationWarning)   # sympy modular-int sort notice
from itertools import product

def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gadd(a, b): return (a[0]+b[0], a[1]+b[1])
def gsub(a, b): return (a[0]-b[0], a[1]-b[1])
Z0 = (0, 0)

# ---- elements of Z[i][zeta_p] stored mod x^p - 1 ----
def pmul(a, b, p):
    c = [Z0]*p
    for i, x in enumerate(a):
        if x == Z0: continue
        for j, y in enumerate(b):
            if y == Z0: continue
            c[(i+j) % p] = gadd(c[(i+j) % p], gmul(x, y))
    return c
def reduce_phi(a, p):            # power-basis coordinates 1..zeta^{p-2}
    t = a[p-1]; return [gsub(x, t) for x in a[:p-1]]
def is_zero(a, p): return all(x == Z0 for x in reduce_phi(a, p))

# ---- residue field F_2[x]/Phi_p (Phi_p irreducible mod 2 for p = 3, 5, 11, 13) ----
def f2_res(a, p):
    """residue mod (1+i) of a in Z[i][zeta_p] (a mod x^p-1) as bitmask poly reduced mod Phi_p."""
    bits = 0
    for k, (x, y) in enumerate(a):
        if (x + y) % 2: bits ^= (1 << k)
    return f2_mod(bits, p)
def f2_mod(a, p):
    m = (1 << p) - 1             # Phi_p mod 2 = 1 + x + ... + x^{p-1}
    for k in range(a.bit_length()-1, p-2, -1):
        if (a >> k) & 1: a ^= (m << (k-(p-1)))
    return a
def f2_mul(a, b, p):
    r = 0
    while b:
        if b & 1: r ^= a
        a <<= 1; b >>= 1
    return f2_mod(r, p)
def f2_inv(a, p):
    q = 2**(p-1)                 # field size
    r, e, base = 1, q-2, a
    while e:
        if e & 1: r = f2_mul(r, base, p)
        base = f2_mul(base, base, p); e >>= 1
    return r

def div2_res(a, p):
    """a in 2 Z[i][zeta_p] ?  then residue of a/2."""
    red = reduce_phi(a, p)
    if any(x % 2 or y % 2 for x, y in red): return None
    half = [(x//2, y//2) for x, y in red] + [Z0]
    return f2_res(half, p)
def divpi_res(a, p):
    """a in (1+i) Z[i][zeta_p] ? then residue of a/(1+i) = a(1-i)/2."""
    return div2_res([gmul(x, (1, -1)) for x in a], p)

def lam_identity(cS, dS, gS, p, u2_over_c_res):
    """S ranges over {empty, {p}} (|I| = 1).  returns (ok, R==1, W==1, details)"""
    c0, c1 = cS; d0, d1 = dS; g0, g1 = gS
    for v in (c0, c1, d0, d1, g0, g1):
        if f2_res(v, p) == 0: return None          # not all 2-units
    NR, DR = pmul(d0, c1, p), pmul(c0, d1, p)     # R = (d0/c0)/(d1/c1)
    NW, DW = pmul(g0, c1, p), pmul(c0, g1, p)
    MR = [gsub(x, y) for x, y in zip(NR, DR)]; MW = [gsub(x, y) for x, y in zip(NW, DW)]
    lr = div2_res(MR, p); lw = divpi_res(MW, p)
    if lr is None or lw is None: return ('membership-fail', lr, lw)
    lamR = f2_mul(lr, f2_inv(f2_res(DR, p), p), p)
    lamW = f2_mul(lw, f2_inv(f2_res(DW, p), p), p)
    return (lamR ^ lamW) == u2_over_c_res, is_zero(MR, p), is_zero(MW, p), (lamR, lamW)

# ---- group ring of C_4 x P, P = C_N (cyclic) ----
def split_h(hvec, N):
    """hvec indexed by m in Z/4N; C_4N = C_4 x C_N via CRT (N odd). returns H_j[z] for j in 0..3, z in Z/N."""
    H = [[0]*N for _ in range(4)]
    for m, s in enumerate(hvec): H[m % 4][m % N] = s
    return H

def evals(H, N, p):
    """c,d,g,b,U,J as lists over Z/N with Gaussian entries; check exact identities; values at chi_empty, chi_p."""
    c = [((H[0][z]+H[1][z]+H[2][z]+H[3][z])//2, 0) for z in range(N)]
    d = [((H[0][z]-H[1][z]+H[2][z]-H[3][z])//2, 0) for z in range(N)]
    g = [((H[0][z]-H[2][z])//2, (H[1][z]-H[3][z])//2) for z in range(N)]
    b = [((H[1][z]+H[3][z])//2, 0) for z in range(N)]
    U = [((H[2][z]-1)//2, (H[3][z]-1)//2) for z in range(N)]
    for z in range(N):
        assert all((H[j][z]) in (1, -1) for j in range(4))
        assert gsub(d[z], c[z]) == (-2*b[z][0], 0)
        rhs = gsub(gsub(gmul((0, 1), b[z]), (1, 0)), gmul((1, -1), U[z]))     # i b - J - (1-i) U
        assert gsub(g[z], c[z]) == gmul((1, 1), rhs)
        assert gadd(gmul((1, 1), (1, 0)), gmul((2, 0), U[z])) == (H[2][z], H[3][z])
    def val(x, S):
        out = [Z0]*p
        for z in range(N):
            k = (z % p) if S else 0                           # chi_p(X^z) = rho^(z mod p)
            out[k] = gadd(out[k], x[z])
        return out
    return [(val(x, False), val(x, True)) for x in (c, d, g)], c

def periodic_ok(hvec):
    n = len(hvec)
    return all(sum(hvec[j]*hvec[(j+t) % n] for j in range(n)) == 0 for t in range(1, n))

def test_random(u, p, trials, rng):
    N = u*u; n = 4*N
    hold = total = rw = memb = 0
    for _ in range(trials):
        hvec = [rng.choice((1, -1)) for _ in range(n)]
        H = split_h(hvec, N)
        (cS, dS, gS), c = evals(H, N, p)
        c_empty = sum(x[0] for x in c)
        res = lam_identity(cS, dS, gS, p, 1)     # u^2/c_empty has residue 1 when c_empty odd (unit)
        if res is None: continue
        if res[0] == 'membership-fail': memb += 1; continue
        total += 1
        if res[0]: hold += 1
        if res[1] and res[2]: rw += 1
    print(f"(2) random sign arrays on C_4 x C_{N} (u={u}): {total} with all values 2-units; "
          f"identity lambda_2(R)+lambda_(1+i)(W) = u^2/c holds in {hold}/{total}; R=W=1 in {rw}; membership failures {memb}")
    return hold == total and memb == 0

def order4():
    hvec = [1, 1, 1, -1]
    H = [[hvec[0]], [hvec[1]], [hvec[2]], [hvec[3]]]
    c = (sum(hvec)) // 2; d = (hvec[0]-hvec[1]+hvec[2]-hvec[3]) // 2
    g = ((hvec[0]-hvec[2])//2, (hvec[1]-hvec[3])//2)
    # R = d/c, W = g/c in Z[i]
    print(f"(1) order 4, row (1,1,1,-1): c={c}, d={d}, g={g[0]}+{g[1]}i  -> R = {d}/{c}, W = g/c = i (order 4, even): "
          f"paper's Remark reproduced; hh*=4: {periodic_ok(hvec)}")

def pba_near_miss():
    # C_4 x C_3 x C_3; subgroups of C_3^2 of order 3: generated by (1,0),(0,1),(1,1),(1,2)
    gens = [(1, 0), (0, 1), (1, 1), (1, 2)]
    E = [(a, b) for a in range(3) for b in range(3)]
    subs = [set(((k*x) % 3, (k*y) % 3) for k in range(3)) for (x, y) in gens]
    D = set()
    for j in range(4):
        part = subs[j] if j < 3 else set(E) - subs[3]
        for e in part: D.add((j, e))
    h = {(j, e): (-1 if (j, e) in D else 1) for j in range(4) for e in E}
    # periodic autocorrelation over C_4 x C_3 x C_3
    ok = True
    for t in [(j, e) for j in range(4) for e in E]:
        if t == (0, (0, 0)): continue
        s = sum(h[(j, e)] * h[((j+t[0]) % 4, ((e[0]+t[1][0]) % 3, (e[1]+t[1][1]) % 3))] for j in range(4) for e in E)
        ok &= (s == 0)
    print(f"(3) near miss: D in C_4 x C_3^2, |D| = {len(D)}, all 35 periodic autocorrelations zero: {ok}"
          f"  (a 12x3 perfect binary array; P = C_3^2 is NOT cyclic)")
    H = [{e: h[(j, e)] for e in E} for j in range(4)]
    c = {e: ((H[0][e]+H[1][e]+H[2][e]+H[3][e])//2, 0) for e in E}
    d = {e: ((H[0][e]-H[1][e]+H[2][e]-H[3][e])//2, 0) for e in E}
    g = {e: ((H[0][e]-H[2][e])//2, (H[1][e]-H[3][e])//2) for e in E}
    p = 3
    def val(x, chi):
        out = [Z0]*p
        for e in E:
            k = (chi[0]*e[0] + chi[1]*e[1]) % 3
            out[k] = gadd(out[k], x[e])
        return out
    def norm_ok(x):   # x x^* = 9 in Z[i][C_3^2]
        for t in E:
            s = Z0
            for e in E:
                y = x[((e[0]-t[0]) % 3, (e[1]-t[1]) % 3)]
                s = gadd(s, gmul(x[e], (y[0], -y[1])))
            if s != ((9, 0) if t == (0, 0) else Z0): return False
        return True
    print(f"    c c* = d d* = g g* = 9: {norm_ok(c) and norm_ok(d) and norm_ok(g)}")
    out = []
    for chi in [(1, 0), (0, 1), (1, 1), (1, 2)]:
        cS = (val(c, (0, 0)), val(c, chi)); dS = (val(d, (0, 0)), val(d, chi)); gS = (val(g, (0, 0)), val(g, chi))
        def rat(x):   # x(chi)/x(1) when x(chi) is rational
            v = reduce_phi(x[1], 3)
            return (v[0], x[0][0]) if v[1] == Z0 else None
        res = lam_identity(cS, dS, gS, 3, 1)
        out.append((chi, rat(cS), rat(dS), rat(gS), res))
    for chi, rc, rd, rg, res in out:
        print(f"    chi with chi(e) = rho^({chi[0]}a+{chi[1]}b): c(chi)/c(1) = {rc[0][0]}/{rc[1][0]}, "
              f"d: {rd[0][0]}/{rd[1][0] if rd else '?'}, g: {rg[0] if rg else 'irr'}/{rg[1] if rg else ''}; "
              f"identity holds: {res[0]}, R=1: {res[1]}, W=1: {res[2]}")
    return ok

# ---------------- general version: any odd u, several primes over 2 ----------------
from sympy import cyclotomic_poly, factor_list, Poly, symbols
warnings.filterwarnings('ignore', category=DeprecationWarning)   # after sympy's own filter
_x = symbols('x')

class Ctx:
    def __init__(self, m):
        self.m = m
        self.phi = [int(c) for c in reversed(Poly(cyclotomic_poly(m, _x), _x).all_coeffs())]  # low->high
        self.deg = len(self.phi) - 1
        facs = factor_list(Poly(cyclotomic_poly(m, _x), _x, modulus=2))[1]
        polys = [sum((int(c) % 2) << k for k, c in enumerate(reversed(f.all_coeffs()))) for f, _ in facs]
        self.f1 = polys[0]; self.others = polys[1:]; self.fdeg = self.f1.bit_length() - 1
        # y = prod_{other factors} f_j(zeta)^2, lifted with 0/1 coefficients
        y = [Z0]*m; y[0] = (1, 0)
        for f in self.others:
            fl = [((f >> k) & 1, 0) for k in range(m)]
            y = self.mul(self.mul(y, fl), fl)
        self.y = y
        self.yres = self.res(y)
        assert self.yres != 0
    def mul(self, a, b):
        m = self.m; c = [Z0]*m
        for i, x in enumerate(a):
            if x == Z0: continue
            for j, z in enumerate(b):
                if z == Z0: continue
                c[(i+j) % m] = gadd(c[(i+j) % m], gmul(x, z))
        return c
    def red(self, a):
        """reduce a (mod x^m - 1 representative) modulo Phi_m -> list of deg coefficients"""
        a = list(a)
        for k in range(len(a)-1, self.deg-1, -1):
            t = a[k]
            if t == Z0: continue
            for j, cj in enumerate(self.phi):          # subtract t * x^{k-deg} * Phi
                idx = k - self.deg + j
                a[idx] = gsub(a[idx], (t[0]*cj, t[1]*cj))
        return a[:self.deg]
    def f2mod(self, a):
        f = self.f1; d = self.fdeg
        for k in range(a.bit_length()-1, d-1, -1):
            if (a >> k) & 1: a ^= (f << (k-d))
        return a
    def f2mul(self, a, b):
        r = 0
        while b:
            if b & 1: r ^= a
            a <<= 1; b >>= 1
        return self.f2mod(r)
    def f2inv(self, a):
        r, e, base = 1, 2**self.fdeg - 2, a
        while e:
            if e & 1: r = self.f2mul(r, base)
            base = self.f2mul(base, base); e >>= 1
        return r
    def res(self, a):
        bits = 0
        for k, (x, y) in enumerate(self.red(a)):
            if (x + y) % 2: bits ^= (1 << k)
        return self.f2mod(bits)
    def half_res(self, a):
        """residue at P1 of a/2, assuming v_P1(a) >= v_P1(2); None if that fails."""
        t = self.red(self.mul(a, self.y))
        if any(x % 2 or z % 2 for x, z in t): return None
        bits = 0
        for k, (x, z) in enumerate(t):
            if ((x//2) + (z//2)) % 2: bits ^= (1 << k)
        return self.f2mul(self.f2mod(bits), self.f2inv(self.yres))

def primes_of(u): return [p for p in range(3, u+1, 2) if u % p == 0 and all(p % q for q in range(2, p))]

def test_random_general(u, trials, rng):
    from itertools import combinations
    from math import prod
    I = primes_of(u); m = prod(I); N = u*u; n = 4*N
    ctx = Ctx(m)
    subsets = [S for r in range(len(I)+1) for S in combinations(I, r)]
    kS = {S: sum(m // p for p in S) % m for S in subsets}
    hold = total = rw = memb = 0
    for _ in range(trials):
        hvec = [rng.choice((1, -1)) for _ in range(n)]
        H = split_h(hvec, N)
        c = [((H[0][z]+H[1][z]+H[2][z]+H[3][z])//2, 0) for z in range(N)]
        d = [((H[0][z]-H[1][z]+H[2][z]-H[3][z])//2, 0) for z in range(N)]
        g = [((H[0][z]-H[2][z])//2, (H[1][z]-H[3][z])//2) for z in range(N)]
        def val(x, S):
            out = [Z0]*m
            for z in range(N):
                k = (z * kS[S]) % m
                out[k] = gadd(out[k], x[z])
            return out
        V = {(name, S): val(x, S) for name, x in (('c', c), ('d', d), ('g', g)) for S in subsets}
        if any(ctx.res(v) == 0 for v in V.values()): continue
        one = [Z0]*m; one[0] = (1, 0)
        NR = DR = NW = DW = one
        for S in subsets:
            if len(S) % 2 == 0:
                NR = ctx.mul(NR, V[('d', S)]); DR = ctx.mul(DR, V[('c', S)])
                NW = ctx.mul(NW, V[('g', S)]); DW = ctx.mul(DW, V[('c', S)])
            else:
                NR = ctx.mul(NR, V[('c', S)]); DR = ctx.mul(DR, V[('d', S)])
                NW = ctx.mul(NW, V[('c', S)]); DW = ctx.mul(DW, V[('g', S)])
        MR = [gsub(a, b) for a, b in zip(NR, DR)]
        MW = [gmul(gsub(a, b), (1, -1)) for a, b in zip(NW, DW)]     # (W-1)/(1+i) = (NW-DW)(1-i)/(2 DW)
        lr = ctx.half_res(MR); lw = ctx.half_res(MW)
        if lr is None or lw is None: memb += 1; continue
        lamR = ctx.f2mul(lr, ctx.f2inv(ctx.res(DR))); lamW = ctx.f2mul(lw, ctx.f2inv(ctx.res(DW)))
        total += 1
        if (lamR ^ lamW) == 1: hold += 1          # u^2 / c_empty has residue 1 (both odd)
        if all(t == Z0 for t in ctx.red(MR)) and all(t == Z0 for t in ctx.red(MW)): rw += 1
    print(f"(2') n = 4u^2 = {n} (u={u}, I={I}, residue field F_2^{ctx.fdeg}, {1+len(ctx.others)} primes over 2): "
          f"{total} random sign arrays with all 3*2^|I| values 2-units; identity holds {hold}/{total}; "
          f"R=W=1: {rw}; membership failures: {memb}")
    return hold == total and memb == 0



if __name__ == '__main__':
    rng = random.Random(7)
    trials = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    order4()
    a = all([test_random(3, 3, trials, rng), test_random(5, 5, trials // 2, rng), test_random(9, 3, trials // 3, rng)])
    a = a and all([test_random_general(u, t, rng) for u, t in ((3, trials), (5, trials // 2), (7, trials // 3),
                                                            (9, trials // 3), (15, trials // 5))])
    b = pba_near_miss()
    print("SUMMARY binary identity:", "consistent" if (a and b) else "INCONSISTENT")

