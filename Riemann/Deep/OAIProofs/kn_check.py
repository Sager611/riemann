"""Numerical (exact rational) test of the lemmas of OpenAI's "The crossing number of complete graphs"
(Sep 23 2026) on actual drawings of K_n (n = 2s+1): straight-line drawings and 1-bend polyline drawings
in general position without adjacent-edge crossings (the paper's "normalized" drawings).
 K1 Lemma 2.1 (cycle pairing): J vanishes on cycle space x cycle space
 K2 Lemma 3.2 (two-triangle interpolation) and Lemma 3.3 (order detection) for random nodes and orders
 K3 Lemma 4.2 (injectivity): the evaluation map E: S -> R^T has rank C(s,2)^2
 K4 Prop 4.3: E^T diag(a) E = 0 (isotropy), hence |T|/2 >= C(s,2)^2
usage: python3 kn_check.py SEED NDRAWINGS"""
import random, sys, itertools
from fractions import Fraction as Fr

def orient(a, b, c):
    v = (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
    return (v > 0) - (v < 0)

def seg_cross(a, b, c, d):
    """+1/-1 sign of proper crossing of a->b with c->d, 0 if none; raises on degeneracy"""
    o1, o2, o3, o4 = orient(a, b, c), orient(a, b, d), orient(c, d, a), orient(c, d, b)
    if 0 in (o1, o2, o3, o4):
        if {a, b} & {c, d}: return 0      # shared endpoint (adjacent spokes / consecutive segments)
        raise ValueError("degenerate")
    if o1 != o2 and o3 != o4:
        dx1, dy1 = b[0]-a[0], b[1]-a[1]; dx2, dy2 = d[0]-c[0], d[1]-c[1]
        det = dx1*dy2 - dy1*dx2
        return 1 if det > 0 else -1
    return 0

def make_drawing(n, rng, bends):
    """vertices + edge polylines; returns (pts, poly dict (i,p)->list of points from i to p)"""
    R = 10**6
    pts = [(rng.randint(-R, R), rng.randint(-R, R)) for _ in range(n)]
    poly = {}
    for i in range(n):
        for p in range(i+1, n):
            path = [pts[i]]
            for _ in range(bends):
                if rng.random() < 0.5:
                    path.append((rng.randint(-2*R, 2*R), rng.randint(-2*R, 2*R)))
            path.append(pts[p])
            poly[(i, p)] = path
    return pts, poly

def crossings(n, poly):
    E = sorted(poly)
    I = {}
    total = 0
    for e, f in itertools.combinations(E, 2):
        ce = poly[e]; cf = poly[f]
        s = 0; cnt = 0
        for a, b in zip(ce, ce[1:]):
            for c, d in zip(cf, cf[1:]):
                x = seg_cross(a, b, c, d)
                if x: s += x; cnt += 1
        if set(e) & set(f):
            if cnt: raise ValueError("adjacent crossing")
        # self-check of edges handled by construction (<=1 bend)
        I[(e, f)] = s; I[(f, e)] = -s
        total += cnt
    return I, total

def angle_key(w, q, cut):
    import math
    a = math.atan2(q[1]-w[1], q[0]-w[0]) - cut
    while a < 0: a += 2*math.pi
    return a

def build_J(n, pts, poly, I, rng):
    """J on oriented basis edges (i<p, oriented i->p): dict (e,f)->Fraction"""
    E = sorted(poly)
    import math
    order = {}
    for w in range(n):
        cut = rng.uniform(0, 2*math.pi)
        inc = []
        for e in E:
            if w in e:
                path = poly[e] if e[0] == w else poly[e][::-1]
                inc.append((angle_key(pts[w], path[1], cut), e))
        inc.sort()
        order[w] = {e: k for k, (_, e) in enumerate(inc)}
    J = {}
    for e in E:
        for f in E:
            if e == f: J[(e, f)] = Fr(0); continue
            v = Fr(I[(e, f)])
            for w in set(e) & set(f):
                se = 1 if e[0] == w else -1; sf = 1 if f[0] == w else -1
                T = 1 if order[w][e] < order[w][f] else -1
                v += Fr(se*sf*T, 2)
            J[(e, f)] = v
    return E, J

def nullspace_int(rows, ncols):
    """basis of kernel of integer matrix (list of rows) over Q"""
    M = [[Fr(x) for x in r] for r in rows]
    piv = []; r = 0
    for c in range(ncols):
        pr = next((k for k in range(r, len(M)) if M[k][c] != 0), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]
        inv = 1/M[r][c]; M[r] = [x*inv for x in M[r]]
        for k in range(len(M)):
            if k != r and M[k][c] != 0:
                fac = M[k][c]; M[k] = [x - fac*y for x, y in zip(M[k], M[r])]
        piv.append(c); r += 1
    free = [c for c in range(ncols) if c not in piv]
    basis = []
    for fc in free:
        v = [Fr(0)]*ncols; v[fc] = Fr(1)
        for k, pc in enumerate(piv): v[pc] = -M[k][fc]
        basis.append(v)
    return basis

def rank(rows):
    M = [list(r) for r in rows]
    if not M: return 0
    ncols = len(M[0]); r = 0
    for c in range(ncols):
        pr = next((k for k in range(r, len(M)) if M[k][c] != 0), None)
        if pr is None: continue
        M[r], M[pr] = M[pr], M[r]
        for k in range(r+1, len(M)):
            if M[k][c] != 0:
                fac = M[k][c]/M[r][c]; M[k] = [x - fac*y for x, y in zip(M[k], M[r])]
        r += 1
        if r == len(M): break
    return r

def K1(n, E, J):
    idx = {e: k for k, e in enumerate(E)}
    bd = [[0]*len(E) for _ in range(n)]
    for e in E:
        bd[e[1]][idx[e]] += 1; bd[e[0]][idx[e]] -= 1
    Z = nullspace_int(bd, len(E))
    assert len(Z) == len(E) - n + 1
    for z in Z:
        for zp in Z:
            val = sum(z[a]*J[(E[a], E[b])]*zp[b] for a in range(len(E)) if z[a] for b in range(len(E)) if zp[b])
            assert val == 0, "cycle pairing fails"
    return len(Z)

def sym_basis(s):
    pairs = [(a, b) for a in range(s-1) for b in range(a, s-1)]
    return [(p1, p2) for p1 in pairs for p2 in pairs]

def mono_sym(a, b, x, u):
    return x**a * u**b + (x**b * u**a if a != b else 0)

def K34(n, E, J, I, tnodes):
    s = (n-1)//2
    w = []
    for i in range(n):
        d = Fr(1)
        for j in range(n):
            if j != i: d *= (tnodes[i]-tnodes[j])
        w.append(1/d)
    basis = sym_basis(s)
    T = []
    for e in E:
        for f in E:
            if e != f and not (set(e) & set(f)) and I[(e, f)] != 0:
                T.append((e, f))
    rows = []
    for (e, f) in T:
        x, u, y, v = (tnodes[e[0]], tnodes[e[1]], tnodes[f[0]], tnodes[f[1]])
        rows.append([mono_sym(a, b, x, u)*mono_sym(c, d, y, v) for ((a, b), (c, d)) in basis])
    rk = rank(rows)
    assert rk == len(basis), ("injectivity fails", rk, len(basis))
    # isotropy: sum_alpha a_alpha E_alpha(P) E_alpha(Q) = 0 for all basis P,Q
    acoef = []
    for (e, f) in T:
        i, p = e; k, q = f
        x, u, y, vv = tnodes[i], tnodes[p], tnodes[k], tnodes[q]
        Delta = (x-u)*(y-vv)*(x-y)*(x-vv)*(u-y)*(u-vv)
        acoef.append(J[(e, f)]*w[i]*w[p]*w[k]*w[q]*Delta)
    nb = len(basis)
    for P in range(nb):
        for Q in range(P, nb):
            val = sum(acoef[r]*rows[r][P]*rows[r][Q] for r in range(len(T)))
            assert val == 0, "isotropy fails"
    return len(T), rk

def K2(smax, rng):
    # Lemma 3.2 and 3.3 with random nodes and random orders
    for s in range(2, smax+1):
        for trial in range(3):
            nodes = rng.sample(range(-50, 50), 2*s)
            perm = list(range(2*s)); rng.shuffle(perm)       # prescribed order
            u = [Fr(x) for x in nodes]
            G = []
            for p in range(2*s):
                d = Fr(1)
                for q in range(2*s):
                    if q != p: d *= (u[p]-u[q])
                G.append(d)
            # two-triangle pattern in the prescribed order (positions in perm)
            Ls = [(perm[a], perm[b]) for a in range(s) for b in range(a)] + [(perm[a], perm[b]) for a in range(s, 2*s) for b in range(a+1, 2*s)]
            for (da, db) in ((s-1, s-2), (s-2, s-1)):
                rows = [[u[p]**i * u[q]**j for i in range(da+1) for j in range(db+1)] for (p, q) in Ls]
                assert rank(rows) == s*(s-1) == len(rows[0]), ("triangles", s)
            # order pairing M on P_{s-1,s-2} x P_{s-2,s-1}
            pos = {perm[k]: k for k in range(2*s)}
            def S(p, q):
                if p == q: return Fr(0)
                return Fr(1, 2) if pos[p] < pos[q] else Fr(-1, 2)
            A = [(i, j) for i in range(s) for j in range(s-1)]
            B = [(i, j) for i in range(s-1) for j in range(s)]
            Mmat = [[sum(S(p, q)*(u[p]-u[q])/(G[p]*G[q]) * u[p]**(a[0]+b[0]) * u[q]**(a[1]+b[1]) for p in range(2*s) for q in range(2*s)) for b in B] for a in A]
            assert rank(Mmat) == len(A), ("detector", s)
    print(f"K2: Lemmas 3.2/3.3 ok for s<={smax} (random nodes and orders)", flush=True)

def make_drawing2(n, rng, bendprob, scale):
    R = 10**6
    pts = [(rng.randint(-R, R), rng.randint(-R, R)) for _ in range(n)]
    poly = {}
    for i in range(n):
        for p in range(i+1, n):
            path = [pts[i]]
            if rng.random() < bendprob:
                (x1, y1), (x2, y2) = pts[i], pts[p]
                L = max(abs(x2-x1), abs(y2-y1))
                mx, my = (x1+x2)//2, (y1+y2)//2
                path.append((mx + int(rng.uniform(-scale, scale)*L), my + int(rng.uniform(-scale, scale)*L)))
            path.append(pts[p])
            poly[(i, p)] = path
    return pts, poly

def tcount(poly, I):
    return sum(1 for (e, f), v in I.items() if e < f and v != 0 and not (set(e) & set(f)))

def climb(n, rng, iters, scale):
    """hill-climb on vertex positions and bends to minimise |T|/2 (pairs with nonzero algebraic crossing)"""
    while True:
        try:
            pts, poly = make_drawing(n, rng, 0)
            I, tot = crossings(n, poly); break
        except ValueError:
            pass
    best = tcount(poly, I)
    R = 10**6
    for it in range(iters):
        pts2 = list(pts); poly2 = {e: list(v) for e, v in poly.items()}
        if rng.random() < 0.3:
            v = rng.randrange(n)
            nv = (pts[v][0] + rng.randint(-R//5, R//5), pts[v][1] + rng.randint(-R//5, R//5))
            pts2[v] = nv
            for e in poly2:
                if e[0] == v: poly2[e][0] = nv
                if e[1] == v: poly2[e][-1] = nv
        else:
            e = rng.choice(sorted(poly2))
            (x1, y1), (x2, y2) = poly2[e][0], poly2[e][-1]
            L = max(abs(x2-x1), abs(y2-y1))
            if rng.random() < 0.2:
                poly2[e] = [poly2[e][0], poly2[e][-1]]
            else:
                mx, my = (x1+x2)//2, (y1+y2)//2
                poly2[e] = [poly2[e][0], (mx + int(rng.uniform(-scale, scale)*L), my + int(rng.uniform(-scale, scale)*L)), poly2[e][-1]]
        try:
            I2, tot2 = crossings(n, poly2)
        except ValueError:
            continue
        c2 = tcount(poly2, I2)
        if c2 <= best:
            pts, poly, I, best = pts2, poly2, I2, c2
    return pts, poly, I, best


if __name__ == "__main__" and not (len(sys.argv) > 3 and sys.argv[3] == "climb"):
    seed = int(sys.argv[1]); nd = int(sys.argv[2])
    rng = random.Random(seed)
    K2(5, rng)
    stats = {}
    done = 0
    while done < nd:
        n = rng.choice([5, 7, 9])
        bends = rng.choice([0, 1, 1])
        try:
            pts, poly = make_drawing2(n, rng, 0.6*bends, 0.35) if bends else make_drawing(n, rng, 0)
            I, total = crossings(n, poly)
        except ValueError:
            continue
        E, J = build_J(n, pts, poly, I, rng)
        nz = K1(n, E, J)
        tn = [Fr(x) for x in rng.sample(range(-30, 30), n)]
        nT, rk = K34(n, E, J, I, tn)
        s = (n-1)//2; Z = (s*(s-1)//2)**2
        assert nT//2 >= Z
        key = (n, bends)
        st = stats.setdefault(key, [0, 10**9, 10**9])
        st[0] += 1; st[1] = min(st[1], nT//2); st[2] = min(st[2], total)
        done += 1
    for k, v in sorted(stats.items()):
        print(f"n={k[0]} bends<={k[1]}: {v[0]} drawings, K1 K3 K4 ok; min |T|/2 = {v[1]}, min crossings = {v[2]}, Z(n) = {(((k[0]-1)//2)*((k[0]-1)//2-1)//2)**2}")

if __name__ == "__main__" and len(sys.argv) > 3 and sys.argv[3] == "climb":
    rng = random.Random(int(sys.argv[1]))
    for n, iters in ((5, 300), (7, 2500), (9, 2500)):
        res = []
        for rep in range(int(sys.argv[2])):
            pts, poly, I, best = climb(n, rng, iters, 0.8)
            E, J = build_J(n, pts, poly, I, rng)
            K1(n, E, J)
            tn = [Fr(x) for x in rng.sample(range(-30, 30), n)]
            nT, rk = K34(n, E, J, I, tn)
            nb = sum(1 for v in poly.values() if len(v) > 2)
            res.append((best, nb))
        s = (n-1)//2
        print(f"climb n={n}: min |T|/2 over {len(res)} runs = {min(r[0] for r in res)} (Z={(s*(s-1)//2)**2}); runs (|T|/2, #bent edges): {res}; K1 K3 K4 ok on all", flush=True)
