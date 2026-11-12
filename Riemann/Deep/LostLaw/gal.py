# gal.py (lostlaw): Galerkin bottom of the even Weil form Q_c and of deleted / model forms
#   Q' = Q - sum_{k in D} 2<f, cos(g_k x)>^2 + sum_{a in A} 2<f, cos(a x)>^2      (even sector)
# Matrix: ChekAudit galerkin_fast.build (closed forms, Arb balls; cosine basis phi_j = N_j cos(2 pi j x/L), j = 0..N).
# Method: block approx-LU solves give x (ground state, subspace iteration), Z = R^-1 b, Z2 = R^-2 b; the deleted ground state
# (A - lam')^-1 e0 lies in span{x, Z, Z2} up to relative lam'/E2; Rayleigh-Ritz there. Rigorous Arb Rayleigh quotients
# (upper bounds) for chosen sets, with zero balls certified by a sign change of Hardy's Z.
# Usage: gal.py c N bits ymax_singles sets.json [rig]
import sys, os, json, math, time, random, resource
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import arb, acb, arb_mat, ctx
import mpmath as mp
import galerkin_fast as gf

CAP = 280e6
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
def guard(tag):
    if rss() > CAP:
        print(json.dumps({'abort': tag, 'rss': rss()}), flush=True); sys.exit(3)

def g2(y):
    y = abs(y); return math.log((2 + math.sqrt(4 - y*y))/y) if y < 2 else 0.0

def Zfun(t):
    s_ = acb(arb(1)/2, t)
    th = acb(arb(1)/4, t/2).lgamma().imag - t/2*arb.pi().log()
    return (s_.zeta()*acb(th.cos(), th.sin())).real

def zero_ball(gs, rexp):
    r = arb(10)**(-rexp); g = arb(gs)
    zl = Zfun(g - r); zr = Zfun(g + r)
    return g.union(g - r).union(g + r), bool((zl < 0 and zr > 0) or (zl > 0 and zr < 0))

def bcol(L, N, g, mid=True):
    pi = arb.pi(); h = L/2; out = []
    for j in range(N+1):
        om = 2*pi*j/L; Nj = (1/L).sqrt() if j == 0 else (2/L).sqrt()
        v = Nj*(((om - g)*h).sin()/(om - g) + ((om + g)*h).sin()/(om + g))
        out.append(v.mid() if mid else v)
    return out

def dot(a, b): return sum((x*y for x, y in zip(a, b)), arb(0))
def midv(a): return [x.mid() for x in a]
def col(Mt, j): return [Mt[i, j] for i in range(Mt.nrows())]
def tomat(cols):
    n = len(cols[0]); A = arb_mat(n, len(cols))
    for j, c in enumerate(cols):
        for i in range(n): A[i, j] = c[i]
    return A
def matvec(M, v):
    r = M*tomat([v]); return [r[i, 0] for i in range(M.nrows())]

def gs(cols):
    Q = []
    for c in cols:
        v = list(c)
        for _ in range(2):
            for q in Q:
                p = dot(q, v).mid(); v = [(vi - p*qi).mid() for vi, qi in zip(v, q)]
        nr = dot(v, v).sqrt().mid()
        if nr == 0: continue
        Q.append([(vi/nr).mid() for vi in v])
    return Q

def ritz(M, Q, dels, adds):
    # H = Q^T (M - 2 sum b_d b_d^T + 2 sum b_a b_a^T) Q (midpoints), smallest eigenpair by mpmath
    MQ = M*tomat(Q); m = len(Q)
    H = mp.matrix(m, m)
    pd = [[dot(q, b).mid() for q in Q] for b in dels]; pa = [[dot(q, b).mid() for q in Q] for b in adds]
    for i in range(m):
        for j in range(i, m):
            h = dot(Q[i], [MQ[r, j] for r in range(MQ.nrows())]).mid()
            for p in pd: h -= 2*p[i]*p[j]
            for p in pa: h += 2*p[i]*p[j]
            H[i, j] = H[j, i] = mp.mpf(h.str(int(ctx.prec*0.30) + 5, radius=False))
    E, V = mp.eigsy(H)
    i0 = min(range(m), key=lambda i: E[i])
    w = [arb(mp.nstr(V[r, i0], int(ctx.prec*0.30) + 5)) for r in range(m)]
    v = [sum((w[r]*Q[r][i] for r in range(m)), arb(0)).mid() for i in range(len(Q[0]))]
    return E[i0], v, sorted(E)

def rayleigh(M, v, dels, adds):
    Mv = matvec(M, v)
    num = dot(v, Mv)
    for b in dels: num -= 2*dot(b, v)**2
    for b in adds: num += 2*dot(b, v)**2
    return num/dot(v, v)

def main():
    c = float(sys.argv[1]); c = int(c) if c.is_integer() else c
    N = int(sys.argv[2]); ctx.prec = int(sys.argv[3]); ymax = float(sys.argv[4])
    sets = json.loads(open(sys.argv[5]).read()) if len(sys.argv) > 5 and sys.argv[5] != '-' else []
    rig = len(sys.argv) > 6 and sys.argv[6] == 'rig'
    mp.mp.prec = ctx.prec; dps = int(ctx.prec*0.30) + 5
    t0 = time.time()
    zs = [l.strip() for l in open(os.path.join(HERE, 'zeros_dps179.txt')) if l.strip()]
    if os.environ.get('ZFILE'):
        for l in open(os.environ['ZFILE']):
            p_ = l.split()
            if len(p_) == 2: zs[int(p_[0]) - 1] = p_[1]
    L = arb(c).log(); twopic = 2*math.pi*c
    ks = [k for k in range(1, len(zs) + 1) if float(zs[k-1])/twopic <= ymax]
    if os.environ.get('KLIST'): ks = [k for k in ks if k in set(int(t) for t in os.environ['KLIST'].split(','))]
    need = set(ks)
    for s in sets:
        need |= set(s.get('del', []))
    # model points a = g_k + delta
    apts = []
    for s in sets:
        for (k, d) in s.get('add', []): apts.append((k, d))
    apts = sorted(set((k, d) for k, d in apts))
    M = gf.build(c, N); guard('build')
    t1 = time.time()
    print(json.dumps({'c': c, 'N': N, 'prec': ctx.prec, 'nzeros': len(ks), 't_build': round(t1 - t0, 1), 'rss': rss()}), flush=True)
    kl = sorted(need)
    B = {k: bcol(L, N, arb(zs[k-1])) for k in kl}
    A = {(k, d): bcol(L, N, arb(zs[k-1]) + arb(d)) for (k, d) in apts}
    random.seed(1)
    X0 = [[arb(random.uniform(-1, 1)) for _ in range(N+1)] for _ in range(3)]
    rhs = X0 + [B[k] for k in kl] + [A[a] for a in apts]
    Y = M.solve(tomat(rhs), algorithm="approx"); guard('solve1')
    X1 = gs([midv(col(Y, j)) for j in range(3)])
    Z = {k: midv(col(Y, 3 + i)) for i, k in enumerate(kl)}
    ZA = {a: midv(col(Y, 3 + len(kl) + i)) for i, a in enumerate(apts)}
    del Y
    rhs2 = X1 + [Z[k] for k in kl] + [ZA[a] for a in apts]
    Y = M.solve(tomat(rhs2), algorithm="approx"); guard('solve2')
    X2 = gs([midv(col(Y, j)) for j in range(3)])
    Z2 = {k: midv(col(Y, 3 + i)) for i, k in enumerate(kl)}
    ZA2 = {a: midv(col(Y, 3 + len(kl) + i)) for i, a in enumerate(apts)}
    del Y
    Y = M.solve(tomat(X2), algorithm="approx"); guard('solve3')
    X3 = gs([midv(col(Y, j)) for j in range(3)]); del Y
    E1, x, Es = ritz(M, X3, [], [])
    x = [xi.mid() for xi in x]
    nx = dot(x, x).sqrt().mid(); x = [(xi/nx).mid() for xi in x]
    Mx = [t.mid() for t in matvec(M, x)]
    E1a = dot(x, Mx).mid()
    t2 = time.time()
    print(json.dumps({'E1': mp.nstr(E1, 10), 'E2': mp.nstr(Es[1], 8), 'E3': mp.nstr(Es[2], 8), 't_solves': round(t2 - t1, 1), 'rss': rss()}), flush=True)
    lnE1 = float(mp.log(E1))
    # singles: 2D Ritz on span{x, Z_k}
    for k in ks:
        b = B[k]; z = Z[k]
        theta = 2*dot(b, z).mid()
        bx = dot(b, x).mid()
        w0 = 2*bx*bx/E1a
        pz = dot(x, z).mid(); u = [(zi - pz*xi).mid() for zi, xi in zip(z, x)]
        nu = dot(u, u).sqrt().mid(); q = [(ui/nu).mid() for ui in u]
        Mq = matvec(M, q)
        bq = dot(b, q).mid()
        h11 = E1a - 2*bx*bx; h12 = (dot(x, Mq) - 2*bx*bq).mid(); h22 = (dot(q, Mq) - 2*bq*bq).mid()
        H = mp.matrix([[mp.mpf(h11.str(dps, radius=False)), mp.mpf(h12.str(dps, radius=False))],
                       [mp.mpf(h12.str(dps, radius=False)), mp.mpf(h22.str(dps, radius=False))]])
        E, V = mp.eigsy(H); lp = min(E[0], E[1])
        y = float(zs[k-1])/twopic
        omt = mp.mpf((1 - theta).str(dps, radius=False))
        print(json.dumps({'k': k, 'gamma': round(float(zs[k-1]), 6), 'y': round(y, 5), 'lam1': mp.nstr(lp, 8),
                          'logratio': round(lnE1 - float(mp.log(lp)), 5), 'pred4g2': round(4*g2(y), 5),
                          'ln_1mtheta': round(float(mp.log(omt)), 3) if omt > 0 else None,
                          'ln_w0': round(float(mp.log(abs(mp.mpf(w0.str(dps, radius=False))))), 3) if w0 != 0 else None}), flush=True)
    guard('singles')
    t3 = time.time()
    # sets: Ritz on span{x, Z_D, Z2_D, ZA, ZA2}
    for s in sets:
        D = s.get('del', []); Ad = [tuple(a) for a in s.get('add', [])]
        cols = [x] + [Z[k] for k in D] + [Z2[k] for k in D] + [ZA[a] for a in Ad] + [ZA2[a] for a in Ad]
        Q = gs(cols)
        lp, v, _ = ritz(M, Q, [B[k] for k in D], [A[a] for a in Ad])
        pred = sum(4*g2(float(zs[k-1])/twopic) for k in D)
        rec = {'set': s['name'], 'del': D, 'add': s.get('add', []), 'lam1': mp.nstr(lp, 10),
               'logratio': round(lnE1 - float(mp.log(lp)), 5) if lp > 0 else None, 'sum4g2_del': round(pred, 5)}
        if rig:
            rexp = int(ctx.prec*0.30) - 5
            bb = []; okall = True
            for k in D:
                gb, ok = zero_ball(zs[k-1], min(rexp, 170)); okall &= ok; bb.append(bcol(L, N, gb, mid=False))
            ab = [bcol(L, N, arb(zs[a[0]-1]) + arb(a[1]), mid=False) for a in Ad]
            Xr = rayleigh(M, v, bb, ab)
            rec.update({'rig_upper': Xr.upper().str(8, radius=False), 'zero_balls_certified': okall})
        print(json.dumps(rec), flush=True)
        guard('set')
    if rig:
        Xf = rayleigh(M, x, [], [])
        print(json.dumps({'full_rig_upper': Xf.upper().str(10, radius=False)}), flush=True)
    print(json.dumps({'t_singles': round(t3 - t2, 1), 't_sets': round(time.time() - t3, 1), 'rss': rss(), 'wall': round(time.time() - t0, 1)}), flush=True)

if __name__ == '__main__':
    main()
