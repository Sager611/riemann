#!/usr/bin/env python3
# code_audit.py -- independent test of cert.py's precomputed tables and cell bounds (CODE-AUDIT.md).
# (1) kernel moment tables: >= 36 random stored entries T[p][i] (floats rounded up) against p int_{lam_i}^oo G(l) l^{p-1} dl
#     recomputed at 256 bits from the DEFINITION of the envelope G (piecewise formulas, Arb's rigorous integrator acb.integral,
#     exact C2/l^2 tail), not from cert.py's G/S_up code; lambda-grid: stored G(lam_i) against G(r^i) at 256 bits; p-grid: balls contain j/200.
# (2) KLN rows: log A(sigma) = c0 + c1 sigma against KLN's formula at 256 bits at random sigma; c1 > 0 and q1 >= 0 (A increasing, B decreasing).
# (3) >= 10 random cells of the certified scan (theta = 0.856, pipeline A1; plus 2 cells of pipeline A0, theta = 0.857): the cell bound
#     recomputed WITHOUT tables or moment majorants: x-powers at the left end, kernel at delta(Lb) (the worst end, as in the certificate),
#     S = 14 log 14 G(14 delta) + int G(delta t)(1 + log t) dt and F(sigma) = N(sigma,T0) G(delta T0) + int G(delta t) dN(sigma,t) with the
#     exact KLN shape, t-integrals by acb.integral (rigorous), far tails bounded above in closed form, sigma-integral by 2-point Gauss-Legendre
#     per sigma-cell for sigma >= 0.8 and the midpoint rule below (< 1% of the sum; not rigorous in sigma; relative error ~1e-4).
#     Cells at 128 bits (acb.integral tolerance 2^-40); tables at 256 bits. It must be <= the certified cell bound.
# usage: uv run --no-project --with python-flint --with mpmath --with sympy python3 Riemann/Deep/Ingham4/code_audit.py   (log: logs/code_audit.log)
import sys, os, math, random, time, gc, resource
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import cert
from flint import arb, acb, ctx, fmpq
OUT = []
def log(*a):
    s = ' '.join(str(x) for x in a); print(s, flush=True); OUT.append(s)
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20
random.seed(20261106)
T_START = time.time()
RT = None
# ------------------------------------------------------------------ the envelope G from its definition (independent of cert.S_up / Kernel.G)
def S1v(): return arb('0.2173')
def ycrit():
    """y_c in (2.5, 3): sin(y_c)/y_c = s1 (bisection at the current precision; sin y/y decreases on (0, pi))."""
    lo_, hi_ = arb('2.5'), arb('3.0')
    for _ in range(ctx.prec + 10):
        mid = (lo_ + hi_)/2
        if (mid.sin()/mid) > S1v(): lo_ = mid
        else: hi_ = mid
    return arb(((lo_ + hi_)/2).mid())
class Env:
    """G(l) = f_a(l) f_b(l), f_c(l) = S(c l/2) + 2 dmax min(c, 1/l), S(y) = sin y/y (y <= y_c), s1 (y_c <= y <= 1/s1), 1/y (y >= 1/s1)."""
    def __init__(self, a, dmax):
        self.a = arb(a); self.b = 1 - self.a; self.dm = arb(dmax); s1 = S1v(); yc = ycrit()
        self.cs = [self.a, self.b] if self.a != self.b else [self.a]
        self.pw2 = self.a == self.b
        bps = []
        for c in self.cs:
            for z in (2*yc/c, 2/(c*s1), 1/c):
                if all(abs(float(z.mid()) - float(w.mid())) > 1e-12 for w in bps): bps.append(z)
        self.bps = sorted(bps, key=lambda z: float(z.mid()))
        self.yc, self.s1 = yc, s1
        self.lend = self.bps[-1]                                  # beyond: G = C2/l^2 exactly
        self.C2 = (2/self.a + 2*self.dm)*(2/self.b + 2*self.dm)
    def factor(self, c, x, piece_mid):
        """f_c on the analytic piece containing the real point piece_mid (formulas chosen by position, not by comparison of x)."""
        y = c*piece_mid/2
        if y <= self.yc: S = (c*x/2).sin()/(c*x/2)
        elif y <= 1/self.s1: S = acb(self.s1)
        else: S = 2/(c*x)
        corr = acb(c) if piece_mid <= 1/c else 1/x
        return S + 2*self.dm*corr
    def G(self, x, piece_mid):
        if self.pw2:
            f = self.factor(self.a, x, piece_mid); return f*f
        return self.factor(self.a, x, piece_mid)*self.factor(self.b, x, piece_mid)
    def Greal(self, l):
        return self.G(acb(l), l).real
    def pieces(self, lo_, hi_):
        pts = [lo_] + [b for b in self.bps if b > lo_ and b < hi_] + [hi_]
        return [(pts[i], pts[i + 1]) for i in range(len(pts) - 1)]
def integrate(env, h, lo_, hi_, tolbits=70):
    """int_lo^hi G(l) h(l) dl over the analytic pieces (acb.integral, rigorous enclosure)."""
    tot = arb(0)
    for (u, v) in env.pieces(lo_, hi_):
        m = (u + v)/2
        f = lambda x, analytic, m=m: env.G(x, m)*h(x, analytic)
        tot += acb.integral(f, u, v, rel_tol=arb(2)**(-tolbits)).real
    return tot
# ------------------------------------------------------------------ (1) moment tables, lambda-grid, p-grid
def test_tables():
    n_mom = n_g = 0; worst_mom = worst_g = 0.0
    for kind in ('tri', 'trap'):
        for dm in ('1e-9', '0.025', '0.15'):
            ctx.prec = 128
            K = cert.Kernel(2, arb(dm), cert.PGRID) if kind == 'tri' else cert.TrapKernel(arb(cert.A_TRAP), arb(dm), cert.PGRID)
            assert list(K.T.keys()) == [str(q) for q in cert.PGRID]
            ctx.prec = 256
            env = Env('0.5' if kind == 'tri' else cert.A_TRAP, dm)
            r = arb(1) + arb(1)/500
            for _ in range(6):
                j = random.randint(1, 300); i = random.randint(0, K.istore)
                p = arb(fmpq(j, 200)); li = r**i
                hpow = lambda x, analytic, p=p: (x.log(analytic=analytic)*(p - 1)).exp()
                if li < env.lend:
                    V = p*integrate(env, hpow, li, env.lend) + p*env.C2*env.lend**(p - 2)/(2 - p)
                else:
                    V = p*env.C2*li**(p - 2)/(2 - p)
                stored = arb(list(K.T.values())[j - 1][i])          # tables are keyed in p-grid order (j/200, j = 1..300)
                assert V.upper() <= stored, (kind, dm, j, i, V, stored)
                n_mom += 1; ratio = float((stored/V).mid()); worst_mom = max(worst_mom, 1/ratio)
                log('  moment %s dmax=%s p=%d/200 i=%d (lambda=%.4f): definition %s <= stored %.10g  (stored/definition = %.6f)' % (
                    kind, dm, j, i, float(li.mid()), V.str(12, radius=False), float(stored.mid()), ratio))
                gi = random.randint(0, K.N)
                Gdef = env.Greal(r**gi)
                assert Gdef.upper() <= K.Gv[gi].upper(), (kind, dm, gi)
                n_g += 1
            del K; gc.collect()
            log('  kernel %s dmax=%s done; rss %.0f MB' % (kind, dm, rss()))
    ctx.prec = 256
    for j in range(1, 301):
        assert cert.PGRID[j - 1].contains(arb(fmpq(j, 200)))
    log('TABLES: %d moment entries and %d lambda-grid G values: every 256-bit value from the definition <= the stored upper value; all 300 p-grid balls contain j/200' % (n_mom, n_g))
# ------------------------------------------------------------------ (2) KLN rows
def test_kln():
    """the 128-bit KLN rows of both pipelines against a 256-bit re-evaluation of KLN's formula (cert.py's module re-executed at 256 bits,
       every constant recomputed) at random sigma; and c1 > 0, q1 >= 0 (A increasing, B decreasing in sigma) in every row."""
    n = 0
    src = open(cert.__file__).read().replace('ctx.prec = 128', 'ctx.prec = 256')
    for pipe in ('A0', 'A1'):
        ctx.prec = 128
        P = cert.build_kln() if pipe == 'A0' else cert.build_kln_H0(cert.H0_RH)
        for p in P:
            assert p.c1 > 0 and p.q1 >= 0, (pipe, p.name)
        ns = {'__name__': 'cert256'}; exec(compile(src, 'cert256', 'exec'), ns)          # sets ctx.prec = 256
        if pipe == 'A1':
            ns['H0'] = arb(cert.H0_RH); ns['HH'] = ns['H0'] - 1; ns['LH0'] = ns['H0'].log(); ns['llH0'] = ns['LH0'].log()
        rows = {r[0]: r for r in cert.KLN_ROWS}
        for p in random.sample(P, 5):
            q = ns['KP'](rows[p.name])                        # all of KLN's constants for this row, at 256 bits
            for _ in range(2):
                s = arb(random.uniform(float(p.smin.mid()) + 0.01, 0.99))
                v256 = ns['logC1'](q.al, q.d, q.de, q.k, ns['HH'], s) - (2*ns['PI']*q.d).log()
                v128 = p.logA(s)
                assert v128.overlaps(v256) and v128.overlaps(q.logA(s)) and p.Bc(s).overlaps(q.Bc(s)), (pipe, p.name)
                n += 1
                log('  KLN %s row %s sigma=%.4f: log A (128-bit affine row) %s overlaps KLN formula at 256 bits %s; B overlaps' % (pipe, p.name, float(s.mid()), v128.str(12, radius=False), v256.str(16, radius=False)))
        ctx.prec = 128
    log('KLN: %d random (row, sigma) pairs in both pipelines: the 128-bit affine rows overlap the 256-bit formula values; c1 > 0 and q1 >= 0 in all 39 rows' % n)
# ------------------------------------------------------------------ (3) cells
def direct_cell(S, dfun, La, Lb, a_str, dmax):
    """the cell bound without tables (see header), at 128 bits with acb.integral tolerance 2^-40."""
    ctx.prec = 128
    La, Lb = arb(La), arb(Lb); T0 = arb('3e12'); B = arb('0.8733819666132076'); sig0 = arb('0.515')
    db = dfun(Lb); da = dfun(La)
    env = Env(a_str, dmax)
    lam0 = db*T0
    # S: int_{14}^oo t log t d(-G(delta t)) = 14 log 14 G(14 delta) + int_{14 delta}^oo G(l)(1 + log(l/delta)) dl/delta
    hS = lambda x, analytic: (1 + (x/db).log(analytic=analytic))/db
    l14 = 14*db; lT = (env.lend.max(l14))*1000
    Sv = 14*arb(14).log()*env.Greal(l14) + integrate(env, hS, l14, lT, 40) + env.C2/db*(2 + (lT/db).log())/lT
    Z1 = 2*((sig0 - 1)*La).exp()*Sv
    # F(sigma) at 2 Gauss nodes per sigma-cell
    gn = [(-1/arb(3).sqrt(), arb(1)), (1/arb(3).sqrt(), arb(1))]
    lT = (env.lend.max(lam0))*1000; tT = lT/db; ltT = tT.log()
    Glam0 = env.Greal(lam0)
    Z2 = arb(0)
    g1 = [(arb(0), arb(2))]                                   # midpoint rule where the contribution is < 1% (sigma < 0.8)
    for (s, t, P, cs) in S.cells:
        s, t = arb(s.mid()), arb(t.mid())
        for (xi, w) in (gn if s >= arb('0.8') else g1):
            sig = (s + t)/2 + (t - s)/2*xi
            A = (P.c0 + P.c1*sig).exp(); Bc = P.q0 - P.q1*sig; a = arb(8)/3*(1 - sig); e = 5 - 2*sig
            def hN(x, analytic, A=A, Bc=Bc, a=a, e=e):
                tt = x/db; lt = tt.log(analytic=analytic)
                return (A*((a - 1)*lt).exp()*((e - 1)*lt.log(analytic=analytic)).exp()*(a*lt + e) + 2*Bc*lt/tt)/db
            lt0 = T0.log()
            N0 = A*lt0**e*T0**a + Bc*lt0**2
            tail = env.C2/db**2*2*(A*ltT**e*tT**(a - 2)/(2 - a - e/ltT) + Bc*ltT**2/tT**2/(2 - 2/ltT))
            F = N0*Glam0 + integrate(env, hN, lam0, lT, 40) + tail
            Z2 += 2*La*w*(t - s)/2*((sig - 1)*La).exp()*F
    J = int(math.floor(float(((Lb + da)/arb(2).log()).upper())))
    PP = 2*(Lb + da)*((da/2).exp()*(-La/2).exp()/4 + J*da.exp()*(-2*La/3).exp()/9 + arb(J).log()/(db*La.exp()))
    Et = 2*(-3*La).exp()
    return Z1 + Z2 + PP + Et
def test_cells():
    res = []
    for pipe, th, a_str in (('A1', '0.856', cert.A_TRAP), ('A0', '0.857', '0.5')):
        ctx.prec = 128
        cert.main_setup(('trap',) if pipe == 'A1' else ('tri',), ('1e-9', '0.025'))
        S = cert.build(False, pipe=pipe)
        num = int(round(float(th)*1000)); dfun = cert.delta_theta(arb(num)/1000)
        cells = []; L = 42.0
        while L < 300 - 1e-12:
            Lb = min(L + cert.theta_width(L), 300); cells.append((L, Lb)); L = Lb
        crit = [c for c in cells if 165 <= c[0] < 215]; other = [c for c in cells if not (165 <= c[0] < 215)]
        pick = (random.sample(crit, 6) + random.sample(other, 4)) if pipe == 'A1' else random.sample(crit, 2)
        worstL = 191.82 if pipe == 'A1' else 191.68                  # plus the cell where the certificate's bound is largest
        pick += [min(cells, key=lambda c_: abs(c_[0] - worstL))]
        for (La, Lb) in sorted(pick):
            ctx.prec = 128
            Rc = cert.cell_bound(S, dfun, La, Lb, 2)
            da = float(dfun(arb(La)).upper())
            dmax = [K.dmax for K in S.kernels[2] if arb(da) < K.dmax][0]
            t0 = time.time()
            Rd = direct_cell(S, dfun, La, Lb, a_str, dmax)
            ok = Rd.upper() <= Rc
            res.append(ok)
            log('  cell %s theta=%s [%.2f, %.2f]: direct (no tables) %s <= certified %s : %s  (certified/direct = %.4f, %.0fs, rss %.0f MB)' % (
                pipe, th, La, Lb, Rd.str(8, radius=False), Rc.str(8, radius=False), ok, float((Rc/Rd).mid()), time.time() - t0, rss()))
            assert ok
        del S; cert.KERNELS = cert.KERNELS_T = None; gc.collect()
    log('CELLS: %d cells (10 + 2 random, plus the worst cell of each pipeline): every table-free direct value <= the certified cell bound' % len(res))
if __name__ == '__main__':
    cert.main_setup((), ())            # PGRID and the input checks only (no kernels)
    log('code_audit.py: python-flint', __import__('flint').__version__, '; tables built at 128 bits as in cert.py and tested at 256 bits; KLN rows re-evaluated at 256 bits; cells recomputed without tables')
    test_tables(); test_kln(); test_cells()
    log('ALL TESTS PASSED; peak rss %.0f MB; %.0fs' % (rss(), time.time() - T_START))
    open(os.path.join(HERE, 'logs', 'code_audit.log'), 'w').write('\n'.join(OUT) + '\n')
