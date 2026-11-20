# Task (6): the worst cell of the theta = 0.856 scan, [La, Lb] = [191.82000000001713, 191.84000000001714] (cert: 0.9711), recomputed
# without cert.py: KLN constants from kln_indep's TeX formulas (H0 = 3e12, H = H0 - 1, 19 admissible rows), the kernel envelope G_a
# (a = 1/10, dmax = 1e-9, the cell's kernel) written from its definition, F(sigma) = min over rows of the Stieltjes integral
# int_{[T0,oo)} N_row(sigma,t) d(-G(delta t)) = N_row(sigma,T0) G(delta T0) + int_{T0}^oo G(delta t) dN_row/dt dt by Gauss-Legendre on
# the smooth pieces of G (+ tanh-sinh on [92, oo)), the sigma-integral by Gauss-Legendre per 0.0025-cell (no moment tables, no p-grid,
# no exponent freezing, no per-cell row choice). mpmath, 256 bits. Not rigorous in the quadratures; convergence shown by two node counts.
import sys, time
from mpmath import mp, mpf, sqrt, log, exp, pi, zeta, quad, inf, euler, findroot, sin, floor
from mpmath.calculus.quadrature import GaussLegendre
mp.prec = 256
t_start = time.time()
NG = int(sys.argv[1]) if len(sys.argv) > 1 else 5          # GL degree: 3*2^(NG-1) nodes per piece
NS = int(sys.argv[2]) if len(sys.argv) > 2 else 4          # GL nodes per sigma-cell (via degree)
exec(open(__file__.replace('worstcell.py', 'kln_core.py')).read())   # run(H0, eta, rows) -> per-row dict with c0, c1, q0, q1, k, smin, ok
H0 = mpf('3e12'); T0 = mpf('3e12')
rows = [r for r in kln_rows(H0, mpf('0.2561882')) if r['ok']]
print('rows admissible at H0 = 3e12:', [r['name'] for r in rows])
# kernel
s1 = mpf('0.2173'); a = mpf(1)/10; b = 1 - a; dm = mpf('1e-9')
def S(y):
    if y == 0: return mpf(1)
    return max(sin(y)/y, s1) if y <= pi else min(s1, 1/y)
def G(l):
    if l == 0: return (1 + 2*dm*a)*(1 + 2*dm*b)
    return (S(a*l/2) + 2*dm*min(a, 1/l))*(S(b*l/2) + 2*dm*min(b, 1/l))
C2 = (2/a + 2*dm)*(2/b + 2*dm)
ys = findroot(lambda y: sin(y)/y - s1, mpf('2.6'))
BRK = sorted(set([2*ys/a, 2*ys/b, 2*pi/a, 2*pi/b, 2/(a*s1), 2/(b*s1), 1/a, 1/b]))
YEND = 2/(a*s1)                                  # beyond: G = C2/l^2 exactly
gl = GaussLegendre(mp)
def nodes(u, v, deg):
    nd = gl.calc_nodes(deg, mp.prec)
    h = (v - u)/2; c = (u + v)/2
    return [(c + h*x, h*w) for (x, w) in nd]
QEDGES = [mpf(0), mpf(2), mpf(6), mpf(14), mpf(30), mpf(62), mpf(126)]
def tailq(h):
    """int_{YEND}^oo (C2/lam^2) h(lam) dlam = int_0^oo C2 e^{-q} h(YEND e^q)/YEND dq, GL on [0,126] (beyond: < 1e-38 relative here)."""
    tot = mpf(0)
    for u, v in zip(QEDGES, QEDGES[1:]):
        tot += sum(w*C2*exp(-q)*h(YEND*exp(q))/YEND for (q, w) in nodes(u, v, NG))
    return tot
th = mpf(856)/1000
def delta(L): return log(1 + exp(-(1 - th)*L))
sig0 = mpf('0.515'); Bh = mpf('0.8733819666132076')
def shape_logs(t, k):
    return log(log(k*t)), log(log(t)), log(t)
def R_bound(La, Lb, dK, xpow_L, pp_mode):
    """Z + PP + 2x^-3 with x-powers at xpow_L, the kernel at delta = dK, PP per pp_mode ('cell' (La, Lb) or 'point')."""
    D = -log(dK); lam0 = dK*T0
    assert dK*mpf(1) < dm
    # S = int_{[14,oo)} t log t d(-G(dK t)) = 14 log 14 G(14 dK) + int_14^oo G(dK t)(1 + log t) dt  (lam = dK t)
    pts = [14*dK] + [x for x in BRK if x > 14*dK] 
    Sv = 14*log(14)*G(14*dK)
    for u, v in zip(pts, pts[1:]):
        Sv += sum(w*G(exp(q))*(1 + log(exp(q)/dK))*exp(q)/dK for (q, w) in nodes(log(u), log(v), NG))
    Sv += tailq(lambda lam: (1 + log(lam/dK))/dK)
    # F(sigma): pieces of [lam0, oo)
    fin = [lam0] + [x for x in BRK if x > lam0]
    ND = []
    for u, v in zip(fin, fin[1:]): ND += [(l, w, G(l)) for (l, w) in nodes(u, v, NG)]
    G0 = G(lam0)
    ks = sorted(set([r['k'] for r in rows]))
    def Fsig(s):
        aa = mpf(8)*(1 - s)/3
        best = None
        # row-independent integrals: I_k = int_{[T0,oo)} (log kt)^{2s}(log t)^{5-4s} t^aa d(-G(dK t)), IB = int log^2 t d(-G)
        Ik = {}
        for k in ks:
            def f(t):
                lk, ll, lt = log(k*t), log(t), log(t)
                val = exp(2*s*log(lk) + (5 - 4*s)*log(ll) + aa*lt)
                return val, val*(2*s/(t*lk) + (5 - 4*s)/(t*ll) + aa/t)
            tot = f(T0)[0]*G0
            for (l, w, g) in ND: tot += w*g*f(l/dK)[1]/dK
            tot += tailq(lambda lam: f(lam/dK)[1]/dK)
            Ik[k] = tot
        IB = log(T0)**2*G0
        for (l, w, g) in ND: IB += w*g*2*log(l/dK)/(l/dK)/dK
        IB += tailq(lambda lam: 2*log(lam/dK)/(lam/dK)/dK)
        for r in rows:
            if not (s > r['smin']): continue
            v = exp(r['c0'] + r['c1']*s)*Ik[r['k']] + (r['q0'] - r['q1']*s)*IB
            if best is None or v < best[0]: best = (v, r['name'])
        return best
    # sigma-integral: 0.025-cells on [0.515, 0.80], 0.0025-cells on [0.80, B]
    edges = [sig0 + mpf(j)*mpf('0.025') for j in range(0, 12)] + [mpf('0.80') + mpf(j)*mpf('0.0025') for j in range(1, 30)]
    edges = [e for e in edges if e < Bh] + [Bh]
    edges = sorted(set(edges))
    integ = mpf(0); used = {}
    for u, v in zip(edges, edges[1:]):
        for (s, w) in nodes(u, v, NS):
            Fv, nm = Fsig(s)
            used[nm] = used.get(nm, 0) + 1
            integ += w*exp((s - 1)*xpow_L)*Fv
    Z = 2*(exp((sig0 - 1)*xpow_L)*Sv + xpow_L*integ)
    if pp_mode == 'cell':
        da = delta(La); db = delta(Lb); top = Lb + da
        J = int(floor(top/log(2)))
        PP = 2*top*(exp(da/2)*exp(-La/2)/4 + J*exp(da)*exp(-2*La/3)/9 + log(J)/(db*exp(La)))
    else:
        d0 = delta(La); top = La + d0
        J = int(floor(top/log(2)))
        PP = 2*top*(exp(d0/2)*exp(-La/2)/4 + J*exp(d0)*exp(-2*La/3)/9 + log(J)/(d0*exp(La)))
    return Z + PP + 2*exp(-3*xpow_L), dict(Z=Z, PP=PP, S_part=2*exp((sig0 - 1)*xpow_L)*Sv, lam0=lam0, D=D, rows_used=used)
La = mpf(191.82000000001713); Lb = mpf(191.84000000001714)
Rc, dc = R_bound(La, Lb, delta(Lb), La, 'cell')
print('NG=%d NS=%d  CELL [%s, %s]: x-powers at La, kernel at delta(Lb) = %s (lambda0 = %s): R_direct = %s  (Z = %s, PP = %s, S-part = %s); rows used %s' % (
    NG, NS, mp.nstr(La, 17), mp.nstr(Lb, 17), mp.nstr(delta(Lb), 10), mp.nstr(dc['lam0'], 8), mp.nstr(Rc, 12), mp.nstr(dc['Z'], 12), mp.nstr(dc['PP'], 6), mp.nstr(dc['S_part'], 4), dc['rows_used']))
Rp, dp = R_bound(La, La, delta(La), La, 'point')
print('NG=%d NS=%d  POINT x = e^La: delta = delta(La): R_direct = %s (Z = %s)' % (NG, NS, mp.nstr(Rp, 12), mp.nstr(dp['Z'], 12)))
print('time %.0fs' % (time.time() - t_start))
