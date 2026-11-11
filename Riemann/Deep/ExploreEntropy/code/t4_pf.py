# t4_pf.py (explore32): Perron-Frobenius structure. Q0 = Weil form minus its pole part = Dirichlet form (jump kernel k(|x|) + prime atoms) + const.
# (A) zeta Q0: ground state simple, even, positive?  (B) full Q (CCM's QW_lambda): ground state sign structure, overlap with Q0's Perron vector.
# (C) L(s, chi_5) control and Davenport-Heilbronn (no pole; atoms of both signs): ground state sign structure.
import sys, math, json, time, numpy as np, mpmath as mp
from scipy.optimize import brentq
from flint import arb, ctx, fmpq
import wmat
part = sys.argv[1]
def fx(v, L, sector, xs):
    if sector == 'even':
        B = [np.full_like(xs, 1/np.sqrt(L))] + [np.sqrt(2/L)*np.cos(2*np.pi*j*xs/L) for j in range(1, len(v))]
    else:
        B = [np.sqrt(2/L)*np.sin(2*np.pi*j*xs/L) for j in range(1, len(v) + 1)]
    return np.array(B).T @ np.array(v)
def signinfo(f, xs, a, sector):
    m = (np.abs(xs) <= 0.98*a) & ((xs > 0.02*a) if sector == 'odd' else True)
    g = f[m]; g = g/np.max(np.abs(g))
    if g[np.argmax(np.abs(g))] < 0: g = -g
    sc = int(np.sum(np.diff(np.sign(g[np.abs(g) > 1e-9])) != 0))
    return float(g.min()), sc
def block_float(F, N, sector, pole):
    F.pole = pole
    T, idx = F.block(N, sector)
    n = len(idx); M = np.array([[float(T[i, j].mid()) for j in range(n)] for i in range(n)])
    return M
if part == 'x0':
    k = lambda x: math.exp(x/2)/(math.exp(x) - math.exp(-x))
    x0 = brentq(lambda x: 2*math.cosh(x/2) - k(x), 0.05, 2)
    print('full-Q kernel 2cosh(x/2) - k(x) changes sign at |x - y| = x0 = %.6f: PF window L < x0, i.e. c < %.6f, half-width < %.6f' % (x0, math.exp(x0), x0/2))
if part == 'A':
    ctx.prec = 128
    for c in [2, 5, 10, 20, 47]:
        F = wmat.Form(fmpq(c), 'zeta'); L = math.log(c); a = L/2; N = 80
        xs = np.linspace(-a, a, 4001); out = {'c': c, 'N': N}
        for sec in ['even', 'odd']:
            M = block_float(F, N, sec, False)
            E, V = np.linalg.eigh(M)
            out[sec] = [round(float(e), 6) for e in E[:3]]
            if sec == 'even':
                f = fx(V[:, 0], L, sec, xs); mn, sc = signinfo(f, xs, a, sec)
                ch = np.cosh(xs/2); ov = np.trapezoid(f*ch, xs)
                out['perron_min_over_max'] = round(mn, 6); out['perron_sign_changes'] = sc; out['<perron,cosh>'] = round(float(abs(ov)), 6)
        print(json.dumps(out), flush=True)
if part in ('B', 'C'):
    kinds = [('zeta', [2, 3, 5, 10, 20])] if part == 'B' else [('chi5', [5, 10, 20, 31]), ('dh', [5, 10, 20, 31])]
    for kind, cs in kinds:
        for c in cs:
            cond = 1 if kind == 'zeta' else 5
            bits = int(4*math.pi*c/cond*1.4427) + 128
            ctx.prec = bits; mp.mp.prec = bits
            N = int(math.ceil(0.6*c*math.log(c))) + 20
            F = wmat.Form(fmpq(c), kind); L = math.log(c); a = L/2
            xs = np.linspace(-a, a, 4001); out = {'kind': kind, 'c': c, 'N': N, 'bits': bits}
            t0 = time.time()
            for sec in ['even', 'odd']:
                T, idx = F.block(N, sec)
                M = wmat.to_mp(T, mp.mp)
                E, Q = mp.mp.eigsy(M)
                order = sorted(range(len(idx)), key=lambda i: E[i])
                out[sec] = [mp.mp.nstr(E[i], 5) for i in order[:3]]
                v = [float(Q[r, order[0]]) for r in range(len(idx))]
                f = fx(v, L, sec, xs); mn, sc = signinfo(f, xs, a, sec)
                out[sec + '_gs_min_over_max'] = round(mn, 5); out[sec + '_gs_sign_changes'] = sc
                if kind == 'zeta' and sec == 'even':
                    try:
                        F2 = wmat.Form(fmpq(c), 'zeta'); ctx.prec = 128
                        M0 = block_float(F2, N, 'even', False); E0, V0 = np.linalg.eigh(M0)
                        out['overlap_with_Q0_perron'] = round(float(abs(np.dot(V0[:, 0], v)/np.linalg.norm(v))), 6)
                    finally:
                        ctx.prec = bits
            out['sec'] = round(time.time() - t0, 1)
            print(json.dumps(out), flush=True)

if part == 'P':
    fr = [0, 0.25, 0.5, 0.7, 0.8, 0.9, 0.95, 0.98, 0.995, 1.0]
    for c in [float(x) for x in sys.argv[2:]]:
        bits = int(4*math.pi*c*1.4427) + 128
        ctx.prec = bits; mp.mp.prec = bits
        N = int(math.ceil(0.6*c*math.log(c))) + 20
        if len(sys.argv) > 1 and 'NMUL' in __import__('os').environ: N = int(N*float(__import__('os').environ['NMUL']))
        F = wmat.Form(fmpq(int(round(c*64)), 64), 'zeta'); L = math.log(c); a = L/2
        xs = np.linspace(-a, a, 20001)
        for sec in ['even', 'odd']:
            T, idx = F.block(N, sec)
            E, Q = mp.mp.eigsy(wmat.to_mp(T, mp.mp))
            i0 = min(range(len(idx)), key=lambda i: E[i])
            v = [float(Q[r, i0]) for r in range(len(idx))]
            f = fx(v, L, sec, xs); f = f/np.max(np.abs(f))
            if f[np.argmax(np.abs(f))] < 0: f = -f
            prof = [round(float(np.interp(t*a, xs, f)), 6) for t in fr]
            inner = np.abs(xs) <= 0.9*a
            if sec == 'odd': inner &= xs > 0.05*a
            print(json.dumps({'c': c, 'N': N, 'sec': sec, 'eig': mp.mp.nstr(E[i0], 5), 'profile_at_x/a': dict(zip(fr, prof)), 'min_|x|<=0.9a': float(f[inner].min()), 'argmin_x/a': round(float(xs[inner][np.argmin(f[inner])]/a), 4)}), flush=True)
