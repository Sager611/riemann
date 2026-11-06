# test_matrix.py: checks of weil_matrix.py against independent numerics (mpmath quadrature), and the smallest
# eigenvalues against weilzeta's certified lower bounds. Usage: test_matrix.py entries c | eigs c N prec
import sys, time, mpmath as mp
from flint import arb, acb, arb_mat, ctx, fmpq
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from weil_matrix import WeilMat

def f2m(x): return mp.mpf(x.mid().str(60, radius=False))

def check_entries(cq, N, K0=200, pairs=None):
    ctx.prec = 128; mp.mp.dps = 15
    W = WeilMat(cq); L = f2m(W.L)
    phi = lambda t: -mp.log(mp.pi) + mp.re(mp.digamma(mp.mpf(1)/4 + 1j*t/2))
    out = []
    for sector in ('even', 'odd'):
        T, idx = W.block(N, sector)
        for a, i in enumerate(idx):
            for b, j in enumerate(idx):
                if b < a: continue
                if pairs is not None and (sector, i, j) not in pairs: continue
                Ni = (1/mp.sqrt(L) if i == 0 else mp.sqrt(2/L)); Nj = (1/mp.sqrt(L) if j == 0 else mp.sqrt(2/L))
                wi, wj = 2*mp.pi*i/L, 2*mp.pi*j/L
                if sector == 'even':
                    bi = lambda y: Ni*mp.cos(wi*y); bj = lambda y: Nj*mp.cos(wj*y); v = lambda y: mp.cosh(y/2); ps = 2
                    Fi = lambda t: Ni*(-1)**i*2*t*mp.sin(t*L/2)/(t*t - wi*wi)
                    Fj = lambda t: Nj*(-1)**j*2*t*mp.sin(t*L/2)/(t*t - wj*wj)
                else:
                    bi = lambda y: Ni*mp.sin(wi*y); bj = lambda y: Nj*mp.sin(wj*y); v = lambda y: mp.sinh(y/2); ps = -2
                    Fi = lambda t: Ni*(-1)**i*2*wi*mp.sin(t*L/2)/(t*t - wi*wi)
                    Fj = lambda t: Nj*(-1)**j*2*wj*mp.sin(t*L/2)/(t*t - wj*wj)
                # spot check of the closed-form transform against quadrature at one t
                t0 = mp.mpf('3.3')
                qi = mp.quad(lambda y: bi(y)*(mp.cos(t0*y) if sector == 'even' else mp.sin(t0*y)), [-L/2, 0, L/2])
                assert abs(qi - Fi(t0)) < 1e-10, (qi, Fi(t0))
                pole = ps*mp.quad(lambda y: bi(y)*v(y), [-L/2, L/2])*mp.quad(lambda y: bj(y)*v(y), [-L/2, L/2])
                Cq = lambda x: mp.quad(lambda y: bi(y)*bj(y - x), [x - L/2, L/2])
                pr = 0
                for n, ln, wt in W.atoms: pr += f2m(wt)*Cq(mp.log(n))
                T0 = K0*2*mp.pi/L
                pts = [k*mp.pi/L + mp.mpf('1e-9') for k in range(0, 2*K0 + 1)]
                ar = 2*mp.quad(lambda t: phi(t)*Fi(t)*Fj(t), pts)/(2*mp.pi)
                if sector == 'even':
                    Kc = Ni*Nj*(-1)**(i + j)*4
                    ar += 2*mp.quad(lambda t: phi(t)*Kc/2/t**2, [T0, mp.inf])/(2*mp.pi)
                direct = pole + ar - 2*pr
                out.append((sector, i, j, mp.nstr(direct, 10), mp.nstr(f2m(T[a, b]), 10), mp.nstr(direct - f2m(T[a, b]), 3)))
                print('entry', *out[-1], flush=True)
    return out

def eigs(cq, N, prec):
    ctx.prec = prec
    W = WeilMat(cq); res = {}
    for sector in ('even', 'odd'):
        T, idx = W.block(N, sector)
        mp.mp.dps = int(prec*0.3)
        A = mp.matrix([[mp.mpf(T[a, b].mid().str(int(prec*0.3), radius=False)) for b in range(len(idx))] for a in range(len(idx))])
        E = mp.eigsy(A, eigvals_only=True)
        res[sector] = sorted(E)[:4]
    return res

if __name__ == '__main__':
    what = sys.argv[1]
    if what == 'entries':
        check_entries(fmpq(int(sys.argv[2])), int(sys.argv[3]) if len(sys.argv) > 3 else 2)
    elif what == 'entries2':   # non-integer c, larger indices
        p, q = sys.argv[2].split('/')
        check_entries(fmpq(int(p), int(q)), 5, K0=150, pairs={('even', 0, 5), ('even', 3, 4), ('even', 5, 5), ('odd', 1, 5), ('odd', 4, 4), ('odd', 2, 5)})
    else:
        c, N, prec = int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); t0 = time.time()
        e = eigs(fmpq(c), N, prec)
        print('c', c, 'N', N, 'even', [mp.nstr(x, 8) for x in e['even']], 'odd', [mp.nstr(x, 8) for x in e['odd']], round(time.time() - t0, 1), 's', flush=True)
