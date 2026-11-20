# hdx_shortcert.py (explore64): rigorous grid argument for short links, window c >= 10.
# For c' in [g_j, g_{j+1}]:  e_k(c') <= e_k(g_j) < L_{j+1} <= Lam*(10, g_{j+1}) <= Lam*(10, c') <= Lam*(c, c')   (domain monotonicity of e_k;
# Lam* nonincreasing in c' and nondecreasing in c).  e_k(g_j) < L_{j+1} is certified by >= k negative pivots of A_N - L_{j+1} I (arb).
import sys, math
sys.path.insert(0, 'Riemann/Deep/ExploreHDX')
from flint import arb, fmpq, ctx
from hdx_ceiling import optimize, certify
import hdx_count

def lam(cq, clq, N=3000):
    u1, v1, t1, F1 = optimize(cq, clq, N, longonly=True)
    best = certify(cq, clq, u1, [t for t in t1 if t[2] == 1], N).lower() if u1 is not None else 0
    u2, v2, t2, F2 = optimize(cq, clq, N)
    if u2 is not None: best = max(float(best), float(certify(cq, clq, u2, t2, N).lower()))
    if u1 is not None: best = max(float(best), float(certify(cq, clq, u1, t2, N).lower()))
    return float(best)

if __name__ == '__main__':
    k = int(sys.argv[1]) if len(sys.argv) > 1 else 2
    grid = sys.argv[2].split(',') if len(sys.argv) > 2 else ['11/10', '6/5', '27/20', '3/2', '7/4', '2', '5/2', '3162277/1000000']
    g = [fmpq(*map(int, x.split('/'))) if '/' in x else fmpq(int(x)) for x in grid]
    L = [lam(fmpq(10), x) for x in g]
    print(f'# k = {k}: certify e_{k}(g_j) < L_(j+1) = arb lower bound of Lam*(10, g_(j+1)); N = 16 modes per sector')
    ok = True
    for j in range(len(g) - 1):
        th = fmpq(int(math.floor(L[j + 1]*10**6)), 10**6)
        ctx.prec = 200
        res = hdx_count.counts(g[j], 16, [th])
        ev, od = res['even'][0], res['odd'][0]
        neg = (ev[0] if ev else 0) + (od[0] if od else 0)
        ok &= neg >= k
        print(f"[{float(g[j]):.4f}, {float(g[j+1]):.4f}]  Lam*(10, right end) >= {L[j+1]:.4f};  #Galerkin eigenvalues of Q_left below it: {neg}  -> e_{k} < Lam* on the interval: {neg >= k}", flush=True)
    print('ALL OK' if ok else 'NOT ALL')
