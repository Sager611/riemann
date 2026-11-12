# theta_rig.py (lostlaw): rigorous upper bound for 1 - theta(g_k), theta = 2<e0, A^-1 e0> = sup_u 2<u,e0>^2/Q(u), via the
# Galerkin test vector u = R^-1 b (any u gives theta >= 2(b.u)^2/Q(u)); Arb balls for Q's matrix (galerkin_fast closed forms) and
# for b (zero ball certified by a sign change of Hardy's Z). Also w0 = 2(b.x)^2/E1 (midpoint, ground state x) and 4 g2(y).
# Usage: theta_rig.py c N bits k1,k2,...
import sys, os, json, math, random
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from flint import arb, arb_mat, ctx
import galerkin_fast as gf
from gal import zero_ball, bcol, dot, midv, col, tomat, matvec, gs, ritz, g2

def main():
    c = int(sys.argv[1]); N = int(sys.argv[2]); ctx.prec = int(sys.argv[3]); ks = [int(t) for t in sys.argv[4].split(',')]
    zs = [l.strip() for l in open(os.path.join(HERE, 'zeros_dps179.txt')) if l.strip()]
    L = arb(c).log(); M = gf.build(c, N)
    B = {k: bcol(L, N, arb(zs[k-1])) for k in ks}
    random.seed(1)
    X0 = [[arb(random.uniform(-1, 1)) for _ in range(N+1)] for _ in range(2)]
    Y = M.solve(tomat(X0 + [B[k] for k in ks]), algorithm="approx")
    X1 = gs([midv(col(Y, j)) for j in range(2)]); Z = {k: midv(col(Y, 2 + i)) for i, k in enumerate(ks)}
    Y = M.solve(tomat(X1), algorithm="approx"); X2 = gs([midv(col(Y, j)) for j in range(2)])
    E1, x, _ = ritz(M, X2, [], [])
    E1a = dot(x, [t.mid() for t in matvec(M, x)]).mid()
    rexp = min(170, int(ctx.prec*0.30) - 10)
    for k in ks:
        gb, ok = zero_ball(zs[k-1], rexp)
        bb = bcol(L, N, gb, mid=False)
        u = Z[k]
        th = 2*dot(bb, u)**2/dot(u, matvec(M, u))
        om = 1 - th
        w0 = 2*dot(B[k], x).mid()**2/E1a
        y = float(zs[k-1])/(2*math.pi*c)
        print(json.dumps({'c': c, 'N': N, 'k': k, 'y': round(y, 5), 'zero_ball_certified': ok,
                          'one_minus_theta_upper': om.upper().str(6, radius=False), 'ln_upper': float(om.upper().log().mid()) if om.upper() > 0 else None,
                          'ln_w0': float(w0.log().mid()), 'm4g2': -4*g2(y)}), flush=True)

if __name__ == '__main__':
    main()
