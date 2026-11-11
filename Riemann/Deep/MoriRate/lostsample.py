# lostsample.py (morirate): bottom of the Weil form after deleting one zero (a "lost sample") from the zero side.
# Q = Galerkin matrix of the zeta Weil form on [-L/2, L/2], L = ln c (c = mu), CCM Fourier basis (explore19 wmat.py, closed forms).
# Under RH (true to these heights) Q(f) = sum_{gamma>0} 2|F(gamma)|^2; deleting gamma_k gives Q - 2 b b^T, b_j = F_j(gamma_k).
# Smallest eigenvalue x* of Q - 2bb^T: root of 2 sum_i c_i^2/(E_i - x) = 1 on [0, E_1), c = V^T b (secular equation).
# Prediction (potential theory, Mori Thm 4 dual certificate): log(E_1/x*) ~ 2 g2(y), y = gamma/(2 pi c), g2(y) = log((2+sqrt(4-y^2))/y) (y<2), 0 (y>=2).
# Usage: lostsample.py c [ymax] [NMUL]
import sys, os, json, math, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ExploreAntenna', 'code'))
from flint import arb, acb, arb_mat, ctx, fmpq
import mpmath as mp
import wmat

def g2(y):
    y = abs(y)
    return math.log((2 + math.sqrt(4 - y*y))/y) if y < 2 else 0.0

def zeros_upto(T, dps):
    fn = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'zeros_dps%d.txt' % dps)
    zs = []
    if os.path.exists(fn):
        zs = [l.strip() for l in open(fn) if l.strip()]
    mp.mp.dps = dps
    while not zs or mp.mpf(zs[-1]) < T:
        z = mp.zetazero(len(zs) + 1).imag
        zs.append(mp.nstr(z, dps))
        with open(fn, 'a') as f: f.write(zs[-1] + '\n')
    return [z for z in zs if mp.mpf(z) <= T]

def main():
    c = float(sys.argv[1]); ymax = float(sys.argv[2]) if len(sys.argv) > 2 else 3.0
    nmul = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
    sectors = (sys.argv[4],) if len(sys.argv) > 4 else ('even', 'odd')
    bits = int(4*math.pi*c*1.4427) + 200
    ctx.prec = bits; mp.mp.prec = bits
    dps = int(bits*0.302) + 10
    t0 = time.time()
    zs = zeros_upto(ymax*2*math.pi*c, dps)
    mp.mp.prec = bits
    cq = fmpq(int(round(c*64)), 64)
    F = wmat.Form(cq, 'zeta')
    N = max(int(24*nmul), int(math.ceil(nmul*2.7*c*math.log(c))) + 8)
    print(json.dumps({'c': float(cq), 'N': N, 'bits': bits, 'nzeros': len(zs), 't_zeros': round(time.time()-t0, 1)}), flush=True)
    for sector in sectors:
        t1 = time.time()
        T, idx = F.block(N, sector)
        M = wmat.to_mp(T, mp.mp)
        E, V = mp.mp.eigsy(M)
        n = len(idx)
        Ev = [E[i] for i in range(n)]
        order = sorted(range(n), key=lambda i: Ev[i])
        E1 = Ev[order[0]]
        print(json.dumps({'sector': sector, 'E1': mp.nstr(E1, 8), 'E2': mp.nstr(Ev[order[1]], 8), 'E3': mp.nstr(Ev[order[2]], 8), 't_eig': round(time.time()-t1, 1)}), flush=True)
        ft = F.ft_even if sector == 'even' else F.ft_odd
        for k, zstr in enumerate(zs, 1):
            g = acb(arb(zstr))
            b = []
            for j in idx:
                v = ft(j, g)
                b.append(mp.mp.mpf((v.real if sector == 'even' else v.imag).mid().str(dps, radius=False)))
            cv = [mp.mp.fsum(V[r, i]*b[r] for r in range(n)) for i in range(n)]
            c2 = [2*cv[i]**2 for i in range(n)]
            def phi(x): return mp.mp.fsum(c2[i]/(Ev[i] - x) for i in range(n))
            phi0 = phi(0)
            # bisection in log scale on (0, E1): phi increasing
            if phi0 >= 1:
                xs = mp.mp.mpf(0)
            else:
                lo, hi = mp.mp.mpf(0), E1
                for _ in range(200):
                    mid = (lo + hi)/2
                    if phi(mid) < 1: lo = mid
                    else: hi = mid
                    if hi - lo < E1*mp.mp.mpf(2)**(-40): break
                xs = (lo + hi)/2
            y = float(mp.mp.mpf(zstr))/(2*math.pi*c)
            lr = float(mp.mp.log(E1/xs)) if xs > 0 else float('inf')
            # weight of the ground state in b: share of 2c_1^2/E1 in phi(0)
            w1 = float(c2[order[0]]/E1/phi0) if phi0 > 0 else 0.0
            print(json.dumps({'sector': sector, 'k': k, 'gamma': float(mp.mp.mpf(zstr)), 'y': round(y, 5), 'phi0': float(phi0), 'logratio': round(lr, 5), 'pred2g2': round(2*g2(y), 5), 'w1': round(w1, 4)}), flush=True)
        print(json.dumps({'sector': sector, 't_sector': round(time.time()-t1, 1)}), flush=True)

if __name__ == '__main__':
    main()
