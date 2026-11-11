# kstar.py (morirate): crossover distance delta* of a close pair. Model: zero gamma_{k+1} moved to gamma_k + delta. For small delta
# lambda(delta) = lambda_del + K delta^2, K = 2 (F0'(gamma_k))^2 with F0 the ground state of Q - 2 b_{k+1} b_{k+1}^T (normalised);
# delta* = sqrt((lambda_orig - lambda_del)/K). One eigendecomposition per c. Usage: kstar.py ytarget c1 c2 ...
import sys, os, json, math, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ExploreAntenna', 'code'))
from flint import arb, acb, arb_mat, ctx, fmpq
import mpmath as mp
import wmat

def main():
    ytar = float(sys.argv[1])
    here = os.path.dirname(os.path.abspath(__file__))
    for c in [float(a) for a in sys.argv[2:]]:
        bits = int(4*math.pi*c*1.4427) + 190
        ctx.prec = bits; mp.mp.prec = bits; dps = int(bits*0.302) + 10
        zs = []
        for f in sorted([f for f in os.listdir(here) if f.startswith('zeros_dps')], key=lambda f: -int(f[9:-4])):
            if int(f[9:-4]) >= dps: zs = [l.strip() for l in open(os.path.join(here, f)) if l.strip()]; break
        mp.mp.prec = bits
        F = wmat.Form(fmpq(int(round(c*64)), 64), 'zeta')
        N = max(24, int(math.ceil(2.7*c*math.log(c))) + 8)
        T, idx = F.block(N, 'even'); M = wmat.to_mp(T, mp.mp); n = len(idx)
        # pick k with y_{k+1} closest to ytar
        k = min(range(1, len(zs)-1), key=lambda k: abs(float(mp.mp.mpf(zs[k]))/(2*math.pi*c) - ytar))
        def bvec(z):
            g = acb(arb(z)) if isinstance(z, str) else z
            return [mp.mp.mpf(F.ft_even(j, g).real.mid().str(dps, radius=False)) for j in idx]
        b1 = bvec(zs[k])
        E0 = mp.mp.eigsy(M, eigvals_only=True); lo = min(E0[i] for i in range(n))
        A = M.copy()
        for i in range(n):
            for j in range(n): A[i, j] -= 2*b1[i]*b1[j]
        E, V = mp.mp.eigsy(A)
        i0 = min(range(n), key=lambda i: E[i]); ld = E[i0]
        v0 = [V[r, i0] for r in range(n)]
        # derivative of F_j at gamma_k by complex step
        h = mp.mp.mpf(2)**(-bits//2)
        g = mp.mp.mpf(zs[k-1])
        dv = [mp.mp.mpf((F.ft_even(j, acb(arb(mp.nstr(g, dps)), arb(mp.nstr(h, 30)))).imag/arb(mp.nstr(h, 30))).mid().str(dps, radius=False)) for j in idx]
        Fp = mp.mp.fsum(v0[r]*dv[r] for r in range(n))
        F0 = mp.mp.fsum(v0[r]*bvec(zs[k-1])[r] for r in range(n))
        K = 2*Fp**2
        ds = mp.mp.sqrt((lo - ld)/K)
        print(json.dumps({'c': c, 'k': k, 'gk': float(g), 'y': float(mp.mp.mpf(zs[k]))/(2*math.pi*c), 'lam_orig': mp.nstr(lo, 6), 'lam_del': mp.nstr(ld, 6),
                          'K': mp.nstr(K, 6), 'F0_at_gk': mp.nstr(F0, 4), 'delta_star': mp.nstr(ds, 6), 'log_delta_star': float(mp.mp.log(ds)), 'log_sqrt_lam': float(mp.mp.log(lo)/2)}), flush=True)

if __name__ == '__main__':
    main()
