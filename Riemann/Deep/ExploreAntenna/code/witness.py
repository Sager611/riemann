# witness.py (explore19): rigorous negativity witness for the Davenport-Heilbronn Weil form on window c.
# v = Galerkin ground vector (even block), rounded to dyadics; Q(f_v) = v^T T v enclosed in Arb (T: closed-form Arb entries).
# Also prints F_v(t) sign changes on a real grid (zeros of the pattern) near the DH off-line zero 85.699 and the values at the DH pair.
import sys, math, json
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, ctx, fmpq, fmpz
import mpmath as mp
import wmat

def main(c, kind='dh', Nfac=0.6, sector='even'):
    bits = int(4*math.pi*c/5*1.4427) + 160
    ctx.prec = bits; mp.mp.prec = bits
    N = int(math.ceil(Nfac*c*math.log(c))) + 20
    cq = fmpq(int(round(c*64)), 64)
    F = wmat.Form(cq, kind)
    T, idx = F.block(N, sector)
    M = wmat.to_mp(T, mp.mp)
    E, Q = mp.mp.eigsy(M)
    i0 = min(range(len(idx)), key=lambda i: E[i])
    v = [Q[r, i0] for r in range(len(idx))]
    s = max(abs(x) for x in v)
    vq = [fmpq(int(mp.mp.nint(x/s*mp.mp.mpf(2)**80)), 2**80) for x in v]
    va = [arb(x) for x in vq]
    num = arb(0)
    for i in range(len(idx)):
        acc = arb(0)
        for j in range(len(idx)): acc += T[i, j]*va[j]
        num += va[i]*acc
    den = sum((x*x for x in va), arb(0))
    rq = num/den
    out = {'kind': kind, 'c': float(cq), 'N': N, 'bits': bits, 'eig_min_mid': mp.mp.nstr(E[i0], 8), 'RQ_ball': rq.str(10, radius=True),
           'RQ_upper_negative': bool(rq < 0)}
    # pattern F_v(t) = sum v_j F_j(t), real t grid: sign changes in [70, 100]
    def Fv(z):
        return sum((va[k]*F.ft_even(idx[k], z) for k in range(len(idx))), acb(0))
    grid = [70 + 0.05*i for i in range(601)]
    vals = [Fv(acb(arb(t))).real.mid() for t in grid]
    zc = [round(grid[i] + 0.025, 3) for i in range(len(grid)-1) if (vals[i] > 0) != (vals[i+1] > 0)]
    out['pattern_real_zeros_70_100'] = zc
    if kind == 'dh':
        z = acb(arb('85.699348'), arb('0.308517'))
        w = Fv(z); out['F_at_DH_pair'] = w.str(6, radius=False)
        out['Re_F2_at_pair_over_norm'] = (w*w).real.mid().str(6) if False else str((w*w).real.mid())[:30]
    print(json.dumps(out), flush=True)

if __name__ == '__main__':
    c = float(sys.argv[1]); kind = sys.argv[2] if len(sys.argv) > 2 else 'dh'
    Nfac = float(sys.argv[3]) if len(sys.argv) > 3 else 0.6
    main(c, kind, Nfac)
