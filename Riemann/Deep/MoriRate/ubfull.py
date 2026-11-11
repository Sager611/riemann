# ubfull.py (morirate): rigorous (Arb) upper bound for lambda_min of the full even Weil form at window c: Rayleigh quotient of the
# Galerkin ground state (wmat.py closed forms, Arb matrix) in a larger basis. Usage: ubfull.py c NMUL
import sys, os, json, math
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ExploreAntenna', 'code'))
from flint import arb, acb, arb_mat, ctx, fmpq
import mpmath as mp
import wmat

def main():
    c = float(sys.argv[1]); nmul = float(sys.argv[2])
    bits = int(4*math.pi*c*1.4427) + 260
    ctx.prec = bits; mp.mp.prec = bits; dps = int(bits*0.302) + 10
    F = wmat.Form(fmpq(int(round(c*64)), 64), 'zeta')
    N = max(24, int(math.ceil(nmul*2.7*c*math.log(c))) + 8)
    T, idx = F.block(N, 'even'); n = len(idx)
    M = wmat.to_mp(T, mp.mp)
    E, V = mp.mp.eigsy(M)
    i0 = min(range(n), key=lambda i: E[i])
    v = [arb(mp.nstr(V[r, i0], dps)) for r in range(n)]
    num = arb(0)
    for i in range(n):
        row = arb(0)
        for j in range(n): row += T[i, j]*v[j]
        num += v[i]*row
    den = sum((x*x for x in v), arb(0))
    X = num/den
    print(json.dumps({'c': c, 'N': N, 'rayleigh_full_even': X.str(10), 'upper': X.upper().str(8), 'galerkin_min': mp.nstr(E[i0], 10)}), flush=True)

if __name__ == '__main__':
    main()
