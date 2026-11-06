# dh.py (explore19): Galerkin spectrum of the Davenport-Heilbronn Weil form on window c (no pole, conductor 5, Gamma((s+1)/2)).
# Usage: dh.py c1 c2 ...   prints the 4 smallest eigenvalues of the even and odd blocks (mpmath eigsy on Arb midpoints).
import sys, time, math, json
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import arb, acb, arb_mat, ctx, fmpq
import mpmath as mp
import wmat

KIND = 'dh'
def run(c, Nfac=0.6, Nadd=20, extra_bits=128):
    bits = int(4*math.pi*c/5*1.4427) + extra_bits
    ctx.prec = bits; mp.mp.prec = bits
    N = int(math.ceil(Nfac*c*math.log(c))) + Nadd
    F = wmat.Form(fmpq(int(round(c*64)), 64), KIND)
    out = {'kind': KIND, 'c': c, 'N': N, 'bits': bits}
    for sec in ('even', 'odd'):
        t0 = time.time()
        T, idx = F.block(N, sec)
        M = wmat.to_mp(T, mp.mp)
        E = sorted(mp.mp.eigsy(M, eigvals_only=True))
        out[sec] = [mp.mp.nstr(e, 6) for e in E[:5]]
        out['s_' + sec] = round(time.time() - t0, 1)
    print(json.dumps(out), flush=True)

if __name__ == '__main__':
    args = sys.argv[1:]
    Nfac = 0.6
    if args and args[0].startswith('kind='):
        KIND = args[0][5:]; args = args[1:]
    if args and args[0].startswith('Nfac='):
        Nfac = float(args[0][5:]); args = args[1:]
    for a in args: run(float(a), Nfac=Nfac)
