# count.py (explore19): near-null counts (eigenvalues < 1e-2 and < 1e-6) per parity vs the Landau excess c - 7/8.
import sys, os, math, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from flint import ctx, fmpq
import mpmath as mp
import wmat
for c in [float(a) for a in sys.argv[1:]]:
    bits = int(4*math.pi*c*1.4427) + 96; ctx.prec = bits; mp.mp.prec = bits
    N = max(24, int(math.ceil(2.7*c*math.log(c))) + 8)
    F = wmat.Form(fmpq(int(round(c*64)), 64), 'zeta')
    out = {'c': c, 'N': N, 'landau_per_parity': round(c - 7/8, 3)}
    for sec in ('even', 'odd'):
        T, idx = F.block(N, sec)
        E = sorted(mp.mp.eigsy(wmat.to_mp(T, mp.mp), eigvals_only=True))
        out[sec + '_lt_1e-2'] = sum(1 for e in E if e < 1e-2)
        out[sec + '_lt_1e-6'] = sum(1 for e in E if e < 1e-6)
        out[sec + '_first'] = [mp.mp.nstr(e, 3) for e in E[:3]]
    print(json.dumps(out), flush=True)
