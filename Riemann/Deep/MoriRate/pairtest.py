# pairtest.py (morirate): (a) a pair at distance delta: zero gamma_{k+1} moved to gamma_k + delta (delta -> 0 = double zero);
# (b) several zeros deleted at once (additivity of the lost-sample law). Galerkin matrix as in lostsample.py (explore19 wmat.py).
# Usage: pairtest.py c sector mode   (mode = pair:k  or  sets)
import sys, os, json, math, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ExploreAntenna', 'code'))
from flint import arb, acb, arb_mat, ctx, fmpq
import mpmath as mp
import wmat

def g2(y):
    y = abs(y); return math.log((2 + math.sqrt(4 - y*y))/y) if y < 2 else 0.0

def main():
    c = float(sys.argv[1]); sector = sys.argv[2]; mode = sys.argv[3]
    extra = int(os.environ.get('EXTRA_BITS', '260'))
    bits = int(4*math.pi*c*1.4427) + extra
    ctx.prec = bits; mp.mp.prec = bits
    dps = int(bits*0.302) + 10
    zfile = [f for f in os.listdir(os.path.dirname(os.path.abspath(__file__)) or '.') if f.startswith('zeros_dps')]
    zs = []
    for f in sorted(zfile, key=lambda f: -int(f[9:-4])):
        if int(f[9:-4]) >= dps:
            zs = [l.strip() for l in open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f)) if l.strip()]; break
    if not zs:
        mp.mp.dps = dps; zs = [mp.nstr(mp.zetazero(i).imag, dps) for i in range(1, 60)]
    mp.mp.prec = bits
    cq = fmpq(int(round(c*64)), 64)
    F = wmat.Form(cq, 'zeta')
    N = max(24, int(math.ceil(2.7*c*math.log(c))) + 8)
    T, idx = F.block(N, sector)
    M = wmat.to_mp(T, mp.mp)
    ft = F.ft_even if sector == 'even' else F.ft_odd
    def bvec(zstr):
        g = acb(arb(zstr)); out = []
        for j in idx:
            v = ft(j, g); out.append(mp.mp.mpf((v.real if sector == 'even' else v.imag).mid().str(dps, radius=False)))
        return out
    n = len(idx)
    def lmin(mods):
        A = M.copy()
        for sgn, b in mods:
            for i in range(n):
                for j in range(n): A[i, j] += 2*sgn*b[i]*b[j]
        E = mp.mp.eigsy(A, eigvals_only=True)
        return min(E[i] for i in range(n))
    t0 = time.time()
    E0 = lmin([])
    print(json.dumps({'c': c, 'sector': sector, 'N': N, 'bits': bits, 'E1': mp.nstr(E0, 10), 't': round(time.time()-t0, 1)}), flush=True)
    if mode.startswith('pair'):
        k = int(mode.split(':')[1])
        gk, gk1 = mp.mp.mpf(zs[k-1]), mp.mp.mpf(zs[k])
        bk1 = bvec(zs[k])
        ldel = lmin([(-1, bk1)])
        print(json.dumps({'k': k, 'gk': float(gk), 'gk1': float(gk1), 'y': float(gk1)/(2*math.pi*c), 'deleted': mp.nstr(ldel, 10), 'logratio_del': float(mp.mp.log(E0/ldel)), 'pred4g2': 4*g2(float(gk1)/(2*math.pi*c))}), flush=True)
        gap = gk1 - gk
        deltas = [gap, gap/2, gap/4, mp.mp.mpf('0.5'), mp.mp.mpf('0.1')] + [mp.mp.mpf(10)**(-e) for e in range(2, int(4*math.pi*c/2/2.3026) + 8, 2)] + [mp.mp.mpf(0)]
        for d in deltas:
            bd = bvec(mp.nstr(gk + d, dps))
            l = lmin([(-1, bk1), (1, bd)])
            print(json.dumps({'delta': mp.nstr(d, 6), 'lmin': mp.nstr(l, 10), 'ratio_to_orig': float(l/E0), 'ratio_to_del': float(l/ldel), 'excess_over_del': mp.nstr(l - ldel, 6)}), flush=True)
    else:
        sets = [[1], [2], [1, 2], [1, 2, 3], [1, 2, 3, 4, 5], [2, 4, 6, 8, 10], list(range(1, 11)), [5], [5, 6], [10], [10, 11], [15], [15, 16, 17], list(range(2, 21, 2))]
        cache = {}
        for S in sets:
            mods = []
            pred = 0.0
            for k in S:
                if k not in cache: cache[k] = bvec(zs[k-1])
                mods.append((-1, cache[k])); pred += 4*g2(float(mp.mp.mpf(zs[k-1]))/(2*math.pi*c))
            l = lmin(mods)
            print(json.dumps({'set': S, 'lmin': mp.nstr(l, 8), 'logratio': round(float(mp.mp.log(E0/l)), 4), 'sum4g2': round(pred, 4)}), flush=True)

if __name__ == '__main__':
    main()
