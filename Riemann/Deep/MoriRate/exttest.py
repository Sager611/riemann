# exttest.py (morirate): exterior test. Delete every s-th zero with y = gamma/(2 pi c) in [y1, y2] (the deleted labels become
# "merged" into neighbours: lost samples). Single deletions above y = 2 cost ~0 (lost-sample law); a whole band of them makes the distinct
# points sub-Nyquist ((1/2pi)log(t/2pi)/s' < (1/2pi)log c) and should cost a lot: (EXT)-type information is necessary, not only (CENT).
# Usage: exttest.py c NMUL y1 y2 stride sector
import sys, os, json, math, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'ExploreAntenna', 'code'))
from flint import arb, acb, arb_mat, ctx, fmpq
import mpmath as mp
import wmat

def g2(y):
    y = abs(y); return math.log((2 + math.sqrt(4 - y*y))/y) if y < 2 else 0.0

def main():
    c = float(sys.argv[1]); nmul = float(sys.argv[2]); y1 = float(sys.argv[3]); y2 = float(sys.argv[4]); stride = int(sys.argv[5]); sector = sys.argv[6]
    bits = int(4*math.pi*c*1.4427) + int(os.environ.get("EXTRA_BITS", "190"))
    ctx.prec = bits; mp.mp.prec = bits
    dps = int(bits*0.302) + 10
    here = os.path.dirname(os.path.abspath(__file__))
    zs = []
    for f in sorted([f for f in os.listdir(here) if f.startswith('zeros_dps')], key=lambda f: -int(f[9:-4])):
        if int(f[9:-4]) >= dps:
            zs = [l.strip() for l in open(os.path.join(here, f)) if l.strip()]; break
    T2 = y2*2*math.pi*c
    if not zs or float(mp.mpf(zs[-1])) < T2:
        raise SystemExit('need zeros to %g at dps %d' % (T2, dps))
    mp.mp.prec = bits
    F = wmat.Form(fmpq(int(round(c*64)), 64), 'zeta')
    N = max(24, int(math.ceil(nmul*2.7*c*math.log(c))) + 8)
    t0 = time.time()
    T, idx = F.block(N, sector)
    M = wmat.to_mp(T, mp.mp)
    n = len(idx)
    ft = F.ft_even if sector == 'even' else F.ft_odd
    def bvec(zstr):
        g = acb(arb(zstr)); out = []
        for j in idx:
            v = ft(j, g); out.append(mp.mp.mpf((v.real if sector == 'even' else v.imag).mid().str(dps, radius=False)))
        return out
    E, V = mp.mp.eigsy(M)
    Ev = [E[i] for i in range(n)]; E1 = min(Ev)
    print(json.dumps({'c': c, 'N': N, 'wmax': float(2*math.pi*N/math.log(c)), 'sector': sector, 'E1': mp.nstr(E1, 8), 't': round(time.time()-t0, 1)}), flush=True)
    sel = []
    ks = [k for k, z in enumerate(zs, 1) if y1*2*math.pi*c <= float(mp.mpf(z)) <= T2]
    sel = ks[::stride] if stride > 1 else ks
    single = []
    for k in sel:
        b = bvec(zs[k-1])
        cv = [mp.mp.fsum(V[r, i]*b[r] for r in range(n)) for i in range(n)]
        c2 = [2*cv[i]**2 for i in range(n)]
        def phi(x): return mp.mp.fsum(c2[i]/(Ev[i] - x) for i in range(n))
        lo, hi = mp.mp.mpf(0), E1
        if phi(lo) >= 1: xs = lo
        else:
            for _ in range(200):
                mid = (lo + hi)/2
                if phi(mid) < 1: lo = mid
                else: hi = mid
                if hi - lo < E1*mp.mp.mpf(2)**(-30): break
            xs = (lo + hi)/2
        y = float(mp.mp.mpf(zs[k-1]))/(2*math.pi*c)
        lr = float(mp.mp.log(E1/xs)) if xs > 0 else float('inf')
        single.append(lr)
        print(json.dumps({'k': k, 'y': round(y, 4), 'single_logratio': round(lr, 4), 'pred4g2': round(4*g2(y), 4), 'b': None}), flush=True)
    A = M.copy()
    merge = os.environ.get('MERGE', '0') == '1'
    for k in sel:
        b = bvec(zs[k-1])
        for i in range(n):
            for j in range(n): A[i, j] -= 2*b[i]*b[j]
        if merge and k > 1:
            # the deleted label joins its lower neighbour: double zero at gamma_{k-1} (weight 2 -> 4)
            bp = bvec(zs[k-2])
            for i in range(n):
                for j in range(n): A[i, j] += 2*bp[i]*bp[j]
    E2 = mp.mp.eigsy(A, eigvals_only=True)
    l = min(E2[i] for i in range(n))
    print(json.dumps({'merge': merge, 'set_size': len(sel), 'band_y': [y1, y2], 'stride': stride, 'lmin_set': mp.nstr(l, 8), 'logratio_set': round(float(mp.mp.log(E1/l)), 4), 'sum_single': round(sum(single), 4), 't': round(time.time()-t0, 1)}), flush=True)

if __name__ == '__main__':
    main()
