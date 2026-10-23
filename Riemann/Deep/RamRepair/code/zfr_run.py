# zfr_run.py -- certified iteration (serial; outward conversions via rig.up/dn) for the per-height zero-free region (own version of Padhi's Thm 4.6).
# Table lam(L) on L in {28.75, 29.00, ..., 46.00}; lam(l) := lam(smallest grid L >= l); lam(l) := 1/5.558691 for l > 46
# (Mossinghoff-Trudgian-Yang 2024, Res. Number Theory 10, paper 11).  Step j assumes the step j-1 table for all heights,
# and certifies M1 > 0 on (l, lam0) in [max(L-0.25, log H0), L] x [lam_{j-1}(L), lam_j(L)] by ball bisection.
import json, math, sys, time, os
import zfr as Z
from flint import arb
from rig import up, dn
MTY = dn(1/arb('5.558691'))            # double <= 1/5.558691
Ls = [28.75 + 0.25*i for i in range(46)] + [40.0 + 0.5*i for i in range(1, 13)]     # 28.75..40 step 1/4, then to 46 step 1/2
PREV = {L: (Ls[i-1] if i > 0 else None) for i, L in enumerate(Ls)}
LOGH0 = dn(Z.LOGH0)                     # double <= log H0
def lookup(tab, l):
    for L in Ls:
        if L >= l: return tab[L]
    return MTY
def ltop(L):     # upper bound of log(K e^L + T0)
    return up((Z.K*arb(L).exp() + Z.T0).log())
def froot(l, A, lo, hi):
    Lt = arb(l); Aa = arb(A)
    if not Z.margin1(Lt, arb(lo), Aa).mid() > 0: return None
    for _ in range(30):
        mid = (lo + hi)/2
        if Z.margin1(Lt, arb(mid), Aa).mid() > 0: lo = mid
        else: hi = mid
    return lo
def certify(l0, l1, x0, x1, A, maxn=20000):
    Aa = arb(A); stack = [(l0, l1, x0, x1)]; n = 0
    while stack:
        a, b, c, d = stack.pop(); n += 1
        if n > maxn: return False, n
        m = Z.margin1(Z.ball(a, b), Z.ball(c, d), Aa)
        if m > 0: continue
        if b - a < 1e-7 and d - c < 1e-9: return False, n
        # split the dimension with the larger estimated effect (l: dM/dl ~ 0.1, lam: dM/dlam ~ 3)
        if (b - a)*0.1 > (d - c)*3: mid = (a + b)/2; stack += [(a, mid, c, d), (mid, b, c, d)]
        else: mid = (c + d)/2; stack += [(a, b, c, mid), (a, b, mid, d)]
    return True, n
def job(args):
    L, prev, A, gap = args
    l0 = PREV[L] if PREV[L] is not None else LOGH0      # cell (previous grid point, L]
    l0 = max(l0, LOGH0)
    if L <= LOGH0: return L, prev, 0, 'below H0'
    r1 = froot(L, A, prev, 0.30); r0 = froot(l0, A, prev, 0.30)
    if r1 is None or r0 is None: return L, prev, 0, 'no float root above prev'
    cand = min(r0, r1) - gap
    tries = 0; ntot = 0
    while cand > prev:
        ok, n = certify(l0, L, prev, cand, A); ntot += n; tries += 1
        if ok: return L, cand, ntot, 'ok after %d tries' % tries
        cand -= max(gap, 1e-4)
    return L, prev, ntot, 'kept prev'
if __name__ == "__main__":
    nsteps = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    gap = float(sys.argv[2]) if len(sys.argv) > 2 else 1e-4
    out = sys.argv[3] if len(sys.argv) > 3 else 'out/zfr_table.json'
    tab = {L: MTY for L in Ls}
    hist = []
    t0 = time.time()
    if True:                                 # serial (one process)
        for j in range(1, nsteps + 1):
            jobs = [(L, tab[L], lookup(tab, ltop(L)), gap) for L in Ls]
            res = sorted(map(job, jobs))
            new = {L: v for L, v, n, msg in res}
            # enforce lam non-increasing in l (reduce later values); keeps new >= old
            run = 1.0
            for L in Ls:
                run = min(run, new[L]); new[L] = run
            assert all(new[L] >= tab[L] for L in Ls)
            ch = max(new[L] - tab[L] for L in Ls)
            tab = new
            hist.append({'step': j, 'max_change': ch, 'boxes': sum(n for _, _, n, _ in res), 'msgs': sorted(set(m for *_, m in res))})
            print('step %d  max change %.2e  %.0fs' % (j, ch, time.time() - t0), [(L, round(tab[L], 5)) for L in Ls[::4]], flush=True)
            json.dump({'table': {repr(L): repr(v) for L, v in tab.items()}, 'above': repr(MTY), 'hist': hist}, open(out, 'w'), indent=0)
            if ch < 2e-5: break
