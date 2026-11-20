# p4_zetaprime.py (explore57): pointwise check of (P4): phi'(t) = -theta'(t) - Re zeta''/zeta'(1/2+it)
#   = (1/2) log 2 + sum_{rho'} (beta'-1/2)/((beta'-1/2)^2 + (t-gamma')^2) + E(t),  E(t) -> 0,
# with the nontrivial zeros rho' of zeta' (both signs of gamma').  Zeros found by Newton from a grid of starts, deduplicated;
# count checked by the argument principle on [1/2+eps, 5] x [tlo, thi] (contour integral of zeta''/zeta').
# usage: python3 p4_zetaprime.py TZ tlo thi npts
import sys, time, math
import mpmath as mp
mp.mp.dps = 18
TZ = float(sys.argv[1]); tlo = float(sys.argv[2]); thi = float(sys.argv[3]); npts = int(sys.argv[4])
t0 = time.time()
d1 = lambda s: mp.zeta(s, derivative=1)
d2 = lambda s: mp.zeta(s, derivative=2)
roots = []
starts = [mp.mpc(sg, tt) for tt in [x*0.5 for x in range(8, int(2*TZ) + 1)] for sg in (0.75, 1.5)]
for s in starts:
    z = s; ok = False
    for it in range(30):
        f1 = d1(z); f2 = d2(z)
        if f2 == 0: break
        st = f1/f2
        if abs(st) > 0.6: st = 0.6*st/abs(st)
        z = z - st
        if z.real < 0.3 or z.real > 6 or z.imag < 1 or z.imag > TZ + 3: break
        if abs(st) < mp.mpf(10)**-14: ok = True; break
    if ok and abs(d1(z)) < 1e-10 and all(abs(z - r) > 1e-6 for r in roots):
        roots.append(z)
roots.sort(key=lambda r: float(r.imag))
print(f'zeta\' zeros found with 1 < gamma\' < {TZ}: {len(roots)}; min beta\' = {mp.nstr(min(r.real for r in roots), 6)}; first {mp.nstr(roots[0], 8)}; elapsed {time.time()-t0:.1f}', flush=True)
# argument-principle count on [0.5+1e-3, 5] x [a, b]
def count(a, b):
    pts = [mp.mpc(0.501, a), mp.mpc(5, a), mp.mpc(5, b), mp.mpc(0.501, b), mp.mpc(0.501, a)]
    tot = 0
    for p, q in zip(pts[:-1], pts[1:]):
        n = max(2, int(abs(q - p)/2))
        sub = [p + (q - p)*k/n for k in range(n + 1)]
        tot += mp.quad(lambda s: d2(s)/d1(s), sub)
    return tot/(2j*mp.pi)
for (a, b) in [(1.5, TZ)]:
    c = count(a, b); nf = sum(1 for r in roots if a < r.imag < b)
    print(f'argument principle [{a},{b}]: {mp.nstr(c.real, 6)} (imag {mp.nstr(c.imag, 2)}) vs found {nf}; elapsed {time.time()-t0:.1f}', flush=True)
H = 0.5*math.log(2)
res = []
for k in range(npts):
    t = tlo + (thi - tlo)*k/(npts - 1)
    s = mp.mpc(0.5, t)
    php = -mp.siegeltheta(t, derivative=1) - mp.re(d2(s)/d1(s))
    P = sum((r.real - 0.5)/((r.real - 0.5)**2 + (t - r.imag)**2) + (r.real - 0.5)/((r.real - 0.5)**2 + (t + r.imag)**2) for r in roots)
    res.append((t, float(php), float(P), float(php - H - P)))
for t, a, b, e in res:
    print(f't={t:8.3f} phi\'={a:9.5f} poisson={b:9.5f} E=phi\'-log2/2-poisson={e:+.5f}', flush=True)
print('elapsed', round(time.time() - t0, 1))
