# gl2audit spot check: for sample certificate cells, sigma, mu from the (untrusted) proposer ceil.build, Lemma D value D(R,lam)
# by OWN mpmath Gauss quadrature, own M(R,lam) (myceil, Chebyshev Ritz K=36): must have target <= D <= M at every sample point.
import sys, os, json, random
ARGS = sys.argv[:]; cert = json.load(open(ARGS[1])); xn, xd = cert['xi'].split('/')
sys.argv = ['m', xn, xd, '36']; import myceil as MY
sys.argv = ['c', xn, xd, cert['kappa_up']]; sys.path.insert(0, os.environ['GL2CEIL']); import ceil as PROP
mp = MY.mp; XI = MY.XI; KUP = mp.mpf(cert['kappa_up']); KP = cert['KP']; HL = cert['HL']
gx, gw = MY.gl(24)
def peval(c, x):
    v = mp.mpf(0)
    for cj in reversed(c): v = v*x + cj
    return v
def pder(c): return [j*c[j] for j in range(1, len(c))] or [0.0]
def Dval(sc, mc, R, lam):
    s = []
    for k in range(KP):
        c = [mp.mpf(float(t)) for t in sc[k]]
        if k: c[0] = peval(s[-1], mp.mpf(1))
        s.append(c)
    m = [[mp.mpf(float(t)) for t in mc[j]] for j in range(KP//2)]
    K2 = S1 = M = K1 = mp.mpf(0)
    for k in range(KP):
        ds = pder(s[k])
        for x, w in zip(gx, gw):
            y = (k + x)/KP; wt = w/KP
            sv = peval(s[k], x); dsv = peval(ds, x)*KP
            mv = peval(m[k-KP//2], x) if k >= KP//2 else peval(m[KP//2-1-k], 1-x)
            K2 += wt*sv*sv; M += wt*mv; K1 += wt*(dsv + mv*mp.exp(-R*y))**2
    s0 = peval(s[0], mp.mpf(0)); s1 = peval(s[-1], mp.mpf(1))
    C, B = MY.CB(lam)
    return mp.mpf(1)/2 + M - 2*s0 - (2*s1 + M*mp.exp(-R))**2/2 - (XI/C)*K1 - K2/(XI*B)
cells = [c for c in cert['cells'] if c[4] == 'dual']
random.seed(int(ARGS[2]) if len(ARGS) > 2 else 7)
pick = sorted(cells, key=lambda c: abs(float(eval(c[0])) - 9.8618) + abs(c[2]*HL - 0.7656))[:1] + random.sample(cells, int(ARGS[3]) if len(ARGS) > 3 else 4)
worst = None
for Rcs, ds, ia, ib, how in pick:
    Rc = mp.mpf(eval(Rcs)); d = mp.mpf(eval(ds)); la, lb = mp.mpf(ia*HL), mp.mpf(ib*HL)
    sc, mc, beta, J = PROP.build(float(Rc), (ia+ib)/2*HL)
    for R, lam in [(Rc-d, la), (Rc+d, lb), (Rc, (la+lb)/2), (Rc-d, lb), (Rc+d, la)]:
        D = Dval(sc, mc, R, lam); Mm = MY.Mlam(R, lam); t = mp.exp((1-KUP)*R)
        ok = t <= D <= Mm*(1+mp.mpf(10)**-12)
        r1, r2 = (D-t)/t, (Mm-D)/Mm
        worst = r1 if worst is None or r1 < worst else worst
        print(f"R={mp.nstr(R,10)} lam={mp.nstr(lam,8)}: (D-t)/t={mp.nstr(r1,3)} (M-D)/M={mp.nstr(r2,3)} {'ok' if ok else 'VIOLATION'}", flush=True)
print('min (D-t)/t over samples', mp.nstr(worst, 4))
