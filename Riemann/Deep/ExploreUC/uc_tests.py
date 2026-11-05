# uc_tests.py - explore2 tests T1-T3 from saved Cauchy jets of H_0 near Lehmer's pair (z = 2*gamma ~ 14010.163)
import mpmath as mp, time, pickle
mp.mp.dps = 80
T0 = time.time()
I = mp.mpc(0, 1)
d = pickle.load(open('jets.pkl', 'rb'))
X0 = mp.mpf(d['X0']); S = mp.mpf(d['S']); raw = [mp.mpc(v) for v in d['raw']]
Z0 = [mp.mpf(z) for z in d['Z0']]
RHO = mp.mpf(7); M = 384; NJ = 150; NREL = 120; NMAX = 70; KMAX = 25
zs = [X0 + RHO*mp.expj(2*mp.pi*m/M) for m in range(M)]
vals = [v/S for v in raw]
W = [mp.expj(-2*mp.pi*m/M) for m in range(M)]
def log(*a): print(*a, flush=True)
def jets(vs):
    c = []
    for j in range(NJ+1):
        acc = mp.mpc(0)
        for m in range(M): acc += vs[m]*W[(j*m) % M]
        c.append(mp.re(acc/M/RHO**j))
    return c
A = [[mp.factorial(n+2*k)/(mp.factorial(n)*mp.factorial(k)) for k in range(KMAX+1)] for n in range(NMAX+1)]
def heat(c, t):
    out = []
    for n in range(NMAX+1):
        acc = mp.mpf(0); tk = mp.mpf(1)
        for k in range(KMAX+1):
            if n+2*k > NREL: break
            acc += A[n][k]*tk*c[n+2*k]; tk *= -t
        out.append(acc)
    return out
def pv(c, h):
    acc = 0
    for x in reversed(c): acc = acc*h + x
    return acc
def der(c): return [(n+1)*c[n+1] for n in range(len(c)-1)]
BIN = [[mp.binomial(m, n) for n in range(NMAX+1)] for m in range(NMAX+1)]
def shift(c, h):
    dd = len(c); q = []
    for n in range(dd):
        acc = mp.mpf(0); hp = mp.mpf(1)
        for m in range(n, dd):
            acc += c[m]*BIN[m][n]*hp; hp *= h
        q.append(acc)
    return q
def mul(p, q):
    r = [mp.mpf(0)]*(len(p)+len(q)-1)
    for i, a in enumerate(p):
        for j, b in enumerate(q): r[i+j] += a*b
    return r
DF = [mp.mpf(1)]
for j in range(1, 200): DF.append(DF[-1]*(2*j-1))   # (2j-1)!!
def gmean(p, v):   # E p(Y), Y ~ N(0, v)
    return sum(p[2*j]*v**j*DF[j] for j in range(len(p)//2 + (len(p) % 2)) if 2*j < len(p))
cH = [mp.mpf(x) for x in d['cH']]
def plant(a, b, Y):
    w = X0 + (a+b)/2 + I*Y
    vs = [vals[m]*(zs[m]-w)*(zs[m]-mp.conj(w))/((zs[m]-X0-a)*(zs[m]-X0-b)) for m in range(M)]
    return jets(vs)
def dz(c, h0, t0):   # double zero of P_t near (h0, t0)
    f = lambda h, t: (pv(heat(c, t), h), pv(der(heat(c, t)), h))
    h, t = mp.findroot(f, (mp.mpf(h0), mp.mpf(t0)))
    return h, t
def freq(c, hb, tb, s):
    q = shift(heat(c, tb + s), hb)
    q1 = der(q)
    v = 2*s
    return s*gmean(mul(q1, q1), v)/gmean(mul(q, q), v)
# ---------- models
a_L, b_L = Z0[4], Z0[5]          # Lehmer pair
a_N, b_N = Z0[2], Z0[3]          # normal-gap pair (gap ~2.70)
log("Lehmer pair h = %s, %s gap %s; NP pair h = %s, %s gap %s" % (mp.nstr(a_L,8), mp.nstr(b_L,8), mp.nstr(b_L-a_L,6), mp.nstr(a_N,6), mp.nstr(b_N,6), mp.nstr(b_N-a_N,4)))
hc, tc = dz(cH, (a_L+b_L)/2, -(b_L-a_L)**2/8)
log("Lehmer collision (forward heat): h_c = %s, t_c = %s  (isolated-pair guess -g^2/8 = %s)" % (mp.nstr(hc,8), mp.nstr(tc,8), mp.nstr(-(b_L-a_L)**2/8,6)))
models = {}
for (nm, a, b, Y) in [("LP0.05", a_L, b_L, mp.mpf('0.05')), ("LP0.15", a_L, b_L, mp.mpf('0.15')), ("LP0.40", a_L, b_L, mp.mpf('0.40')),
                      ("NP0.30", a_N, b_N, mp.mpf('0.30')), ("NP0.75", a_N, b_N, mp.mpf('0.75'))]:
    c = plant(a, b, Y)
    m = (a+b)/2
    try:
        hs, ts = dz(c, m, Y**2/2*mp.mpf('0.6'))
    except Exception as e:
        hs, ts = None, None
    models[nm] = (c, a, b, Y, hs, ts)
    log("plant %s: off-line pair %s +- %si lands (creation) at h* = %s, t* = %s; de Bruijn bound Y^2/2 = %s" % (nm, mp.nstr(m,6), mp.nstr(Y,3), mp.nstr(hs,8) if hs is not None else '-', mp.nstr(ts,8) if ts is not None else '-', mp.nstr(Y**2/2,5)))
log("t=%.1fs" % (time.time()-T0))
# ---------- T1: doubling index and Jensen counts around base h = 1.7 (no zero nearby)
hb = mp.mpf('1.7')
def Mr(c, r, n=720):
    return max(abs(pv(c, hb + r*mp.expj(2*mp.pi*k/n))) for k in range(n))
def count(zlist, r): return sum(1 for z in zlist if abs(z - hb) < r)
c0 = heat(cH, 0)
log("T1 doubling index DI(r)=ln(M(2r)/M(r)), Jensen bound J(r)=ln(M(2r)/|P(hb)|)/ln 2, at hb=1.7")
for nm, c in [("H", c0), ("LP0.15", heat(models["LP0.15"][0], 0)), ("NP0.75", heat(models["NP0.75"][0], 0))]:
    if nm == "H": zl = list(Z0)
    else:
        _, a, b, Y, _, _ = models[nm]; mm = (a+b)/2
        zl = [z for z in Z0 if abs(z-a) > 1e-6 and abs(z-b) > 1e-6] + [mm + I*Y, mm - I*Y]
    row = []
    for r in [mp.mpf('0.5'), mp.mpf(1), mp.mpf(2), mp.mpf(3)]:
        m1 = Mr(c, r); m2 = Mr(c, 2*r)
        row.append("r=%s DI=%s J=%s n=%d nreal=%d" % (mp.nstr(r,2), mp.nstr(mp.log(m2/m1),5), mp.nstr(mp.log(m2/abs(pv(c, hb)))/mp.log(2),4), count(zl, r), sum(1 for z in zl if abs(z-hb) < r and abs(mp.im(z)) < 1e-9)))
    log("  %-7s " % nm + " | ".join(row))
log("t=%.1fs" % (time.time()-T0))
# ---------- T2: Poon frequency N(s) = s int u_x^2 G_s / int u^2 G_s, base (hb, tb), data at time tb + s
sgrid = [mp.mpf(x) for x in ['1e-6', '1e-5', '1e-4', '1e-3', '0.01', '0.03', '0.06', '0.1', '0.1164', '0.14', '0.2']]
bases = [("H regular zero h=3.316 t=0", cH, Z0[6], mp.mpf(0)),
         ("H regular zero h=-5.457 t=0", cH, Z0[1], mp.mpf(0)),
         ("H Lehmer collision", cH, hc, tc)]
for nm in ["LP0.15", "NP0.75"]:
    c, a, b, Y, hs, ts = models[nm]
    if hs is None: continue
    bases.append(("F_%s creation" % nm, c, hs, ts))
    ct = heat(c, ts)
    # regular real zero of F at time t* near h = 3.316
    zr = mp.findroot(lambda x: pv(ct, x), Z0[6])
    bases.append(("F_%s regular zero at t*" % nm, c, zr, ts))
log("T2 Poon frequency; s grid: " + " ".join(mp.nstr(s, 3) for s in sgrid))
for (nm, c, h_, t_) in bases:
    Ns = [freq(c, h_, t_, s) for s in sgrid]
    mono = all(Ns[i] <= Ns[i+1] + mp.mpf('1e-30') for i in range(len(Ns)-1))
    log("  %-34s " % nm + " ".join(mp.nstr(x, 4) for x in Ns) + ("  monotone" if mono else "  NOT MONOTONE"))
k = mp.log(X0/(4*mp.pi))/4
log("cos model at this height: k = %s; N_model(s) = s k^2 coth(2 k^2 s): " % mp.nstr(k,5) + " ".join(mp.nstr(s*k*k*mp.coth(2*k*k*s), 4) for s in sgrid))
log("t=%.1fs" % (time.time()-T0))
# ---------- T3: windowed Laguerre W(v) = E L1(m + Y), Y ~ N(0,v), L1 = P'^2 - P P''
vgrid = [mp.mpf(10)**(mp.mpf(e)/4) for e in range(-20, 5)]
def wlag(c, m, v):
    q = shift(c, m); q1 = der(q); q2 = der(q1)
    a1 = mul(q1, q1); a2 = mul(q, q2)
    n = max(len(a1), len(a2)); a1 += [0]*(n-len(a1)); a2 += [0]*(n-len(a2))
    L = [x - y for x, y in zip(a1, a2)]
    return gmean(L, v)
log("T3 windowed Laguerre at the pair centre: pointwise L1(m) and smallest window sd where W >= 0 for all larger windows")
for nm, c, m in [("H (Lehmer mid)", cH, (a_L+b_L)/2), ("H (NP mid)", cH, (a_N+b_N)/2)] + [("F_"+k_, models[k_][0], (models[k_][1]+models[k_][2])/2) for k_ in models]:
    c0_ = heat(c, 0)
    Ws = [wlag(c0_, m, v) for v in vgrid]
    neg = [i for i, w in enumerate(Ws) if w < 0]
    sdstar = mp.sqrt(vgrid[max(neg)+1]) if neg and max(neg)+1 < len(vgrid) else (mp.mpf(0) if not neg else None)
    if neg and max(neg)+1 < len(vgrid):
        lo, hi = mp.log(vgrid[max(neg)]), mp.log(vgrid[max(neg)+1])   # bisect the last sign change in log-variance
        for _ in range(40):
            mid = (lo+hi)/2
            if wlag(c0_, m, mp.exp(mid)) < 0: lo = mid
            else: hi = mid
        sdx = mp.sqrt(mp.exp(hi))
    else:
        sdx = None
    Yp = models[nm[2:]][3] if nm.startswith("F_") else None
    log("  %-16s L1(m)=%s  last sign change of W at window sd = %s%s ; W(sd=0.1,0.3,1)=%s" % (nm, mp.nstr(wlag(c0_, m, mp.mpf('1e-30')), 4), mp.nstr(sdx, 4) if sdx is not None else 'none (W>0)', (" = %s Y" % mp.nstr(sdx/Yp, 3)) if (sdx is not None and Yp) else "", " ".join(mp.nstr(wlag(c0_, m, v), 3) for v in [mp.mpf('0.01'), mp.mpf('0.09'), mp.mpf(1)])))
log("done t=%.1fs" % (time.time()-T0))
