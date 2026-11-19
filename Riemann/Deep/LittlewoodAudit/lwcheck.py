# Independent re-implementation (written from PROOF.md, not from cert.py) of Prop. 10 (F), Prop. 11 (Fup) and Lemma T,
# in mpmath interval arithmetic at 256 bits, with its own float optimiser and ADAPTIVE cells (independent covering).
# Run: uv run --no-project --with mpmath python3 lwcheck.py [B|78] > log
import sys, math, time
from mpmath import iv
iv.prec = 256
MODE = sys.argv[1] if len(sys.argv) > 1 else 'B'
def I(n, d=1): return iv.mpf(n)/iv.mpf(d)
if MODE == '78': Bi, Bf = I(7, 8), 0.875
else: Bi, Bf = I(87338196661320752, 10**17), 0.87338196661320752
bi, bf = 1 - Bi, 1 - Bf
G, PI, LN2, HALF = iv.euler, iv.pi, iv.log(2), I(1, 2)
C0 = iv.log(PI**2/6) - G
EUP = I(-133258227573, 10**11)            # >= E = -1.332582275733221 (Dusart D:695)
Gf = 0.5772156649015329; C0f = math.log(math.pi**2/6) - Gf; EUPf = -1.33258227573

def taubar(l):
    return iv.exp(-l/2) + ((l/LN2)*(iv.exp(l/3) + 3)/2 + 8)*iv.exp(-l)
def taubarf(l):
    return math.exp(-l/2) + (l/math.log(2))*(math.exp(-2*l/3) + 3*math.exp(-l))/2 + 8*math.exp(-l)

def kap(L, h, L1, h1, B=None, b=None):
    B = Bi if B is None else B; b = bi if b is None else b
    g0 = lambda be: (iv.exp(-(1+h)*L*(1-be))/(1+h) + iv.exp(-L*(1-be)))/(h*L*L)
    g1 = lambda be: (iv.exp(-(1+h1)*L1*(1-be)) + iv.exp(-L1*(1-be)))/(h1*L1)
    return g0(B)*B/b + g0(HALF), g1(B)*B/b + g1(HALF)
def kapf(L, h, L1, h1):
    g0 = lambda be: (math.exp(-(1+h)*L*(1-be))/(1+h) + math.exp(-L*(1-be)))/(h*L*L)
    g1 = lambda be: (math.exp(-(1+h1)*L1*(1-be)) + math.exp(-L1*(1-be)))/(h1*L1)
    return g0(Bf)*Bf/bf + g0(0.5), g1(Bf)*Bf/bf + g1(0.5)

def F(ell, P):          # Prop. 10: log|L(1,chi)| >= F; P exact rationals (L, h, L1, h1)
    L, h, L1, h1 = P
    k0, k1 = kap(L, h, L1, h1)
    if not (k1.b < 1): return None
    kp = k0/(1 - k1); x = iv.exp(L); x1 = iv.exp(L1)
    Eh = ((1+h)*iv.log(1+h) - h)/h
    eD = 1/(10*L**2) + I(4, 15)/L**3
    T0 = iv.exp(-L)/(L*(1 - iv.exp(-2*L))); T1 = iv.exp(-L1)/(1 - iv.exp(-2*L1))
    W = (L1*(1 + h1/2) + EUP + I(1, 5)/L1 + I(1, 5)/L1**2 + L1*taubar(L1) + (1+h1)*L1/(2*(x1 - 1))
         + (ell - iv.log(PI))/2 - G/2 + T1)
    return C0 - iv.log(L) - Eh - (eD + 1/(x - 1) + taubar(L)/2 + T0) - kp*W
def Ff(ell, P):
    L, h, L1, h1 = P
    if L < 10 or L1 < 8 or not (0 < h <= 2) or not (0 < h1 <= 2): return -1e9
    k0, k1 = kapf(L, h, L1, h1)
    if k1 >= 1: return -1e9
    kp = k0/(1 - k1); Eh = ((1+h)*math.log1p(h) - h)/h
    eD = 1/(10*L*L) + 4/(15*L**3); ex, ex1 = math.exp(-L), math.exp(-L1)
    W = (L1*(1 + h1/2) + EUPf + 0.2/L1 + 0.2/L1**2 + L1*taubarf(L1) + (1+h1)*L1*ex1/(2*(1 - ex1))
         + (ell - math.log(math.pi))/2 - Gf/2 + ex1/(1 - ex1*ex1))
    return C0f - math.log(L) - Eh - (eD + ex/(1 - ex) + taubarf(L)/2 + ex/L) - kp*W
def Fup(ell, P):        # Prop. 11: log|L(1,chi)| <= Fup
    L, h, L1, h1 = P
    k0, k1 = kap(L, h, L1, h1)
    if not (k1.b < 1): return None
    kp = k0/(1 - k1); x = iv.exp(L)
    Eh = ((1+h)*iv.log(1+h) - h)/h
    eD = 1/(10*L**2) + I(4, 15)/L**3
    T1 = iv.exp(-L1)/(1 - iv.exp(-2*L1))
    Wup = L1*(1 + h1/2) - G + I(1, 5)/L1 + I(1, 5)/L1**2 + (ell - iv.log(PI))/2 - G/2 + T1
    return G + iv.log(L) + Eh + eD + 1/(2*(x - 1)) + kp*Wup
def Fupf(ell, P):       # returns -Fup (to maximise)
    L, h, L1, h1 = P
    if L < 10 or L1 < 8 or not (0 < h <= 2) or not (0 < h1 <= 2): return -1e9
    k0, k1 = kapf(L, h, L1, h1)
    if k1 >= 1: return -1e9
    kp = k0/(1 - k1); Eh = ((1+h)*math.log1p(h) - h)/h; ex1 = math.exp(-L1)
    Wup = L1*(1 + h1/2) - Gf + 0.2/L1 + 0.2/L1**2 + (ell - math.log(math.pi))/2 - Gf/2 + ex1/(1 - ex1*ex1)
    return -(Gf + math.log(L) + Eh + 1/(10*L*L) + 4/(15*L**3) + math.exp(-L)/2 + kp*Wup)

def optimise(ell, start=None, obj=Ff):
    if start is None:
        best = None
        for L in [10*1.25**i for i in range(30)]:
            for h in (0.01, 0.03, 0.1, 0.3, 0.6, 1.0, 1.5, 2.0):
                for r in (0.2, 0.35, 0.5, 0.7, 1.0):
                    for h1 in (0.03, 0.1, 0.3, 1.0, 2.0):
                        P = (L, h, max(8.0, r*L), h1); v = obj(ell, P)
                        if best is None or v > best[0]: best = (v, P)
        P = list(best[1])
    else: P = list(start)
    st = [P[0]*0.05, 0.05, P[2]*0.05, 0.05]; cur = obj(ell, P)
    while max(st[1], st[3]) > 2e-5 or max(st[0], st[2]) > 2e-4:
        moved = False
        for i in range(4):
            for s in (1, -1):
                Q = list(P); Q[i] += s*st[i]
                if i in (1, 3): Q[i] = min(2.0, Q[i])
                v = obj(ell, Q)
                if v > cur: cur, P, moved = v, Q, True
        if not moved: st = [t/2 for t in st]
    r = (max(10000, round(P[0]*1000)), min(2000, max(1, round(P[1]*1000))), max(8000, round(P[2]*1000)), min(2000, max(1, round(P[3]*1000))))
    return r, tuple(t/1000 for t in r)
def Pq(r): return tuple(I(t, 1000) for t in r)
def flo(x): return float(x.a)

def cover(c, lo_num, hi_num, den, upper=False, wmax=None):
    """Adaptive independent covering of [lo_num/den, hi_num/den] in Lambda = loglog q.
    lower: F(e^{Lam_j}) > log(c/Lam_i); upper: Fup(e^{Lam_j}) < log(c Lam_i). Cell widths halve on failure (min 1/den)."""
    i = lo_num; w = 1; Pst = None; n = 0; worst = (9, None); wmax = wmax or 64
    cI = I(*c) if isinstance(c, tuple) else c
    while i < hi_num:
        w = min(w, hi_num - i)
        j = i + w; Li, Lj = I(i, den), I(j, den); ell = iv.exp(Lj)
        r, Pst = optimise(float(ell.a), Pst, Fupf if upper else Ff)
        if upper:
            v = Fup(ell, Pq(r)); tg = iv.log(cI*Li); ok = v is not None and v.b < tg.a; m = flo(tg - v) if ok else None
        else:
            v = F(ell, Pq(r)); tg = iv.log(cI/Li); ok = v is not None and v.a > tg.b; m = flo(v - tg) if ok else None
        if ok:
            n += 1
            if m < worst[0]: worst = (m, (i/den, j/den) + Pst)
            i = j; w = min(2*w, wmax)
        else:
            if w == 1: return False, (i/den, j/den, Pst), n
            w //= 2
    return True, worst, n


# ------------------------------------------------------------------ main
import random
random.seed(20261119)
t0 = time.time()
def say(*a): print(*a, flush=True)
say('MODE', MODE, ' limit e^C0 b =', flo(iv.exp(C0)*bi), ' e^gamma/b =', flo(iv.exp(G)/bi))
if MODE == '78': CA, CB, LB, CB2, CUP = (105, 1000), (11, 100), 3000, (114, 1000), (15, 1)
else:            CA, CB, LB, CB2, CUP = (11, 100), (115, 1000), 3130, (116, 1000), (15, 1)
# (0) Lemma T(ii) brackets at L = 8, h = 2
L8 = I(8); T18 = iv.exp(-L8)/(1 - iv.exp(-2*L8))
om = EUP + I(1, 5)/L8 + I(1, 5)/L8**2 + L8*taubar(L8) + 3*L8/(2*(iv.exp(L8) - 1)) - iv.log(PI)/2 - G/2 + T18
omu = -G + I(1, 5)/L8 + I(1, 5)/L8**2 - iv.log(PI)/2 - G/2 + T18
say('(0) omega(8)|h=2 <=', float(om.b), om.b < 0, ' upper bracket <=', float(omu.b), omu.b < 0)
# (1) tightest cells and small-q value at the cert's own parameters (B mode only)
if MODE == 'B':
    for lab, lo_, hi_, P, c, up in (('A', 3000, 3005, (12863, 1016, 8000, 1834), CA, False), ('B', 3130, 3135, (13133, 1003, 8000, 1879), CB, False),
                                    ('B2', 99990, 100000, (760908, 25, 274644, 88), CB2, False), ('G', 3000, 3005, (12764, 1038, 8000, 1845), CUP, True)):
        ell = iv.exp(I(hi_, 1000)); Li = I(lo_, 1000)
        if up: v = Fup(ell, Pq(P)); m = iv.log(I(*c)*Li) - v
        else:  v = F(ell, Pq(P)); m = v - iv.log(I(*c)/Li)
        say(f'(1) tightest {lab} [{lo_/1000},{hi_/1000}] P={P}: margin in [{float(m.a):.15f}, {float(m.b):.15f}]')
    v = F(iv.exp(I(3)), Pq((12854, 1016, 8000, 1832))); say(f'(1) small q at cert P: F(e^3) in [{float(v.a):.9f},{float(v.b):.9f}]  log(0.11/3) = {float(iv.log(I(11,300)).a):.9f}')
# (2) small q with own parameters
ell3 = iv.exp(I(3)); r, P = optimise(float(ell3.a)); v = F(ell3, Pq(r)); tg = iv.log(I(*CA)/3)
say(f'(2) small q: own P={P}: F(e^3) >= {flo(v):.6f} > log(cA/3) = {float(tg.b):.6f}: {v.a > tg.b}; |L(1)| >= {math.exp(flo(v)):.5f}')
# (3) adaptive independent coverings
for lab, c, lo_, hi_, up in (('A', CA, 3000, 40000, False), ('B', CB, LB, 40000, False), ('B2', CB2, 4000, 100000, False), ('G', CUP, 3000, 40000, True)):
    ok, worst, n = cover(c, lo_, hi_, 1000, upper=up, wmax=256)
    say(f'(3) cover {lab} c={c} on [{lo_/1000}, {hi_/1000}]: {"PASS" if ok else "FAIL"}, {n} cells, worst {worst}  t={time.time()-t0:.0f}s')
# (4) tails, Lemma T and T'
def tailT(LamN, h):
    assert (LamN/bi).a >= 10
    L = LamN/bi; ell = iv.exp(LamN); x = iv.exp(L)
    k0, k1 = kap(L, h, L, h); kp = k0/(1 - k1)
    Eh = ((1+h)*iv.log(1+h) - h)/h; eD = 1/(10*L**2) + I(4, 15)/L**3; T0 = iv.exp(-L)/(L*(1 - iv.exp(-2*L)))
    Err = eD + 1/(x - 1) + taubar(L)/2 + T0 + kp*(L*(1 + h/2) + ell/2)
    Eru = eD + 1/(2*(x - 1)) + kp*(L*(1 + h/2) + ell/2)
    return k1, iv.exp(C0 - Eh - Err)*bi, iv.exp(G)/bi*iv.exp(Eh + Eru)
for LN, h in ((I(40), I(15, 1000)), (I(100), I(6, 1000)), (I(1000), I(1, 1000)), (I(10000), I(1, 10000))):
    k1, c, cu = tailT(LN, h)
    say(f'(4) tail Lam>={float(LN.a):g} h={float(h.a):g}: kappa1 <= {float(k1.b):.3e} <1:{k1.b < 1}  c >= {float(c.a):.6f}  Cup <= {float(cu.b):.4f}')
# (5) class-number thresholds, independent cells in ell = log d from just below log(d0)
def cls(H, d0, cA):
    den = 10**4; lo_ = math.floor(math.log(d0)*den) - 1
    e = lo_
    while not (math.log(e/den) >= 3 and (iv.exp(I(e, den)/2)*I(*cA)/(PI*iv.log(I(e, den)))).a > H): e += 10
    i, w, Pst, n = lo_, 16, None, 0
    while i < e:
        w = min(w, e - i); j = i + w
        r, Pst = optimise(float(I(j, den).a), Pst)
        v = F(I(j, den), Pq(r)); ok = v is not None and (iv.exp(I(i, den)/2)*iv.exp(v)/PI).a > H
        if ok: i, n, w = j, n + 1, min(2*w, 512)
        elif w == 1: return f'FAIL at ell={i/den}'
        else: w //= 2
    return f'PASS from ell={lo_/den} <= log({d0:g}) = {math.log(d0):.5f} to {e/den} ({n} cells), claim (A) beyond'
for H, d0 in ((100, 6.973e7), (1000, 7.292e9), (2000, 2.957e10), (10**4, 7.627e11), (10**5, 7.976e13)) if MODE == 'B' else ((2000, 3.047e10), (10**5, 8.137e13)):
    say(f'(5) h(-d) > {H} for d >= {d0:g}: {cls(H, d0, CA)}')
# (6) 20 random cells of the cert's own partitions (0.005 cells), own parameters
bad = 0
for lab, c, lo_, hi_, st in (('A', CA, 3000, 40000, 5), ('B', CB, LB, 40000, 5), ('B2', CB2, 4000, 100000, 10)):
    for _ in range(20):
        i = random.randrange(lo_, hi_, st); j = i + st
        r, P = optimise(float(iv.exp(I(j, 1000)).a)); v = F(iv.exp(I(j, 1000)), Pq(r))
        if not (v is not None and v.a > iv.log(I(*c)/I(i, 1000)).b): bad += 1; say('  random cell FAIL', lab, i, j, P)
say(f'(6) random cert cells: {60 - bad}/60 pass  t={time.time()-t0:.0f}s')

