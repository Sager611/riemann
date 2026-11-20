# cz_arb.py -- own Arb certificate (stc2) of C_zeta >= sum_p M_p (Casana-Giner, Zenodo 23866074, (4.7), Lemmas 4.2-4.4,
# Prop. 4.5), all six kernels in one pass over the primes (segmented sieve, streamed; no prime list held).
# Usage: cz_arb.py N0 [names]   (the paper: N0 = 1e8, P0 = 1000)
# (1) p <= P0: M_p <= max_x phi_p(x) + sum_{p^k > N0} H_{p,k} (Lemma 4.2(a)), phi_p = sum_{p^k<=N0} lam_k (a_k cos kx + b_k sin kx);
#     max phi_p by branch and bound on arcs [c-r,c+r]: phi(c) + |phi'(c)| r + Lam2 r^2/2, Lam2 = sum k^2 H_{p,k} >= |phi''|;
#     also the Euclidean bound sum_k H_{p,k}; the min of the two is used (as (4.11)).
# (2) P0 < p <= N0: sum_{p^k <= N0} H_{p,k}, H_{p,k} = log p p^(-k s1) h(k log p), one Arb ball per prime power.
# (3) n > N0: Lemma 4.3 with G = g, g(x) = h(log x) x^(-s1); g is decreasing on [N0, oo): nu(L) = s1 (A^2+B^2) - (A A' + B B') > 0
#     (A = pi a, B = pi b) by interval bisection on [log N0, LM].
#     int_{N0}^oo G dx <= int_{L0}^{LM} h e^{(1-s1)L} dL (acb.integral, rigorous) + int_{LM}^oo gt (closed form)
#     + gt(LM) e^LM, with G = max(g, gt(LM)) on [N0, e^LM] and G = gt beyond (non-increasing, >= g); LM <= 400, nu > 0 on [L0, LM].
import json, sys, time, math
from fractions import Fraction as Fr
from math import factorial
import numpy as np
from flint import arb, acb, arb_poly, ctx
ctx.prec = 64
KJ = json.load(open('Riemann/Deep/STC2/kernels.json'))
N0 = int(float(sys.argv[1])); names = sys.argv[2].split(',') if len(sys.argv) > 2 else list(KJ)
P0 = 1000; PI = arb.pi(); t0 = time.time()
T2 = {'L5': 0.70360002, 'L20': 1.46807534, 'L60': 2.10426949, 'L150': 2.25578346, 'C400': 4.46374163, 'D6000': 4.58951035}
KS = []
for nm in names:
    k = KJ[nm]; dq = Fr(k['d']).limit_denominator(10**6); d = arb(dq.numerator)/dq.denominator
    B = [(X, int(j)) for X, j in k['B']]; C = [arb(c) for c in k['c']]; DEG = max(j for _, j in B) - 1
    pa = [arb(0)]*(DEG+1); pb = [arb(0)]*(DEG+1)
    for (X, j), c in zip(B, C):
        if X == 'R': pa[j-1] -= c/factorial(j-1)
        else: pb[j-1] -= c/factorial(j-1)
    KS.append(dict(name=nm, dq=dq, d=d, s1=d + arb(1)/2, pa=pa, pb=pb, Pa=arb_poly(pa), Pb=arb_poly(pb),
                   Pa1=arb_poly(pa).derivative(), Pb1=arb_poly(pb).derivative(), DEG=DEG, small=arb(0), small_e=arb(0),
                   small_first=arb(0), med=arb(0), gridf=0.0))
DS = sorted(set(K['dq'] for K in KS))
def hpi(K, L, ed):   # pi*h(L) given ed = e^{-dL}
    A = K['Pa'](L); Bv = (ed - 2)/L + K['Pb'](L)
    return (A*A + Bv*Bv).sqrt()
# ---- (1) small primes
def smallprimes(n):
    s = np.ones(n+1, bool); s[:2] = False
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0].tolist()
psi = arb(0); nprimes = 0
TWO_PI = 2*PI
def bb_max(al, be):
    # rigorous upper bound of max_x sum_k (al_k cos kx + be_k sin kx); al, be lists of (k, arb)
    Lam2 = sum((kk*kk*(a*a + b*b).sqrt() for kk, a, b in zip(range(1, 10**6), al, be)), arb(0)).upper()
    ks = list(range(1, len(al)+1))
    def ev(c):
        e = acb(0, c).exp(); z = acb(1); v = arb(0); dv = arb(0)
        for kk, a, b in zip(ks, al, be):
            z = z*e; v += a*z.real + b*z.imag; dv += kk*(b*z.real - a*z.imag)
        return v, dv
    n0 = 128; r0 = PI/n0; arcs = [(TWO_PI*(2*i+1)/(2*n0), r0, 0) for i in range(n0)]
    vals = []; best = None
    for (c, r, dep) in arcs:
        v, dv = ev(c); vals.append((c, r, dep, v, dv)); lb = v.lower(); best = lb if best is None or lb > best else best
    tol = arb(10)**-13; out = None; stack = vals
    while stack:
        c, r, dep, v, dv = stack.pop()
        up = (v + abs(dv)*r + Lam2*r*r/2).upper()
        if up <= best + tol: continue
        if dep > 45:
            out = up if out is None or up > out else out; continue
        for cc in (c - r/2, c + r/2):
            v2, dv2 = ev(cc); lb = v2.lower(); best = lb if lb > best else best; stack.append((cc, r/2, dep+1, v2, dv2))
    res = (best + tol).upper()
    return res if out is None or out < res else out, Lam2
for p in smallprimes(P0):
    lp = arb(p).log(); kmax = 0; n = p
    while n <= N0: kmax += 1; n *= p
    psi += lp*kmax; nprimes += 1
    for K in KS:
        al = []; be = []; Hs = arb(0); Hk1 = arb(0)
        for kk in range(1, kmax+1):
            L = kk*lp; ed = (-K['d']*L).exp(); w = lp*(-K['s1']*L).exp()
            A = K['Pa'](L)/PI; Bv = ((ed - 2)/L + K['Pb'](L))/PI
            al.append(w*A); be.append(w*Bv); Hk = w*(A*A + Bv*Bv).sqrt(); Hs += Hk; Hk1 += kk*Hk
        mx, Lam2 = bb_max(al, be)
        eu = Hs.upper(); K['small'] += (mx if mx < eu else eu); K['small_e'] += Hs
        K['small_first'] += Hk1   # sum k H (for the paper's first-order correction pi/M * sum k H)
        # float grid max at M = 4096 (estimate of the paper's (4.11) term, not used in the bound)
        xs = np.arange(4096)*2*math.pi/4096; f = np.zeros(4096)
        for kk, a, b in zip(range(1, kmax+1), al, be): f += float(a.mid())*np.cos(kk*xs) + float(b.mid())*np.sin(kk*xs)
        K['gridf'] += min(float(f.max()) + math.pi/4096*float(Hk1.mid()), float(Hs.mid()))
print('small primes done: %d primes, %.1f s' % (nprimes, time.time()-t0)); sys.stdout.flush()
# ---- (2) medium primes, segmented sieve
base = np.array(smallprimes(int(math.isqrt(N0)) + 1), dtype=np.int64)
SEG = 2_000_000; lo = P0 + 1
while lo <= N0:
    hi = min(lo + SEG, N0 + 1); s = np.ones(hi - lo, bool)
    for q in base.tolist():
        if q*q >= hi: break
        st = max(q*q, ((lo + q - 1)//q)*q); s[st - lo::q] = False
    for p in (np.nonzero(s)[0] + lo).tolist():
        nprimes += 1; L = arb(p).log(); sq = 1/arb(p).sqrt()
        kmax = 2 if p*p <= N0 else 1
        psi += L*kmax
        eds = {dq: (-(arb(dq.numerator)/dq.denominator)*L).exp() for dq in DS}
        for K in KS:
            ed = eds[K['dq']]; K['med'] += L*ed*sq*hpi(K, L, ed)
            if kmax == 2:
                L2 = 2*L; ed2 = ed*ed; K['med'] += L*ed2*sq*sq*hpi(K, L2, ed2)
    lo = hi
    if (lo // SEG) % 10 == 0: print('  sieved to %d, primes %d, %.0f s' % (lo, nprimes, time.time()-t0)); sys.stdout.flush()
for K in KS: K['med'] = K['med']/PI
print('medium primes done: %d primes total, psi(N0) = %s, %.1f s' % (nprimes, psi.str(15), time.time()-t0)); sys.stdout.flush()
# ---- (3) tail
L0 = arb(N0).log(); L1 = arb(800)
def inc(kk, X, al):   # int_X^oo L^k e^{-al L} dL
    return (-al*X).exp()*sum((arb(factorial(kk))/factorial(i)*X**i/al**(kk-i+1) for i in range(kk+1)), arb(0))
res = {}
for K in KS:
    d = K['d']; s1 = K['s1']; DEG = K['DEG']; pa, pb = K['pa'], K['pb']
    q = [abs(pa[i]) + abs(pb[i]) for i in range(DEG+1)]
    # nu > 0 (g decreasing) on [L0, LM] by bisection (natural interval extension); LM = last certified edge <= 400
    def nu(Lb):
        ed = (-d*Lb).exp(); A = K['Pa'](Lb); A1 = K['Pa1'](Lb); Bv = (ed - 2)/Lb + K['Pb'](Lb)
        B1 = -d*ed/Lb - (ed - 2)/(Lb*Lb) + K['Pb1'](Lb)
        return s1*(A*A + Bv*Bv) - (A*A1 + Bv*B1)
    def okint(a, b, dep=0):
        if nu(a.union(b)) > 0: return True, 1
        if dep >= 30: return False, 1
        m = (a + b)/2; o1, n1 = okint(a, m, dep+1)
        if not o1: return False, n1
        o2, n2 = okint(m, b, dep+1); return o2, n1 + n2
    edges = [arb(float(L0.lower()))]
    while float(edges[-1].mid()) < 400: edges.append(arb(min(float(edges[-1].mid())*1.01, 400.0)))
    nint = 0; bad = -1
    for i in range(len(edges)-1):
        ok, nn = okint(edges[i], edges[i+1]); nint += nn
        if not ok: bad = i
    LS = edges[bad+1]; LM = edges[-1]; LB = float(LS.mid())
    assert float(LS.mid()) < 60, 'monotonicity of g not certified on [%s, 400]' % LS
    # G = max(g, gt(LM)) on [e^LS, e^LM] (g decreasing there), G = gt beyond; on [N0, e^LS] (if LS > L0) the step majorant
    # G_i = max(upper bounds of g on the steps j >= i, g(e^LS)) (the paper's construction), steps of relative width 1e-5
    assert LM > DEG/s1
    gLS = hpi(K, LS, (-d*LS).exp())/PI*(-s1*LS).exp(); Istep = arb(0); GN0 = None
    if bad >= 0:
        st = [float(L0.lower())]
        while st[-1] < float(LS.mid()): st.append(min(st[-1]*(1 + 1e-5), float(LS.mid())))
        st = [arb(x) for x in st[:-1]] + [LS]; run = gLS.upper(); gs = []
        for i in range(len(st)-2, -1, -1):
            Lb = st[i].union(st[i+1]); gu = (hpi(K, Lb, (-d*Lb).exp())/PI*(-s1*st[i]).exp()).upper()
            run = run if run > gu else gu; gs.append((st[i], st[i+1], run))
        Istep = sum((G*(b_.exp() - a_.exp()) for a_, b_, G in gs), arb(0)); GN0 = gs[-1][2]
    def fz(z, analytic):
        A = sum((acb(pa[i])*z**i for i in range(len(pa))), acb(0))
        Bv = (((-acb(d)*z).exp() - 2)/z + sum((acb(pb[i])*z**i for i in range(len(pb))), acb(0)))
        return (A*A + Bv*Bv).sqrt(analytic=analytic)*((1 - acb(s1))*z).exp()/PI
    Imain = acb.integral(fz, acb(LS if bad >= 0 else L0), acb(LM))
    assert Imain.is_finite()
    al = s1 - 1
    Ifar = (sum((q[i]*inc(i, LM, al) for i in range(DEG+1)), arb(0)) + 3*inc(0, LM, al))/PI
    Ifar += (sum((q[i]*LM**i for i in range(DEG+1)), arb(0)) + 3)/PI*((1 - s1)*LM).exp()   # int_{N0}^{e^LM} gt(LM) dx <= gt(LM) e^LM
    gN0 = hpi(K, L0, (-d*L0).exp())/PI*(-s1*L0).exp() if GN0 is None else GN0
    tail = gN0*(arb('1.03883')*N0 - psi) + arb('1.03883')*(Istep + Imain.real + Ifar)
    Cz = K['small'] + K['med'] + tail
    est_paper = K['gridf'] + float(K['med'].mid()) + float(tail.mid())
    res[K['name']] = dict(small=K['small'].upper().str(12), med=K['med'].str(12), tail=tail.str(12), Imain=Imain.real.str(12),
                          LS=LB, steps=bool(bad >= 0), Istep=Istep.str(8), nu_intervals=nint, Cz_upper=float(Cz.upper()), printed=T2.get(K['name']),
                          margin=T2.get(K['name'], 0) - float(Cz.upper()), est_paper_4_11_float=est_paper)
    print('%-5s Cz <= %.10f (small %.10f med %.10f tail %.10f [Imain %.8f Ifar %.2e gN0-term %.3e], g decreasing from L = %.2f, nu-intervals %d) printed %.8f margin %.3e | float est. of the paper (4.11) with own tail: %.8f' % (
        K['name'], float(Cz.upper()), float(K['small'].upper()), float(K['med'].upper()), float(tail.upper()), float(Imain.real.mid()), float(Ifar.upper()),
        float((gN0*(arb('1.03883')*N0 - psi)).upper()), LB, nint, T2.get(K['name'], float('nan')), T2.get(K['name'], 0) - float(Cz.upper()), est_paper)); sys.stdout.flush()
json.dump(res, open('Riemann/Deep/STC2/cz_arb_%s.json' % ('%.0e' % N0).replace('+', ''), 'w'), indent=1)
print('done %.0f s' % (time.time()-t0))

