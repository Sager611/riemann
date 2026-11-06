# Certificate: |zeta(1/2+it)| <= A t^(27/164) for all t >= 3 (see REPORT.md).
# uv run --no-project --with python-flint --with numpy python3 certify.py A [Lmax]
# Ball arithmetic (python-flint arb, 128 bits). Every final comparison is arb "certainly <".
import sys, time
import numpy as np
from flint import arb, ctx
ctx.prec = 128
PI = arb.pi()
def fr(p, q): return arb(p)/arb(q)
def up(x):  # upper end of a ball as a float (rounded up) for reporting only
    return float(x.mid() + x.rad()) * (1 + 1e-15)

# ---------- Yang's explicit kth derivative test constants (2302.13444v1 Lemma 1.4 = Yang 2301.03165 Lem 2.4/2.5)
def A3B3(eta, h):
    lam0 = (1/eta + 32*eta.sqrt()*h/(15*PI.sqrt()))**(-3)
    d3 = (fr(1, 2) + (1 + fr(3, 8)*PI.sqrt()*eta**fr(3, 2)).sqrt()/2).sqrt()
    l03 = lam0**fr(1, 3)
    A3 = (1/(eta*h) + 32/(15*PI.sqrt())*(eta + l03).sqrt() + (eta + l03)*l03/3).sqrt()*d3
    B3 = arb(32).sqrt()/(arb(3).sqrt()*PI**fr(1, 4)*eta**fr(1, 4))*d3
    return A3, B3
def AkBk(k, eta, h):
    A_, B_ = A3B3(eta, h)
    for j in range(3, k):
        J = 2**(j - 1)
        d = (1 + 2/arb(2337)**(1 - fr(2, J))*(9*PI/1024*eta)**fr(1, J)).sqrt()
        A_ = d*(h**(-fr(1, J)) + arb(2)**fr(19, 12)*(J - 1)/arb((2*J - 1)*(4*J - 3)).sqrt()*A_.sqrt())
        B_ = d*arb(2)**fr(3, 2)*(J - 1)/arb((2*J - 3)*(4*J - 5)).sqrt()*B_.sqrt()
    return A_, B_

LAMC = 53632*arb(2).sqrt()*PI**4/729            # lambda_5 = LAMC a^9/(t^4 r^4)   (TeX l. 714-718)
HC5 = 76545*arb(2).sqrt()/107264                # h_5 = HC5 h^9
B1 = (1/PI)*LAMC**fr(1, 30)*HC5**fr(1, 8)       # b_1 (TeX l. 735)
B2 = PI**(-fr(7, 8))*(1/LAMC)**fr(1, 30)        # b_2 (TeX l. 738)
K3 = fr(1800, 2911)*B1*PI.sqrt()                # 0.47500...
K4 = fr(28800, 54481)*B2*PI.sqrt()              # 0.25303...
E6 = 6 + 5/PI*arb(2).log()
LOG2PI = (2*PI).log()

# monomial exponents: (p, e) -> a^p t^e with a = e^{alpha_j} t^{1/2} (s-coordinates): t-exponent p/2+e
P_ABA = {'T1': (fr(5, 41), -fr(5, 82)), 'T3': (-fr(17, 164), fr(7, 164)), 'T4': (fr(87, 82), -fr(57, 82)),
         'T5a': (-fr(46, 123), -fr(7, 82)), 'T5b': (-fr(5, 123), -fr(7, 82)), 'T7a': (fr(51, 41), -fr(65, 82)),
         'T7b': (fr(36, 41), -fr(49, 82)), 'F1': (-fr(5, 41), -fr(11, 41))}
P_4TH = {'D1': (fr(3, 14), fr(1, 14) - fr(27, 164)), 'D2': (fr(15, 28), -fr(1, 14) - fr(27, 164))}
P_TRV = (fr(1, 2), -fr(27, 164))

THETAS = [0.03*1.07**k for k in range(81)]     # theta grid (q = floor(theta a^{36/41} t^{-11/41})); floats are exact balls
ETA1 = [0.75*1.1**k for k in range(23)]         # eta_1 grid (Lemma 1.4, k = 5)
ETA2 = [0.5 + k/16 for k in range(17)]          # eta_2 grid (Lemma 1.4, k = 4)
GROUP = int(__import__('os').environ.get('SW_GROUP', '8'))                                       # blocks per parameter group (float search on the middle block)
TH_ARB = {th: dict(m1=1/arb(th), p1130=arb(th)**fr(11, 30), p61120=arb(th)**fr(61, 120), mh=1/arb(th).sqrt(),
                   p13=arb(th)**fr(1, 3), m1930=arb(th)**(-fr(19, 30)), m1130=arb(th)**(-fr(11, 30)), p1=arb(th)) for th in THETAS}

def region_consts(H):
    """H: arb, upper bound for b/a of every block in the region (1 < H <= 2)."""
    assert H > 1 and H < 2
    amin = arb(2336)/(2*H*H*(H - 1))
    kap = (amin/(amin - 1))**fr(9, 2)            # h_5 inflation for nu in (alpha-1, alpha]
    H5p = HC5*H**9*kap
    aba = {}
    for eta in ETA1:
        A5, B5 = AkBk(5, arb(eta), H5p)
        aba[eta] = dict(c2=K3*A5*kap**fr(1, 8)*H**fr(21, 8)*(H - 1), c3=K4*B5*H**fr(3, 2)*(H - 1)**fr(7, 8),
                        c7a=fr(1800, 451)*PI.sqrt()*H**fr(3, 2)*A5*H5p**fr(1, 8)*LAMC**fr(1, 30),
                        c7b=fr(1800, 931)*PI.sqrt()*H**fr(3, 2)*B5/LAMC**fr(1, 30))
    c4 = 8*fr(4686, 1000)*PI.sqrt()*H**fr(3, 2)/3
    c5a = fr(9, 14)*fr(893, 1000)*(5/PI)*((H - 1)/PI)**fr(1, 3)
    c5b = fr(9, 14)*(6/PI)*H**3*(H - 1)
    fth = {}
    for eta in ETA2:
        A4, B4 = AkBk(4, arb(eta), H**4)
        fth[eta] = (A4*(3/PI)**fr(1, 14)*H**fr(5, 7)*(H - 1), B4*(PI/3)**fr(1, 14)*H**fr(2, 7)*(H - 1)**fr(3, 4))
    return dict(H=H, aba=aba, c4=c4, c5a=c5a, c5b=c5b, fth=fth, trv=2*(H.sqrt() - 1))

def floats_of(rc):
    f = lambda x: float(x.mid())
    return dict(H=f(rc['H']), aba={e: {k: f(v) for k, v in d.items()} for e, d in rc['aba'].items()}, c4=f(rc['c4']),
                c5a=f(rc['c5a']), c5b=f(rc['c5b']), fth={e: (f(a), f(b)) for e, (a, b) in rc['fth'].items()}, trv=f(rc['trv']))

def choose(alph, L0, L1, rf):
    """float model: per block (alpha_j array) choose method/eta/theta. returns list of (meth, eta, theta)."""
    E6f = float(E6.mid()); H = rf['H']
    la0 = alph + L0/2; la1 = alph + L1/2           # log a at L0, L1
    th = np.array(THETAS)[None, :]; lth = np.log(th)
    best = np.full(alph.shape, np.inf); arg = [None]*len(alph)
    # ABA at L0
    la = la0[:, None]; L = L0
    lQ = lth + 36/41*la - 11/41*L
    ok = lQ > np.log(2.0)
    F1 = (H - 1) + np.exp(lth - 5/41*la - 11/41*L)
    g = 1/(1 - np.exp(-lQ))
    for eta in ETA1:
        d = rf['aba'][eta]
        F2 = ((H - 1)/th*g + d['c2']*th**(11/30))*np.exp(5/41*la - 5/82*L) + d['c3']*th**(61/120)*np.exp(-17/164*la + 7/164*L) \
            + rf['c4']/np.sqrt(th)*np.sqrt(g)*np.exp(87/82*la - 57/82*L) \
            + th**(1/3)*(rf['c5a']*np.exp(-46/123*la - 7/82*L) + rf['c5b']*np.exp(-5/123*la - 7/82*L)) + E6f*np.exp(-27/82*L) \
            + d['c7a']*th**(-19/30)*g**(19/30)*np.exp(51/41*la - 65/82*L) + d['c7b']*th**(-11/30)*g**(11/30)*np.exp(36/41*la - 49/82*L)
        v = np.where(ok, np.sqrt(F1*F2), np.inf)
        i = np.argmin(v, axis=1); m = v[np.arange(len(i)), i]
        for b in np.nonzero(m < best)[0]:
            arg[b] = ('aba', eta, THETAS[i[b]])
        best = np.minimum(best, m)
    # 4th at L1
    for eta in ETA2:
        c1, c2 = rf['fth'][eta]
        v = c1*np.exp(3/14*la1 + (1/14 - 27/164)*L1) + c2*np.exp(15/28*la1 + (-1/14 - 27/164)*L1)
        for b in np.nonzero(v < best)[0]:
            arg[b] = ('4th', eta, None)
        best = np.minimum(best, v)
    v = rf['trv']*np.exp(0.5*la1 - 27/164*L1)
    for b in np.nonzero(v < best)[0]:
        arg[b] = ('trv', None, None)
    best = np.minimum(best, v)
    return arg, best

def interval_bound(L0, L1, Dhi=fr(12, 1000), Dlo=fr(5, 100), ghi=0.385, glo=0.245, verbose=False):
    """Upper bound (arb) for |zeta(1/2+it)| t^{-27/164}, t in [e^L0, e^L1]. L0, L1 dyadic floats."""
    Lb = arb((L0 + L1)/2, (L1 - L0)/2)              # contains [L0, L1]
    aL0, aL1 = arb(L0), arb(L1)
    Lf = {}
    def LF(c):
        k = (c.mid().str(30))
        if k not in Lf: Lf[k] = (c*Lb).exp()
        return Lf[k]
    # block system: ell_j(L) = -j*Dhi - log(2pi)/2 + L/2 (j<=J), then Dlo steps up to J2
    J = int(np.ceil(((0.5 - ghi)*L1 - 0.5*float(LOG2PI.mid()))/float(Dhi.mid())))
    J2 = J + int(np.ceil((ghi - glo)*L1/float(Dlo.mid())))
    alpha = lambda j: -(arb(min(j, J))*Dhi + arb(max(j - J, 0))*Dlo) - LOG2PI/2
    ell_lo = alpha(J2) + aL0/2
    assert ell_lo > 10
    ff = lambda p: (1 - (-ell_lo).exp())**p      # floor factor for negative p (uniform over blocks)
    Hhi = Dhi.exp()/(1 - (-(alpha(J) + aL0/2)).exp())
    Hlo = Dlo.exp()/(1 - (-ell_lo).exp())
    rhi, rlo = region_consts(Hhi), region_consts(Hlo)
    tot = arb(0)
    parts = {'aba': arb(0), '4th': arb(0), 'trv': arb(0)}
    for (j0, j1, rc) in ((1, J, rhi), (J + 1, J2, rlo)):
        rf = floats_of(rc)
        js = list(range(j0, j1 + 1))
        mids = [js[min(i + GROUP//2, len(js) - 1)] for i in range(0, len(js), GROUP)]
        alph = np.array([float(alpha(j).mid()) for j in mids])
        argm, _ = choose(alph, L0, L1, rf)
        arg = [argm[i//GROUP] for i in range(len(js))]
        H = rc['H']
        for idx, j in enumerate(range(j0, j1 + 1)):
            meth, eta, th = arg[idx]
            al = alpha(j)
            if meth == 'aba':
                d = rc['aba'][eta]; tp = TH_ARB[th]
                M = {k: (p*al).exp()*LF(p/2 + e) for k, (p, e) in P_ABA.items()}
                ellj0 = al + aL0/2
                Qmin = tp['p1']*(fr(36, 41)*((ellj0.exp() - 1).log()) - fr(11, 41)*aL0).exp()
                assert Qmin > 2, (L0, j, th)
                g = 1/(1 - 1/Qmin)
                F2 = ((H - 1)*tp['m1']*g + d['c2']*tp['p1130'])*M['T1'] + d['c3']*tp['p61120']*M['T3']*ff(-fr(17, 164)) \
                    + rc['c4']*tp['mh']*g.sqrt()*M['T4'] + tp['p13']*(rc['c5a']*M['T5a']*ff(-fr(46, 123)) + rc['c5b']*M['T5b']*ff(-fr(5, 123))) \
                    + E6*LF(-fr(27, 82)) + d['c7a']*tp['m1930']*g**fr(19, 30)*M['T7a'] + d['c7b']*tp['m1130']*g**fr(11, 30)*M['T7b']
                F1 = (H - 1) + tp['p1']*M['F1']*ff(-fr(5, 41))
                assert F1 > 0 and F2 > 0
                v = (F1*F2).sqrt()
            elif meth == '4th':
                c1, c2 = rc['fth'][eta]
                v = sum(c*(p*al).exp()*LF(p/2 + e) for c, (p, e) in zip((c1, c2), P_4TH.values()))
            else:
                p, e = P_TRV
                v = rc['trv']*(p*al).exp()*LF(p/2 + e)
            parts[meth] += v
            tot += v
    # remainder n <= a_{J2}(t) <= exp(ell_{J2}(L1)) : trivial, increasing in t -> at L1
    rem = 2*((alpha(J2) + aL1/2).exp() + fr(1, 2)).sqrt()*(-fr(27, 164)*aL1).exp()
    rs = (fr(148, 100)*(-aL0/4).exp() + fr(127, 1000)*(-3*aL0/4).exp())*(-fr(27, 164)*aL0).exp()
    B = 2*(tot + rem) + rs
    if verbose:
        print('  L=[%g,%g] J=%d J2=%d Hhi=%.6f  aba %.4f 4th %.4f trv %.4f rem %.2e  -> %.4f' % (L0, L1, J, J2, up(Hhi), up(2*parts['aba']), up(2*parts['4th']), up(2*parts['trv']), up(2*rem), up(B)))
    return B

# ---------- tail t >= e^LT: Patel-Yang Lemmas 3.2-3.5 structure (fixed geometric blocks), with the corrections of REPORT.md
def tail_bound(LT, h, h2, eta1, eta2, th1, th2, th3, verbose=False):
    T = arb(LT).exp(); aLT = arb(LT)
    h, h2, eta1, eta2, th1, th2, th3 = [arb(x) for x in (h, h2, eta1, eta2, th1, th2, th3)]
    C0 = 2*(th3*(1 + 1/(2*th3*T**fr(27, 82)))).sqrt()                       # S_1 (Lemma 3.2, theta_3 restored)
    h3 = h2/(1 - h2/(th3*T**fr(27, 82)))                                     # ratio bound of S_2 blocks (theta_3, not theta_2)
    A4, B4 = AkBk(4, eta2, h3**4)                                            # Lemma 1.4 with h = h3^4 (ratio of |f''''|)
    D3 = (3/PI)**fr(1, 14)*A4*h3**fr(5, 7)*(h3 - 1)*th2**fr(3, 14)/(h2**fr(3, 14) - 1)
    D4 = (PI/3)**fr(1, 14)*B4*h3**fr(2, 7)*(h3 - 1)**fr(3, 4)*th2**fr(15, 28)/(h2**fr(15, 28) - 1)
    eps = h/(th2*T**fr(7, 17))
    h0 = h/(1 - eps)                                                         # ratio bound of S_3 blocks
    assert h0 > 1 and h0 < 2
    th2m = th2/h0                                                            # theta_2 used in Lemma 3.4: th2m t^{7/17} < a_K
    assert th2m > 1
    q0 = th1*th2m**fr(7, 17)*(fr(65, 697)*aLT).exp()
    assert q0 > 2
    al = (h0 - 1 + th1/(th2m**fr(5, 41)*(fr(222, 697)*aLT).exp())).sqrt()
    amin = arb(2336)/(2*h0*h0*(h0 - 1))
    kap = (amin/(amin - 1))**fr(9, 2)                                        # h_5 inflation (nu in (alpha-1, alpha])
    H5 = HC5*h0**9*kap
    A5, B5 = AkBk(5, eta1, H5)
    C1 = al*((h0 - 1)/th1/(1 - 1/q0) + K3*th1**fr(11, 30)*A5*kap**fr(1, 8)*h0**fr(21, 8)*(h0 - 1)).sqrt()
    # N <= beta-alpha+1 correction terms (vii-a), (vii-b), in Lemma 3.4 form (sqrt split, q >= theta_1 a^{36/41} t^{-11/41} - 1)
    E7a = al*(fr(1800, 451)*PI.sqrt()*h0**fr(3, 2)*A5*H5**fr(1, 8)*LAMC**fr(1, 30)*th1**(-fr(19, 30))*(1 - 1/q0)**(-fr(19, 30))).sqrt()
    E7b = al*(fr(1800, 931)*PI.sqrt()*h0**fr(3, 2)*B5/LAMC**fr(1, 30)*th1**(-fr(11, 30))*(1 - 1/q0)**(-fr(11, 30))).sqrt()
    C2 = al*(K4*th1**fr(61, 120)*B5*h0**fr(3, 2)*(h0 - 1)**fr(7, 8)).sqrt()
    E1 = al*(8*fr(4686, 1000)*PI.sqrt()/3).sqrt()*h0**fr(3, 4)*(th1 - th2m**(-fr(36, 41))*(-fr(65, 697)*aLT).exp())**(-fr(1, 4))
    E2 = al*(fr(9, 14)*th1**fr(1, 3)*(fr(893, 1000)*5/PI*((h0 - 1)/PI)**fr(1, 3)/(th2m**fr(1, 3)*(fr(7, 51)*aLT).exp()) + 6/PI*h0**3*(h0 - 1))).sqrt()
    E3 = al*E6.sqrt()
    mu3 = lambda a: (2*PI)**(-a/2)/(h**a - 1)
    mu4 = lambda a: (1 - eps)**(-a)*(h/th2)**a/(1 - h**(-a))
    Kb = (fr(3, 34)*aLT - (th2*(2*PI).sqrt()).log())/h.log() + 1      # log-affine envelope of K(t), at T
    assert Kb > fr(3, 34)/(fr(27, 164)*h.log())                          # => Kb(t) t^{-27/164} decreasing for t >= T
    terms = dict(S1=C0, S2a=D3*((fr(19, 119) - fr(27, 164))*aLT).exp(), S2b=D4*((fr(71, 476) - fr(27, 164))*aLT).exp(),
                 C1=C1*mu3(fr(5, 82)), C2=C2*mu4(fr(17, 328)), E1=E1*mu3(fr(87, 164))*(-fr(27, 328)*aLT).exp(),
                 E2=E2*mu4(fr(5, 246))*(-fr(427, 8364)*aLT).exp(), E3=E3*Kb*(-fr(27, 164)*aLT).exp(),
                 E7a=E7a*mu3(fr(51, 82))*(-fr(14, 164)*aLT).exp(), E7b=E7b*mu3(fr(18, 41))*(-fr(13, 164)*aLT).exp())
    B = 2*sum(terms.values())
    if verbose: print('  tail t >= e^%g:' % LT, {k: round(up(2*v), 5) for k, v in terms.items()}, '->', up(B))
    return B

TAIL = (2000, 1.006973, 1.000094, 1.12054, 0.87152, 1.51415, 5.82204e14, 5.40246e-73)

def main(Atgt, Lmax=2000.0, L_hpy=75.0):
    t0 = time.time(); A = arb(Atgt)
    print('target A =', Atgt)
    # [3, e^75]: Hiary-Patel-Yang 2207.02366 Thm 1.1: |zeta| <= 0.618 t^{1/6} log t (t >= 3); ratio increasing in t
    hpy = fr(618, 1000)*(arb(L_hpy)/492).exp()*arb(L_hpy)
    print('part H: t in [3, e^%g]: 0.618 e^{L/492} L <= %.5f' % (L_hpy, up(hpy))); assert hpy < A
    worst = (0, None); n = 0
    L = L_hpy; sched = []
    while L < Lmax:
        d = 25.0 if L < 300 else (8.0 if L < 500 else (2.0 if L < 950 else (8.0 if L < 1300 else 50.0)))
        sched.append((L, min(L + d, Lmax))); L = min(L + d, Lmax)
    for (a, b) in sched:
        todo = [(a, b)]
        while todo:
            L0, L1 = todo.pop()
            B = interval_bound(L0, L1)
            n += 1
            if not (B < A):
                assert L1 - L0 > 0.2, ('FAIL', L0, L1, up(B))
                m = (L0 + L1)/2; todo += [(m, L1), (L0, m)]; continue
            if up(B) > worst[0]: worst = (up(B), (L0, L1))
        if int(a) % 100 < int(b) % 100 or b >= Lmax:
            pass
        else:
            print('  intervals up to L=%g ok (%d evaluated), worst so far %.5f on %s, %.0f s' % (b, n, worst[0], worst[1], time.time() - t0)); sys.stdout.flush()
    print('part C: all intervals in [%g, %g] certified (%d evaluations); max bound %.5f on L in %s' % (L_hpy, Lmax, n, worst[0], worst[1]))
    Bt = tail_bound(*TAIL, verbose=True)
    print('part T: t >= e^%g: %.5f' % (TAIL[0], up(Bt))); assert Bt < A
    assert TAIL[0] <= Lmax
    print('ALL PARTS OK: |zeta(1/2+it)| <= %s t^(27/164) for t >= 3   (%.0f s)' % (Atgt, time.time() - t0))

if __name__ == '__main__':
    if len(sys.argv) > 1 and ',' in sys.argv[1]:
        t0 = time.time()
        for a in sys.argv[1:]:
            L0, L1 = [float(x) for x in a.split(',')]
            B = interval_bound(L0, L1, verbose=True)
            print('   time %.1f s' % (time.time() - t0)); sys.stdout.flush()
    elif len(sys.argv) > 1 and sys.argv[1] == 'tail':
        tail_bound(*TAIL, verbose=True)
    else:
        main(sys.argv[1] if len(sys.argv) > 1 else '55.5', float(sys.argv[2]) if len(sys.argv) > 2 else 2000.0)
