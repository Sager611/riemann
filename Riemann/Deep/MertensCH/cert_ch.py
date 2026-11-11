# cert_ch.py -- |M(x)| <= c x^theta for ALL x >= 1: Chirre-Helfgott black box (x <= X0), hybrid contour
# (RH to Platt-Trudgian height on Re s = sigma_L for |t| <= T', half-plane certificate of Mertens78 on
# Re s = sigma_H for |t| >= T') on [X0, X1], monotone final range [X1, oo).  See DERIVATION.md.
# Float (numpy) code only chooses parameters; every reported constant is re-evaluated in Arb (python-flint).
# Run: uv run --no-project --with python-flint --with numpy python3 cert_ch.py [thetas...]
import json, math, sys, os, time
import numpy as np
from flint import arb, ctx
HERE = os.path.dirname(os.path.abspath(__file__))
M78 = os.path.join(HERE, '..', 'Mertens78') if os.path.exists(os.path.join(HERE, '..', 'Mertens78')) else \
      '/Users/adriansager/playground/riemann/Riemann/Deep/Mertens78'
sys.path.insert(0, M78)
from cert import Cert, J0_direct          # Mertens78 checker (B = 7/8 certificates, Arb)
ctx.prec = 120
PI = arb.pi(); LN10 = math.log(10)
# ---------------- fixed inputs (DERIVATION.md section 0) ----------------
T_PT = arb(3000175332800)                  # Platt-Trudgian: zeros with 0 < Im rho <= T_PT are on Re s = 1/2
TPS = os.environ.get('MCH_TP', '3000000000000')     # T' (contour switch height); MCH_TP=1000 gives the routine baseline
TP = arb(int(TPS))
T1 = arb(1000)                             # J_L0 computed on [0, T1]
SIGL = arb("0.85")                          # sigma_L
S0P = arb("1.02"); ETAP = arb("0.1")      # RH-strip Harnack: centre line sigma_0', Rademacher eta'
BETA = 3/(2*PI**2)                          # squarefree unsmoothing coefficient (Lemma H0)
SQ = arb(6)                                 # 2 * 3 (|R(y)| <= 3 sqrt y), Lemma H0
CH_A = 3/(PI*arb(10)**10); CH_B = arb("11.39")   # Chirre-Helfgott Cor. 1.5 (black box), x >= 1
assert TP + (S0P - arb(1)/2) < T_PT and TP >= T1
# ---------------- RH-strip Harnack constants (Lemma H1) ----------------
KAP = (arb(1)/2 + ETAP)/2
AP = (3*(1 + ETAP).zeta()).log() + KAP*(arb(2)/3).log()
assert AP > 0
LZ = (S0P.zeta()/(2*S0P).zeta()).log(); assert LZ > 0
def LamL(s): return (2*S0P - s - arb(1)/2)/(s - arb(1)/2)
def logCL(s, s_in_R): return (LamL(s) - 1)*AP + LamL(s)*LZ + (s_in_R + LamL(s)*((1 + KAP)*(1 + S0P) + 2))/T1
LAML = LamL(SIGL); AL = (LAML - 1)*KAP
CL = logCL(SIGL, SIGL).exp()               # |1/zeta(sigma_L+it)| <= CL t^AL, T1 <= t <= T'
def CLstar(sH): return logCL(SIGL, sH).exp()   # uniform on sigma in [sigma_L, sigma_H] at t = T'
# ---------------- Mertens78 certificates (half-plane Re s > 7/8) ----------------
CFG = json.load(open(os.path.join(M78, 'all.json')))
M78ROWS = json.load(open(os.path.join(M78, 'results.json')))   # Mertens78 record rows (comparison only)
def load(cf):
    if cf.get('harnack'): d = dict(cf)
    else:
        d = json.load(open(os.path.join(M78, cf['cert']))); d['eta'] = cf['eta']; d['T1'] = '1000'
    return Cert(d)
def main():
    t0 = time.time()
    certs = []
    BPMODE = '--bp' in sys.argv
    pairs = [(cf['name'], load(cf)) for cf in CFG]
    if BPMODE:                                # add the B' = 0.87338... translated certificates (bprime.py)
        import bprime
        pairs += [(cf['name'] + '_Bp', bprime.make(bprime.translate(cf))) for cf in CFG]
    for nm, c in pairs:
        cf = dict(name=nm)
        lC = c.tail(TP)                       # log C_inf with T_1 := T' (Corollary A4)
        assert c.a > 0 and c.a < 1 and c.s1 > SIGL
        assert float(c.d["sigma1"]) >= (0.888 if nm.endswith('_Bp') else 0.89)   # B2 contour shift (bprime.py asserts a' < 2)
        certs.append(dict(name=cf['name'], c=c, a=c.a, s=c.s1, lC=lC, CS=CLstar(c.s1),
                          af=float(c.a.mid()), sf=float(c.s1.mid()), lCf=float(lC.upper()),
                          lCSf=float(CLstar(c.s1).log().upper()), thstar=float(((c.s1 + c.a)/(1 + c.a)).upper())))
    jl = os.path.join(HERE, 'logs', 'JL0.json')
    if os.path.exists(jl): JL0 = arb(json.load(open(jl))['JL0_upper'])
    else:
        JL0, n, mx = J0_direct("0.85", "1000"); JL0 = arb(JL0.upper())
        os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
        json.dump(dict(sigma_L=str(SIGL), T1=1000, JL0_upper=str(JL0.upper()), n=n, max_inv_zeta=mx), open(jl, 'w'))
    RHP = JL0 + 2*CL*TP**AL/AL                # (X^sigma_L / 2 pi) * RHP bounds the Re s = sigma_L piece
    print("RH strip: Lam_L", LAML.str(5), "a_L", AL.str(5), "log C_L", CL.log().str(5), "J_L0", JL0.str(5),
          "log RHP", RHP.log().str(5), "secs", round(time.time() - t0, 1), flush=True)
    consts = dict(RHPf=float(RHP.upper()), lTP=math.log(int(TPS)))
    thetas = [float(t) for t in sys.argv[1:] if not t.startswith('--')] or [0.9409, 0.94, 0.9382, 0.937, 0.935, 0.9325, 0.93, 0.9296, 0.9277,
              0.925, 0.92, 0.9192, 0.915, 0.9121, 0.91, 0.9092, 0.905, 0.9024, 0.9, 0.8992]
    outf = os.path.join(HERE, 'logs', 'cert_ch_results%s%s.json' % ('_bp' if BPMODE else '', '' if TPS == '3000000000000' else '_T' + TPS))
    out = json.load(open(outf)) if ('--append' in sys.argv and os.path.exists(outf)) else []
    out = [o for o in out if all(abs(o['theta'] - t) > 1e-9 for t in thetas)]
    for th in thetas:
        plan = make_plan(th, certs, consts)
        if plan is None: print(th, "no plan"); continue
        res = certify(th, plan, certs, RHP)
        res['theta'] = th; out.append(res)
        json.dump(out, open(outf, 'w'), indent=1)
        rec = [r['allx']['Cprime_log10'] for r in M78ROWS if r['theta'] <= th + 1e-12]
        res['record_row_log10'] = min(rec) if rec else None
        print("theta %.4f  c = %s (log10 %.4f; float plan %.4f)  record row log10 %s  X0 1e%.2f  X1 1e%.2f  ranges %d  final %s  [%.0fs]" %
              (th, res['c'], res['log10c'], plan['log10c'], ("%.3f" % min(rec)) if rec else "-", plan['L0'], plan['L1'],
               len(plan['ranges']), plan['fin'][0], time.time() - t0), flush=True)
    json.dump(out, open(outf, 'w'), indent=1)

# ---------------- float planning (untrusted) ----------------
def G(la, ltau, lT):                          # log of int_T^oo t^a min(1/t, tau/t^2) dt
    small = ltau + (la - 1)*lT - np.log(1 - la)
    with np.errstate(over='ignore', invalid='ignore'):
        big = np.log(np.maximum((np.exp(la*ltau) - np.exp(la*lT))/la, 0) + np.exp(la*ltau)/(1 - la))
    big = np.where(la*ltau > 600, la*ltau + np.log(1/la + 1/(1 - la)), big)
    return np.where(ltau <= lT, small, big)
def lse(*v):
    m = np.maximum.reduce(v); return m + np.log(sum(np.exp(x - m) for x in v))
def rng_bound(th, ce, consts, Lk, Lk1, lch):
    """log of the range constant on [10^Lk, 10^Lk1] (vectorised), mirrors certify()."""
    lx = Lk*LN10; lx1 = Lk1*LN10; s = ce['sf']; a = ce['af']; sl = float(SIGL.mid())
    lr = np.log1p(np.exp(lch + (th - 1)*lx)/2)            # log(1 + (c_h/2) x_k^(theta-1))
    la_ = np.log(BETA_F) + lch
    lb = np.log(6.0) + (0.5 - th)*lx + lr/2
    ld = (sl - th)*lx + sl*lr + np.log(consts['RHPf']/(2*math.pi))
    ldec = (s - th)*(lx if s <= th else lx1) + s*lr
    ltau1 = np.logaddexp(math.log(2) + (1 - th)*lx1 - lch, 0.0)                          # no overflow cap
    lT = consts['lTP']
    le = np.log(max(s - sl, 1e-300)) + ce['lCSf'] + float(AL.upper())*lT - math.log(math.pi) + ldec + np.minimum(-lT, ltau1 - 2*lT)
    lc = ce['lCf'] - math.log(math.pi) + ldec + G(a, ltau1, lT)
    v = lse(la_, lb, ld, le, lc)
    return np.where(lch <= (1 - th)*lx, v, np.inf)
def fin_bound(th, ce, consts, Lk, lch):
    lx = Lk*LN10; s = ce['sf']; a = ce['af']; sl = float(SIGL.mid())
    lr = np.log1p(np.exp(lch + (th - 1)*lx)/2)
    la_ = np.log(BETA_F) + lch
    lb = np.log(6.0) + (0.5 - th)*lx + lr/2
    ld = (sl - th)*lx + sl*lr + np.log(consts['RHPf']/(2*math.pi))
    ldec = (s - th)*lx + s*lr
    lT = consts['lTP']
    le = np.log(s - sl) + ce['lCSf'] + float(AL.upper())*lT - math.log(math.pi) + ldec - lT
    lq = (s + a*(1 - th) - th)*lx + s*lr + a*np.logaddexp(math.log(2) - lch, (th - 1)*lx)   # no underflow
    lc = ce['lCf'] - math.log(math.pi) + math.log(1/a + 1/(1 - a)) + lq
    v = lse(la_, lb, ld, le, lc)
    return np.where(lch <= (1 - th)*lx, v, np.inf)
def golden(f, lo, hi, it=45):
    g = (math.sqrt(5) - 1)/2
    a, b = lo.copy(), hi.copy()
    for _ in range(it):
        c = b - g*(b - a); d = a + g*(b - a); fc = f(c); fd = f(d)
        m = fc < fd
        b = np.where(m, d, b); a = np.where(m, a, c)
    x = (a + b)/2; return x, f(x)
STEP = 0.05; DEN = 10**6; NUM = 50000
def make_plan(th, certs, consts):
    ext = np.arange(3000.0, 100000.0, 1.0) if th < 0.91 else np.zeros(0)       # long ranges only where needed (theta < 0.91)
    Lh = np.concatenate([np.arange(round(60/STEP), round(3000/STEP))*STEP, ext])  # range grid (log10 x)
    Lg = np.exp(np.arange(math.log(100000.0 if th < 0.91 else 3000.0), math.log(10**6.3), 0.0005))
    Lall = np.concatenate([Lh, Lg])
    Lnext = np.append(Lall[1:], Lall[-1] + 1.0)
    # CH part
    lxa = Lall*LN10
    chp = np.logaddexp(math.log(float(CH_A.upper())) + (1 - th)*lxa, math.log(11.4) + (0.5 - th)*lxa)
    chp = np.maximum(chp, math.log(11.4))
    # best hybrid range constant at each grid point of Lh (range [L, L+STEP])
    cr = np.full(len(Lall), np.inf); cri = np.zeros(len(Lall), int); crh = np.zeros(len(Lall))
    for i, ce in enumerate(certs):
        lo = np.full(len(Lh), 0.0); hi = (1 - th)*Lh*LN10
        x, v = golden(lambda q: rng_bound(th, ce, consts, Lh, Lnext[:len(Lh)], q), lo, hi)
        m = v < cr[:len(Lh)]
        cr[:len(Lh)] = np.where(m, v, cr[:len(Lh)]); cri[:len(Lh)] = np.where(m, i, cri[:len(Lh)]); crh[:len(Lh)] = np.where(m, x, crh[:len(Lh)])
    cf = np.full(len(Lall), np.inf); cfi = np.zeros(len(Lall), int); cfh = np.zeros(len(Lall))
    for i, ce in enumerate(certs):
        if ce['sf'] + ce['af']*(1 - th) - th > -1e-9: continue
        lo = np.full(len(Lall), 0.0); hi = (1 - th)*Lall*LN10
        x, v = golden(lambda q: fin_bound(th, ce, consts, Lall, q), lo, hi)
        m = v < cf
        cf = np.where(m, v, cf); cfi = np.where(m, i, cfi); cfh = np.where(m, x, cfh)
    if not np.isfinite(cf).any(): return None
    def feas(lc):
        ok = np.nonzero(chp <= lc)[0]
        if len(ok) == 0: return None
        i0 = ok[-1]
        if i0 >= len(Lh): # beyond the hybrid grid: final must start at i0
            return (i0, i0) if cf[i0] <= lc else None
        bad = np.nonzero(cr[i0:] > lc)[0]; kb = i0 + (bad[0] if len(bad) else len(Lall) - 1 - i0)
        js = np.nonzero(cf[i0:kb + 1] <= lc)[0]
        return (i0, i0 + js[0]) if len(js) else None
    lo, hi = 0.0, 2000.0
    if feas(hi) is None: return None
    for _ in range(60):
        mid = (lo + hi)/2
        if feas(mid) is None: lo = mid
        else: hi = mid
    i0, j = feas(hi)
    ranges = [(int(round(Lall[k]*DEN)), int(round(Lall[k + 1]*DEN)), int(cri[k]), float(crh[k])) for k in range(i0, j)]
    return dict(log10c=hi/LN10, L0=float(Lall[i0]), L1=float(Lall[j]), ranges=ranges,
                fin=(certs[cfi[j]]['name'], int(cfi[j]), float(cfh[j]), float(Lall[j])))
BETA_F = 3/(2*math.pi**2)
# ---------------- Arb certification ----------------
def Garb(a, tau, T):
    if tau <= T: return tau*T**(a - 1)/(1 - a)
    if tau >= T: return (tau**a - T**a)/a + tau**a/(1 - a)
    return tau**a*(1/a + 1/(1 - a))          # >= G for every tau > 0 (DERIVATION.md, Lemma H3)
def L2x(Lnum, Lden): return (arb(Lnum)/Lden*arb(10).log()).exp()   # 10^(Lnum/Lden), exact rational exponent
def dyadic(lv):  # exact short decimal for exp(lv) (parameters need not be optimal, only exact); no overflow
    l10 = lv/LN10; e = math.floor(l10); return arb("%.6fe%d" % (10**(l10 - e), e))
def certify(th, plan, certs, RHP):
    TH = arb("%.4f" % th) if abs(th*1e4 - round(th*1e4)) < 1e-9 else arb(repr(th))
    lines = []
    # [1, X0]: CH black box; f(x) = CH(x)/x^theta is convex in log x -> max at the endpoints
    X0 = L2x(int(round(plan['L0']*DEN)), DEN)
    cCH = CH_A*X0**(1 - TH) + CH_B*X0**(arb(1)/2 - TH)
    cCH1 = CH_A + CH_B
    pieces = [cCH1.upper(), cCH.upper()]
    lines.append(("CH", str(cCH.upper())))
    # hybrid ranges [x_k, x_k+1]
    rmax = arb(0)
    for (nlo, nhi, ci, lch) in plan['ranges']:
        ce = certs[ci]; s = ce['s']; a = ce['a']; Cinf = ce['lC'].exp()
        xk = L2x(nlo, DEN); xk1 = L2x(nhi, DEN)
        ch = dyadic(lch)
        assert ch <= xk**(1 - TH), "h <= x fails"
        assert xk >= 2*arb(10)**8
        r = 1 + (ch/2)*xk**(TH - 1)
        Ta = BETA*ch
        Tb = SQ*xk**(arb(1)/2 - TH)*r.sqrt()
        Td = xk**(SIGL - TH)*r**SIGL*RHP/(2*PI)
        dec = (xk**(s - TH) if s <= TH else xk1**(s - TH))*r**s
        tau1 = 2*xk1**(1 - TH)/ch + 1
        Te = (s - SIGL)*ce['CS']*TP**AL/PI*dec*(tau1/TP**2 if tau1 <= TP else 1/TP)
        Tc = Cinf/PI*dec*Garb(a, tau1, TP)
        ck = Ta + Tb + Td + Te + Tc
        if ck.upper() > rmax.upper(): rmax = ck          # exact comparison of upper endpoints (no float overflow)
    lines.append(("ranges", str(rmax.upper()) if plan['ranges'] else "none"))
    if plan['ranges']: pieces.append(rmax.upper())
    # final range [X1, oo)
    name, ci, lch, L1 = plan['fin']; ce = certs[ci]; s = ce['s']; a = ce['a']; Cinf = ce['lC'].exp()
    X1 = L2x(int(round(L1*DEN)), DEN)
    n0, n1 = int(round(plan['L0']*DEN)), int(round(L1*DEN))      # exact exponents: X0 = 10^(n0/DEN), X1 = 10^(n1/DEN)
    assert TH < 1 and SIGL <= TH                                   # Lemmas H4, H5 need sigma_L <= theta < 1
    assert n1 >= n0 and (n1 == plan['ranges'][-1][1] if plan['ranges'] else n1 == n0)   # no gap between X0 and X1
    assert all(lo_ < hi_ for lo_, hi_, _, _ in plan['ranges'])
    assert all(plan['ranges'][k + 1][0] == plan['ranges'][k][1] for k in range(len(plan['ranges']) - 1))   # contiguous
    if plan['ranges']: assert plan['ranges'][0][0] == n0
    ch = dyadic(lch)
    assert ch <= X1**(1 - TH) and X1 >= 2*arb(10)**8
    expo = s + a*(1 - TH) - TH; assert expo <= 0, "theta below theta* of final certificate"
    assert SIGL <= TH and s <= TH
    r = 1 + (ch/2)*X1**(TH - 1)
    Ta = BETA*ch; Tb = SQ*X1**(arb(1)/2 - TH)*r.sqrt(); Td = X1**(SIGL - TH)*r**SIGL*RHP/(2*PI)
    dec = X1**(s - TH)*r**s
    Te = (s - SIGL)*ce['CS']*TP**AL/PI*dec/TP
    Tc = Cinf/PI*(1/a + 1/(1 - a))*X1**expo*r**s*(2/ch + X1**(TH - 1))**a
    cfin = Ta + Tb + Td + Te + Tc
    lines.append(("final", name, str(cfin.upper())))
    pieces.append(cfin.upper())
    worst = pieces[0]
    for q in pieces:
        if q > worst: worst = q
    # reported constant: 3 significant digits, rounded up, then re-checked against every piece in Arb
    lw = float((arb(worst).upper().log()/arb(10).log()).upper()); e = math.floor(lw)
    m3 = math.ceil(10**(lw - e + 2)*(1 + 1e-12))         # no float overflow for c > 1e308
    crep = "%de%d" % (m3, e - 2); C = arb(crep)
    for q in [cCH, cCH1, rmax, cfin]: assert C >= q.upper(), "rounding check"
    lc = math.log10(m3) + e - 2
    return dict(log10c=lc, c=crep, pieces=lines, L0=plan['L0'], L1=plan['L1'], nranges=len(plan['ranges']),
                final=plan['fin'][0], plan_log10c=plan['log10c'])
if __name__ == '__main__':
    main()
