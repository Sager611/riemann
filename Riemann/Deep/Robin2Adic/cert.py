#!/usr/bin/env python3
# cert.py -- Robin2Adic: Arb upper bound for E* = sup_{x >= e^{108.57}} E(x) (E = Nicolas' quantity, Gracia Llorente Prop 2)
# and the 2-adic range K (Robin for every n > 5040 with v_2(n) <= K). See DERIVATION.md (sections "Ours").
# One command (repo root):
#   uv run --no-project --with python-flint python3 Riemann/Deep/Robin2Adic/cert.py
# Imports Riemann/Deep/Ingham4/cert.py as a library (unchanged): KLN rows with RH(3e12), CS rows with b4 raised.
import sys, os, time, json, importlib.util, resource
from flint import arb, ctx
ctx.prec = 128
sys.dont_write_bytecode = True                             # no __pycache__ next to the imported files
HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location('ingham4_cert', os.path.join(HERE, '..', 'Ingham4', 'cert.py'))
I4 = importlib.util.module_from_spec(spec); spec.loader.exec_module(I4)
up, lo = I4.up, I4.lo                      # exact arb >= (<=) every point of a ball
T_START = time.time(); LOG = []
def say(*a):
    s = ' '.join(str(v) for v in a); print(s, flush=True); LOG.append(s)
def s6(x, n=7): return x.str(n, radius=False)
def emin(a, b): return a if (a <= b) else (b if (b <= a) else None)   # exact (zero-radius) arbs only
QUICK = 'quick' in sys.argv
# ------------------------------------------------------------------ inputs (DERIVATION.md "Inputs")
PI = arb.pi(); EG = arb.const_euler(); LN2 = arb(2).log(); LN10 = arb(10).log()
B2 = (arb(777) - arb(1081).sqrt())/852                        # Richey (QuasiRiemannTracker PR 129): no zero with Re s > B2
TSTAR = B2 + arb('2.2e-254')                                  # Bebel (Lean): no zero with Re s > T*
B_BEBEL = arb('0.87338196661320752')                          # Bebel's decimal statement
NB = 2**44
_bu = int(float(up(B_BEBEL*NB).mid())) + 2
B_UP = arb(_bu)/NB                                            # exact dyadic, >= every threshold above
assert B_UP.rad() == 0 and B_UP > B_BEBEL and B_BEBEL > TSTAR and B_UP - B_BEBEL < arb('1e-12')
H0 = arb(3*10**12); h = H0.log()                              # Platt-Trudgian: RH for 0 < Im rho <= 3e12
GAMMA1 = arb('14.13')                                         # first zero ordinate 14.1347... > 14.13
S_LOW = 2 + EG - (4*PI).log()                                 # sum_rho 1/(rho(1-rho)) = 2 + gamma - log 4 pi; all terms have Re > 0
M0_100 = arb('2.4530e-12')                                    # BKLNW (Math. Comp. 2021) Table ThetaEpsilon, log X0 = 100: |theta(x)-x| < M0 x
M0_50 = arb('1.1199e-9')                                      # same table, log X0 = 50 (upper constant M0)
RS = arb('1.01624')                                           # Rosser-Schoenfeld 1962 Thm 9: theta(x) < 1.01624 x (x > 0)
HSW = (arb('0.1038'), arb('0.2573'), arb('9.3675'))           # Hasanalizade-Shen-Wong: |N(T) - T/2pi log(T/2pi e) - 7/8| <= a log T + b loglog T + c
Y_R = arb('46.79')                                            # RobinCH: Robin for 5041 <= n <= 10^(10^46.79)
L1 = arb(10857)/100                                           # claimed: every n > 10^(10^46.79) has x = p_k >= e^{L1}
I2 = (h*h/2 + h/2 + arb(1)/4)/H0**2                           # int_{H0}^oo t^-3 log^2 t dt
I1 = (h/2 + arb(1)/4)/H0**2                                   # int_{H0}^oo t^-3 log t dt
a_, b_, c_ = HSW
TAU = up((h - (2*PI).log() + 1)/(PI*H0) + 8*((a_ + b_)*(h/2 + arb(1)/4) + c_/2)/H0**2)   # >= sum_{|lam|>H0} lam^-2 (both signs)
def threshold_check():
    Lam = lo(arb(10)**Y_R*LN10)                               # log n > Lam for n > 10^(10^46.79)
    c1 = RS*2*arb(100).exp() + 100 + LN2 < Lam                # x <= e^100 impossible (Bertrand + RS)
    c2 = (1 + M0_100)*L1.exp() + L1 + LN2 < Lam               # x < e^{L1} impossible (BKLNW + Bertrand)
    say('threshold: log Lam = %s; x <= e^100 excluded: %s; x < e^%s excluded: %s' % (s6(Lam.log(), 8), c1, s6(L1, 6), c2))
    # Lemma 3's c = 0.998/(xL): (T+U) log(T+U) < xL/0.998 with T+U < x(1+M0) + log(2x); the ratio decreases in x, so check at L1
    x = L1.exp(); TU = x*(1 + M0_100) + L1 + LN2
    c3 = TU*TU.log() < x*L1/arb('0.998')
    say('Lemma 3 slot condition (T+U)log(T+U) < xL/0.998 at x = e^L1: %s (ratio %s)' % (c3, s6(TU*TU.log()/(x*L1), 15)))
    assert c1 and c2 and c3
    return Lam
# ------------------------------------------------------------------ density envelopes W(s) >= sum_{|lam|>H0, beta>s} lam^-2
def Ig(q, e): return (q*h).gamma_upper(e + 1)/q**(e + 1)     # int_{H0}^oo t^{-q-1} (log t)^e dt
KLN_P = None; CS_R = None
def W_kln(s):
    best = None
    for p in KLN_P:
        if not (s > p.smin): continue
        a = 8*(1 - s)/3; q = 2 - a; e = 5 - 2*s                # (log kt)^{2s}(log t)^{5-4s} <= (log t)^{5-2s} since 1e9/H0 <= k <= 1
        v = up(4*(p.logA(s).exp()*Ig(q, e) + p.Bc(s).max(arb(0))*I2))
        best = v if best is None else emin(best, v)
    return best
def W_cs(s, raw=False):
    best = None
    for (s1, s2, B1r, B1c, Bb2, Bb3) in CS_R:
        if not (s >= s1 and s <= s2): continue
        a = 3*(1 - s)/(2 - s); e = (7 - 5*s)/(2 - s); q = 2 - a
        v = up(4*((B1r if raw else B1c)*Ig(q, e) + Bb2*I2 + Bb3*I1))
        best = v if best is None else emin(best, v)
    return best
def grid():
    """exact dyadic s_0 = 1/2 < ... < s_J = B_UP: step 2^-9 below 0.8, 2^-16 from 0.8 (2^-12 / 2^-10 in quick mode)"""
    f1, f2 = (2**10, 2**12) if QUICK else (2**9, 2**16)
    g = [arb(j)/f1 for j in range(f1//2, int(0.8*f1))] + [arb(j)/f2 for j in range(int(0.8*f2), int(float(B_UP.mid())*f2) + 1)]
    g = [s for s in g if s < B_UP] + [B_UP]
    assert all(s.rad() == 0 for s in g) and all(g[i] < g[i+1] for i in range(len(g) - 1))
    return g
def envelope(G, Wfun):
    """W~_j: exact, nonincreasing, W~_0 = tau/2 (s_0 = 1/2), W~_j = min(tau/2, W(s_i), i <= j); returns c_j >= 0 (Abel form)."""
    Wt = [TAU/2]
    for s in G[1:-1]:
        v = Wfun(s)
        m = Wt[-1] if v is None else emin(Wt[-1], v)
        Wt.append(m)
    c = [TAU - 2*Wt[0]] + [2*(Wt[j-1] - Wt[j]) for j in range(1, len(Wt))] + [2*Wt[-1]]
    assert len(c) == len(G) and all(cj >= 0 for cj in c)
    return c, Wt
def Tsum(G, c, L):
    x = L.exp()
    return sum((cj*(x**(s - 1) + x**(-s))/2 for cj, s in zip(c, G)), arb(0))
# ------------------------------------------------------------------ the majorant of E(x), x = e^L (DERIVATION.md "Ours")
def Gdd_sup(L, ds):
    d = arb(0, ds)                                            # ball [-ds, ds]
    return up((L + 1 + (1 + d).log())/((1 + d)**2*(L + (1 + d).log())**2))
def majorant(L, G, c, detail=False):
    x = L.exp(); f = 1/L + 1/L**2
    eps_hi = 2*(L + 2)/(L*(L + 1)*H0); eps_lo = 2*(L + 2)/(L*(L + 1)*GAMMA1)
    cpsi = 1 + M0_50 + RS*(-L/6).exp()*L/LN2                  # psi(y)-theta(y) <= cpsi sqrt(y) for y >= x
    t_low = f*(-L/2).exp()*(2*cpsi + S_LOW*(1 + eps_lo))       # psi-theta and zeros with |lam| <= H0
    T = Tsum(G, c, L)
    t_high = f*(1 + eps_hi)*T                                  # zeros with |lam| > H0
    t_pi = (2*PI).log()/(x*L)
    t_G = M0_100**2/2*Gdd_sup(L, M0_100)                       # G(d,L) <= d^2/2 sup G_dd, |d| < M0_100
    tot = t_low + t_high + t_pi + t_G
    if detail: return tot, dict(low=t_low, high=t_high, pi=t_pi, G=t_G, T=T)
    return tot
def kval(E):
    lg = -(up(E).log()/LN2)                                    # -log2 E* (lower end)
    K = int(float(lo(lg).mid()) // 1) - 1                     # floor(-log2 E*) - 1
    assert arb(2)**(-(K + 1)) > up(E)                          # 2^{-K-1} > E*
    return lg, K
# ------------------------------------------------------------------ Gracia Llorente's (20) as printed (reproduction)
GL_TAB = [('9.387','8.290','147.0'),('10.70','7.997','141.9'),('12.11','7.704','136.7'),('13.66','7.411','131.5'),
          ('15.35','7.118','126.3'),('17.21','6.826','121.1'),('19.31','6.533','115.9'),('21.77','6.240','110.7')]
def gl_terminal(L, smax=None):
    """GL eq. (20) with GL's Table 1, J = 2048 grid s_j = 5/8 + j/8192 (s_J = 7/8), 28^{p-3}, tau = 3e-12, 2.134, 2/(xL), G0/L;
       smax: optional cap on Re rho (W = 0 for s > smax), as in Harvest6/k_est.py."""
    J = 2048; tau = arb('3e-12'); HG = arb(3*10**12); hg = HG.log()
    S = [arb(5)/8 + arb(j)/(4*J) for j in range(J + 1)]
    Wt = []
    for j, s in enumerate(S[:J]):
        i = min(7, j//256)                                     # GL's row index min(7, int((s - 5/8)*32)); s - 5/8 = j/8192
        b1, b2, b3 = [arb(v) for v in GL_TAB[i]]
        al = 3*(1 - s)/(2 - s); q = 2 - al; p = (7 - 5*s)/(2 - s)
        I3 = arb(28)**(p - 3)*HG**(-q)*(hg**3/q + 3*hg**2/q**2 + 6*hg/q**3 + 6/q**4)
        v = up(4*(b1*I3 + b2*(hg**2/2 + hg/2 + arb(1)/4)/HG**2 + b3*(hg/2 + arb(1)/4)/HG**2))
        v = emin(up(tau/2), v)
        if smax is not None and s > smax: v = arb(0)
        Wt.append(v if not Wt else emin(Wt[-1], v))
    c = [up(tau) - 2*Wt[0]] + [2*(Wt[j-1] - Wt[j]) for j in range(1, J)] + [2*Wt[J-1]]
    assert all(cj >= 0 for cj in c)
    x = L.exp(); f = 1/L + 1/L**2
    TJ = sum((cj*(x**(s - 1) + x**(-s))/2 for cj, s in zip(c, S)), arb(0))
    G0 = arb(13)/25*arb('2.794e-12')**2
    return (1 + arb('1e-12'))*f*TJ + arb('2.134')*f*(-L/2).exp() + 2*(-L).exp()/L + G0/L
# ------------------------------------------------------------------ main
def rss(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20     # macOS: bytes
if __name__ == '__main__':
    res = {}
    say('Robin2Adic cert.py, prec %d bits%s' % (ctx.prec, ' (QUICK grid)' if QUICK else ''))
    say('B2 = %s, T* = B2 + 2.2e-254, Bebel 0.87338196661320752, used B_UP = %s (exact dyadic)' % (s6(B2, 20), s6(B_UP, 20)))
    say('tau (HSW) = %s; S_low = 2 + gamma - log 4pi = %s' % (s6(TAU), s6(S_LOW)))
    # (1) reproduction of GL's terminal bound at L = 106
    v = gl_terminal(arb(106)); res['GL_L106'] = s6(v, 10)
    say('(1) GL (20) as printed at L = 106: %s  (GL: < 1.227035e-23)' % s6(v, 10))
    assert v < arb('1.227035e-23') and v > arb('1.2270e-23')
    for lab, sm in [('7/8', None), ('B', B_UP)]:
        v = gl_terminal(L1, sm); res['GL_L1_' + lab] = s6(v, 8)
        say('    GL (20) at L1 = 108.57 with Re rho <= %s: %s = 2^-%s' % (lab, s6(v, 8), s6(-(v.log()/LN2), 6)))
    say('    elapsed %.0f s' % (time.time() - T_START))
    # (2) threshold: n > 10^(10^46.79) => x >= e^{108.57}
    threshold_check()
    # (3) density rows
    KLN_P = I4.build_kln_H0('3e12'); KLN_P1 = I4.build_kln_H0('3e12', de_override='1')
    for p in KLN_P + KLN_P1: assert p.k <= 1 and p.k*H0 >= arb(10)**9
    CS_R = []
    for row in I4.CS_ROWS:
        r = I4.CSRow(row); CS_R.append((r.s1, r.s2, arb(row[2]), r.B1, arb(row[3]), arb(row[4])))
    assert I4.T0 == H0
    G = grid(); say('(3) grid: %d points s_0 = 1/2 ... s_J = B_UP' % len(G))
    variants = [('CS_b4_0.54', lambda s: W_cs(s, raw=True)), ('CS_b4_corrected', W_cs), ('KLN_RH3e12', None), ('KLN_RH3e12_delta1', None),
                ('min(CScorr,KLN)', None)]
    out = {}
    for name, Wf in variants:
        t = time.time()
        if name == 'KLN_RH3e12': KLN_P_use = KLN_P
        if name == 'KLN_RH3e12_delta1': KLN_P_use = KLN_P1
        if name.startswith('KLN'):
            globals()['KLN_P'] = KLN_P_use; Wf = W_kln
        if name.startswith('min'):
            globals()['KLN_P'] = KLN_P
            def Wf(s):
                a, b = W_cs(s), W_kln(s)
                return a if b is None else (b if a is None else emin(a, b))
        c, Wt = envelope(G, Wf)
        E, d = majorant(L1, G, c, detail=True)
        lg, K = kval(E)
        # share of the zero sum from s in [0.85, B]: Abel-inverse form T = g_0 tau + sum_j (g_j - g_{j-1}) 2 W~_{j-1}
        x = L1.exp(); g = [(x**(s - 1) + x**(-s))/2 for s in G]
        part = sum((2*Wt[j-1]*(g[j] - g[j-1]) for j in range(1, len(G)) if G[j-1] >= arb('0.85')), arb(0))
        mono = [majorant(arb(Lv), G, c) for Lv in ('108.57', '109', '110', '115', '130', '200', '500', '2000')]
        dec = all(mono[i+1] < mono[i] for i in range(len(mono) - 1))
        out[name] = dict(E=s6(E, 8), log2=s6(lg, 6), K=K, inv=s6(1/up(E), 6), low=s6(d['low'], 4), high=s6(d['high'], 6), G=s6(d['G'], 4),
                         pi=s6(d['pi'], 3), share085=s6(part/d['T'], 3), W_at_085=s6(Wf(arb(int(0.85*2**16))/2**16), 5),
                         W_at_B=s6(Wt[-1], 5), decreasing_sample=dec, sec=round(time.time() - t, 1))
        lE = (1/up(E)).log(); slots = ' '.join('%d:%d' % (p, int(float(lo(lE/arb(p).log()).mid()) // 1)) for p in (2, 3, 5, 7, 11, 13))
        out[name]['slot_min_exponent'] = slots
        say('    slot form: a counterexample has p_j^{a_j+1} > 1/E*, so a_j >= (p:min a_j) %s; a_j >= 2 while p_j <= %s' % (slots, s6((1/up(E)).sqrt(), 4)))
        say('(4) %s: E* <= %s = 2^-%s -> K = %d; 1/E* >= %s; terms: zeros>H0 %s, low %s, G %s, 2pi %s; share s>=0.85: %s; W(0.85) %s, W(B-) %s; decreasing at 8 L: %s; %.0f s' % (
            name, out[name]['E'], out[name]['log2'], K, out[name]['inv'], out[name]['high'], out[name]['low'], out[name]['G'], out[name]['pi'],
            out[name]['share085'], out[name]['W_at_085'], out[name]['W_at_B'], dec, time.time() - t))
        assert dec
    res['variants'] = out; res['rss_MB'] = round(rss(), 1); res['sec'] = round(time.time() - T_START, 1)
    say('peak RSS %.1f MB, %.0f s' % (rss(), time.time() - T_START))
    os.makedirs(os.path.join(HERE, 'logs'), exist_ok=True)
    json.dump(dict(res=res, log=LOG), open(os.path.join(HERE, 'logs', 'cert_result%s.json' % ('_quick' if QUICK else '')), 'w'), indent=1)

