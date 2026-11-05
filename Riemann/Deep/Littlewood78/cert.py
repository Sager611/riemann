# Certificate for the explicit Littlewood-type bound under the zero-free half-plane Re s > B.
# Run: uv run --no-project --with python-flint python3 cert.py cert_output.txt   (single-threaded, ~3 min, < 100 MB)
# Fallback with OpenAI's 7/8:  LW78_B=7/8 LW78_CA=0.105 LW78_CB=0.11 LW78_CB2=0.114 uv run ... cert.py cert_output_7_8.txt
# Every claim printed as PASS is checked with arb ball arithmetic (python-flint); floats only choose parameters.
import math, sys
from flint import arb, ctx
ctx.prec = 256

import os
B_STR = os.environ.get('LW78_B', '0.87338196661320752')   # default: Lean-checked half-plane (Bebel, 31 Oct 2026)
B = arb(7)/8 if B_STR == '7/8' else arb(B_STR)              # ball contains the exact value
b = 1 - B
CA = arb(os.environ.get('LW78_CA', '0.11')); CB = arb(os.environ.get('LW78_CB', '0.115'))
CB2 = arb(os.environ.get('LW78_CB2', '0.116')); CUP = int(os.environ.get('LW78_CUP', '15'))
EG = arb.const_euler(); PI = arb.pi(); LOG2 = arb(2).log()
C0 = (PI**2/6).log() - EG             # log((pi^2/6) e^{-gamma})
bf = 1 - (0.875 if B_STR == '7/8' else float(B_STR)); EGf = 0.5772156649015329; C0f = math.log(math.pi**2/6) - EGf

def Q(num, den=1000):                 # exact rational num/den as an arb ball
    return arb(num)/den

E_UP = arb('-1.33258227573')         # upper bound for Dusart's E = -1.3325822757332...
Ef_UP = -1.33258227573

def kappas(L, h, L1, h1, exp=lambda t: t.exp()):
    """kappa0 (for G, weight (L,h)) and kappa1 (for V at s=1, weight (L1,h1)), with the pairing rho <-> 1-conj(rho):
    kappa = g(B)*B/b + g(1/2), g increasing in beta."""
    def g0(beta):   # |G(rho)| <= g0(beta)/|1-rho|^2
        u = 1 - beta
        return (exp(-(1+h)*L*u)/(1+h) + exp(-L*u))/(h*L*L)
    def g1(beta):   # |V(rho-1)| <= g1(beta)/|1-rho|^2
        u = 1 - beta
        return (exp(-(1+h1)*L1*u) + exp(-L1*u))/(h1*L1)
    return g0(B)*B/b + g0(arb(1)/2), g1(B)*B/b + g1(arb(1)/2)

def tau_bar(Lx):
    """Lemma 9(d): tau(e^Lx) <= tau_bar(Lx) = e^{-Lx/2} + ((Lx/log 2)(e^{Lx/3}+3)/2 + 8) e^{-Lx}.
    tau itself jumps up at every prime; tau_bar is decreasing for Lx >= 3/2 and Lx*tau_bar(Lx) for Lx >= 3 (PROOF.md)."""
    xx = Lx.exp(); K = Lx/LOG2
    return (-Lx/2).exp() + (K*((Lx/3).exp()+3)/2 + 8)/xx

def F(ell, P, parts=False):
    """Master inequality: log|L(1,chi)| >= F(ell, P) for every primitive chi mod q, ell = log q, q >= 3,
    P = (L, h, L1, h1): L >= 10, L1 >= 8, 0 < h, h1 <= 2, kappa1 < 1.  Returns an arb ball (use its lower end)."""
    L, h, L1, h1 = P
    assert (L >= 10) and (L1 >= 8) and (h > 0) and (h <= 2) and (h1 > 0) and (h1 <= 2)
    x = L.exp(); x1 = L1.exp()
    k0, k1 = kappas(L, h, L1, h1)
    if not (k1 < 1): return None
    kp = k0/(1-k1)
    T1 = (-L1).exp()/(1 - (-2*L1).exp())           # trivial zeros, formula at s=1 with (L1,h1)
    T0 = (-L).exp()/(L*(1 - (-2*L).exp()))         # trivial zeros, G with (L,h)
    Eh = ((1+h)*(1+h).log() - h)/h
    epsD = 1/(10*L*L) + arb(4)/(15*L*L*L)
    main_err = epsD + 1/(2*(x-1)) + tau_bar(L)/2 + 1/(2*(x-1)) + T0      # Mertens, p>=x squares, p<x tails, Lemma 9(c), trivial
    W = (L1*(1+h1/2) + E_UP + arb('0.2')/L1 + arb('0.2')/(L1*L1) + L1*tau_bar(L1) + (1+h1)*L1/(2*(x1-1))
         + (ell - PI.log())/2 - EG/2 + T1)
    val = C0 - L.log() - Eh - main_err - kp*W
    if parts:
        return val, dict(logL=L.log(), Eh=Eh, kpW=kp*W, small=main_err, k0=k0, k1=k1, kp=kp)
    return val

def Ff(ell, P):                       # float copy of F, used only to choose parameters
    L, h, L1, h1 = P
    if L < 10 or L1 < 8 or h <= 0 or h > 2 or h1 <= 0 or h1 > 2: return -1e9
    Bf = 1 - bf
    g0 = lambda be: (math.exp(-(1+h)*L*(1-be))/(1+h) + math.exp(-L*(1-be)))/(h*L*L)
    g1 = lambda be: (math.exp(-(1+h1)*L1*(1-be)) + math.exp(-L1*(1-be)))/(h1*L1)
    k0 = g0(Bf)*Bf/bf + g0(0.5); k1 = g1(Bf)*Bf/bf + g1(0.5)
    if k1 >= 1: return -1e9
    kp = k0/(1-k1)
    T1 = math.exp(-L1)/(1-math.exp(-2*L1)); T0 = math.exp(-L)/(L*(1-math.exp(-2*L)))
    tau = lambda Lx: math.exp(-Lx/2) + (Lx/math.log(2)*(math.exp(min(Lx/3,700))+3)/2 + 8)*math.exp(-Lx)
    Eh = ((1+h)*math.log1p(h) - h)/h
    epsD = 1/(10*L*L) + 4/(15*L**3)
    x = math.exp(min(L, 700)); x1 = math.exp(min(L1, 700))
    main_err = epsD + 1/(x-1) + tau(L)/2 + T0
    W = L1*(1+h1/2) + Ef_UP + 0.2/L1 + 0.2/L1**2 + L1*tau(L1) + (1+h1)*L1/(2*(x1-1)) + (ell - math.log(math.pi))/2 - EGf/2 + T1
    return C0f - math.log(L) - Eh - main_err - kp*W

def opt(ell, P0=None, obj=None):
    """maximise obj (default Ff) over P=(L,h,L1,h1); returns exact rationals as integer numerators over 1000."""
    Ff_ = obj if obj else Ff
    if P0 is None:
        best = max(((Ff_(ell, (L, h, L, h)), L, h) for L in [10+2*i for i in range(400)] for h in [0.04*j for j in range(1, 51)]))
        P = [best[1], best[2], best[1], best[2]]
    else:
        P = list(P0)
    steps = [2.0, 0.02, 2.0, 0.02]
    cur = Ff_(ell, P)
    while max(steps[0], steps[2]) > 1e-3 or max(steps[1], steps[3]) > 1e-4:
        moved = False
        for i in range(4):
            for sgn in (1, -1):
                Q2 = list(P); Q2[i] += sgn*steps[i]
                v = Ff_(ell, Q2)
                if v > cur: cur, P, moved = v, Q2, True
        if not moved: steps = [t/2 for t in steps]
    r = [max(10000, round(P[0]*1000)), min(2000, max(1, round(P[1]*1000))), max(8000, round(P[2]*1000)), min(2000, max(1, round(P[3]*1000)))]
    return tuple(r)

def Pq(r):  # integer numerators -> arb tuple
    return tuple(Q(t) for t in r)

def Fup(ell, P):
    """Upper-bound master inequality (PROOF.md, Prop. 11): log|L(1,chi)| <= Fup(ell, P), same parameter ranges."""
    L, h, L1, h1 = P
    assert (L >= 10) and (L1 >= 8) and (h > 0) and (h <= 2) and (h1 > 0) and (h1 <= 2)
    x = L.exp(); x1 = L1.exp()
    k0, k1 = kappas(L, h, L1, h1)
    if not (k1 < 1): return None
    kp = k0/(1-k1)
    T1 = (-L1).exp()/(1 - (-2*L1).exp())
    Eh = ((1+h)*(1+h).log() - h)/h
    epsD = 1/(10*L*L) + arb(4)/(15*L*L*L)
    Wup = L1*(1+h1/2) - EG + arb('0.2')/L1 + arb('0.2')/(L1*L1) + (ell - PI.log())/2 - EG/2 + T1
    return EG + L.log() + Eh + epsD + 1/(2*(x-1)) + kp*Wup

def Fupf(ell, P):                     # float copy, objective = -Fup (to maximise)
    L, h, L1, h1 = P
    if L < 10 or L1 < 8 or h <= 0 or h > 2 or h1 <= 0 or h1 > 2: return -1e9
    Bf = 1 - bf
    g0 = lambda be: (math.exp(-(1+h)*L*(1-be))/(1+h) + math.exp(-L*(1-be)))/(h*L*L)
    g1 = lambda be: (math.exp(-(1+h1)*L1*(1-be)) + math.exp(-L1*(1-be)))/(h1*L1)
    k0 = g0(Bf)*Bf/bf + g0(0.5); k1 = g1(Bf)*Bf/bf + g1(0.5)
    if k1 >= 1: return -1e9
    kp = k0/(1-k1)
    T1 = math.exp(-L1)/(1-math.exp(-2*L1))
    Eh = ((1+h)*math.log1p(h) - h)/h
    x = math.exp(min(L, 700))
    Wup = L1*(1+h1/2) - EGf + 0.2/L1 + 0.2/L1**2 + (ell - math.log(math.pi))/2 - EGf/2 + T1
    return -(EGf + math.log(L) + Eh + 1/(10*L*L) + 4/(15*L**3) + 1/(2*(x-1)) + kp*Wup)

out = []
def say(*a):
    s = ' '.join(str(t) for t in a); print(s, flush=True); out.append(s)

def lower(a): return float(a.lower()) if hasattr(a, 'lower') else float(a)

say("B =", B, " b = 1-B =", b)
say("limit constant (pi^2/6) e^-gamma (1-B) =", (C0.exp()*b))

# ---------- (A) uniform claim: |L(1,chi)| >= cA / max(LamA, loglog q) for all q >= 3 ----------
def check_claim(c, Lam_lo, Lam_hi, step_num, den=1000, label="", quiet=False):
    """cells [Lam_i, Lam_{i+1}] with Lam_i = Lam_lo + i*step_num/den; checks F(exp(Lam_{i+1})) > log(c/Lam_i)."""
    lo = round(Lam_lo*den); hi = round(Lam_hi*den)
    i = lo; worst = (1e9, None); Lh = None; ncell = 0
    while i < hi:
        j = min(i + step_num, hi)
        Lam_i = Q(i, den); Lam_j = Q(j, den)
        ell_j = Lam_j.exp()
        r = opt(float(ell_j.mid()), Lh)
        Lh = tuple(t/1000 for t in r)
        val = F(ell_j, Pq(r))
        target = (arb(c)/Lam_i).log() if isinstance(c, (int,)) else (c/Lam_i).log()
        if val is None or not (val > target):
            if not quiet: say(f"  FAIL {label} cell [{i/den},{j/den}] P={Lh}")
            return False, worst, ncell
        m = lower(val - target)
        if m < worst[0]: worst = (m, (i/den, j/den) + Lh)
        ncell += 1; i = j
    return True, worst, ncell

def tail_terms(Lam, h):
    """Lemma T with L = L1 = Lam/b, h1 = h: returns (kappa_1, Err(Lam), Err_up(Lam)); Err uses the decreasing majorant tau_bar."""
    L = Lam/b; ell = Lam.exp(); x = L.exp()
    k0, k1 = kappas(L, h, L, h)
    kp = k0/(1-k1)
    T0 = (-L).exp()/(L*(1 - (-2*L).exp()))
    epsD = 1/(10*L*L) + arb(4)/(15*L*L*L)
    Err = epsD + 1/(x-1) + tau_bar(L)/2 + T0 + kp*(L*(1+h/2) + ell/2)
    Err_up = epsD + 1/(2*(x-1)) + kp*(L*(1+h/2) + ell/2)
    return k1, Err, Err_up

def tail(Lam_N, h, show=True):
    """Lemma T (PROOF.md): if Lam_N/b >= 10, 0 < h <= 2 and kappa_1(Lam_N) < 1, then for all Lam >= Lam_N
    log|L(1,chi)| >= C0 - log(Lam/b) - E_h - Err(Lam_N)  (kappa_1, Err are decreasing in Lam).  Returns the constant c."""
    assert Lam_N/b >= 10 and h > 0 and h <= 2
    k1, Err, _ = tail_terms(Lam_N, h)
    ok = k1 < 1
    Eh = ((1+h)*(1+h).log() - h)/h
    c = (C0 - Eh - Err).exp()*b if ok else arb(0)
    if show: say(f"  tail Lam>={float(Lam_N.mid())}, h={float(h.mid())}: kappa_1(Lam_N) <= {k1.upper().str(4, radius=False)} < 1: {'PASS' if ok else 'FAIL'};"
                 f" E_h={float(Eh.mid()):.6f} Err={float(Err.upper()):.3e}  => |L(1,chi)| >= {float(c.lower()):.6f}/loglog q")
    return c

say("\n(0) side conditions of the tail lemma (PROOF.md, Lemma T (ii)): brackets at L = 8, h = 2 (non-increasing in L >= 8, increasing in h)")
L8 = arb(8); T1_8 = (-L8).exp()/(1 - (-2*L8).exp())
br_lo = E_UP + arb('0.2')/L8 + arb('0.2')/(L8*L8) + L8*tau_bar(L8) + 3*L8/(2*(L8.exp()-1)) - PI.log()/2 - EG/2 + T1_8
br_up = -EG + arb('0.2')/L8 + arb('0.2')/(L8*L8) - PI.log()/2 - EG/2 + T1_8
say(f"  lower bound: E + 0.2/L + 0.2/L^2 + L tau_bar + (1+h)L/(2(e^L-1)) - log(pi)/2 - gamma/2 + T_1 <= {float(br_lo.upper()):.6f} < 0: {'PASS' if br_lo < 0 else 'FAIL'}")
say(f"  upper bound: -gamma + 0.2/L + 0.2/L^2 - log(pi)/2 - gamma/2 + T_1 <= {float(br_up.upper()):.6f} < 0: {'PASS' if br_up < 0 else 'FAIL'}")
say("  monotonicity spot check (not needed for the proof): kappa_1 and Err at Lam = 40*{1, 1.25, 1.5, 2, 3, 5, 10}, h = 0.015")
prev = None; mono = True
for f in ('1', '1.25', '1.5', '2', '3', '5', '10'):
    k1s, Errs, Errus = tail_terms(40*arb(f), Q(15))
    if prev is not None and not (Errs.upper() < prev[1].lower() and k1s.upper() < prev[0].lower()): mono = False
    prev = (k1s, Errs)
    say(f"    Lam = {40*float(f):6.1f}: kappa_1 = {float(k1s.mid()):.3e}, Err = {float(Errs.mid()):.6e}, Err_up = {float(Errus.mid()):.6e}")
say(f"    strictly decreasing on the sample: {mono}")

say(f"\n(A) |L(1,chi)| >= {float(CA.mid())}/max(3, loglog q), every primitive chi mod q >= 3")
ell3 = arb(3).exp()
r3 = opt(float(ell3.mid())); v = F(ell3, Pq(r3))
okA1 = v > (CA/3).log()
say(f"  small q (loglog q <= 3): F(e^3; P={tuple(t/1000 for t in r3)}) = {float(v.lower()):.6f} > log(cA/3) = {float((CA/3).log().mid()):.6f}: {'PASS' if okA1 else 'FAIL'}  (|L(1,chi)| >= {math.exp(float(v.lower())):.5f})")
okA2, worstA, nA = check_claim(CA, 3.0, 40.0, 5, label="A")
say(f"  grid loglog q in [3, 40], {nA} cells: {'PASS' if okA2 else 'FAIL'}; tightest cell {worstA}")
cT = tail(arb(40), Q(15))
okA3 = cT > CA
say(f"  tail loglog q >= 40: {'PASS' if okA3 else 'FAIL'}")

say(f"\n(B) |L(1,chi)| >= {float(CB.mid())}/loglog q for loglog q >= LamB")
# find LamB: smallest multiple of 0.01 where the optimised bound at the cell passes
for LamB_num in range(300, 500, 1):
    LamB = LamB_num/100
    ok, worst, n = check_claim(CB, LamB, LamB + 0.02, 5, quiet=True)   # search for the start
    if ok: break
okB, worstB, nB = check_claim(CB, LamB, 40.0, 5, label="B")
say(f"  LamB = {LamB}: grid [LamB, 40], {nB} cells: {'PASS' if okB else 'FAIL'}; tightest cell {worstB}")
okB3 = cT > CB
say(f"  tail loglog q >= 40 (constant above): {'PASS' if okB3 else 'FAIL'}")

say(f"\n(B2) |L(1,chi)| >= {float(CB2.mid())}/loglog q for loglog q >= LamB2")
for LamB2_num in range(400, 3000, 5):
    LamB2 = LamB2_num/100
    ok, worst, n = check_claim(CB2, LamB2, LamB2 + 0.02, 5, quiet=True)  # search for the start
    if ok: break
okB2, worstB2, nB2 = check_claim(CB2, LamB2, 100.0, 10, label="B2")
cT100 = tail(arb(100), Q(6))
say(f"  LamB2 = {LamB2}: grid [LamB2, 100], {nB2} cells: {'PASS' if okB2 else 'FAIL'}; tightest cell {worstB2}; tail >= 100: {'PASS' if cT100 > CB2 else 'FAIL'}")

say(f"\n(C) tail lemma constants (approach (pi^2/6)e^-gamma(1-B) = {float((C0.exp()*b).mid()):.7f})")
for LN, hh in ((40, 15), (100, 6), (1000, 1)):
    tail(arb(LN), Q(hh))
tail(arb(10000), Q(1, 10000))

say("\n(D) pointwise table: certified lower bound for |L(1,chi)| at q = exp(exp(Lam)) (optimised L, h)")
for Lam in [1, 2, 2.5, 3, 3.2, 3.5, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40]:
    LamA = Q(round(Lam*1000)); ell = LamA.exp()
    r = opt(float(ell.mid())); v, P = F(ell, Pq(r), parts=True)
    lb = math.exp(float(v.lower()))
    say(f"  loglog q={Lam:5}: L={r[0]/1000:8.3f} h={r[1]/1000:5.3f} L1={r[2]/1000:8.3f} h1={r[3]/1000:5.3f} F={float(v.lower()):.9f} |L(1,chi)| >= {lb:.6f}  (x loglog q = {lb*Lam:.5f})  logL={float(P['logL'].mid()):.3f} E_h={float(P['Eh'].mid()):.4f} kappa'W={float(P['kpW'].mid()):.4f}")

say("\n(E) class numbers: h(-d) = sqrt(d) L(1,chi_-d)/pi for fundamental -d < -4")
def class_threshold(H, ell_guess):
    # find smallest ell0 (multiple of 0.01) such that cells [ell_i, ell_i+0.01] from ell0 to ell_end certify h > H,
    # where beyond ell_end claim (A) gives h(-d) >= cA sqrt(d)/(pi loglog d) > H (increasing in d once loglog d >= 3).
    ell_end = None
    for e in range(int(ell_guess*100), 100000):
        ee = Q(e, 100)
        if math.log(e/100) >= 3 and (ee/2).exp()*CA/(PI*(ee.log())) > H:   # claim (A), loglog d >= 3
            ell_end = e; break
    def cells_ok(s):
        Lh = None; i = s
        while i < ell_end:
            lo = Q(i, 100); hi = Q(i+1, 100)
            r = opt(float(hi.mid()), Lh); Lh = tuple(t/1000 for t in r)
            v = F(hi, Pq(r))
            if v is None or not ((lo/2).exp()*v.exp()/PI > H): return False
            i += 1
        return True
    s = int(ell_guess*100)
    if cells_ok(s):
        while s > 0 and cells_ok(s-1): s -= 1
    else:
        while not cells_ok(s): s += 1
    return s/100, ell_end/100
for H in (100, 1000, 2000, 10**4, 10**5):
    guess = 2*math.log(H*math.pi/0.035)
    e0, e1 = class_threshold(H, guess)
    v0 = math.exp(e0); k0 = math.floor(math.log10(v0)); up = math.ceil(v0/10**k0*1000)/1000   # rounded up
    say(f"  h(-d) > {H} for every d >= e^{e0} (so for d >= {up}e{k0})  (cells of 0.01 in log d up to {e1}, claim (A) beyond)")

say("  (2610.21511 remark, same half-plane, via its Prop. L1: h(-d) > 2000 for log d >= 27.2, h(-d) > 10^5 for log d >= 35.4)")
say("\n(F) comparison with h(-d) >= 0.0001828446364 sqrt(d)/log d (2610.21511 Prop. L1): ratio of our bound (A)")
for d_log in (2, 5, 10, 20, 35, 100):
    ours = float(CA.mid())/(math.pi*max(3, math.log(d_log)))
    theirs = 0.0001828446364/d_log
    say(f"  log d = {d_log}: ours/sqrt(d) >= {ours:.5f}, theirs = {theirs:.3e}, ratio {ours/theirs:.0f}")

say("\n(G) side result, upper bound: |L(1,chi)| <= Cup * loglog q")
def check_up(C, Lam_lo, Lam_hi, step_num, den=1000):
    lo = round(Lam_lo*den); hi = round(Lam_hi*den); i = lo; worst = (1e9, None); Lh = None; n = 0
    while i < hi:
        j = min(i + step_num, hi)
        Lam_i = Q(i, den); ell_j = Q(j, den).exp()
        r = opt(float(ell_j.mid()), Lh, Fupf); Lh = tuple(t/1000 for t in r)
        v = Fup(ell_j, Pq(r))
        target = (C*Lam_i).log()
        if v is None or not (v < target): return False, worst, n
        m = lower(target - v)
        if m < worst[0]: worst = (m, (i/den, j/den) + Lh)
        n += 1; i = j
    return True, worst, n
def tail_up(Lam_N, h):
    """Upper-bound tail (PROOF.md, after Prop. 11): same hypotheses as Lemma T."""
    assert Lam_N/b >= 10 and h > 0 and h <= 2
    k1, _, Err = tail_terms(Lam_N, h)
    ok = k1 < 1
    Eh = ((1+h)*(1+h).log() - h)/h
    C = EG.exp()/b*(Eh + Err).exp() if ok else arb(10)**9
    say(f"  tail Lam>={float(Lam_N.mid())}, h={float(h.mid())}: kappa_1(Lam_N) <= {k1.upper().str(4, radius=False)} < 1: {'PASS' if ok else 'FAIL'}; |L(1,chi)| <= {float(C.upper()):.4f} loglog q")
    return C
say("  limit constant e^gamma/(1-B) =", EG.exp()/b)
for Cup, Lam_lo in ((CUP + 1, 3.0), (CUP, None)):
    if Lam_lo is None:
        for k in range(300, 4000, 5):
            ok, _, _ = check_up(arb(Cup), k/100, k/100 + 0.02, 5)
            if ok: Lam_lo = k/100; break
    okU, worstU, nU = check_up(arb(Cup), Lam_lo, 40.0, 5)
    CT = tail_up(arb(40), Q(15))
    say(f"  |L(1,chi)| <= {Cup} loglog q for loglog q >= {Lam_lo}: grid [{Lam_lo}, 40] {nU} cells {'PASS' if okU else 'FAIL'} (tightest {worstU}); tail >= 40: {'PASS' if CT < Cup else 'FAIL'}")
for LN, hh in ((100, 6), (1000, 1)):
    tail_up(arb(LN), Q(hh))
for Lam in (3, 4, 5, 10, 20):
    ell = Q(Lam*1000).exp(); r = opt(float(ell.mid()), None, Fupf); v = Fup(ell, Pq(r))
    say(f"  loglog q={Lam}: |L(1,chi)| <= {math.exp(float(v.upper())):.3f}  (= {math.exp(float(v.upper()))/Lam:.3f} loglog q; (1/2)log q = {0.5*math.exp(Lam):.1f})")

open(sys.argv[1] if len(sys.argv) > 1 else 'cert_output.txt', 'w').write('\n'.join(out) + '\n')
