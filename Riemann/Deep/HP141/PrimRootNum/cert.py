#!/usr/bin/env python3
"""
cert.py -- Arb certificate: least prime primitive roots from the Lean-checked zero-free half-plane.
Run (one command, one thread):
    cd Riemann/Deep/PrimRootNum && uv run --no-project --with python-flint --with mpmath --with numpy python3 cert.py all
Sections:  core | sieve | omega | tail | small | all     (DERIVATION.md section numbers in comments)
Inputs: those of Riemann/Deep/ExplProg/cert.py (imported by path: half-plane B, kernel 3.1, Lemmas EF, P, Z, Zzeta, K1),
        sieve.py (exact sieve main terms / costs), small.py (exact g*(p), p <= 10^7).
Every inequality used as a proof step is a certified sign of an arb ball (pos()); floats only choose parameters.
"""
import sys, os, math, json, time, importlib.util
from fractions import Fraction as Fr
sys.dont_write_bytecode = True        # do not leave __pycache__ in ExplProg/ or here

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location("explprog_cert", os.path.join(HERE, "..", "ExplProg", "cert.py"))
EP = importlib.util.module_from_spec(_spec); _spec.loader.exec_module(EP)
from flint import arb, ctx
ctx.prec = 200
sys.setrecursionlimit(100000)
import sieve as SV

A, pos, amax = EP.A, EP.pos, EP.amax
B, b, INV_b = EP.B, EP.b, EP.INV_b
KER = EP.KER                      # kernel 3.1 of ExplProg: ell = 7, c = 21/10
LOG9 = arb(9).log()
ONE = arb(2) ** 64
E1 = arb(1).exp()

def ZoverL(co, nn, Lm):
    """sup_{L >= Lm} Z(max(L, log 9))/L  (Lemma Z uses L = max(log conductor, log 9))"""
    if Lm >= LOG9:
        return EP.Z_over_L_majorant(co, nn, Lm)
    return amax(EP.Zof(co, nn, LOG9) / Lm, EP.Z_over_L_majorant(co, nn, LOG9))

def sieve_margin(kappa, Lm, Tlo, C):
    """DERIVATION sec 4, inequality (4.3) divided by x^B L, at x = (kappa L)^(1/b), with every L-dependent term replaced
    by its value at L = Lm (non-increasing in L, sec 4.4).  Positive => g*(p) <= (kappa log p)^(1/b) for every prime p with
    log p >= Lm and omega(p-1) = w, where (Tlo, C) are the main term (lower bound) and cost of the shifted sieve for w."""
    kappa = A(kappa)
    X = (kappa * Lm).log() / b
    Y = X - KER.ell
    if not (Y > 1):
        return None
    (n0, k1, k2), (m0, k1L, k2L) = KER.ks(Y)
    co = EP.zsum_coeffs(n0, k1, k2); coL = EP.zsum_coeffs(m0, k1L, k2L)
    ZoL = ZoverL(co, n0, Lm); ZLoL = ZoverL(coL, m0, Lm)
    xB = (B * X).exp(); xhB = ((A(Fr(1, 2)) - B) * X).exp()
    z1 = EP.SUM_INV_G2 * k2 / 2; z1L = EP.SUM_INV_G2 * k2L / 2
    triv = KER.W0 / ((1 - (-2 * Y).exp()) * xB * Lm)
    Xo = KER.W1 * kappa - (z1 + xhB * z1L) / Lm - triv - KER.winf * (7 / Lm + 1) / xB      # X_0/(x^B L) lower bound
    Ro = ZoL / 2 + xhB * ZLoL / 2 + triv                                                     # R/(x^B L) upper bound
    return Tlo * Xo - (C - 1) * Ro

def best_kappa(Lm, Tlo, C, digits=3):
    """smallest kappa = m 10^e (m with 'digits' digits) with certified positive margin (float bisection, Arb check)"""
    lo, hi = 1e-6, 1e60
    def ok(k):
        m = sieve_margin(Fr(k).limit_denominator(10 ** 12) if k < 1e12 else Fr(int(k)), Lm, Tlo, C)
        return m is not None and m > 0
    while not ok(hi):
        hi *= 10
    for _ in range(200):
        mid = math.sqrt(lo * hi)
        if ok(mid):
            hi = mid
        else:
            lo = mid
        if hi / lo < 1 + 1e-6:
            break
    e = math.floor(math.log10(hi)) - (digits - 1)
    m = math.ceil(hi / 10 ** e - 1e-12)
    k = Fr(m) * Fr(10) ** e
    mg = sieve_margin(k, Lm, Tlo, C)
    pos(mg, "sieve margin at kappa = %s" % k)
    return k, mg


# ---------------------------------------------------------------- sec 3: sieve tables (exact, recomputed here)
def load_choice(path=os.path.join(HERE, "logs", "sieve_best.json")):
    d = json.load(open(path))
    return d["Wmax"], [row[0] for row in d["best"]]

def sieve_tables(Wmax, ychoice):
    """recompute T_w (exact lower-bound numerator, unit 2^-64) and C_w for the chosen level y_w of every w <= Wmax"""
    P = SV.primes_first(Wmax)
    need = {}
    for w, y in enumerate(ychoice, start=1):
        need.setdefault(y, []).append(w)
    T = [None] * Wmax; C = [None] * Wmax
    t0 = time.time()
    for y in sorted(need):
        Tq, Cq = SV.enumerate_weights(y, P)
        Tc, Cc = SV.cumulative(Tq, Cq)
        for w in need[y]:
            T[w - 1] = Tc[w - 1]; C[w - 1] = Cc[w - 1]
    print("[sieve] %d levels recomputed for w <= %d in %.0f s" % (len(need), Wmax, time.time() - t0), flush=True)
    return P, T, C

def L_of_w(P, w):
    """log(P_w + 1) <= log p for every prime p with omega(p-1) = w (p - 1 >= p_1 ... p_w); at least log 3"""
    prod = 1
    for q in P[:w]:
        prod *= q
    return arb(prod + 1).log() if prod + 1 >= 3 else arb(3).log()

def kappa_fast(Lm, Tlo, C, digits=3):
    """fixed-point estimate of the least kappa, rounded up to 'digits' significant digits, then certified"""
    k = float(((C - 1) * 0.5 * 2.3 / (Tlo * 2.42)).upper()) + 1.0
    for _ in range(6):
        kA = Fr(k).limit_denominator(10 ** 9)
        X = (A(kA) * Lm).log() / b; Y = X - KER.ell
        if not (Y > 1):
            k *= 4; continue
        (n0, k1, k2), (m0, k1L, k2L) = KER.ks(Y)
        co = EP.zsum_coeffs(n0, k1, k2); coL = EP.zsum_coeffs(m0, k1L, k2L)
        ZoL = ZoverL(co, n0, Lm); ZLoL = ZoverL(coL, m0, Lm)
        xB = (B * X).exp(); xhB = ((A(Fr(1, 2)) - B) * X).exp()
        z1 = EP.SUM_INV_G2 * k2 / 2; z1L = EP.SUM_INV_G2 * k2L / 2
        triv = KER.W0 / ((1 - (-2 * Y).exp()) * xB * Lm)
        sub = (z1 + xhB * z1L) / Lm + triv + KER.winf * (7 / Lm + 1) / xB
        Ro = ZoL / 2 + xhB * ZLoL / 2 + triv
        knew = float((((C - 1) * Ro + Tlo * sub) / (Tlo * KER.W1)).upper())
        k = max(knew, 1e-9)
    e = math.floor(math.log10(k)) - (digits - 1)
    m = math.ceil(k / 10 ** e * (1 + 1e-9))
    for _ in range(50):
        kk = Fr(m) * Fr(10) ** e
        mg = sieve_margin(kk, Lm, Tlo, C)
        if mg is not None and mg > 0:
            return kk, mg
        m += 1
    raise SystemExit("kappa_fast failed")

def omega_section(Wmax=None, Lsmall=None, out="logs/kappa_table.txt", tables=None):
    """sec 4.5: for every w <= Wmax: kappa_w with g*(p) <= (kappa_w log p)^{1/b} for all p with omega(p-1) = w,
    log p >= max(L_w, Lsmall).  Returns rows (w, p_w, Lm, kappa, c_w) with c_w^3 >= kappa_w/Lm^2 (exponent 3/b form)."""
    Wm, ych = load_choice()
    if Wmax is None:
        Wmax = Wm
    if tables is None:
        tables = sieve_tables(Wmax, ych[:Wmax])
    P, T, C = tables
    rows = []; t0 = time.time()
    prod = 1
    f = open(os.path.join(HERE, out), "w")
    f.write("# w p_w y_w T_w(lower) C_w Lm kappa_w c_w   [g*(p) <= (kappa_w log p)^(1/(1-B)) <= (c_w log p)^(3/(1-B)) if omega(p-1)=w, log p >= Lm]\n")
    for w in range(1, Wmax + 1):
        prod *= P[w - 1]
        Lw = arb(prod + 1).log()
        Lm = Lw if Lsmall is None else amax(Lw, A(Lsmall))
        Tlo = arb(T[w - 1]) / ONE
        pos(Tlo, "T_w > 0 (w = %d)" % w)
        kap, mg = kappa_fast(Lm, Tlo, C[w - 1])
        cw3 = A(kap) / (Lm * Lm)
        cw = float(cw3.upper()) ** (1.0 / 3.0) * (1 + 1e-12)
        pos(A(Fr(cw).limit_denominator(10 ** 12)) ** 3 - cw3, "c_w^3 >= kappa_w/Lm^2")
        rows.append((w, P[w - 1], Lm, kap, cw))
        f.write("%d %d %d %.8g %d %.6f %s %.6f\n" % (w, P[w - 1], ych[w - 1], float(Tlo.lower()), C[w - 1], float(Lm.lower()), str(float(kap)), cw))
        if w in (1, 2, 3, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000):
            print("[omega] w=%5d p_w=%6d L>=%.4g kappa=%.4g c_w=%.4f  (%.0f s)" % (w, P[w - 1], float(Lm.lower()), float(kap), cw, time.time() - t0), flush=True)
    f.close()
    return rows


# ---------------------------------------------------------------- sec 5: omega(p-1) > W1 (Buchstab + Selberg; tail.py)
def tail_section(s=Fr(23, 10), delta=Fr(1, 10), W1=10000, w0=50, Q2=100000, uA=60, out="logs/tail_output.txt"):
    import tail as TL
    t0 = time.time()
    T = TL.Tail(s, w0=w0, Q2=Q2)
    IJK = T.integrals(n=4000)
    I, J, K3, P1 = IJK
    P = SV.primes_first(W1 + 1)
    Z1 = P[-1]
    u1 = arb(Z1).log()
    prod = 1
    for q in P:
        prod *= q
    L1 = arb(prod + 1).log()                        # log p >= log(P_{W1+1} + 1) when omega(p-1) >= W1 + 1
    lines = []
    def say(m):
        print(m, flush=True); lines.append(m)
    say("[tail] s = %s, delta = %s, w0 = %d, Q2 = %d, W1 = %d: Z1 = p_{W1+1} = %d, u1 = log Z1 = %s, L1 = %s" % (s, delta, w0, Q2, W1, Z1, u1.str(8), L1.str(8)))
    say("[tail] integrals: I = %s, J = %s, K3 = %s, Psi(1) = %s, acap = %s, c_sigma <= %s" % (I.str(6), J.str(6), K3.str(6), P1.str(6), T.acap, T.c_sig.str(8)))
    pos(1 - (1 + 1 / (2 * T.LQ ** 2)) / (1 - arb(1) / Q2) * I, "(1+eps_V) I < 1")
    # cells in u = log z
    EMG = (-arb.const_euler()).exp()
    tau = None; worst_u = None
    u = u1; n = 0
    while u < A(uA):
        h = A(Fr(1, 10)) if u < 20 else (A(Fr(1, 4)) if u < 40 else A(Fr(1, 2)))
        uh = u + h
        pa = T.partA(u, uh)
        Vlow = EMG * (1 - 1 / (2 * u * u)) / uh                  # R3: V+(z) > e^-g (1 - 1/(2 log^2 z))/log z, z >= 286
        Arel = pa / Vlow
        rr = T.B_rel_restricted(uh, IJK)
        if rr is None:
            raise SystemExit("restricted grid failed")
        Brel = T.B_rel(u, IJK, rr)
        tc = arb((1 - Arel - Brel).lower())            # exact lower endpoint: comparisons below are a total order
        if tau is None or tc < tau:
            tau = tc; worst_u = u
        u = uh; n += 1
    # u >= uA: Part A by Rankin decay with one alpha, Part B by the full-range bound at uA
    al = Fr(12, 100)
    AR = arb(0)
    for q in T.PA:
        lq = arb(q).log()
        lx = (T.s * A(uA) - lq) / 2
        r = (-A(al) * lx).exp() * T.Pia[al][q]
        pos(1 - r, "Rankin r < 1 at uA (q = %d)" % q)
        AR += T.V[q] / q * r / (1 - r)
    kap = A(al) * T.s / 2
    pos(kap * A(uA) - 1, "kappa uA >= 1 (u e^{-kappa u} decreasing)")
    Arel_inf = AR / EMG * A(uA) / (1 - 1 / (2 * A(uA) ** 2))        # e^gamma u / (1 - 1/(2u^2)) >= 1/V+(z)
    Brel_inf = T.B_rel(A(uA), IJK)
    tau_inf = arb((1 - Arel_inf - Brel_inf).lower())
    say("[tail] %d cells u in [%s, %d]: min tau = %s (at u = %s); u >= %d: A_rel <= %s, B_rel <= %s, tau >= %s" %
        (n, u1.str(6), uA, tau.str(6), worst_u.str(5), uA, Arel_inf.str(4), Brel_inf.str(6), tau_inf.str(6)))
    tauT = tau if tau < tau_inf else tau_inf
    tauT = arb(tauT.lower())
    pos(tauT, "tau > 0: T >= tau V+(z) for every z >= Z1")
    KA = T.K_A(u1); KB = T.K_B()
    say("[tail] cost: C - 1 <= 2^%d + c_sigma^2 (K_A + K_B) z^s, K_A <= %s, K_B <= %s" % (T.pi_w0, KA.str(6), KB.str(6)))
    beta = 1 + A(s) + A(delta)
    kz = 1 - 1 / u1
    ell1 = (L1 / kz).log()
    # m_delta <= L^delta / log(L/kz) for every L >= L1: the function decreases up to L* = kz e^{1/delta}, then increases
    if A(delta) * ell1 >= 1:
        mdel = L1 ** A(delta) / ell1
    else:
        mdel = A(delta) * E1 * kz ** A(delta)
    say("[tail] m_delta = %s (lower bound of L^delta/log(L/(1-1/u1)) on L >= L1)" % mdel.str(6))
    def margin(cT):
        cT = A(cT)
        X = beta * (cT * L1).log() / b
        Y = X - KER.ell
        (n0, k1, k2), (m0, k1L, k2L) = KER.ks(Y)
        co = EP.zsum_coeffs(n0, k1, k2); coL = EP.zsum_coeffs(m0, k1L, k2L)
        ZoL = ZoverL(co, n0, L1); ZLoL = ZoverL(coL, m0, L1)
        xB = (B * X).exp(); xhB = ((A(Fr(1, 2)) - B) * X).exp()
        z1 = EP.SUM_INV_G2 * k2 / 2; z1L = EP.SUM_INV_G2 * k2L / 2
        triv = KER.W0 / ((1 - (-2 * Y).exp()) * xB * L1)
        sub = (z1 + xhB * z1L) / L1 + triv + KER.winf * (7 / L1 + 1) / xB
        Ro = ZoL / 2 + xhB * ZLoL / 2 + triv
        taup = tauT * EMG * (1 - 1 / (2 * u1 * u1))
        lhs = taup * (KER.W1 * cT ** beta * mdel - sub / (L1 ** A(s) * ell1))
        rhs = Ro * (arb(2) ** T.pi_w0 / L1 ** A(s) + T.c_sig ** 2 * (KA + KB) / (1 - 1 / u1) ** A(s))
        return lhs - rhs
    lo_, hi_ = 0.01, 100.0
    for _ in range(60):
        mid = math.sqrt(lo_ * hi_)
        if margin(Fr(mid).limit_denominator(10 ** 9)) > 0:
            hi_ = mid
        else:
            lo_ = mid
    cT = Fr(math.ceil(hi_ * 1000 * (1 + 1e-9)), 1000)
    mg = margin(cT)
    pos(mg, "tail final inequality at c_T = %s" % cT)
    say("[tail] every prime p with omega(p-1) > %d: g*(p) <= (%s log p)^((1+s+delta)/(1-B)) = (%.4g log p)^%s  (margin %s; %.0f s)" %
        (W1, cT, float(cT), (beta * INV_b).str(6), mg.str(4), time.time() - t0))
    with open(os.path.join(HERE, out), "a") as f:
        f.write("\n".join(lines) + "\n")
    return cT, beta * INV_b, tauT


# ---------------------------------------------------------------- sec 7: corollaries and the uniform statement
def uniform_c(rowsA, cT, betaT):
    """sec 7.1: c with g*(p) <= (c log p)^{betaT} for every prime p >= 3: max of c_T and, for w <= W1,
    (kappa_w / L_w^{b betaT - 1})^{1/(b betaT)}  (kappa_w L <= c^{b betaT} L^{b betaT} for L >= L_w)."""
    e = betaT * b                     # = 1 + s + delta
    worst = A(cT); ww = None
    for (w, pw, Lm, kap, cw) in rowsA:
        need = (A(kap) / Lm ** (e - 1)) ** (1 / e)
        if need > worst:
            worst = need; ww = w
    c = Fr(math.ceil(float(worst.upper()) * 1000 * (1 + 1e-12)), 1000)
    for (w, pw, Lm, kap, cw) in rowsA:
        pos(A(c) ** e * Lm ** (e - 1) - A(kap), "uniform c covers w = %d" % w)
    pos(A(c) - A(cT), "uniform c >= c_T")
    return c, ww

def corollary_thresholds(tables, rowsA, cT, betaT, wmax_scan=120):
    """sec 7.2: least Lambda such that g*(p) < p (resp. < sqrt(p) - 2) for every prime p with log p >= Lambda.
    For w <= wmax_scan: smallest grid point Lm >= L_w (step 1/4) with certified kappa(Lm) and (1/b) log(kappa Lm) < Lm
    (resp. < Lm/2 - 1/1000; log(x + 2) <= log x + 2/x and x > e^8). Larger w: checked at L_w with the table kappa."""
    P, T, C = tables
    res = {}
    for name, fac, slack in (("p", A(1), A(0)), ("sqrt", A(Fr(1, 2)), A(Fr(1, 1000)))):
        Lam = arb(0); wlam = None
        prod = 1
        for w in range(1, len(rowsA) + 1):
            prod *= P[w - 1]
            Lw = arb(prod + 1).log()
            kap = A(rowsA[w - 1][3])
            def ok(k, Lm):
                # (1/b) log(k L) + slack < fac L at L = Lm, and the difference decreases for L >= Lm (fac b Lm > 1)
                return ((k * Lm).log() / b + slack < fac * Lm) and (fac * b * Lm > 1)
            if ok(kap, Lw):
                continue
            if w > wmax_scan:
                raise SystemExit("corollary: w = %d needs a scan" % w)
            Tlo = arb(T[w - 1]) / ONE
            Lm = Lw
            while True:
                Lm = Lm + A(Fr(1, 4))
                k2, mg = kappa_fast(Lm, Tlo, C[w - 1])
                if ok(A(k2), Lm):
                    break
            Lm_up = arb(Lm.upper())                  # exact upper endpoint: total order
            if Lm_up > Lam:
                Lam = Lm_up; wlam = w
        res[name] = (Lam, wlam)
    # tail (w > W1): log x = betaT log(cT L) < L/2 - 1/1000 at L >= 1e5 (and increasing margin)
    L1 = A(100000)
    pos(L1 / 2 - A(Fr(1, 1000)) - betaT * (A(cT) * L1).log(), "tail bound < sqrt(p) - 2 for log p >= 1e5")
    pos(L1 - 2 * betaT, "monotone: L/2 - betaT log(cT L) increasing for L >= 2 betaT")
    return res

def small_section(N, c, E, Pmin=5):
    """sec 6: exact g*(p) for 3 <= p <= N (small.py); check g*(p) <= (c log p)^E for Pmin <= p <= N (integer/float check
    with a 1e-9 relative safety margin on log g*)."""
    import numpy as np, small as SM
    spf = SM.spf_table(N); pr = SM.primes_upto(N)
    plist = [int(x) for x in pr[:2000]]
    rec = 0; records = []; bad = 0; worst = 0.0; worstp = None; cnt = 0
    lc = math.log(float(c))
    for p in pr:
        p = int(p)
        if p < 3:
            continue
        g = SM.gstar(p, SM.distinct_primes(p - 1, spf), plist)
        cnt += 1
        if g > rec:
            rec = g; records.append((p, g))
        if p >= Pmin:
            rhs = E * math.log(float(c) * math.log(p))
            r = math.log(g) / rhs
            if r > worst:
                worst = r; worstp = (p, g)
            if math.log(g) > rhs * (1 - 1e-9):
                bad += 1
    if bad:
        raise SystemExit("CERT FAIL: small primes: %d violations" % bad)
    # rigorous form: (c log p)^E is increasing in p, so max g* < (c log Pmin)^E proves the bound for all Pmin <= p <= N
    E_ball = 3 * INV_b
    vmin = (A(c) * arb(Pmin).log()) ** E_ball
    pos(vmin - rec, "max g*(p) (p <= N) < (c log %d)^(3/(1-B))" % Pmin)
    print("[small] proof: max g*(p) = %d < (%s log %d)^(3/(1-B)) = %s, increasing in p: bound holds for %d <= p <= %d" % (rec, c, Pmin, vmin.str(6), Pmin, N), flush=True)
    print("[small] %d primes 3 <= p <= %d: max g*(p) = %d; records %s" % (cnt, N, rec, ", ".join("%d:%d" % t for t in records)), flush=True)
    print("[small] g*(p) <= (%s log p)^%.6f for 5 <= p <= %d; worst log g*/log bound = %.4f at p = %d (g* = %d)" % (c, E, N, worst, worstp[0], worstp[1]), flush=True)
    return records


if __name__ == "__main__":
    mode = sys.argv[1] if len(sys.argv) > 1 else "all"
    t0 = time.time()
    if mode == "probe":
        Wmax = int(sys.argv[2])
        Ls = float(sys.argv[3]) if len(sys.argv) > 3 else None
        rows = omega_section(Wmax, Ls, out="logs/kappa_probe_%s.txt" % ("A" if Ls is None else "B"))
        cmax = max(r[4] for r in rows); wmax = max(rows, key=lambda r: r[4])[0]
        print("[omega] max c_w = %.5f at w = %d" % (cmax, wmax))
        for lo, hi in [(1, 5), (6, 10), (11, 20), (21, 50), (51, 100), (101, 300), (301, 1000), (1001, 3000), (3001, 10000)]:
            sub = [r for r in rows if lo <= r[0] <= hi]
            if sub:
                m = max(sub, key=lambda r: r[4])
                print("   w in [%d,%d]: max c_w = %.5f (w = %d)" % (lo, hi, m[4], m[0]))
    if mode == "all":
        EP.core()
        Wm, ych = load_choice()
        tables = sieve_tables(Wm, ych)
        rowsA = omega_section(Wm, None, out="logs/kappa_table.txt", tables=tables)
        cA = max(r[4] for r in rowsA)
        print("[omega] Theorem B (all p, omega(p-1) <= %d): g*(p) <= (%.4f log p)^(3/(1-B)), worst w = %d" % (Wm, cA, max(rowsA, key=lambda r: r[4])[0]), flush=True)
        L7 = arb(10 ** 7).log()
        # version B: log p >= log(1e7) (p > 1e7); p <= 1e7 by the exact table (sec 6)
        P_, T_, C_ = tables
        cB = 0.0; wB = None; prod = 1
        for w in range(1, Wm + 1):
            prod *= P_[w - 1]
            Lw = arb(prod + 1).log()
            Lm = Lw if Lw > L7 else L7
            Tlo = arb(T_[w - 1]) / ONE
            kap, mg = kappa_fast(Lm, Tlo, C_[w - 1])
            cw = float((A(kap) / (Lm * Lm)).upper()) ** (1 / 3) * (1 + 1e-12)
            pos(A(Fr(cw).limit_denominator(10 ** 12)) ** 3 * Lm * Lm - A(kap), "c_w (p > 1e7)")
            if cw > cB:
                cB = cw; wB = w
        cBr = Fr(math.ceil(cB * 1000), 1000)
        print("[omega] Theorem B' (p > 1e7, omega(p-1) <= %d): g*(p) <= (%s log p)^(3/(1-B)), worst w = %d" % (Wm, cBr, wB), flush=True)
        small_section(10 ** 7, cBr, float((3 * INV_b).lower()))
        cT, betaT, tauT = tail_section(Fr(23, 10), Fr(1, 10), W1=Wm, out="logs/tail_output.txt")
        cU, wU = uniform_c(rowsA, cT, betaT)
        print("[final] every prime p >= 3: g*(p) <= (%s log p)^%s   (binding: %s)" % (cU, betaT.str(8), "tail" if wU is None else "w = %d" % wU), flush=True)
        cor = corollary_thresholds(tables, rowsA, cT, betaT)
        for name, (Lam, wl) in cor.items():
            print("[cor] g*(p) < %s for every prime p with log p >= %s (binding w = %s); sufficient: p >= 10^%.3f" % ("p" if name == "p" else "sqrt(p) - 2", Lam.str(6), wl, math.ceil(float(Lam.upper()) / math.log(10) * 1000) / 1000), flush=True)
    if mode == "tailprobe":
        s_ = Fr(sys.argv[2]); d_ = Fr(sys.argv[3]); Q2_ = int(sys.argv[4]) if len(sys.argv) > 4 else 100000
        tail_section(s_, d_, Q2=Q2_, out="logs/tail_probe.txt")
    print("runtime %.1f s" % (time.time() - t0))
