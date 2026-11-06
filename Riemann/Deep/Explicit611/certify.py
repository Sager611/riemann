# Rigorous certificate for |zeta(1/2+it)| <= C t^{1/6} log t  (refinement of arXiv 2602.05614).
# Interval arithmetic: python-flint arb balls (every comparison is "certainly <"; rounding handled by arb).
# Floats (numpy) are used ONLY to choose free parameters (K, c_r, phi); every inequality is re-checked in arb.
# Run: uv run --no-project --with python-flint --with numpy python3 certify.py C [parts]
import sys, math, time
import numpy as np
from flint import arb, acb, ctx
ctx.prec = 96
PI = arb.pi(); SQPI = PI.sqrt(); ONE = arb(1)
C_STR = sys.argv[1] if len(sys.argv) > 1 else "0.573"
PARTS = sys.argv[2] if len(sys.argv) > 2 else "ABCDEF"
C = arb(C_STR)
def log(msg):
    print(msg, flush=True)
def iv(a, b):  # arb ball containing [a, b] (a, b arb)
    return a.union(b)
def rhs_low(t_lo):  # C t^{1/6} log t, increasing in t >= 1; evaluated at the left end
    return C * t_lo.root(6) * t_lo.log()
# ---------- exact partial sums H(N) = sum_{n<=N} n^{-1/2} ----------
_HC = [arb(0)]
def H(N):
    if N < 200000:
        while len(_HC) <= N:
            n = len(_HC); _HC.append(_HC[-1] + 1/arb(n).sqrt())
        return _HC[N]
    return arb(0.5).zeta() - arb(0.5).zeta(arb(N + 1))   # Hurwitz zeta, rigorous in arb
def rem(t_lo):  # Gabcke remainder bound used in the paper (TeX l. 756-760), valid t >= 200; decreasing in t
    return arb('1.463') * t_lo ** (-ONE / 4) + arb('0.127') * t_lo ** (-3 * ONE / 4)
# ---------- block bound (proof of the paper's Lemma 2.9 = Lemma 2.4 + 2.7 + 2.8, per block, partial length L <= K) ----------
def blk_arb(t, K, r, L, c):
    """upper bound for 2|sum_{n=Kr}^{Kr+L-1} n^{-1/2+it}|, t an arb ball, K,r,L ints, c float in (0,3)."""
    cc = arb(c)
    W = PI * arb(r + 1) ** 3 * arb(K) ** 3 / t
    W3 = W.root(3)
    mu = (1 + ONE / r) * (2 + ONE / r + ONE / (K - 1)) / 2           # Lemma 2.5 with K0 = K
    al = 1 / cc + 32 * mu / (15 * SQPI) * (cc + 1 / W3).sqrt() + mu * cc / (3 * W3) + mu / (3 * W3 ** 2)
    be = 32 / (3 * (PI * cc).sqrt()) + (6 / PI - 1) / (cc * W3 ** 2) + (1 - 6 / PI) / W3
    Q = (L + cc * W3) * (L / (cc * W3) + (K / W3) * (al - 1 / cc) + W3 * be)
    vdc = 2 * (Q / (K * r)).sqrt()
    triv = 2 * arb(L) / arb(K * r).sqrt()
    return vdc if vdc.upper() < triv.upper() else triv
# float twin for optimisation
CG = np.linspace(0.30, 2.99, 270)
def blk_f(t, K, r, L):
    W = math.pi * (r + 1) ** 3 * K ** 3 / t; W3 = W ** (1 / 3)
    mu = 0.5 * (1 + 1 / r) * (2 + 1 / r + 1 / (K - 1)); c = CG
    al = 1 / c + 32 * mu / (15 * math.sqrt(math.pi)) * np.sqrt(c + 1 / W3) + mu * c / (3 * W3) + mu / (3 * W3 ** 2)
    be = 32 / (3 * np.sqrt(math.pi * c)) + (6 / math.pi - 1) / (c * W3 ** 2) + (1 - 6 / math.pi) / W3
    Q = (L + c * W3) * (L / (c * W3) + (K / W3) * (al - 1 / c) + W3 * be)
    i = int(np.argmin(Q))
    return min(2 * L / math.sqrt(K * r), 2 / math.sqrt(K * r) * math.sqrt(Q[i])), float(CG[i])
HF = np.concatenate([[0.0], np.cumsum(1 / np.sqrt(np.arange(1, 200001, dtype=float)))])
# ---------- exact Weyl variant: integer M (1 <= M, M-1 < 3 W^{1/3}, i.e. M = ceil(c W^{1/3}) for some c in (0,3)),
# exact weights sum_{m<=mm} (1-m/M) m^e instead of Lemma 2.8 (a),(c),(d); mm = min(M, L-1) since s'_m(L) = 0 for m >= L.
def blk_exact_f(t, K, r, L):
    W = math.pi * (r + 1) ** 3 * K ** 3 / t; W3 = W ** (1 / 3)
    mu = 0.5 * (1 + 1 / r) * (2 + 1 / r + 1 / (K - 1))
    Mmax = int(math.floor(3 * W3 + 1 - 1e-7))
    m = np.arange(1, Mmax + 1, dtype=float)
    X = 4 * mu * K * np.sqrt(m) / math.sqrt(math.pi * W) + mu * K * m / W + 4 * math.sqrt(W / math.pi) / np.sqrt(m) + 1 - 6 / math.pi
    CX = np.concatenate([[0.0], np.cumsum(X)]); CmX = np.concatenate([[0.0], np.cumsum(m * X)])
    Ms = np.arange(1, Mmax + 1); mm = np.minimum(Ms, L - 1)
    S = 2 * (CX[mm] - CmX[mm] / Ms)
    Q = (L + Ms - 1) / Ms * (L + S)
    i = int(np.argmin(Q))
    return min(2 * L / math.sqrt(K * r), 2 * math.sqrt(Q[i] / (K * r))), int(Ms[i])
_WC = {}
def weights(M, mm):
    key = (M, mm)
    if key not in _WC:
        wh = wo = wmh = wz = arb(0)
        for m in range(1, mm + 1):
            f = arb(M - m) / M; sm = arb(m).sqrt()
            wh += f * sm; wo += f * m; wmh += f / sm; wz += f
        _WC[key] = (wh, wo, wmh, wz)
    return _WC[key]
def blk_exact_arb(t, K, r, L, M):
    W = PI * arb(r + 1) ** 3 * arb(K) ** 3 / t
    W3 = W.root(3)
    assert M >= 1 and arb(M - 1) < 3 * W3 and M <= K - 1
    mu = (1 + ONE / r) * (2 + ONE / r + ONE / (K - 1)) / 2
    mm = min(M, L - 1)
    if mm >= 1:
        wh, wo, wmh, wz = weights(M, mm)
        S = 2 * (4 * mu * K / (PI * W).sqrt() * wh + mu * K / W * wo + 4 * (W / PI).sqrt() * wmh + (1 - 6 / PI) * wz)
    else:
        S = arb(0)
    Q = arb(L + M - 1) / M * (L + S)
    vdc = 2 * (Q / (K * r)).sqrt()
    triv = 2 * arb(L) / arb(K * r).sqrt()
    return vdc if vdc.upper() < triv.upper() else triv
def F_exact_f(n1, K, t):
    R = n1 // K
    tot = 2 * HF[4 * K - 1] + 1.463 * t ** -0.25 + 0.127 * t ** -0.75; info = []
    for r in range(4, R):
        b1_, c = blk_f(t, K, r, K); b2_, M = blk_exact_f(t, K, r, K)
        tot += min(b1_, b2_); info.append((c, M))
    LR = n1 - K * R + 1
    b1_, c = blk_f(t, K, R, LR); b2_, M = blk_exact_f(t, K, R, LR)
    tot += min(b1_, b2_); info.append((c, M))
    return tot, info
def F_f(n1, K, t):
    R = n1 // K
    tot = 2 * HF[4 * K - 1] + 1.463 * t ** -0.25 + 0.127 * t ** -0.75; cs = []
    for r in range(4, R):
        b, c = blk_f(t, K, r, K); tot += b; cs.append(c)
    b, c = blk_f(t, K, R, n1 - K * R + 1); tot += b; cs.append(c)
    return tot, cs
# ---------- Part A: direct evaluation on [3, TA] ----------
TA = 300
def partA():
    t0 = time.time(); a = arb(3); worst = 0.0; n = 0
    step = arb('0.02')
    while a < TA:
        b = a + step
        if b > TA: b = arb(TA)
        s = acb(arb(0.5), iv(a, b))
        z = s.zeta().abs_upper()
        rl = rhs_low(a)
        if z < rl:
            worst = max(worst, float((z / (a.root(6) * a.log())).upper())); a = b; n += 1
            step = min(step * 2, arb('0.05'))
        else:
            step = step / 2
            if step < arb('1e-6'):
                log("PART A FAIL at t=%s" % a); return False
    log("PART A ok: |zeta(1/2+it)| < %s t^(1/6) log t on [3,%d], %d balls, max ratio upper %.5f (%.0fs)" % (C_STR, TA, n, worst, time.time() - t0))
    return True
# ---------- Part B: trivial Riemann-Siegel bound (paper Lemma 2.2 style, exact sum) ----------
def partB(nB):
    t0 = time.time(); worst = 0.0; worst_n = None
    for n1 in range(int(math.floor(math.sqrt(TA / (2 * math.pi)))), nB + 1):
        tl = max(arb(TA), 2 * PI * n1 ** 2); th = 2 * PI * (n1 + 1) ** 2
        if not (tl < th): continue
        Fv = 2 * H(n1) + rem(tl)
        rl = rhs_low(tl)
        if not (Fv < rl):
            log("PART B FAIL n1=%d" % n1); return False
        q = float((Fv / (tl.root(6) * tl.log())).upper())
        if q > worst: worst, worst_n = q, n1
    log("PART B ok: trivial RS bound for n1 in [%d,%d] (t in [%d, %.6g)), max ratio %.5f at n1=%s (%.0fs)" % (int(math.floor(math.sqrt(TA / (2 * math.pi)))), nB, TA, 2 * math.pi * (nB + 1) ** 2, worst, worst_n, time.time() - t0))
    return True
# ---------- Part C: refined van der Corput, one n1 interval at a time ----------
def candidates(n1, tlo, thi):
    Kmin = math.ceil(thi ** (1 / 3) + 1e-9); Kmax = math.floor(tlo ** 0.35)
    cs = {Kmin, Kmin + 1, Kmin + 2}
    for Rp in range(max(4, n1 // max(Kmax, 1)), n1 // Kmin + 2):
        k = n1 // Rp
        for kk in (k - 1, k, k + 1):
            cs.add(kk)
    return sorted(k for k in cs if Kmin <= k <= Kmax and n1 // k >= 4)
NEX = int(sys.argv[6]) if len(sys.argv) > 6 else 0     # use the exact-Weyl variant for n1 <= NEX
EXACT_TAIL = NEX > 0                                   # exact-Weyl blocks r < r1 in the tail bound (Parts E, F)
def amin(a, b):
    return a if a.upper() < b.upper() else b
def partC(nC0, nC1, log_every=2000):
    t0 = time.time(); worst = 0.0; worst_n = None; fails = 0
    for n1 in range(nC0, nC1 + 1):
        tlo = 2 * math.pi * n1 ** 2; thi = 2 * math.pi * (n1 + 1) ** 2
        ex = n1 <= NEX
        best = None
        for K in candidates(n1, tlo, thi):
            v, cs = (F_exact_f if ex else F_f)(n1, K, tlo)
            if best is None or v < best[0]: best = (v, K, cs)
        if best is None:
            log("PART C FAIL n1=%d: no admissible K" % n1); return False
        _, K, cs = best
        tl = 2 * PI * n1 ** 2; th = 2 * PI * (n1 + 1) ** 2; t = iv(tl, th)
        # hypotheses: K >= t^{1/3} on the interval, K - 1 < t^{0.35} (so K = ceil(t^phi), phi in [1/3,0.35]), t >= 2012, R >= 4
        assert arb(K) ** 3 > th and arb(K - 1) < tl ** arb('0.35') and tl > 2012
        R = n1 // K; assert R >= 4
        tot = 2 * H(4 * K - 1) + rem(tl)
        for i, r in enumerate(range(4, R + 1)):
            L = K if r < R else n1 - K * R + 1
            if ex:
                c, M = cs[i]
                tot += amin(blk_arb(t, K, r, L, c), blk_exact_arb(t, K, r, L, M))
            else:
                tot += blk_arb(t, K, r, L, cs[i])
        tot = amin(tot, 2 * H(n1) + rem(tl))                 # trivial Riemann-Siegel bound as a fallback
        rl = rhs_low(tl)
        if not (tot < rl):
            log("PART C FAIL n1=%d K=%d bound=%s rhs=%s" % (n1, K, tot, rl)); fails += 1
            if fails > 5: return False
            continue
        q = float((tot / (tl.root(6) * tl.log())).upper())
        if q > worst: worst, worst_n = q, (n1, K)
        if (n1 - nC0) % log_every == 0:
            log("  C n1=%d K=%d R=%d ratio %.5f (%.0fs)" % (n1, K, R, q, time.time() - t0))
    if fails: return False
    log("PART C ok: refined vdC for n1 in [%d,%d] (t in [%.6g, %.6g)), max ratio %.5f at (n1,K)=%s (%.0fs)" % (nC0, nC1, 2 * math.pi * nC0 ** 2, 2 * math.pi * (nC1 + 1) ** 2, worst, worst_n, time.time() - t0))
    return True
# ---------- Part D: subregions [ta, tb], phi fixed, all blocks full (paper style, per-block mu, W, c) ----------
def subregion_f(ta, tb, phi):
    """float estimate of the uniform bound on [ta,tb] with K = ceil(t^phi)."""
    Kup = math.ceil(tb ** phi); Klo = math.ceil(ta ** phi)
    Rmax = int(math.floor(tb ** (0.5 - phi) / math.sqrt(2 * math.pi)))
    tot = 2 * HF[min(4 * Kup - 1, 200000)] if 4 * Kup - 1 <= 200000 else 4 * math.sqrt(4 * Kup) - 2.92
    cs = []
    for r in range(4, Rmax + 1):
        # worst case over t,K approximated at (ta, Klo) for alpha and (tb, Kup) for prefactor: just use mid
        b, c = blk_f(math.sqrt(ta * tb), round(math.sqrt(Klo * Kup)), r, round(math.sqrt(Klo * Kup)))
        tot += b * (tb / ta) ** (1 / 6); cs.append(c)
    return tot, cs, Rmax
def blk_unif_arb(t, Kiv, r, c):
    """uniform full-block bound: 2 t^{1/6}/(pi^{1/6} sqrt(r(r+1))) sqrt((1+c W^{1/3}/K)(alpha+beta W^{2/3}/K)), t, Kiv arb balls"""
    cc = arb(c)
    Wm3 = t.root(3) / (PI.root(3) * (r + 1) * Kiv)        # W^{-1/3}
    W3K = PI.root(3) * (r + 1) / t.root(3)                 # W^{1/3}/K
    W23K = PI.root(3) ** 2 * arb(r + 1) ** 2 * Kiv / t.root(3) ** 2   # W^{2/3}/K
    Km1 = Kiv - 1
    mu = (1 + ONE / r) * (2 + ONE / r + 1 / Km1) / 2
    al = 1 / cc + 32 * mu / (15 * SQPI) * (cc + Wm3).sqrt() + mu * cc * Wm3 / 3 + mu * Wm3 ** 2 / 3
    be = 32 / (3 * (PI * cc).sqrt()) + (6 / PI - 1) * Wm3 ** 2 / cc + (1 - 6 / PI) * Wm3
    return 2 * t.root(6) / (PI.root(6) * arb(r * (r + 1)).sqrt()) * ((1 + cc * W3K) * (al + be * W23K)).sqrt()
def partD(tD0, tD1, rho0=1.02):
    t0 = time.time(); worst = 0.0; nsub = 0
    ta = tD0
    while ta < tD1:
        rho = rho0
        while True:
            tb = min(ta * rho, tD1)
            # choose phi by a small float search
            best = None
            for phi in np.linspace(1 / 3 + 1e-6, 0.35 - 1e-6, 25):
                v, cs, Rm = subregion_f(ta, tb, float(phi))
                if best is None or v < best[0]: best = (v, float(phi), cs, Rm)
            _, phi, cs, Rmax = best
            ta_a = arb(repr(ta)); tb_a = arb(repr(tb)); t = iv(ta_a, tb_a)
            ph = arb(repr(phi))
            assert ph > ONE / 3 and ph < arb('0.35')
            Kiv = iv(ta_a ** ph, tb_a ** ph + 1)                 # K = ceil(t^phi) in [ta^phi, tb^phi + 1)
            Rbound = tb_a ** (ONE / 2 - ph) / (2 * PI).sqrt()   # R <= t^{1/2-phi}/sqrt(2 pi)
            assert arb(Rmax + 1) > Rbound
            Kup = math.ceil(float((tb_a ** ph).upper()))
            tot = 2 * H(4 * Kup - 1) + rem(ta_a)
            for i, r in enumerate(range(4, Rmax + 1)):
                tot += blk_unif_arb(t, Kiv, r, cs[i])
            rl = rhs_low(ta_a)
            if tot < rl:
                q = float((tot / (ta_a.root(6) * ta_a.log())).upper()); worst = max(worst, q); nsub += 1
                break
            rho = 1 + (rho - 1) / 2
            if rho < 1 + 1e-6:
                log("PART D FAIL at ta=%.6g" % ta); return False
        ta = tb
    log("PART D ok: %d subregions covering [%.6g, %.6g], max ratio %.5f (%.0fs)" % (nsub, tD0, tD1, worst, time.time() - t0))
    return True
# ---------- Part E/F: uniform bound for all t >= T (phi = 1/3), per-block mu_r, no sqrt split ----------
def alpha_hat(r, c, T3):   # T3 = T^{1/3} (arb); bound for alpha_r valid for all t >= T, K = ceil(t^{1/3})
    cc = arb(c)
    mu = (1 + ONE / r) * (2 + ONE / r + 1 / (T3 - 1)) / 2
    w = 1 / (PI.root(3) * (r + 1))                          # W^{-1/3} <= 1/(pi^{1/3}(r+1))
    return 1 / cc + 32 * mu / (15 * SQPI) * (cc + w).sqrt() + mu * cc * w / 3 + mu * w ** 2 / 3
def beta_hat(r, c, kappa):  # bound for beta_r: W^{-2/3} <= 1/(pi^{2/3}(r+1)^2), W^{1/3} <= pi^{1/3}(r+1)kappa
    cc = arb(c)
    return 32 / (3 * (PI * cc).sqrt()) + (6 / PI - 1) / (cc * PI.root(3) ** 2 * arb(r + 1) ** 2) + (1 - 6 / PI) / (PI.root(3) * (r + 1) * kappa)
def tail_coeffs(T, r1, c_tail, c_list, M_list=None):
    """|zeta(1/2+it)| <= b1 t^{1/6} log t + b2 t^{1/6} + b3 for all t >= T (T arb). c_list[r-4] for r in [4, r1)."""
    T3 = T.root(3); T6 = T.root(6)
    assert T6 / (2 * PI).sqrt() > r1 + 1
    kappa = 1 + 1 / T3                                        # K/t^{1/3} <= 1 + t^{-1/3}
    e1 = 1 / ((2 * PI).sqrt() * T6) + 1 / T3                  # (r+1)/t^{1/3} <= e1 for r <= R
    S1 = arb(0)
    kiv = iv(ONE, kappa); uiv = iv(arb(0), 1 / T3)              # k = K/t^{1/3} in [1, kappa], u = t^{-1/3} in (0, T^{-1/3}]
    for r in range(4, r1):
        c = c_list[r - 4]; cc = arb(c)
        epsc = cc * PI.root(3) * e1
        term = ((1 + epsc) * (alpha_hat(r, c, T3) + beta_hat(r, c, kappa) * PI.root(3) ** 2 * kappa * arb(r + 1) ** 2 / T3)).sqrt() / arb(r * (r + 1)).sqrt()
        if M_list is not None:                                  # exact-Weyl block (integer M, exact weights), W = pi (r+1)^3 k^3
            M = M_list[r - 4]
            assert M >= 1 and arb(M - 1) < 3 * PI.root(3) * (r + 1) and arb(M) < T3 - 1
            mu = (1 + ONE / r) * (2 + ONE / r + 1 / (T3 - 1)) / 2
            wh, wo, wmh, wz = weights(M, M)
            r32 = arb(r + 1) ** 3
            S1c = 2 * (4 * mu / (PI * r32.sqrt() * kiv.sqrt())) * wh + 2 * (mu / (PI * r32 * kiv ** 2)) * wo
            S0c = 8 * (r32 * kiv ** 3).sqrt() * wmh + 2 * (1 - 6 / PI) * wz   # 2 sum (1-m/M)(4 sqrt(W/pi) m^{-1/2} + 1 - 6/pi)
            ex = PI.root(6) * ((kiv + (M - 1) * uiv) * (kiv + S1c + S0c * uiv) / (M * kiv * r)).sqrt()   # = (pi^{1/6}/2) B_r / t^{1/6}
            term = amin(term, ex)
        S1 += term
    cc = arb(c_tail)
    epsc = cc * PI.root(3) * e1
    A = alpha_hat(r1, c_tail, T3); a = A.sqrt()
    Bhat = 32 / (3 * (PI * cc).sqrt()) + (6 / PI - 1) / (cc * PI.root(3) ** 2 * arb(r1 + 1) ** 2)
    Bq = Bhat * PI.root(3) ** 2 * kappa
    Bq1 = Bq * (1 + ONE / r1)
    s_r1 = (A / (r1 * (r1 + 1)) + Bq * (r1 + 1) / (r1 * T3)).sqrt()
    bstar = Bq1.sqrt() / (2 * PI).sqrt()
    GR = (A + bstar ** 2).sqrt() - a * (a / bstar).asinh()
    mG = -a + a * (2 * a / (r1 * Bq1.sqrt())).log() + arb(r1) ** 2 * Bq1 / (4 * a * T3)
    sq = (1 + epsc).sqrt()
    b1 = sq * a / (3 * PI.root(6))
    b2 = 2 / PI.root(6) * (S1 + sq * (s_r1 + GR + mG)) + 8 * (1 + 3 / (8 * T3))
    N0 = int(math.floor(float((4 * T3).lower())))             # N0 <= 4 T^{1/3} <= 4K, so 4K-1 >= N0-1
    assert arb(N0) <= 4 * T3
    b3 = 2 * H(N0 - 1) - 4 * arb(N0 - 1).sqrt() + rem(T)
    return b1, b2, b3
def tail_ratio(T, b1, b2, b3):
    """sup_{t>=T} (b1 t^{1/6}log t + b2 t^{1/6} + b3)/(t^{1/6} log t) <= value returned (checks monotonicity)."""
    L = T.log(); T6 = T.root(6)
    assert b2 > 0
    if not (b3 >= 0):
        assert b2 > (-b3) * (L / 6 + 1) / T6     # derivative negative for all t >= T (RHS decreasing in t for log t > 6)
        assert L > 6
    return b1 + b2 / L + b3 / (T6 * L)
def tail_opt(Tf, r1):
    """float optimisation of c_r (r<r1) and c_tail for the tail bound at T = Tf."""
    T3 = Tf ** (1 / 3); kappa = 1 + 1 / T3; e1 = 1 / (math.sqrt(2 * math.pi) * Tf ** (1 / 6)) + 1 / T3
    p3 = math.pi ** (1 / 3)
    def ah(r, c):
        mu = 0.5 * (1 + 1 / r) * (2 + 1 / r + 1 / (T3 - 1)); w = 1 / (p3 * (r + 1))
        return 1 / c + 32 * mu / (15 * math.sqrt(math.pi)) * np.sqrt(c + w) + mu * c * w / 3 + mu * w ** 2 / 3
    def bh(r, c):
        return 32 / (3 * np.sqrt(math.pi * c)) + (6 / math.pi - 1) / (c * p3 ** 2 * (r + 1) ** 2) + (1 - 6 / math.pi) / (p3 * (r + 1) * kappa)
    cl = []
    for r in range(4, r1):
        v = (1 + CG * p3 * e1) * (ah(r, CG) + bh(r, CG) * p3 ** 2 * kappa * (r + 1) ** 2 / T3)
        cl.append(float(CG[int(np.argmin(v))]))
    return cl
def tail_optM(Tf, r1):
    """float choice of integer M_r (r < r1) for the exact-Weyl tail blocks, at k = 1, u = T^{-1/3}."""
    T3 = Tf ** (1 / 3); u = 1 / T3; Ml = []
    for r in range(4, r1):
        mu = 0.5 * (1 + 1 / r) * (2 + 1 / r + 1 / (T3 - 1)); r32 = (r + 1) ** 1.5
        Mmax = int(math.floor(3 * math.pi ** (1 / 3) * (r + 1) + 1 - 1e-7))
        best = None
        for M in range(1, Mmax + 1):
            m = np.arange(1, M + 1, dtype=float); w = (M - m) / M
            S1c = 2 * 4 * mu / (math.pi * r32) * np.sum(w * np.sqrt(m)) + 2 * mu / (math.pi * r32 ** 2) * np.sum(w * m)
            S0c = 8 * r32 * np.sum(w / np.sqrt(m)) + 2 * (1 - 6 / math.pi) * np.sum(w)
            v = (1 + (M - 1) * u) * (1 + S1c + S0c * u) / M
            if best is None or v < best[0]: best = (v, M)
        Ml.append(best[1])
    return Ml
def best_tail(T_str, r1s=(5, 6, 8, 10, 12, 16, 20, 40, 80, 160, 320), cts=None):
    T = arb(T_str); Tf = float(T_str)
    best = None
    cts = cts if cts is not None else [1.0 + 0.02 * i for i in range(60)]
    for r1 in r1s:
        if Tf ** (1 / 6) / math.sqrt(2 * math.pi) <= r1 + 2: continue
        cl = tail_opt(Tf, r1)
        Ml = tail_optM(Tf, r1) if EXACT_TAIL else None
        for ct in cts:
            b1, b2, b3 = tail_coeffs(T, r1, ct, cl, Ml)
            q = tail_ratio(T, b1, b2, b3)
            if best is None or q.upper() < best[0].upper(): best = (q, r1, ct, b1, b2, b3)
    return best
def partE(T_str):
    t0 = time.time()
    q, r1, ct, b1, b2, b3 = best_tail(T_str)
    ok = q < C
    log("PART E %s: for t >= %s: |zeta| <= %.6f t^(1/6) log t + %.6f t^(1/6) + %.6f  (upper ends), sup ratio <= %.6f (r1=%d, c_tail=%.2f) (%.0fs)" % ("ok" if ok else "FAIL", T_str, float(b1.upper()), float(b2.upper()), float(b3.upper()), float(q.upper()), r1, ct, time.time() - t0))
    return ok
def partF(T_str="8.97e17"):
    t0 = time.time()
    q, r1, ct, b1, b2, b3 = best_tail(T_str)
    log("PART F: for t >= %s: |zeta| <= %.6f t^(1/6) log t + %.6f t^(1/6) + %.6f (upper ends); sup ratio <= %.6f (r1=%d, c_tail=%.2f) (%.0fs)" % (T_str, float(b1.upper()), float(b2.upper()), float(b3.upper()), float(q.upper()), r1, ct, time.time() - t0))
    return q
if __name__ == "__main__":
    log("C = %s, parts %s, prec %d" % (C_STR, PARTS, ctx.prec))
    nB = int(sys.argv[3]) if len(sys.argv) > 3 else 560
    nC1 = int(sys.argv[4]) if len(sys.argv) > 4 else 39894
    TD1 = sys.argv[5] if len(sys.argv) > 5 else "1e13"
    ok = True
    if "A" in PARTS: ok &= partA()
    if "B" in PARTS: ok &= partB(nB)
    if "C" in PARTS: ok &= partC(nB + 1, nC1)
    if "D" in PARTS: ok &= partD(2 * math.pi * (nC1 + 1) ** 2, float(TD1))
    if "E" in PARTS: ok &= partE(TD1)
    if "F" in PARTS: partF("8.97e17")
    log("ALL PARTS OK" if ok else "SOME PART FAILED")
