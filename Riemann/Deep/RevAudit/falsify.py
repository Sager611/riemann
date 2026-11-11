# RevAudit falsification tests for arXiv 2602.05614v1 (Revers), Lemmas 2.4-2.9 and the Riemann-Siegel input.
# Exact exponential sums in Arb balls (python-flint), 128 bits; a case FAILS only if the LHS lower end exceeds the RHS upper end;
# a case is SUSPECT if the balls overlap (never happened = reported). Single thread, fixed seed.
# Run: uv run --no-project --with python-flint --with mpmath --with numpy python3 falsify.py [tests] [ncases] [nscan]
# (logs/falsify.log: tests 28,24,25,29,29x,27,27x,26,P,b,bind,RS, ncases 3000, nscan 80; about 3 min, < 60 MB)
import sys, math, random, time
import numpy as np
from flint import arb, acb, ctx
ctx.prec = 128
PI = arb.pi(); SQPI = PI.sqrt(); ONE = arb(1)
TESTS = sys.argv[1] if len(sys.argv) > 1 else "28,24,29,29x,27,27x,26,bind,P,b,RS"
NC = int(sys.argv[2]) if len(sys.argv) > 2 else 1500
NSCAN = int(sys.argv[3]) if len(sys.argv) > 3 else 40
rng = random.Random(20261111)
def log(*a): print(*a, flush=True)
def mu_of(r, K0):  # Lemma 2.5: mu = 1/2 lambda^{2/3}(1 + 1/((1-1/K0) lambda^{1/3})) = 1/2 (1+1/r)(2+1/r+1/(K0-1))
    return (1 + ONE / r) * (2 + ONE / r + 1 / (K0 - 1)) / 2
def Wof(t, K, r): return PI * arb(r + 1) ** 3 * arb(K) ** 3 / t
def E_list(t, A, L):  # e(f(n)) = (A+n)^{it}, n = 0..L-1, A = Kr
    return [acb(0, t * arb(A + n).log()).exp() for n in range(L)]
def absq(z): return (z.real ** 2 + z.imag ** 2)
class Stat:
    def __init__(s, name): s.name, s.n, s.fail, s.susp, s.worst, s.arg = name, 0, 0, 0, -1.0, None
    def add(s, lhs, rhs, info):  # lhs, rhs arb; ratio = lhs/rhs
        s.n += 1
        if lhs.lower() > rhs.upper():
            s.fail += 1
            if s.fail <= 5: log("  FAIL", s.name, info, lhs, rhs)
        elif not (lhs < rhs): s.susp += 1
        q = float((lhs / rhs).mid()) if rhs > 0 else float('inf')
        if q > s.worst: s.worst, s.arg = q, info
    def rep(s): log("%s: %d cases, %d FAIL, %d overlap; max LHS/RHS = %.5f at %s" % (s.name, s.n, s.fail, s.susp, s.worst, s.arg))
# ---------------- Lemma 2.8 (a)-(d), all M <= MMAX, cumulative Arb sums ----------------
def test28(MMAX=1000000):
    t0 = time.time(); Ah = Ao = At = arb(0); worst = [(-9.0, 0)] * 4; bad = 0
    for M in range(1, MMAX + 1):
        sm = arb(M).sqrt(); Ah += 1 / sm; Ao += sm; At += M * sm
        a = Ah - Ao / M; c = Ao - At / M                     # sum (1-m/M) m^{-1/2}, sum (1-m/M) m^{1/2}
        ra = a - 4 * sm / 3; rc = c - 4 * M * sm / 15        # must be < 0
        if not (ra < 0 and rc < 0):
            bad += 1
            if bad <= 5: log("  2.8 FAIL/overlap M=%d" % M, ra, rc)
        if M % 1000 == 1 or M <= 2000:
            worst[0] = max(worst[0], (float(ra.mid()), M)); worst[2] = max(worst[2], (float(rc.mid()), M))
    for M in range(1, 3001):   # (b), (d) exact integers (times M): sum(M-m) = M(M-1)/2 ; sum (M-m) m = M(M^2-1)/6 <= M^3/6
        b2 = sum(2 * (M - m) for m in range(1, M + 1)); d6 = sum(6 * (M - m) * m for m in range(1, M + 1))
        if b2 != M * (M - 1) or d6 != M * (M * M - 1) or d6 > M ** 3:
            bad += 1
            if bad <= 5: log("  2.8(b,d) FAIL M=%d" % M)
    log("2.8: M<=%d (a),(c) in Arb, M<=3000 (b),(d) exact: %d failures; max (LHS-RHS) (a) %.5f at M=%d, (c) %.5f at M=%d (%.0fs)"
        % (MMAX, bad, worst[0][0], worst[0][1], worst[2][0], worst[2][1], time.time() - t0))
# ---------------- Lemma 2.5 in the stronger form K0 = infinity (exact rationals): for 0 <= m <= K-1,
# K^2 [1/(Kr(Kr+m)) - 1/((Kr+K-1-m)(Kr+K-1))] <= 1/r^2 - 1/(r+1)^2  (= 2 mu/(r+1)^3 with K0 = infinity) ----------------
def test25(KMAX=120, RMAX=40):
    from fractions import Fraction as F
    t0 = time.time(); bad = 0; worst = (F(0), None); n = 0
    for K in range(2, KMAX + 1):
        for r in range(1, RMAX + 1):
            A = K * r; rhs = F(1, r * r) - F(1, (r + 1) ** 2)
            for m in range(0, K):
                lhs = K * K * (F(1, A * (A + m)) - F(1, (A + K - 1 - m) * (A + K - 1))); n += 1
                if lhs > rhs:
                    bad += 1
                    if bad <= 5: log("  2.5 FAIL", K, r, m)
                q = lhs / rhs
                if q > worst[0]: worst = (q, (K, r, m))
    log("2.5 (K0 = infinity, exact): %d cases K<=%d r<=%d all m: %d failures, max LHS/RHS = %.6f at (K,r,m)=%s (%.0fs)"
        % (n, KMAX, RMAX, bad, float(worst[0]), worst[1], time.time() - t0))
# ---------------- Lemma 2.4 (Weyl differencing), random phases, double precision with margin ----------------
def test24(N=3000):
    t0 = time.time(); worst = 0.0; bad = 0
    for i in range(N):
        L = rng.randint(1, 300); M = rng.randint(1, 80); typ = rng.random()
        n = np.arange(L, dtype=float)
        if typ < 0.4: ph = rng.uniform(1e2, 1e7) * np.log(rng.randint(1, 10 ** 6) + n)
        elif typ < 0.8: ph = rng.uniform(-1, 1) * n ** 2 + rng.uniform(-1, 1) * n ** 3 / 100 + rng.uniform(0, 7) * n
        else: ph = np.array([rng.uniform(0, 7) for _ in range(L)])   # arbitrary phases (lemma holds for any f)
        e = np.exp(1j * ph); S = abs(e.sum()) ** 2
        acc = L + 2 * sum((1 - m / M) * abs((e[m:] * np.conj(e[:L - m])).sum()) for m in range(1, min(M, L)))
        rhs = (L + M - 1) / M * acc
        q = S / rhs
        if S > rhs * (1 + 1e-9): bad += 1; log("  2.4 FAIL", L, M, S, rhs)
        worst = max(worst, q)
    log("2.4: %d random cases, %d failures, max LHS/RHS = %.6f (%.0fs)" % (N, bad, worst, time.time() - t0))
# ---------------- parameter sampler (Lemma 2.7/2.9 setting) ----------------
def sample(tlo, thi, ext=False):
    while True:
        u = rng.random()
        t = arb(repr(10 ** rng.uniform(math.log10(tlo), math.log10(thi)))) if u > 0.03 else arb(repr(float(tlo)))
        tf = float(t.mid()); n1 = int(math.floor(math.sqrt(tf / (2 * math.pi))))
        if ext:   # Explicit611 item 8 setting: any integer K with K^3 > t, K-1 < t^0.35, R >= 4
            Kmin = math.ceil(tf ** (1 / 3) + 1e-9); Kmax = min(math.floor(tf ** 0.35) + 1, n1 // 4)
            if Kmax < Kmin: continue
            K = rng.choice([Kmin, Kmax, rng.randint(Kmin, Kmax)]); phi = None
        else:
            phi = rng.choice([1 / 3, 0.35, rng.uniform(1 / 3, 0.35)])
            K = math.ceil(tf ** phi)
        assert arb(K) ** 3 >= t
        R = n1 // K
        if R < 4: continue
        r = rng.choice([4, R, rng.randint(4, R)])
        c = rng.choice([0.01, 2.999, rng.uniform(0.01, 2.999), rng.uniform(0.5, 2.5)])
        return t, tf, phi, K, n1, R, r, c
# ---------------- Lemma 2.9 full inequality (exact sums) ----------------
def rhs29(t, K, r, c, K0):
    W = Wof(t, K, r); W3 = W.root(3); mu = mu_of(r, K0); cc = arb(repr(c))
    al = 1 / cc + 32 * mu / (15 * SQPI) * (cc + 1 / W3).sqrt() + mu * cc / (3 * W3) + mu / (3 * W3 ** 2)
    be = 32 / (3 * (PI * cc).sqrt()) + (6 / PI - 1) / (cc * W3 ** 2) + (1 - 6 / PI) / W3
    assert al > 0 and be > 0
    return (K / W3 + cc) * (al * K + be * W3 ** 2)
def admissible(tf, K, r):  # K = ceil(t^phi) for some phi in [1/3, 0.35] <=> K^3 >= t and K-1 < t^0.35; and r <= R
    n1 = int(math.floor(math.sqrt(tf / (2 * math.pi))))
    return K ** 3 >= tf and K - 1 < tf ** 0.35 and 4 <= r <= n1 // K
def scan_t(tf, K, r, L, nS, m=None, dlf=None):
    """adversarial scan: t_j = t + j*delta (delta = 0.25 Kr/L, the peak width), j < nS; returns list of (t_j arb, value arb)
    value = |sum_{n<L} (Kr+n)^{i t_j}|^2 (m None) or |s'_m(L)| at t_j. E updated by multiplication with exp(i delta log(Kr+n)),
    re-anchored by direct evaluation every 16 steps (rectangular acb balls grow by up to sqrt 2 per multiplication)."""
    A = K * r; dl = arb(repr(dlf if dlf is not None else 0.25 * A / L)); t = arb(repr(tf))
    U = [acb(0, dl * arb(A + n).log()).exp() for n in range(L)]; out = []
    for j in range(nS):
        tj = t + j * dl
        if j % 16 == 0: E = E_list(tj, A, L)
        if m is None: v = absq(sum(E, acb(0)))
        else: v = abs(sum((E[n + m] * E[n].conjugate() for n in range(0, L - m)), acb(0)))
        out.append((tj, v)); E = [E[n] * U[n] for n in range(L)]
    return out
def test29(N, ext=False, deadline=900, nscan=40, nS=250):
    st = Stat("2.9" + (" ext t in [1.5e6,5e6)" if ext else " t in [5e6,1e11]")); t0 = time.time()
    sc = Stat("2.9 adversarial t-scans" + (" ext" if ext else ""))
    for i in range(N):
        if time.time() - t0 > deadline: break
        t, tf, phi, K, n1, R, r, c = sample(1.5e6, 5e6, True) if ext else sample(5e6, 1e11)
        L = rng.choice([1, 2, K, K - 1, rng.randint(1, K), rng.randint(K // 2, K)])
        if not ext: assert tf >= 5e6
        S = sum(E_list(t, K * r, L), acb(0))
        st.add(absq(S), rhs29(t, K, r, c, arb(K)), "t=%.4g K=%d r=%d L=%d c=%.3f" % (tf, K, r, L, c))
        if i < nscan and L >= 8:
            tend = tf + nS * 0.25 * K * r / L
            if not (admissible(tf, K, r) and admissible(tend, K, r)): continue
            for tj, v in scan_t(tf, K, r, L, nS):
                st_rhs = rhs29(tj, K, r, c, arb(K)); sc.add(v, st_rhs, "t=%.6g K=%d r=%d L=%d c=%.3f" % (float(tj.mid()), K, r, L, c))
    st.rep(); sc.rep(); log("   (%.0fs)" % (time.time() - t0))
# ---------------- Lemma 2.7 per-m bound (exact sums) and 2.7(b),(c) ----------------
def rhs27(t, K, r, m, K0):
    W = Wof(t, K, r); mu = mu_of(r, K0)
    return 4 * mu * K / (PI * W).sqrt() * arb(m).sqrt() + mu * K / W * m + 4 * (W / PI).sqrt() / arb(m).sqrt() + 1 - 6 / PI
def test27(N, ext=False, deadline=600, nscan=40):
    st = Stat("2.7" + (" ext t in [1.5e6,5e6)" if ext else " t in [5e6,1e11]")); t0 = time.time(); bcheck = [9e9, 9e9]
    sc = Stat("2.7 adversarial t-scans" + (" ext" if ext else ""))
    for i in range(N):
        if time.time() - t0 > deadline: break
        t, tf, phi, K, n1, R, r, c = sample(1.5e6, 5e6, True) if ext else sample(5e6, 1e11)
        W = Wof(t, K, r); W3 = W.root(3); cc = arb(repr(c))
        M = int(math.ceil(float((cc * W3).mid())))
        if not (arb(M - 1) < cc * W3 and cc * W3 <= M): continue     # M = ceil(c W^{1/3}) certified
        # 2.7(b): K >= c W^{1/3} + 2 ; 2.7(c): Delta = sqrt(M/(pi W)) < 0.154088
        bb = arb(K) - cc * W3 - 2; Dl = (arb(M) / (PI * W)).sqrt()
        assert bb > 0 and Dl < arb('0.154088'), ("2.7(b)/(c) FAIL", tf, K, r, c)
        bcheck[0] = min(bcheck[0], float(bb.mid())); bcheck[1] = min(bcheck[1], float((arb('0.154088') - Dl).mid()))
        m = rng.choice([1, M, rng.randint(1, M)]); L = rng.choice([K, m + 1, rng.randint(m + 1, K)]) if m + 1 <= K else K
        E = E_list(t, K * r, L)
        s = sum((E[n + m] * E[n].conjugate() for n in range(0, L - m)), acb(0))
        st.add(abs(s), rhs27(t, K, r, m, arb(K)), "t=%.4g K=%d r=%d m=%d M=%d L=%d" % (tf, K, r, m, M, L))
        if nscan > 0 and L - m >= 8:
            nscan -= 1; dlf = 0.5 * (K * r) ** 2 / (m * L)     # s'_m(t) has frequencies ~ m/(Kr+n): peak width ~ 2 pi (Kr)^2/(mL)
            tend = tf + 250 * dlf
            if admissible(tf, K, r) and admissible(tend, K, r):
                for tj, v in scan_t(tf, K, r, L, 250, m, dlf):
                    Wj = Wof(tj, K, r); Mj = (cc * Wj.root(3)).ceil() if hasattr(arb, 'ceil') else None
                    if not (arb(m) <= 3 * Wj.root(3)): continue      # m <= M = ceil(c W^{1/3}) still possible for some c < 3
                    sc.add(v, rhs27(tj, K, r, m, arb(K)), "t=%.6g K=%d r=%d m=%d L=%d" % (float(tj.mid()), K, r, m, L))
    st.rep(); sc.rep(); log("   min K-(cW^{1/3}+2) = %.3f, min 0.154088-Delta(M) = %.5f (%.0fs)" % (bcheck[0], bcheck[1], time.time() - t0))
# ---------------- Lemma 2.6(b) with random Delta (independent implementation), incl. m = 1 ----------------
def test26(N, deadline=300):
    st = Stat("2.6(b) random Delta"); t0 = time.time(); skipped = 0
    for i in range(N):
        if time.time() - t0 > deadline: break
        tf = 10 ** rng.uniform(2, 9); t = arb(repr(tf)); K = rng.randint(2, 150); r = rng.randint(2, 40)
        m = rng.choice([1, 2, rng.randint(1, K - 1)]); L = rng.randint(m + 1, K) if m + 1 <= K else None
        if L is None: continue
        D = rng.choice([0.01, 0.49, rng.uniform(0.01, 0.49)]); Dl = arb(repr(D)); K0 = rng.choice([arb(K), arb(K + 1) / 2])
        W = Wof(t, K, r); mu = mu_of(r, K0); A = K * r
        gp = lambda x: t / (2 * PI) * (1 / arb(A + m + x) - 1 / arb(A + x))
        g0, g1 = gp(0), gp(L - 1 - m)
        if g0.floor().upper() != g0.floor().lower() or g1.floor().upper() != g1.floor().lower(): skipped += 1; continue
        e1 = g0 - g0.floor(); e2 = g1 - g1.floor()
        Bc = 1 + 2 * W * Dl / m + 2 / (PI * Dl)
        def p(e):
            if e < Dl: return 1 + 2 / (PI * Dl) + W / m * (Dl - e)
            if e.upper() <= (1 - Dl).lower() and e.lower() >= Dl.upper(): return (1 / Dl + 1 / e) / PI
            if e > 1 - Dl: return arb(0)
            return None
        p1, p2 = p(e1), p(1 - e2)
        if p1 is None or p2 is None: skipped += 1; continue
        rhs = (m * K / W * mu - 1) * Bc + (Bc * e1 + p1) + (Bc * (1 - e2) + p2) - 2 / (PI * Dl)
        E = E_list(t, A, L)
        s = sum((E[n + m] * E[n].conjugate() for n in range(0, L - m)), acb(0))
        st.add(abs(s), rhs, "t=%.4g K=%d r=%d m=%d L=%d D=%.3f" % (tf, K, r, m, L, D))
    st.rep(); log("   skipped (ball straddles integer / eps at a case boundary): %d (%.0fs)" % (skipped, time.time() - t0))
# ---------------- 2.7(d): P(eps) <= 1 + 5/(pi Delta) - 3/pi, grid, then Arb at the worst points ----------------
def testP():
    t0 = time.time()
    D = np.linspace(1e-4, 1 / 6, 3000)[:, None]; e = np.linspace(0, 1, 4001)[None, :]
    B = 1 + 4 / (np.pi * D)
    p = np.where(e < D, 1 + 2 / (np.pi * D) + (D - e) / (np.pi * D ** 2), np.where(e <= 1 - D, (1 / D + 1 / np.maximum(e, 1e-300)) / np.pi, 0.0))
    P = B * e + p; marg = (1 + 5 / (np.pi * D) - 3 / np.pi) - P
    i = np.unravel_index(np.argmin(marg), marg.shape)
    # Arb recheck at the worst grid point and at the case boundaries eps = Delta, 1 - Delta for that Delta
    Dw = arb(repr(float(D[i[0], 0]))); worstA = None
    for ee in (arb(repr(float(e[0, i[1]]))), Dw, 1 - Dw):
        if ee < Dw: pp = 1 + 2 / (PI * Dw) + (Dw - ee) / (PI * Dw ** 2)
        elif ee <= 1 - Dw: pp = (1 / Dw + 1 / ee) / PI
        else: pp = arb(0)
        mg = 1 + 5 / (PI * Dw) - 3 / PI - ((1 + 4 / (PI * Dw)) * ee + pp)
        worstA = mg if worstA is None or mg < worstA else worstA
    log("2.7(d) P-bound: grid 3000 Delta x 4001 eps: min margin %.6f at Delta=%.5f eps=%.5f; Arb recheck min margin %s (%.0fs)"
        % (marg[i], D[i[0], 0], e[0, i[1]], worstA.str(5), time.time() - t0))
# ---------------- 2.7(b) inequality at t = 2012 ----------------
def testb():
    for tt in ('2012', '2011', '1.5e6', '5e6'):
        t = arb(tt); v = t.root(3) - (3 * PI.root(3) * (t.root(6) / (2 * PI).sqrt() + 1) + 2)
        log("2.7(b): t=%s: t^{1/3} - (3 pi^{1/3}(t^{1/6}/sqrt(2pi)+1)+2) = %s" % (tt, v.str(6)))
# ---------------- binding interval of 0.568: n1 = 605, K = 133, R = 4, last block L = 74 ----------------
def testbind(npts=400):
    import mpmath as mp
    t0 = time.time(); n1, K = 605, 133; R = n1 // K; L = n1 - K * R + 1; r = R; A = K * r
    tlo, thi = 2 * math.pi * n1 ** 2, 2 * math.pi * (n1 + 1) ** 2
    worst_blk = worst_m = worst_chain = 0.0; worst_zeta = 0.0
    H = sum(1 / arb(n).sqrt() for n in range(1, 4 * K))
    for j in range(npts):
        tf = tlo + (thi - tlo) * (j + rng.random()) / npts; t = arb(repr(tf))
        W = Wof(t, K, r); W3 = W.root(3); mu = mu_of(r, arb(K))
        E = E_list(t, A, L); blk = sum(E, acb(0))
        actual = 2 * abs(sum((E[n] / arb(A + n).sqrt() for n in range(L)), acb(0)))
        # exact-Weyl bound over integer M with M-1 < 3W^{1/3}, M <= K-1 (Lemma 2.4 + 2.7), and per-m check vs actual |s'_m|
        Mmax = int(math.floor(float((3 * W3).lower()))) ; best = None
        for M in range(1, Mmax + 1):
            S = arb(0)
            for m in range(1, min(M, L - 1) + 1):
                X = rhs27(t, K, r, m, arb(K)); S += 2 * (1 - arb(m) / M) * X
            Q = arb(L + M - 1) / M * (L + S)
            best = Q if best is None or Q.upper() < best.upper() else best
        for m in range(1, min(Mmax, L - 1) + 1):
            s = sum((E[n + m] * E[n].conjugate() for n in range(0, L - m)), acb(0))
            worst_m = max(worst_m, float((abs(s) / rhs27(t, K, r, m, arb(K))).mid()))
        worst_blk = max(worst_blk, float((absq(blk) / best).mid()))
        bound = 2 * (best / (K * r)).sqrt()
        worst_chain = max(worst_chain, float((actual / bound).mid()))
        if j % 40 == 0:
            z = abs(mp.zeta(mp.mpc(0.5, tf))); rhsC = 0.568 * tf ** (1 / 6) * math.log(tf)
            tot = 2 * H + bound + arb('1.463') * t ** (-ONE / 4) + arb('0.127') * t ** (-3 * ONE / 4)
            worst_zeta = max(worst_zeta, float(z) / float(tot.mid()))
            log("  bind t=%.6g: |zeta|=%.4f  chain bound=%.4f  0.568 t^(1/6)log t=%.4f  block actual=%.4f bound=%.4f"
                % (tf, float(z), float(tot.mid()), rhsC, float(actual.mid()), float(bound.mid())))
    log("bind n1=605,K=133,R=4,L=74 (%d t-points): max |block|^2/Q_exact %.4f; max |s'_m|/Lemma2.7 %.4f; max 2|sum n^{-1/2+it}|/bound %.4f; max |zeta|/chain %.4f (%.0fs)"
        % (npts, worst_blk, worst_m, worst_chain, worst_zeta, time.time() - t0))
# ---------------- Gabcke R0 bound and C0 term (Riemann-Siegel), mpmath ----------------
def testRS(N=600):
    import mpmath as mp
    mp.mp.dps = 30; t0 = time.time(); worstR = 0.0; argR = None; worstC = 0.0
    ts = [10 ** rng.uniform(math.log10(200), math.log10(3e4)) for _ in range(N)]
    for Mj in range(6, 70):   # just below / above the jumps t = 2 pi M^2
        for d in (-1e-6, 1e-6, -0.01, 0.01): ts.append(2 * math.pi * Mj * Mj + d)
    for tf in ts:
        if tf < 200: continue
        t = mp.mpf(tf); a = mp.sqrt(t / (2 * mp.pi)); Nn = int(mp.floor(a)); z = 1 - 2 * (a - Nn)
        th = mp.siegeltheta(t)
        main = 2 * mp.fsum(mp.cos(th - t * mp.log(n)) / mp.sqrt(n) for n in range(1, Nn + 1))
        C0 = mp.cos(mp.pi * (z * z / 2 + mp.mpf(3) / 8)) / mp.cos(mp.pi * z)
        c0t = (-1) ** (Nn - 1) * C0 / mp.sqrt(a)
        R0 = mp.siegelz(t) - main - c0t
        q = float(abs(R0) * t ** 0.75); qc = float(abs(c0t) * t ** 0.25)
        if q > worstR: worstR, argR = q, tf
        worstC = max(worstC, qc)
    log("RS (Gabcke Satz 3.2.2(b)): %d t in [200,3e4] incl. jumps 2 pi M^2: max |R0| t^(3/4) = %.5f (bound 0.127) at t=%.4f; max |C0 term| t^(1/4) = %.6f (bound 1.463) (%.0fs)"
        % (len(ts), worstR, argR, worstC, time.time() - t0))
if __name__ == "__main__":
    log("falsify.py tests=%s ncases=%d prec=%d" % (TESTS, NC, ctx.prec)); T = TESTS.split(',')
    if '28' in T: test28()
    if '24' in T: test24()
    if '25' in T: test25()
    if '29' in T: test29(NC, deadline=300, nscan=NSCAN)
    if '29x' in T: test29(NC, ext=True, deadline=300, nscan=NSCAN)
    if '27' in T: test27(NC, deadline=240, nscan=NSCAN)
    if '27x' in T: test27(NC, ext=True, deadline=240, nscan=NSCAN)
    if '26' in T: test26(NC, deadline=120)
    if 'P' in T: testP()
    if 'b' in T: testb()
    if 'bind' in T: testbind()
    if 'RS' in T: testRS()
    log("DONE")
