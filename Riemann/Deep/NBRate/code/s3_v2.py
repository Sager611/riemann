# nbrate stage 3 (lean memory): RIGOROUS constrained distance (V2 = 2611.03924's d_n(1/2)^2: sum c_k/k = 0)
# and unconstrained V1 recheck at selected n.  Same Gram pipeline as s2_rig.py, but limbs stored as int32,
# cotangent sums as int64 limb arrays, lambda certificate in row blocks.
# usage: uv run --no-project --with python-flint --with numpy --with scipy python3 s3_v2.py NMAX
import sys, math, time, gc
import numpy as np
from math import gcd
from fractions import Fraction as Fr
from flint import arb, fmpq, fmpz, ctx
import scipy.linalg as sla
t0 = time.time()
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1280
TARGETS = [n for n in [40, 80, 160, 320, 640, 1280] if n <= NMAX]
P, PC, Q = 120, 140, 60
ctx.prec = 220
u = 2.0 ** -53
def log(*a): print(f'[{time.time()-t0:7.1f}s]', *a, flush=True)
def fix(xa, S):
    man, e = xa.mid().man_exp(); man = int(man); e = int(e) + S
    if e >= 0: return man << e
    q = (abs(man) + (1 << (-e - 1))) >> (-e)
    return q if man >= 0 else -q
def radf(xa):
    man, e = xa.rad().man_exp(); return math.ldexp(float(int(man)) * (1 + 1e-15), int(e)) if int(man) else 0.0
LB = 26; MASK = (1 << LB) - 1
Vl = [None] * (NMAX + 1); Vpos = [None] * (NMAX + 1); evk = [0.0] * (NMAX + 1)
for k in range(2, NMAX + 1):
    h = (k - 1) // 2
    a_arr = np.array([x for x in range(1, k) if gcd(x, k) == 1], dtype=np.int64)
    pos = np.full(k, -1, dtype=np.int32); pos[a_arr] = np.arange(len(a_arr), dtype=np.int32); Vpos[k] = pos
    if h == 0:
        Vl[k] = np.zeros((len(a_arr), 1), dtype=np.int64); continue
    Cm = []; emax = 0.0
    for m in range(1, h + 1):
        c = arb(fmpq(m, k)).cot_pi(); Cm.append(fix(c, PC)); emax = max(emax, radf(c))
    emax += 2.0 ** (-PC - 1); evk[k] = h * emax * (1 + 1e-12)
    nl = (max(Cm).bit_length() + LB - 1) // LB
    climb = np.array([[(x >> (LB * l)) & MASK for l in range(nl)] for x in Cm], dtype=np.float64)
    m_arr = np.arange(1, h + 1, dtype=np.int64)
    R = (2 * (np.outer(a_arr, m_arr) % k) - k).astype(np.float64)
    M = R @ climb
    assert np.abs(M).max() < 2.0 ** 52
    Vl[k] = M.astype(np.int64)
def Sval(a, k):
    row = Vl[k][Vpos[k][a]]
    return sum(int(row[l]) << (LB * l) for l in range(len(row)))
log('cot sums done')
two = arb(2); tPC = two ** (-PC)
LG = (2 * arb.pi()).log() - arb.const_euler()
LOGS = [None] + [arb(n).log() for n in range(1, NMAX + 1)]
PI = arb.pi()
GL = 20; NGL = 7; GMASK = (1 << GL) - 1
Gl = np.zeros((NGL, NMAX, NMAX), dtype=np.int32)
Gf = np.zeros((NMAX, NMAX), dtype=np.float64)
epsG_ulps = 0.0
for bb in range(1, NMAX + 1):
    for aa in range(1, bb + 1):
        if gcd(aa, bb) != 1: continue
        if aa == bb: g = LG
        else:
            vab = arb(fmpz(Sval(aa, bb))) * tPC / bb + arb(0, evk[bb])
            vba = arb(0) if aa == 1 else arb(fmpz(Sval(bb % aa, aa))) * tPC / aa + arb(0, evk[aa])
            g = (LG * (aa + bb) + (bb - aa) * (LOGS[aa] - LOGS[bb]) - PI * (vab + vba)) / (2 * aa * bb)
        gt = fix(g, P); e1 = radf(g) * 2.0 ** P + 0.5
        assert gt >= 0
        for dd in range(1, NMAX // bb + 1):
            x = (gt + dd // 2) // dd
            epsG_ulps = max(epsG_ulps, e1 / dd + 0.5)
            i, j = dd * aa - 1, dd * bb - 1
            for l in range(NGL):
                v = (x >> (GL * l)) & GMASK; Gl[l, i, j] = v; Gl[l, j, i] = v
            f = math.ldexp(float(x), -P); Gf[i, j] = f; Gf[j, i] = f
del Vl, Vpos; gc.collect()
epsG = epsG_ulps * 2.0 ** -P * (1 + 1e-12)
log(f'Gram done, epsG = {epsG:.3e}')
EGa = arb.const_euler()
bt = []; epsb_ulps = 0.0
for k in range(1, NMAX + 1):
    x = (LOGS[k] + 1 - EGa) / k; bt.append(fix(x, P)); epsb_ulps = max(epsb_ulps, radf(x) * 2.0 ** P + 0.5)
epsb = epsb_ulps * 2.0 ** -P * (1 + 1e-12)
bf = np.array([math.ldexp(float(x), -P) for x in bt])
# lambda_min certificate, row blocks
ev = np.linalg.eigvalsh(Gf); lmin_est = ev[0]; lmax_est = ev[-1]; del ev; gc.collect()
mu = 2.0 ** math.floor(math.log2(lmin_est * 0.6))
Lt = np.linalg.cholesky(Gf - mu * np.eye(NMAX))
Z = np.round(Lt * 2.0 ** 50).astype(np.int64); del Lt; gc.collect()
def limb(Zm, a):
    return (((Zm >> (17 * a)) & ((1 << 17) - 1)) if a < 2 else (Zm >> 34)).astype(np.float64)
Y = sum(np.abs(limb(Z, a)) * 2.0 ** (17 * a - 50) for a in range(3))
yn = np.sqrt((Y * Y).sum(axis=1)) * (1 + 1e-12); del Y; gc.collect()
sumU2 = 0.0; BS = 160
for i0 in range(0, NMAX, BS):
    I = slice(i0, min(NMAX, i0 + BS))
    S = np.zeros((I.stop - I.start, NMAX))
    for a in range(3):
        Aa = limb(Z[I], a)
        for b2 in range(3):
            Bb = limb(Z, b2)
            Pm = Aa @ Bb.T
            assert np.abs(Pm).max() < 2.0 ** 52
            S += Pm * 2.0 ** (17 * (a + b2) - 100); del Bb, Pm
    Gb = Gf[I]
    Eh = Gb - S; Eh[np.arange(I.stop - I.start), np.arange(I.start, I.stop)] -= mu
    U = np.abs(Eh) + epsG + 12 * u * (np.abs(Gb) + mu + np.outer(yn[I], yn))
    sumU2 += float((U * U).sum()); del S, Eh, U
normE = math.sqrt(sumU2) * (1 + 1e-10); lam_lb = mu - normE; del Z; gc.collect()
log(f'lambda_min(G_{NMAX}) >= {lam_lb:.6e} (mu {mu:.4e}, ||E||_F <= {normE:.3e}, est {lmin_est:.6e}, lmax est {lmax_est:.4e})')
assert lam_lb > 0
def gc_exact(N, ctint):
    """w = Gt[:N,:N] @ ctint exactly (python ints, scale 2^(P+Q))"""
    sg = [1 if x >= 0 else -1 for x in ctint]
    nl = max(1, (max(abs(x) for x in ctint).bit_length() + 19) // 20)
    Cl = np.array([[sg[i] * ((abs(ctint[i]) >> (20 * m)) & GMASK) for m in range(nl)] for i in range(N)], dtype=np.float64)
    w = np.zeros(N, dtype=object)
    Mt = [Gl[l, :N, :N].astype(np.float64) @ Cl for l in range(NGL)]
    for t in range(NGL + nl - 1):
        Dt = np.zeros(N, dtype=np.int64)
        for l in range(max(0, t - nl + 1), min(NGL, t + 1)):
            assert np.abs(Mt[l][:, t - l]).max() < 2.0 ** 52
            Dt += Mt[l][:, t - l].astype(np.int64)
        w = w + (Dt.astype(object) << (GL * t))
    return [int(x) for x in w]
def Gt_entry(i, j):
    return sum(int(Gl[l, i, j]) << (GL * l) for l in range(NGL))
Carb = 2 + EGa - (4 * arb.pi()).log()
def ball(lo, hi):
    m = (lo + hi) / 2; r = (hi - lo) / 2
    return arb(fmpq(m.numerator, m.denominator)) + arb(0, float(r) * (1 + 1e-9) + 1e-300)
ks = np.arange(1, NMAX + 1)
SC = 1 << (P + Q)
for n in TARGETS:
    # ---- V2: constrained, basis phi_k = rho_{1/k} - rho_1/k, k = 2..n
    G = Gf[:n, :n]
    A = G[1:, 1:] - np.outer(G[1:, 0], 1 / ks[1:n]) - np.outer(1 / ks[1:n], G[0, 1:]) + G[0, 0] * np.outer(1 / ks[1:n], 1 / ks[1:n])
    b2 = np.log(ks[1:n]) / ks[1:n]
    cp = np.linalg.solve(A, b2)
    ct = [0] + [int(x) for x in np.round(cp * 2.0 ** Q)]           # c'_k at scale 2^Q, slot 1 = 0
    c1 = -sum(Fr(ct[k - 1], k) for k in range(2, n + 1)) / (1 << Q)   # exact rational
    w = gc_exact(n, ct)                                              # (G c_hat) scale 2^(P+Q)
    G1 = [Gt_entry(j, 0) for j in range(n)]                          # G_{j1} scale 2^P
    # exact quantities (as Fractions, true value up to entry errors)
    Gc = [Fr(w[j], SC) + Fr(G1[j], 1 << P) * c1 for j in range(n)]
    cvec = [c1] + [Fr(ct[k - 1], 1 << Q) for k in range(2, n + 1)]
    cGc = sum(cvec[j] * Gc[j] for j in range(n))
    bc = sum(Fr(bt[j], 1 << P) * cvec[j] for j in range(n))
    Qv = 1 - 2 * bc + cGc
    S1 = float(sum(abs(x) for x in cvec)) * (1 + 1e-12)
    errQ = (epsG * S1 * S1 + 2 * epsb * S1) * (1 + 1e-9)
    r = [Gc[j] - Fr(bt[j], 1 << P) for j in range(n)]
    er = epsG * S1 + epsb
    rp2 = sum((abs(float(r[k - 1] - r[0] / k)) * (1 + 1e-12) + er * (1 + 1 / k)) ** 2 for k in range(2, n + 1)) * (1 + 1e-9)
    lo = Qv - Fr(errQ) - Fr(rp2 / lam_lb); hi = Qv + Fr(errQ)
    B = ball(lo, hi)
    # sanity: constraint exactly satisfied by construction; double value
    d2dbl = 1 - b2 @ cp
    log(f'V2 n={n}: d^2 in {B.str(16, radius=True)} (width {float(hi-lo):.1e}); d^2 log n = {(B*LOGS[n]).str(10, radius=True)}; minus C = {(B*LOGS[n]-Carb).str(5, radius=True)}; double {d2dbl:.12f}; sum|c| {S1:.3g}')
log('done')
