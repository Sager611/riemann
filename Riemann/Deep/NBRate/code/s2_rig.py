# nbrate stage 2: RIGOROUS d_N^2 (Baez-Duarte, L^2(0,inf), unconstrained = V1) for all N <= NMAX.
# G_{jk} = <rho_{1/j}, rho_{1/k}> by Vasyunin's formula with cotangent sums summed EXACTLY in integers
# (cot(pi m/k) as Arb balls -> fixed point 2^-PC, limb products exact in float64), entries fixed point 2^-P.
# Certificate: d_N^2 = Q(c) - r^T G_N^{-1} r for any c (r = G_N c - b_N), lambda_min(G_N) >= lambda_min(G_NMAX)
# (interlacing) >= mu - ||E||_F with G = Lhat Lhat^T + mu I + E (Lhat dyadic, Lhat Lhat^T exact by limbs).
# usage: uv run --no-project --with python-flint --with numpy --with scipy python3 s2_rig.py NMAX
import sys, math, time, gc
import numpy as np
from math import gcd
from fractions import Fraction as Fr
from flint import arb, fmpq, fmpz, ctx
import scipy.linalg as sla
t0 = time.time()
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 1280
OUT = '/Users/adriansager/playground/riemann/.scratch/nbrate/'
P, PC, Q = 120, 140, 60
ctx.prec = 220
u = 2.0 ** -53
def log(*a):
    print(f'[{time.time()-t0:7.1f}s]', *a, flush=True)
def fix(xa, S):  # nearest integer to mid(xa)*2^S
    man, e = xa.mid().man_exp(); man = int(man); e = int(e) + S
    if e >= 0: return man << e
    q = (abs(man) + (1 << (-e - 1))) >> (-e)
    return q if man >= 0 else -q
def radf(xa):  # float upper bound of radius
    man, e = xa.rad().man_exp(); return math.ldexp(float(int(man)) * (1 + 1e-15), int(e)) if int(man) else 0.0
# ---- Part 1: exact cotangent sums S(a,k): V(a,k) = S/(k 2^PC) +- evk[k]
Vs = [None] * (NMAX + 1); evk = [0.0] * (NMAX + 1)
LB = 26; MASK = (1 << LB) - 1
for k in range(2, NMAX + 1):
    h = (k - 1) // 2
    a_arr = np.array([x for x in range(1, k) if gcd(x, k) == 1], dtype=np.int64)
    if h == 0:
        Vs[k] = {int(a): 0 for a in a_arr}; continue
    Cm = []; emax = 0.0
    for m in range(1, h + 1):
        c = arb(fmpq(m, k)).cot_pi()
        Cm.append(fix(c, PC)); emax = max(emax, radf(c))
    emax = emax + 2.0 ** (-PC - 1)
    evk[k] = h * emax * (1 + 1e-12)   # (1/k) sum |2r-k| e_m <= h*emax
    nl = (max(Cm).bit_length() + LB - 1) // LB
    climb = np.array([[(x >> (LB * l)) & MASK for l in range(nl)] for x in Cm], dtype=np.float64)
    m_arr = np.arange(1, h + 1, dtype=np.int64)
    R = (2 * (np.outer(a_arr, m_arr) % k) - k).astype(np.float64)
    M = R @ climb   # exact: |R|<2^11, limbs<2^26, h<=640 terms -> < 2^47
    assert np.abs(M).max() < 2.0 ** 52
    Mi = M.astype(np.int64)
    d = {}
    for i, a in enumerate(a_arr.tolist()):
        s = 0
        for l in range(nl):
            s += int(Mi[i, l]) << (LB * l)
        d[a] = s
    Vs[k] = d
log('cot sums done; max evk', max(evk))
# ---- Part 2: g(a,b) coprime, Gram fixed point Gt = round(G 2^P), limbs base 2^20 (7 limbs), float copy Gf
two = arb(2); tPC = two ** (-PC); tP = two ** P
LG = (2 * arb.pi()).log() - arb.const_euler()
LOGS = [None] + [arb(n).log() for n in range(1, NMAX + 1)]
PI = arb.pi()
GL = 20; NGL = 7; GMASK = (1 << GL) - 1
Gl = np.zeros((NGL, NMAX, NMAX), dtype=np.float64)
Gf = np.zeros((NMAX, NMAX), dtype=np.float64)
epsG_ulps = 0.0   # max over entries of |G 2^P - Gt|
def put(i, j, x):
    for l in range(NGL):
        v = (x >> (GL * l)) & GMASK
        Gl[l, i, j] = v; Gl[l, j, i] = v
    f = math.ldexp(float(x), -P)
    Gf[i, j] = f; Gf[j, i] = f
npairs = 0
for bb in range(1, NMAX + 1):
    for aa in range(1, bb + 1):
        if gcd(aa, bb) != 1: continue
        if aa == bb:  # (1,1)
            g = LG
        else:
            vab = arb(fmpz(Vs[bb][aa])) * tPC / bb + arb(0, evk[bb])
            if aa == 1:
                vba = arb(0)
            else:
                vba = arb(fmpz(Vs[aa][bb % aa])) * tPC / aa + arb(0, evk[aa])
            g = (LG * (aa + bb) + (bb - aa) * (LOGS[aa] - LOGS[bb]) - PI * (vab + vba)) / (2 * aa * bb)
        gt = fix(g, P)
        e1 = radf(g) * 2.0 ** P + 0.5     # |g 2^P - gt| <= e1 ulps
        assert gt >= 0
        npairs += 1
        for dd in range(1, NMAX // bb + 1):
            x = (gt + dd // 2) // dd       # round(gt/dd)
            epsG_ulps = max(epsG_ulps, e1 / dd + 0.5)
            put(dd * aa - 1, dd * bb - 1, x)
del Vs; gc.collect()
epsG = epsG_ulps * 2.0 ** -P * (1 + 1e-12)
assert Gl[NGL - 1].max() < 2 ** GL
log(f'Gram done: {npairs} coprime pairs, epsG = {epsG:.3e}')
# b_k = (log k + 1 - gamma)/k
EGa = arb.const_euler()
bt = []; epsb_ulps = 0.0
for k in range(1, NMAX + 1):
    x = (LOGS[k] + 1 - EGa) / k
    bt.append(fix(x, P)); epsb_ulps = max(epsb_ulps, radf(x) * 2.0 ** P + 0.5)
epsb = epsb_ulps * 2.0 ** -P * (1 + 1e-12)
bf = np.array([math.ldexp(float(x), -P) for x in bt])
Gd = np.load(OUT + 'G_double.npy')
log('max |Gf - G_double(stage1)| =', float(np.abs(Gf - Gd[:NMAX, :NMAX]).max())); del Gd
# ---- Part 3: lambda_min certificate for G_NMAX
ev = np.linalg.eigvalsh(Gf)
lmin_est, lmax_est = ev[0], ev[-1]; del ev
mu = 2.0 ** math.floor(math.log2(lmin_est * 0.6))
Lt = np.linalg.cholesky(Gf - mu * np.eye(NMAX))
Z = np.round(Lt * 2.0 ** 50).astype(np.int64); del Lt
Z0 = (Z & ((1 << 17) - 1)); Z1 = ((Z >> 17) & ((1 << 17) - 1)); Z2 = (Z >> 34)
Zs = [Z0.astype(np.float64), Z1.astype(np.float64), Z2.astype(np.float64)]; del Z0, Z1, Z2
assert all(np.abs(z).max() < 2.0 ** 18 for z in Zs)
S = np.zeros((NMAX, NMAX))
for a in range(3):
    for b2 in range(3):
        Pm = Zs[a] @ Zs[b2].T          # exact: < 2^36 * 1280 < 2^47
        S += Pm * 2.0 ** (17 * (a + b2) - 100)
        del Pm
Y = sum(np.abs(Zs[a]) * 2.0 ** (17 * a - 50) for a in range(3)); del Zs
yn = np.sqrt((Y * Y).sum(axis=1)) * (1 + 1e-12); del Y
Ehat = Gf - mu * np.eye(NMAX) - S; del S
U = np.abs(Ehat) + epsG + 12 * u * (np.abs(Gf) + mu + np.outer(yn, yn))
normE = math.sqrt(float((U * U).sum())) * (1 + 1e-10); del U, Ehat
lam_lb = mu - normE
log(f'lambda_min est {lmin_est:.6e}, lambda_max est {lmax_est:.6e}, mu = {mu:.4e}, ||E||_F <= {normE:.3e}, CERTIFIED lambda_min(G_N) >= {lam_lb:.6e} for all N <= {NMAX}')
assert lam_lb > 0
gc.collect()
# ---- Part 4: per-N certificates
Lc = np.linalg.cholesky(Gf)
yv = sla.solve_triangular(Lc, bf, lower=True)
EGf = 0.57721566490153286061
Carb = 2 + EGa - (4 * arb.pi()).log()
mu_m = np.ones(NMAX + 1, dtype=np.int64); mu_m[0] = 0
for p in range(2, NMAX + 1):
    if all(p % q for q in range(2, int(p ** 0.5) + 1)):
        mu_m[p::p] *= -1; mu_m[p * p::p * p] = 0
def certify(N, c):
    """rigorous [lo, hi] for Q_N(c_hat) and residual bound; c float array length N; returns (Qlo,Qhi,r2,S1)"""
    ct = [int(x) for x in np.round(c * 2.0 ** Q)]
    sg = [1 if x >= 0 else -1 for x in ct]
    nl = max(1, (max(abs(x) for x in ct).bit_length() + 19) // 20)
    Cl = np.array([[sg[i] * ((abs(ct[i]) >> (20 * m)) & GMASK) for m in range(nl)] for i in range(N)], dtype=np.float64)
    Mt = np.stack([Gl[l, :N, :N] @ Cl for l in range(NGL)])   # (NGL, N, nl) exact (< 2^20*2^20*1280 < 2^51)
    assert np.abs(Mt).max() < 2.0 ** 52
    Mi = Mt.astype(np.int64)
    w = np.zeros(N, dtype=object)
    for t in range(NGL + nl - 1):
        Dt = np.zeros(N, dtype=np.int64)
        for l in range(max(0, t - nl + 1), min(NGL, t + 1)):
            Dt += Mi[l, :, t - l]
        w = w + (Dt.astype(object) << (GL * t))
    ctv = np.array(ct, dtype=object)
    cGc = int((ctv * w).sum())                  # scale 2^(P+2Q)
    bc = int(sum(bt[i] * ct[i] for i in range(N)))   # scale 2^(P+Q)
    Qex = Fr((1 << (P + 2 * Q)) - 2 * (bc << Q) + cGc, 1 << (P + 2 * Q))
    S1 = math.ldexp(float(sum(abs(x) for x in ct)), -Q) * (1 + 1e-12)
    errQ = (epsG * S1 * S1 + 2 * epsb * S1) * (1 + 1e-9)
    # residual r_j = w_j/2^(P+Q) - b_j
    rj = [Fr(int(w[j]) - (bt[j] << Q), 1 << (P + Q)) for j in range(N)]
    er = epsG * S1 + epsb
    r2 = sum((abs(float(x)) * (1 + 1e-12) + er) ** 2 for x in rj) * (1 + 1e-9)
    return Qex, errQ, r2, S1, float(np.abs(c).max())
def ball(lo, hi):
    lo = Fr(lo); hi = Fr(hi)
    m = (lo + hi) / 2; r = (hi - lo) / 2
    return arb(fmpq(m.numerator, m.denominator)) + arb(0, float(r) * (1 + 1e-9) + 1e-300)
res = {}
sel = set([10, 20, 40, 41, 48, 49, 80, 160, 205, 206, 320, 411, 524, 583, 584, 640, 695, 1280])
rows = []
signs = np.zeros(NMAX + 1, dtype=np.int64)
lo_arr = np.zeros(NMAX + 1); hi_arr = np.zeros(NMAX + 1)
bcf_lo = np.zeros(NMAX + 1); bcf_hi = np.zeros(NMAX + 1)
for N in range(1, NMAX + 1):
    c = sla.solve_triangular(Lc[:N, :N], yv[:N], lower=True, trans='T')
    Qex, errQ, r2, S1, cmax = certify(N, c)
    corr = r2 / lam_lb
    lo = Qex - Fr(errQ) - Fr(corr); hi = Qex + Fr(errQ)
    lo_arr[N] = math.nextafter(float(lo), -math.inf); hi_arr[N] = math.nextafter(float(hi), math.inf)
    f = ball(lo, hi) * LOGS[N] - Carb
    if f < 0: signs[N] = -1
    elif f > 0: signs[N] = 1
    # BCF choice c_k = -mu(k)(1 - log k/log N), N >= 2
    if N >= 2:
        cb = -mu_m[1:N + 1] * (1 - np.log(np.arange(1, N + 1)) / math.log(N))
        Qb, eQb, r2b, S1b, _ = certify(N, cb)
        # |Q(c) - Q(c_hat)| <= 2 ||G c_hat - b|| ||delta|| + lmax ||delta||^2, ||delta|| <= sqrt(N) 2^-(Q+1)
        dl = math.sqrt(N) * 2.0 ** (-Q - 1)
        eb = eQb + 2 * math.sqrt(r2b) * dl + 6.0 * dl * dl
        bcf_lo[N] = math.nextafter(float(Qb) - eb, -math.inf); bcf_hi[N] = math.nextafter(float(Qb) + eb, math.inf)
    if N in sel:
        rows.append((N, ball(lo, hi).str(12, radius=True), (ball(lo, hi) * LOGS[N]).str(10, radius=True), f.str(4, radius=True), float(hi - lo), cmax, S1, corr))
        log(f'N={N}: d^2 in {ball(lo,hi).str(14, radius=True)}  d^2 logN - C = {f.str(5, radius=True)}  width {float(hi-lo):.2e} |c|max {cmax:.3g} sum|c| {S1:.3g} resid-corr {corr:.2e}')
np.savez(OUT + 'rig_V1.npz', lo=lo_arr, hi=hi_arr, signs=signs, bcf_lo=bcf_lo, bcf_hi=bcf_hi, lam_lb=lam_lb, epsG=epsG, mu=mu, normE=normE)
neg = [N for N in range(2, NMAX + 1) if signs[N] < 0]
und = [N for N in range(2, NMAX + 1) if signs[N] == 0]
print('UNDECIDED N:', und)
first = neg[0] if neg else None
print('first N>=2 with certified d_N^2 log N < C:', first, '; all N in [2, first) certified > C:', all(signs[N] > 0 for N in range(2, first)) if first else None)
runs = []; N = 2
while N <= NMAX:
    s = signs[N]; M2 = N
    while M2 + 1 <= NMAX and signs[M2 + 1] == s: M2 += 1
    runs.append((int(s), N, M2)); N = M2 + 1
print('sign runs (sign, from, to):', runs)
print('max width hi-lo over N:', float((hi_arr[1:] - lo_arr[1:]).max()))
print('TABLE')
for r in rows: print(r)
log('done')
