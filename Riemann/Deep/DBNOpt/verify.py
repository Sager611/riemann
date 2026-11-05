# /// script
# dependencies = ["python-flint==0.9.0"]
# ///
"""Independent check of a re-optimised Hakimoglu-Brown table (arXiv 2610.12803, Theorem derivative-step, p = 61/20).

Usage (single-threaded, from this directory):
  uv run --no-project --with python-flint==0.9.0 python verify.py table_theta.csv 5458637291332547/6250000000000000
  uv run --no-project --with python-flint==0.9.0 python verify.py table_78.csv 7/8
Runtime (measured 2026-11-05, one core): 124 s for table_theta.csv (474 rows), 128 s for table_78.csv, 9 s for
table_theta_h200.csv; peak RSS 37 MB. The time is the Arb sums to N1 = 65536 per row.

What is checked (everything the theorem needs from the table; the analysis is HB's, audited in ../HBDBN):
  G. table-independent scalars: L/4 > 5.5716 and N >= 69098 for x >= X (exp bound E160bar, pi < 355/113),
     J(2) < 1, log N1 < 12, int_{N1}^oo log y y^(-81/40) dy < 1e-3 (Arb).
  S. self-tests: HB's row 0 scalar table entries and HB's least Bernstein coefficient (row 21) are reproduced.
  R. per row i: the nine printed U-values are strict upper bounds (Arb, 128 bits) of Z_*, psi_*, psi2_*,
     delta_0..2 (t^- = t_i, t^+ = t_i + d_i) and rho_0..2; exact rationals give eps0 < 1/2, C_i <= C^cert,
     G_i >= G^cert; all 400 Bernstein coefficients (16 U-pieces) of Psi_i are > 0 (so Phi_i >= 0 on D >= 0,
     0 <= r <= q_i); k_i > 0 exactly; gate E_20(k_i d_i)(2 + k_i q_{i+1}) >= 2 + k_i q_i exactly
     (E_20 = degree-20 Taylor sum <= exp); chain q_0 >= (2 theta - 1)^2, q_0 < B, q strictly decreasing,
     q_n = 1/25, t_0 = 0, t_{i+1} = t_i + d_i, t_n <= 1/5, all C, G, lambda, mu >= 0.
  Conclusion printed: Lambda <= t_n + q_n/2.
"""
import sys, time, math, resource
from math import comb
from flint import arb, ctx, fmpq, fmpz, fmpq_mpoly_ctx
ctx.prec = 128
T0 = time.time()
FAIL = []
def check(name, cond):
    if not cond: FAIL.append(name); print("FAIL", name, flush=True)
def Q(s):
    if '/' in s: a, b = s.split('/'); return fmpq(int(a), int(b))
    return fmpq(int(s), 1)

# ---------------- G: table-independent scalars ----------------
def E(x, n):  # sum_{j<=n} x^j/j!  (lower bound of exp for x >= 0)
    s = fmpq(0); term = fmpq(1)
    for j in range(n + 1):
        if j > 0: term = term * x / j
        s += term
    return s
def Ebar(x, n):  # upper bound of exp(x) for 0 <= x < n+2
    assert x < n + 2
    return E(x, n) + x ** (n + 1) / math.factorial(n + 1) / (1 - x / (n + 2))
X = fmpq(600000839515, 10)
check("G1 Ebar160(13929/625) < 113X/1420 (L/4 > 5.5716)", Ebar(fmpq(13929, 625), 160) < 113 * X / 1420)
check("G2 69098^2 < 113X/1420 (N >= 69098)", fmpq(69098 ** 2) < 113 * X / 1420)
J2 = sum(arb(n).log() / (n * n) for n in range(2, 17)) + (arb(16).log() + 1) / 16
check("G3 J(2) < 1", J2 < 1)
N1 = 65536
check("G4 log N1 < 12", arb(N1).log() < 12)
v0 = arb(81) / 40
check("G5 int_{N1}^oo log y y^(-81/40) < 1e-3", arb(N1) ** (1 - v0) * (arb(N1).log() / (v0 - 1) + 1 / (v0 - 1) ** 2) < arb(10) ** -3)

# ---------------- Arb scalar quantities ----------------
logs = [None, arb(0)] + [arb(n).log() for n in range(2, N1 + 1)]
sv = bytearray([1]) * (N1 + 1); sv[0] = sv[1] = 0
for i in range(2, int(N1 ** 0.5) + 1):
    if sv[i]: sv[i * i::i] = bytearray(len(sv[i * i::i]))
ppow = []
for l in range(2, N1 + 1):
    if sv[l]:
        m, j = l, 1
        while m <= N1: ppow.append((m, j, l)); m *= l; j += 1
p_a = arb(61) / 20; sigp = arb(81) / 40; Nst = arb(69098); Lst = 2 * Nst.log() + 3; logN1 = logs[N1]
def Tj(j, w):
    return arb(N1) ** (1 - w) * sum((arb(math.factorial(j) // math.factorial(j - k)) * logN1 ** (j - k) / (w - 1) ** (k + 1) for k in range(j + 1)), arb(0))
def scal(tm, tp):
    u = sigp + arb(fmpq(557, 100) * tm); u8 = sigp + arb(fmpq(278, 100) * tm); e = (p_a - 1) / 2 - arb(fmpq(278, 100) * tm)
    tq = arb(tp) / 4
    Z = arb(1); d0 = arb(0); d1 = arb(0); d2 = arb(0)
    for n in range(2, N1 + 1):
        ln = logs[n]; w = (-u * ln).exp(); Z += w
        x = (tq * ln * ln).expm1() * w; d0 += x; x = x * ln; d1 += x; d2 += x * ln
    ps = arb(0); ps2 = arb(0)
    for (n, j, l) in ppow:
        w = (-u * logs[n]).exp(); ps += logs[l] * w; ps2 += j * logs[l] * logs[l] * w
    r0 = arb(101) / 100 * Nst ** (-p_a) * (Nst ** (e + 1) / (e + 1) + Nst ** e)
    return [Z + Tj(0, u), ps + Tj(1, u), ps2 + Tj(2, u), d0 + Tj(0, u8), d1 + Tj(1, u8), d2 + Tj(2, u8), r0, Lst * r0, (Lst * Lst + arb(10) ** -9) * r0]
UDEN = 10 ** 12

# ---------------- exact C^cert, G^cert (HB eqs new-errors, new-C, new-G) ----------------
c1 = fmpq(1, 2) + fmpq(1, 10 ** 10); c2 = fmpq(1, 4) + fmpq(1, 10 ** 10)
def CGcert(U):
    Z, ps, ps2, d0, d1, d2, r0, r1, r2 = [fmpq(x, UDEN) for x in U]
    e0 = d0 + r0 + fmpq(1, 10 ** 6); ep0 = Z * e0
    e1 = c1 * d1 + r1 + fmpq(4, 10 ** 5); ep1 = Z * (e1 + c1 * ps * e0)
    e2 = c2 * d2 + r2 + fmpq(16, 10 ** 4) + fmpq(1, 10 ** 9)
    ep2 = Z * (e2 + 2 * c1 * ps * e1 + c2 * (ps2 + ps * ps) * e0)
    C = fmpq(55716, 10 ** 4) - fmpq(1, 10 ** 6) - c1 * ps - ep1 / (1 - ep0)
    G = c2 * ps2 + ep2 / (1 - ep0) + (ep1 / (1 - ep0)) ** 2 + fmpq(2, 10 ** 10)
    return C, G, ep0

# ---------------- Bernstein certificate (flint multivariate polynomials) ----------------
cDr = fmpq_mpoly_ctx.get(('D', 'r'), 'lex'); D, r = cDr.gens()
cUV = fmpq_mpoly_ctx.get(('U', 'V'), 'lex'); U_, V_ = cUV.gens()
cuv = fmpq_mpoly_ctx.get(('u', 'v'), 'lex'); u_, v_ = cuv.gens()
Bq = fmpq(3721, 400); pq = fmpq(61, 20)
MU = 16
def bern_min(A2, lam, mu):
    w = D + Bq - r; Dl = w * w + 4 * D * r
    Nc = 2 * lam * w * Dl - mu * (2 * w + 4 * Bq) * Dl + 8 * mu * Bq * w * w
    Ph = 2 * Dl * Dl - (D + A2 + 3 * r) * Nc
    Ps = cUV.from_dict({})
    for (i, j), c in Ph.to_dict().items():
        assert i <= 4 and j <= 4
        Ps += c * Bq ** (i - 4) * A2 ** j * U_ ** i * (1 - U_) ** (4 - i) * V_ ** j
    mn = None; cnt = 0
    for l in range(MU):
        a = Ps.compose((u_ + l) / MU, v_).to_dict()
        for aa in range(5):
            for bb in range(5):
                beta = fmpq(0)
                for k in range(aa + 1):
                    for m in range(bb + 1):
                        if (k, m) in a: beta += a[(k, m)] * fmpq(comb(aa, k), comb(4, k)) * fmpq(comb(bb, m), comb(4, m))
                cnt += 1
                if mn is None or beta < mn: mn = beta
    return mn, cnt

# ---------------- S: self-tests against HB's printed values ----------------
hb0 = [1622101766394, 548568036655, 839050132415, 2333875727, 6985223600, 28309848180, 5463325, 138148683, 3493305246]
s0 = scal(fmpq(0), fmpq(1, 200))
check("S1 HB row-0 entries are strict upper bounds within 1e-11", all(arb(Uh) > s * UDEN and arb(Uh) - s * UDEN < 10 for Uh, s in zip(hb0, s0)))
save = MU; MU = 8
mn, _ = bern_min(fmpq(78037, 10 ** 6), fmpq(1910, 1000), fmpq(1318, 1000))
check("S2 HB least Bernstein coefficient (row 21) reproduced", mn == fmpq(623287907570432169922330885087, 479268282493202500000000000000000))
MU = save

# ---------------- R: the table ----------------
path = sys.argv[1]; theta = Q(sys.argv[2])
rows = [l.strip().split(',') for l in open(path) if l.strip() and not l.startswith('#')]
b = fmpq(1, 25); t = fmpq(0); prev_qn = None
mins = dict(Cm=None, Gm=None, gate=None, beta=None); maxeps = fmpq(0); kmin = kmax = None; nbern = 0
for idx, row in enumerate(rows):
    i = int(row[0]); ti = Q(row[1]); di = Q(row[2])
    qi = fmpq(int(row[3]), 10 ** 10); qn = fmpq(int(row[4]), 10 ** 10)
    Ci = fmpq(int(row[5]), 10 ** 8); Gi = fmpq(int(row[6]), 10 ** 8)
    li = fmpq(int(row[7]), 10 ** 6); mi = fmpq(int(row[8]), 10 ** 6)
    Uv = [int(x) for x in row[9:18]]
    check("R%d index" % idx, i == idx)
    check("R%d time chain" % i, ti == t and di > 0)
    if idx == 0:
        check("R0 q0 >= (2 theta - 1)^2", qi >= (2 * theta - 1) ** 2)
        check("R0 q0 < B", qi < Bq)
    else:
        check("R%d q chain" % i, qi == prev_qn)
    check("R%d q decreasing, >= b" % i, qn < qi and qn >= b)
    check("R%d nonneg C,G,lam,mu" % i, Ci >= 0 and Gi >= 0 and li >= 0 and mi >= 0)
    sc = scal(ti, ti + di)
    check("R%d U strict upper bounds" % i, all(arb(Uk) > s * UDEN for Uk, s in zip(Uv, sc)))
    Cc, Gc, ep0 = CGcert(Uv)
    check("R%d eps0 < 1/2" % i, ep0 < fmpq(1, 2)); maxeps = max(maxeps, ep0)
    check("R%d C_i <= C^cert" % i, Ci <= Cc); check("R%d G_i >= G^cert" % i, Gi >= Gc)
    mins['Cm'] = Cc - Ci if mins['Cm'] is None else min(mins['Cm'], Cc - Ci)
    mins['Gm'] = Gi - Gc if mins['Gm'] is None else min(mins['Gm'], Gi - Gc)
    mn, cnt = bern_min(qi, li, mi); nbern += cnt
    check("R%d Bernstein > 0" % i, mn > 0)
    mins['beta'] = mn if mins['beta'] is None else min(mins['beta'], mn)
    k = 4 * (li * Ci / pq - mi * Gi - 2 * li / (Bq - qi) - 2 * mi * (Bq + qi) / (Bq - qi) ** 2)
    check("R%d k > 0" % i, k > 0)
    kmin = k if kmin is None else min(kmin, k); kmax = k if kmax is None else max(kmax, k)
    g = E(k * di, 20) * (2 + k * qn) - (2 + k * qi)
    check("R%d gate" % i, g >= 0)
    mins['gate'] = g if mins['gate'] is None else min(mins['gate'], g)
    t = ti + di; prev_qn = qn
    if idx % 50 == 0: print("row %d ok so far (%d fails)  %.0fs" % (i, len(FAIL), time.time() - T0), flush=True)
check("final q_n = 1/25", prev_qn == b)
check("t_n <= 1/5", t <= fmpq(1, 5))
bound = t + prev_qn / 2
def fl(x):  # fmpq -> float without overflow
    pz, qz = int(x.p), int(x.q); sh = max(0, max(pz.bit_length(), qz.bit_length()) - 900)
    return (pz >> sh) / (qz >> sh) if sh else pz / qz
print("rows: %d   Bernstein coefficients checked: %d" % (len(rows), nbern))
print("max eps0 = %.6f   min(C^cert-C) = %.3e   min(G-G^cert) = %.3e   min gate = %.3e   min Bernstein = %.3e   k in [%.5f, %.5f]" % (
    fl(maxeps), fl(mins['Cm']), fl(mins['Gm']), fl(mins['gate']), fl(mins['beta']), fl(kmin), fl(kmax)))
print("theta =", theta, " q0 =", rows[0][3] + "/10^10", " t_n =", t, " = %.10f" % fl(t))
print("peak RSS %.0f MB, runtime %.1f s" % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20, time.time() - T0))
if FAIL: print("RESULT: FAIL (%d)" % len(FAIL), FAIL[:10])
else: print("RESULT: ALL CHECKS PASS.  Lambda <= t_n + q_n/2 =", bound, "= %.8f" % fl(bound))
