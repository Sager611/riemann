#!/usr/bin/env python3
"""BeyondB2 verification (node N7): moving below B_2 = (777 - sqrt 1081)/852 by signed (cross-row) estimates.

Source: B. Liu, Hexagon 2611.00082v1 (number_field_new_v2-1.tex, "l."), and its input NF3 = arXiv 2610.26835v2
(number_fields_v4.tex, "NF3 l.").  Scalar model: Riemann/Deep/Hexagon2/verify.py (99 checks).

  A  dual lever: signed class row sum.  Hypothesis H(sigma): on every amplitude class B (count U^R) the row sum
     of the correction (NF3 (v:cmp:6),(v:cmp:11); 2611 Lemma sc:reference) has
         |sum_{u in B} c_u N_u| <= U^{-sigma R + eps} (#B)^{3/4} (sum_B |N_u|^4)^{1/4},
     i.e. (3R+1)/4 -> (3R+1)/4 - sigma R.  Then E -> E - h sigma R, S -> S - sigma R, and the exact infimum is
     B(sigma) (closed form below; B(0) = B_2; dB/dsigma(0) = -h_2 R_2/(3+4 delta_2)).
  B  direct lever: the boundary depends on B only through l = 11/3 - 4B (the direct exponent of NF3 Prop
     direct-estimate) and the principal box; a saving Z^{-lam} in that Cauchy pairing gives B_2 - lam exactly.
  C  signed large sieve (the one provable signed tool): exact deficits at the contact; it cannot move B_2.
  D  finite model: exact L-values of all Dirichlet characters mod a prime q (Hurwitz zeta + FFT).
Run: uv run --no-project --with numpy --with sympy --with mpmath python verify.py   (single thread, ~30 s, <250 MB)
"""
import math, time, resource, sys
import numpy as np
import sympy as sp
import mpmath

t0 = time.time()
RES = []
def check(name, ok, detail=""):
    RES.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name + (("  [" + detail + "]") if detail else ""))
    sys.stdout.flush()

# ---------------------------------------------------------------- scalar model (2611 sec 6, l.3195-3460)
B, b, d, x, s, lam = sp.symbols('B b delta x sigma lambda', real=True)
z0 = sp.sqrt(1081)
B2 = (777 - z0) / 852; D2 = (61 - z0) / 66; b2 = -sp.Rational(11, 71) + 595 * z0 / 76751
R2 = (281 - z0) / 396
def geom(Bv, bv):
    l = sp.Rational(11, 3) - 4 * Bv
    return dict(l=l, lx=(1 - bv - l) / 2, ly=(1 + bv - l) / 2, h=(1 + bv + 3 * l) / 2,
                T=(1 + bv - l) - (1 + bv + 3 * l) / 2, C0=-sp.Rational(1, 4) + bv / 6 + 5 * l / 4)
G = geom(B, b)
md = (sp.Rational(5, 6) - d * x) / (sp.Rational(5, 2) - d - 2 * d * x)
Rd = 1 - d * (x + 2 * (1 - x) * md)
g = sp.Rational(1, 4) - b / 6 - 5 * G['l'] / 4 - d / 2 - d * x * G['l'] + 3 * G['h'] * (1 - Rd) / 4 + d * x * G['T'] / 4
gs = g + G['h'] * s * Rd            # active margin under H(sigma)

print("== A. Dual lever: exact infimum B(sigma) under H(sigma)")
g12 = sp.together(gs.subs(x, sp.Rational(1, 2)))
a1 = sp.simplify(sp.diff(g12, B)); be = sp.simplify(sp.diff(g12, b)); a0 = sp.simplify(g12 - B * a1 - b * be)
check("A1 g_sigma(x=1/2) is affine in B and b (no B*b term)", sp.simplify(sp.diff(g12, B, 2)) == 0 and
      sp.simplify(sp.diff(g12, b, 2)) == 0 and sp.simplify(sp.diff(g12, B, b)) == 0)
qd = (33 - 36 * s) * d**2 - (61 - 98 * s) * d + (20 - 60 * s)
check("A2 b-coefficient = [(33-36s)d^2-(61-98s)d+(20-60s)]/(24(4d-5)); its smaller root is delta*(sigma)",
      sp.simplify(be - qd / (24 * (4 * d - 5))) == 0)
disc = 964 * s**2 - 1156 * s + 1081
dstar = ((61 - 98 * s) - sp.sqrt(disc)) / (2 * (33 - 36 * s))
Bs = sp.simplify(-a0.subs(d, dstar) / a1.subs(d, dstar))
Bclosed = (103528 * s**2 + 7316 * s * sp.sqrt(disc) + 51710 * s - 3823 * sp.sqrt(disc) - 51257) / \
          (12 * (8720 * s**2 + 640 * s * sp.sqrt(disc) + 4912 * s - 356 * sp.sqrt(disc) - 5179))
okc = all(abs(sp.N(Bs.subs(s, v) - Bclosed.subs(s, v), 40)) < 1e-35 for v in [0, sp.Rational(1, 1000), sp.Rational(1, 100), sp.Rational(3, 100)])
check("A3 closed form B(sigma) = (103528s^2+7316s r+51710s-3823r-51257)/(12(8720s^2+640s r+4912s-356r-5179)), r=sqrt(964s^2-1156s+1081)", okc)
check("A4 B(0) = B_2 and delta*(0) = delta_2 exactly", sp.simplify(sp.radsimp(Bs.subs(s, 0)) - B2) == 0 and
      sp.simplify(dstar.subs(s, 0) - D2) == 0)
dB0 = sp.diff(Bclosed, s).subs(s, 0)
h2 = (1 + b2 + 3 * (sp.Rational(11, 3) - 4 * B2)) / 2
check("A5 dB/dsigma(0) = -h_2 R_2/(3+4 delta_2) = -0.1078420638", abs(sp.N(dB0 + h2 * R2 / (3 + 4 * D2), 30)) < 1e-25,
      str(sp.N(dB0, 12)))
bstar = -(sp.diff(a0 + B * a1, d)) / sp.diff(be, d)
vals = {}
for v in [sp.Rational(1, 1000), sp.Rational(35, 10000), sp.Rational(1, 100), sp.Rational(296, 10000)]:
    vals[v] = (sp.N(Bclosed.subs(s, v), 20), sp.N(dstar.subs(s, v), 15))
print("     B(0.001) = %s  B(0.0035) = %s  B(0.01) = %s  B(0.0296) = %s" % tuple(str(vals[k][0]) for k in vals))
sig873 = sp.nsolve(Bclosed - sp.Rational(873, 1000), s, 0.003)
sig87 = sp.nsolve(Bclosed - sp.Rational(87, 100), s, 0.03)
check("A6 B(sigma)=0.873 at sigma=0.0035178; =0.87 at 0.0295728", abs(sig873 - 0.00351784198) < 1e-9 and abs(sig87 - 0.0295728437) < 1e-9,
      "%s, %s" % (sp.N(sig873, 10), sp.N(sig87, 10)))

# float model for the full-system check (all entries of Prop sc:rational-certificate / Lemma sc:outside)
P1 = 1124661 / 250000; LST = 0.21
def fmodel(Bv, bv, dd, xx, sv):
    l = 11 / 3 - 4 * Bv; lx = (1 - bv - l) / 2; ly = (1 + bv - l) / 2; h = (1 + bv + 3 * l) / 2; T = 2 * ly - h
    C0 = -0.25 + bv / 6 + 5 * l / 4; th0 = dd * xx / P1; A = 5 / 6 - dd
    mdd = (5 / 6 - dd * xx) / (5 / 2 - dd - 2 * dd * xx); Fd = xx + 2 * (1 - xx) * mdd; Rdd = 1 - dd * Fd
    c = 2 / P1; Dc = 3 - (1 + 2 * c) * xx; Pc = (2 - 2 * c * xx) * (1 - xx); J = A * Dc + dd * Pc
    Rs = 1 - dd + A * dd * Pc / (2 * J)
    E = lambda R: C0 + dd / 2 + dd * xx * (l + (h / 2 - ly) / P1) - 3 * h * (1 - R) / 4 - h * sv * R
    S = lambda R: (3 * R + 1) / 4 + th0 / 2 - 17 / 50 - sv * R
    N = lambda R: (3 * R + 1) / 4 - th0 / 2
    Dsh = lambda R: R + (0.5 - th0) * 3 / 40 - N(R)
    JN = 5 / 2 - 5 / P1
    return dict(l=l, ly=ly, h=h, T=T, C0=C0, th0=th0, Rd=Rdd, Rs=Rs, E=E, S=S, Dsh=Dsh,
                GN=dd * xx * JN / 10 * T, GL=th0 * (0.70 * h - ly), GR=(dd / 2 - th0) * (ly - 0.53 * h))
DOM = ((11 / 50, 13 / 20), (6 / 25, 0.5)); D1 = ((1 / 8, 5 / 6), (0, 6 / 25)); D2r = ((1 / 8, 11 / 50), (6 / 25, 0.5))
D3 = ((13 / 20, 5 / 6), (6 / 25, 0.5)); OMG = ((1 / 8, 5 / 6), (0, 0.5))
ENT = [("G_L-E", lambda m: m["GL"] - m["E"](m["Rd"]), DOM), ("G_R-E", lambda m: m["GR"] - m["E"](m["Rd"]), DOM),
       ("S(R_d)-1/5", lambda m: m["S"](m["Rd"]) - 0.2, DOM), ("1/250-E", lambda m: 1 / 250 - m["E"](m["Rd"]), DOM),
       ("(p1-2dx)L-R_d", lambda m: (P1 - 2 * m["th0"] * P1) * LST - m["Rd"], DOM),
       ("-E-hDsh", lambda m: -m["E"](m["Rd"]) - m["h"] * m["Dsh"](m["Rd"]), DOM),
       ("-E(R_s) D1", lambda m: -m["E"](m["Rs"]), D1), ("-E(R_s) D2", lambda m: -m["E"](m["Rs"]), D2r),
       ("-E(R_s) D3", lambda m: -m["E"](m["Rs"]), D3), ("S(R_s)", lambda m: m["S"](m["Rs"]), OMG)]
def gridmin(f, rect, Bv, bv, sv, n=301):
    dd, xx = np.meshgrid(np.linspace(*rect[0], n), np.linspace(*rect[1], n))
    v = f(fmodel(Bv, bv, dd, xx, sv)); i = np.unravel_index(np.argmin(v), v.shape)
    return float(v[i]), float(dd[i]), float(xx[i])
okall = True; rep = []
for sv in [0.001, 0.0035, 0.005, 0.0075]:
    Bv = float(Bclosed.subs(s, sv)); dv = float(dstar.subs(s, sv)); bv = float(bstar.subs({s: sv, B: Bv, d: dv}))
    act, da, xa = gridmin(lambda m: m["GN"] - m["E"](m["Rd"]), DOM, Bv + 1e-9, bv, sv)
    others = [gridmin(f, r, Bv + 1e-9, bv, sv)[0] for _, f, r in ENT]
    Gm = fmodel(Bv, bv, 0.3, 0.3, sv); l, ly, h, C0 = Gm["l"], Gm["ly"], Gm["h"], Gm["C0"]; m0 = 1 - ly / h
    cnames = ["b", "l-1/6", "1/5-l", "lx-l", "ly-l-11b/6", "p1 l-3/4", "T", "1/2-ly", "B-0.87", "l/h-L", "2/9-l/(h-1/40)", "m0-0.30", "0.47-m0", "ly/h-0.57", "0.59-ly/h", "-Qsm", "small-amp", "floor", "small-rows"]
    consts = [bv, l - 1 / 6, 0.2 - l, (1 - bv - l) / 2 - l, ly - l - 11 * bv / 6, P1 * l - 0.75, 2 * ly - h, 0.5 - ly, Bv - 0.87,
              l / h - LST, 2 / 9 - l / (h - 1 / 40), m0 - 0.30, 0.47 - m0, ly / h - 0.57, 0.59 - ly / h,
              -(l * (0.5 + 2 * l / h) - 11 * h / 54),
              -1 / 250 - C0 - (1 - h) / 16 + h * sv * 7 / 8,          # small amplitude: class count R >= 1-delta >= 7/8
              -1 / 1000 - C0 - (0.5 + 1.5 * l) / 50 + h * sv,          # floor bins: trivial count R = 1
              -1 / 20 - 13 * h / 75 + ly / 2 - 1 / 50]
    ok = act > -1e-7 and abs(da - dv) < 3e-3 and xa == 0.5 and min(others) > 0 and min(consts) > 0
    okall &= ok
    badc = [nm for nm, cv in zip(cnames, consts) if cv <= 0]
    rep.append("s=%g: B=%.10f active min %.1e at (%.4f,%.2f), other entries min %.2e, consts min %.2e %s" % (sv, Bv, act, da, xa, min(others), min(consts), badc))
for r_ in rep: print("     " + r_)
check("A7 for sigma<=0.0075 the KKT contact (delta*,1/2) is the only active margin; all other entries of the system positive", okall)
def floor_on_path(sv):
    Bv = float(Bclosed.subs(s, sv)); dv = float(dstar.subs(s, sv)); bv = float(bstar.subs({s: sv, B: Bv, d: dv}))
    l = 11 / 3 - 4 * Bv; h = (1 + bv + 3 * l) / 2; C0 = -0.25 + bv / 6 + 5 * l / 4
    return -1 / 1000 - C0 - (0.5 + 1.5 * l) / 50 + h * sv
lo, hi = 0.001, 0.01
for _ in range(40):
    mid = (lo + hi) / 2
    if floor_on_path(mid) > 0: lo = mid
    else: hi = mid
sfl = lo
check("A8 along the KKT path the floor requirement (sc:floor-small-rows, with the saving h sigma) first binds at sigma_f = %.6f, B(sigma_f) = %.10f" % (sfl, float(Bclosed.subs(s, sfl))),
      0.0075 < sfl < 0.008)
print("     beyond sigma_f the infimum stays within 1.2e-5 of B(sigma) up to sigma = 0.03 and saturates at 13/15 (l <= 1/5) for sigma >= ~0.06 (infimum_scan.py)")

print("== B. Direct lever (NF3 Prop direct-estimate l.8505-8538: Cauchy for the pairing Q'^{-1/2} sum A~ B)")
L = sp.symbols('ell', real=True)
gl = g.subs(sp.Rational(11, 3) - 4 * B, L)
gl = sp.simplify(g.subs(B, (sp.Rational(11, 3) - L) / 4))
check("B1 the active margin is a function of (ell, b, delta, x) only: B enters only via ell = 11/3 - 4B",
      not gl.has(B))
lx_, ly_ = (1 - b - L) / 2, (1 + b - L) / 2
Bdir = sp.simplify(lx_ / 2 + b / 12 + (4 + b) / 6)
check("B2 direct exponent l_x/2 + b/12 = C(B_low) gives B_low = 11/12 - ell/4 for every b", sp.simplify(Bdir - (sp.Rational(11, 12) - L / 4)) == 0)
print("     => with |I_mod| << Z^{l_x/2+b/12-lam}: B_low = 11/12 - ell/4 - lam, same margins, infimum B_2 - lam exactly")
s0 = sp.symbols('s0')
mono = [s0 + 0.95, s0 + 0.99, s0 + 1.94, s0 + 0.94, s0 + 1.89, 6 * s0 - 3.01, 7 * s0 - 3.01, 6 * s0 - 1.11, 6 * s0 - 0.12,
        6 * s0 - 2.01, 6 * s0 - 1.06, sp.Float(1.94)]
sel = [s0, sp.Float(0.95), sp.Float(0.99), sp.Float(0.94), 5 * s0 + 0.99 - 4]
ok = all(float(m.subs(s0, 0.67)) > 1 for m in mono) and all(float(m.subs(s0, 0.67)) > 0 for m in sel) and float(mono[5].subs(s0, 0.668)) < 1
check("B3 principal box (2611 Lemma principal-box): with (Re w, Re z) = (0.95, 0.165) all 12 decay exponents >1 and 5 selected exponents >0 for Re s >= 0.67; binding monomial R: 6Re s - 3.01", ok)
lam87 = sp.N(B2 - sp.Rational(87, 100), 12); lam1315 = sp.N(B2 - sp.Rational(13, 15), 12)
print("     lam = %s gives 0.87; lam = %s gives 13/15 (R109 floor, itself the l<=1/5 branch of the direct estimate)" % (lam87, lam1315))

print("== C. Signed large sieve in the Moebius variable (the only provable signed tool here)")
a2 = (1 + D2) / 2; l2 = sp.Rational(11, 3) - 4 * B2; g2 = D2 / 2
ly2 = (1 + b2 - l2) / 2; m02 = 1 - ly2 / h2; zz2 = (1 - 2 * m02) / 2
Del1 = sp.nsimplify(sp.radsimp(sp.expand(1 - a2 - h2 * R2 / 2)))
Del2 = sp.radsimp(sp.expand(1 - a2 - g2 * l2 + h2 / 4 - 3 * h2 * R2 / 4 + h2 * g2 * zz2 / 2))
check("C1 LS1 (|sum b_u F_u| <= Z ||b||_2, conjectural diagonal sextic LS) loses Z^{1-a-hR/2} at the contact: 1-a_2-h_2R_2/2 = %s > 0" % sp.N(Del1, 12),
      sp.N(Del1, 30) > 0, str(sp.simplify(Del1)))
check("C2 LS2 (LS + M2 for S_m Q_full, no pointwise prime bound) loses Z^{%s} > 0; M2 length m_0+2l/h = %s <= 1" % (sp.N(Del2, 10), sp.N(m02 + 2 * l2 / h2, 6)),
      sp.N(Del2, 30) > 0 and sp.N(m02 + 2 * l2 / h2, 30) <= 1)
h2f = float(sp.N(h2, 20))
ddg, xxg = np.meshgrid(np.linspace(1 / 8, 5 / 6, 401), np.linspace(0, 0.5, 201))
mdg = (5 / 6 - ddg * xxg) / (5 / 2 - ddg - 2 * ddg * xxg); Rg = 1 - ddg * (xxg + 2 * (1 - xxg) * mdg)
gap = ((1 - ddg) - h2f * Rg) / 2
check("C3 at the B_2 geometry LS1 is worse than the absolute route on all of [1/8,5/6]x[0,1/2]: min (1-a-h_2R_d/2) = %.4f > 0" % gap.min(), gap.min() > 0)
mus = sp.N(-(sp.N(Del1, 20)) / (h2 * R2), 8)
print("     in H(sigma) units the LS route is sigma = %s < 0: the min(absolute, LS) system keeps the contact, infimum = B_2" % mus)

print("== D. Finite model: Dirichlet characters mod prime q, exact L-values (Hurwitz zeta + FFT)")
BERN = [1/6, -1/30, 1/42, -1/30, 5/66, -691/2730, 7/6, -3617/510]
def hurwitz(sv, xv, N=12):
    tot = np.zeros_like(xv)
    for n in range(N): tot += (n + xv) ** (-sv)
    y = N + xv; tot += y ** (1 - sv) / (sv - 1) + 0.5 * y ** (-sv)
    fact = 1.0; poch = sv; ypow = y ** (-sv - 1)
    for k, bb in enumerate(BERN, start=1):
        fact *= (2 * k - 1) * (2 * k) if k > 1 else 2
        tot += bb / fact * poch * ypow; poch *= (sv + 2 * k - 1) * (sv + 2 * k); ypow = ypow / (y * y)
    return tot
def proot(q):
    f, m, p = [], q - 1, 2
    while p * p <= m:
        if m % p == 0:
            f.append(p)
            while m % p == 0: m //= p
        p += 1
    if m > 1: f.append(m)
    return next(gg for gg in range(2, q) if all(pow(gg, (q - 1) // p, q) != 1 for p in f))
def isp(n): return n > 1 and all(n % i for i in range(2, int(n ** 0.5) + 1))
err = max(abs(float(hurwitz(sv, np.array([xv]))[0]) / float(mpmath.zeta(sv, xv)) - 1) for sv in [0.45264, 0.75] for xv in [1e-3, 0.37, 1.0])
q = 1009; gq = proot(q); ak = [pow(gq, k, q) for k in range(q - 1)]
Lc = q ** (-0.75) * (q - 1) * np.fft.ifft(hurwitz(0.75, np.array(ak) / q))
chi = [0] * q
for k in range(q - 1): chi[ak[k]] = mpmath.exp(2j * mpmath.pi * 5 * k / (q - 1))
check("D0 Hurwitz (rel err %.0e) and FFT L-values agree with mpmath.dirichlet (q=1009, chi_5, s=0.75)" % err,
      err < 1e-12 and abs(Lc[5] - complex(mpmath.dirichlet(0.75, chi))) < 1e-10)
rng = np.random.default_rng(20261105)
def rows(q, th=0.04736, sa=0.75, eta=7, zP=0.21):
    gq = proot(q); n = q - 1; ak = np.empty(n, dtype=np.int64); a = 1
    for k in range(n): ak[k] = a; a = a * gq % q
    xv = ak / q
    Lw = q ** (-(0.5 - th)) * n * np.fft.ifft(hurwitz(0.5 - th, xv)); Ls = q ** (-sa) * n * np.fft.ifft(hurwitz(sa, xv))
    j = np.arange(n); c = 1.0 / Ls[(eta - j) % n]
    P = q ** zP; ka = np.empty(q, dtype=np.int64); ka[ak] = np.arange(n); arr = np.zeros(n)
    for p in range(int(P) + 1, int(2 * P) + 1):
        if isp(p) and p != q: arr[ka[p]] += math.exp(-1.0 / (p / P - 0.999))
    Q = n * np.fft.ifft(arr / math.sqrt(P))[(-j) % n]
    msk = (j != 0) & (((eta - j) % n) != 0)
    return Lw[msk], Q[msk], c[msk]
summary = []
for q in [65537, 163841]:
    Lw, Q, c = rows(q); T = Lw * Q * c; K = int(round(len(T) ** 0.6266))
    def rr(idx): return abs(T[idx].sum()) / np.abs(T[idx]).sum()
    fits = {}
    for nm, key in [("top|T|", np.abs(T)), ("top|c||Q|", np.abs(c) * np.abs(Q))]:
        o = np.argsort(-key); Ks = [int(round(len(T) ** R)) for R in (0.4, 0.5, 0.6266, 0.7, 0.8)]
        lr = [math.log(rr(o[:k])) for k in Ks]; sl = -np.polyfit([math.log(k) for k in Ks], lr, 1)[0]
        fits[nm] = (rr(o[:K]), sl)
    rnd = float(np.median([rr(rng.choice(len(T), K, replace=False)) for _ in range(100)]))
    o = np.argsort(-(np.abs(c) * np.abs(Q)))[:K]; ph = np.angle(T[o])
    best = max((o[np.abs(np.angle(np.exp(1j * (ph - p0)))) < math.pi / 4] for p0 in np.linspace(-math.pi, math.pi, 73)), key=len)
    summary.append((q, K, fits, rnd, abs(T.mean()) / np.abs(T).mean(), len(best), rr(best)))
    print("     q=%d K=%d: top|T| r=%.3f slope sigma=%.3f; top|c||Q| r=%.3f slope %.3f; random r=%.3f (family-mean ratio %.3f); adversarial quarter-arc subset (%d rows) r=%.3f" %
          (q, K, fits["top|T|"][0], fits["top|T|"][1], fits["top|c||Q|"][0], fits["top|c||Q|"][1], rnd, summary[-1][4], len(best), rr(best)))
check("D1 magnitude-selected classes cancel at a power rate (fitted sigma in [0.1,0.4] >> 0.0035 needed)",
      all(0.1 < f[1] < 0.4 for sm in summary for f in sm[2].values()))
check("D2 adversarial subsets with the same magnitude data exist: quarter-arc subsets keep r > 0.85 (no power saving)",
      all(sm[6] > 0.85 and sm[5] > sm[1] / 6 for sm in summary))
print("== %d checks, %d failed; time %.1f s; maxrss %.0f MB" % (len(RES), RES.count(False), time.time() - t0,
      resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2**20))
