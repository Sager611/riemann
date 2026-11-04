# arbcheck.py -- independent rigorous re-verification (Arb ball arithmetic via python-flint), written from the TeX.
# Differences from the authors' verifier: direct part via exact power series of P*R (no quadrature);
# transformed part via Arb's rigorous adaptive integration acb.integral (not the Fejer rule + analytic lemma);
# finite inverses via arb_mat.inv (not the supplied W matrices).
# usage: python arbcheck.py scalar | residual | matrix | bernstein
import sys, re, time, json, os
from fractions import Fraction as Fr
from math import comb
from flint import arb, acb, arb_mat, arb_series, acb_series, ctx
ctx.prec = 128
ctx.cap = 64   # series truncation length (python-flint default cap is 10)
import os as _os
_ROOT = _os.path.abspath(_os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', '..'))
SRC = _os.environ.get('CIRC2D_SRC', _os.path.join(_ROOT, '.scratch', 'circ2d', 'src'))
SCR = _os.environ.get('CIRC2D_SCRATCH', _os.path.join(_ROOT, '.scratch', 'circ2d'))
pi = arb.pi(); I = acb(0, 1)
b = arb(3).sqrt() / 2; h = arb(2) / 5; B = arb(4) / 3
A = [0, 1, 3, 4, 7, 9]
def R(x): return arb(Fr(x).numerator) / arb(Fr(x).denominator)   # exact decimal -> ball
OK = []; BAD = []
def check(name, cond, info=''):
    (OK if cond else BAD).append(name)
    print(('PASS ' if cond else 'FAIL ') + name + ('  ' + info if info else ''), flush=True)
# ---- P_l by expanding the product (own) ----
poly = {0: acb(1)}
for a in A:
    f = {0: acb(2), 1: -(-I * pi * a / 6).exp(), -1: -(I * pi * a / 6).exp()}
    new = {}
    for k1, v1 in poly.items():
        for k2, v2 in f.items():
            new[k1 + k2] = new.get(k1 + k2, acb(0)) + v1 * v2
    poly = new
P = {l: poly[l] for l in range(-6, 7)}
ts = [arb(j) / 6 for j in range(7)]
def nodes(nmax): return [n for n in range(1, nmax + 1) if n % 12 in A]
def Pder(s, k):
    return sum((P[l] * (I * pi * ts[abs(l)] * (1 if l >= 0 else -1)) ** k * (I * pi * arb(l) / 6 * s).exp() for l in range(-6, 7)), acb(0)).real
QD = {}
for a in A:
    P2 = Pder(arb(a), 2); P3 = Pder(arb(a), 3)
    QD[a] = (P2 / 2, P3 / (3 * P2))
def Qn(n): return QD[n % 12][0]
def Dn(n): return QD[n % 12][1]
# ---- table ----
tab = {}
for line in open(SRC + '/build/data/coefficients.tex'):
    m = re.findall(r'-?\d+', line)
    if len(m) == 5:
        tab[int(m[0])] = [Fr(int(x), 10 ** 10) for x in m[1:]]
LIST = {1: {0: (Fr(1), Fr(44, 100))}, 2: {0: (Fr(0), Fr(-368, 1000))}}
for n, v in tab.items():
    LIST[1][n] = (v[0], v[1]); LIST[2][n] = (v[2], v[3])
CONST = {1: Fr(-13, 1000), 2: Fr(17, 1000)}
SIG = {1: -1, 2: 1}
def lam(t): return I / (b * (t + I * h))
def zf(t): return -B / (t + I * h) - I * h
# ---- spectral densities on piece j (own derivation, matches the paper's folded formula) ----
def alpha_beta(n, j):
    al = acb(0); be = acb(0)
    for l in range(j, 7):
        e = P[l] * (I * pi * ts[l] * n).exp()
        al += e; be += ts[l] * e
    return al, be
def piece_coefs(L, j):
    out = []
    for n, (c, d) in sorted(L.items()):
        if c == 0 and d == 0: continue
        al, be = alpha_beta(n, j)
        out.append((n, -2 * pi ** 2 * R(c) * be + 2 * I * pi * R(d) * al, 2 * pi ** 2 * R(c) * al))
    return out
def dens_eval(coefs, t):
    w = (-I * pi * t).exp()
    tot = acb(0); wp = acb(1); last = 0
    for n, A0, A1 in coefs:
        wp = wp * w ** (n - last); last = n
        tot += wp * (A0 + A1 * t)
    return tot
FACT = [arb(1)]
for q in range(1, 40): FACT.append(FACT[-1] * q)
NEVAL = [0]
def K_moment(coefs_by_piece, C, m, l):
    """X_l(m) = K^{(l)}(m)/l! rigorously (Re of folded integral + atoms)."""
    tot = acb(0)
    for j in range(1, 7):
        cf = coefs_by_piece[j]
        if not cf: continue
        def f(t, an):
            NEVAL[0] += 1
            z = zf(t)
            return lam(t) * (I * pi * z) ** l / FACT[l] * (I * pi * z * m).exp() * dens_eval(cf, t)
        tot += acb.integral(f, ts[j - 1], ts[j], rel_tol=arb(2) ** -70, abs_tol=arb(2) ** -100)
    if C != 0:
        for j in range(7):
            mass = (1 if j == 0 else 2) * R(C) * P[j]
            z = zf(ts[j])
            tot += mass * lam(ts[j]) * (I * pi * z) ** l / FACT[l] * (I * pi * z * m).exp()
    return tot.real
def coefs_all(L):
    return {j: piece_coefs(L, j) for j in range(1, 7)}
# ---- direct part: exact Taylor coefficients of p(m+u)=P(m+u)R(m+u) at a node m ----
def p_taylor(L, C, m, order):
    N = order + 1
    a0 = m % 12; k = (m - a0) // 12
    # (2 sin(pi u/12)/u)^2 series
    cs = []
    for q in range(N):
        if q % 2 == 0:
            r = q // 2
            cs.append((-1) ** r * (pi / 12) ** (2 * r + 1) / FACT[2 * r + 1])
        else:
            cs.append(arb(0))
    sdu = arb_series(cs, prec=N)
    ser = (2 * sdu) ** 2
    for a in A:
        if a == a0: continue
        x = arb_series([pi * (m - a) / 12, pi / 12], prec=N)
        ser = ser * (2 * x.sin()) ** 2
    # bracket: c_m + d_m u + u^2 (C + sum_{n != m} ...)
    rest = arb_series([R(C)], prec=N)
    for n, (c, d) in L.items():
        if n == m: continue
        inv = arb_series([arb(1)], prec=N) / arb_series([arb(m - n), arb(1)], prec=N)
        rest = rest + R(c) * inv * inv + R(d) * inv
    cm, dm = L.get(m, (Fr(0), Fr(0)))
    brk = arb_series([R(cm), R(dm)], prec=N) + arb_series([arb(0), arb(0), arb(1)], prec=N) * rest
    out = ser * brk
    return [out.coeffs()[q] if q < len(out.coeffs()) else arb(0) for q in range(N)]

def part_scalar():
    print('== scalar envelopes ==')
    bs = [R(x) for x in ['2.887', '2.665', '2.218', '1.804', '1.486', '1.250', '1.073']]
    as_ = [R(x) for x in ['2.934', '2.713', '2.269', '1.859', '1.548', '1.320', '1.151']]
    gs = [R(x) for x in ['2.933', '2.440', '1.567', '.900', '.482', '.224', '.0597']]
    Ms = [None] + [R(x) for x in ['56', '48', '35', '25', '19', '6.284']]
    betas = [None] + [R(x) for x in ['0', '4.8', '3.17', '1.946', '1.217', '.792']]
    for j in range(7):
        t = ts[j]
        check('|lam(t_%d)|<=b*' % j, abs(lam(t)) <= bs[j], str(abs(lam(t))))
        check('|z(t_%d)|<=a*' % j, abs(zf(t)) <= as_[j], str(abs(zf(t))))
        check('Im z(t_%d)>=g' % j, zf(t).imag >= gs[j], str(zf(t).imag))
    for j in range(1, 7):
        dz = lambda t: 2 * h * B * t / (t ** 2 + h ** 2) ** 2
        mn = dz(ts[j - 1]).min(dz(ts[j]))
        check('beta_%d<=min endpoint |Im z\'|' % j, betas[j] <= mn, str(mn))
        mx = arb(0)
        for a in A:
            s1 = sum((P[l] * (I * pi * ts[l] * a).exp() for l in range(j, 7)), acb(0))
            mx = mx.max(2 * pi * abs(s1))
            for u in (ts[j - 1], ts[j]):
                s2 = sum(((ts[l] - u) * P[l] * (I * pi * ts[l] * a).exp() for l in range(j, 7)), acb(0))
                mx = mx.max(2 * pi ** 2 * abs(s2))
        check('density bound M_%d' % j, mx < Ms[j], str(mx))
    ms = [abs(P[0])] + [2 * abs(P[j]) for j in range(1, 7)]
    check('atom masses (5,4,4,4,2,4,2)', all(abs(ms[j] - [5, 4, 4, 4, 2, 4, 2][j]) < arb(10) ** -30 for j in range(7)))
    check('Im z(1) = 26/435', abs(zf(arb(1)).imag - arb(26) / 435) < arb(10) ** -30)
    check('Q_n>1.75, |D_n|<2.12', all(QD[a][0] > R('1.75') and abs(QD[a][1]) < R('2.12') for a in A))
    exact = {0: (pi ** 2 / 3, -7 * pi * arb(3).sqrt() / 18), 1: (2 * pi ** 2 / 3 * (2 - arb(3).sqrt()), pi * (3 + arb(3).sqrt()) / 18),
             3: (2 * pi ** 2 / 3 * (2 - arb(3).sqrt()), -pi * (3 + arb(3).sqrt()) / 18), 4: (pi ** 2 / 3, 7 * pi * arb(3).sqrt() / 18),
             7: (2 * pi ** 2 / 3 * (2 + arb(3).sqrt()), -pi * (3 - arb(3).sqrt()) / 18), 9: (2 * pi ** 2 / 3 * (2 + arb(3).sqrt()), pi * (3 - arb(3).sqrt()) / 18)}
    check('Q_n,D_n table', all(abs(QD[a][0] - exact[a][0]) < arb(10) ** -25 and abs(QD[a][1] - exact[a][1]) < arb(10) ** -25 for a in A))
    def Lj(j, s):
        mn = arb(1) / 6 if betas[j] == 0 or s == 0 else (arb(1) / 6).min(1 / (pi * betas[j] * s))
        return bs[j - 1] * (-pi * gs[j] * s).exp() * mn
    def LP(j, s): return ms[j] * bs[j] * (-pi * gs[j] * s).exp()
    for x, bound, cert in ((0, '27.7', '27.618723'), (4, '.272', '.269094'), (16, '.0145', '.014122'), (100, '4.7e-10', '4.617010e-10')):
        tot = arb(0)
        for n in [n for n in range(x + 1, x + 13) if n % 12 in A]:
            for j in range(1, 7):
                tot += (1 + abs(Dn(n)) + pi * as_[j - 1]) / Qn(n) * Ms[j] * Lj(j, arb(n)) / (1 - (-12 * pi * gs[j]).exp())
        check('row envelope x=%d < %s' % (x, bound), tot < R(bound), str(tot) + ' (paper: <' + cert + ')')
    tot = arb(0)
    for n in [n for n in range(101, 113) if n % 12 in A]:
        for j in range(7):
            tot += (1 + abs(Dn(n)) + pi * as_[j]) / Qn(n) * LP(j, arb(n)) / (1 - (-12 * pi * gs[j]).exp())
    check('atom tail x=100 < 3.2e-8', tot < R('3.2e-8'), str(tot))
    delta = R('3e-5')
    for (sst, nu, bound) in ((0, 0, '.004'), (Fr(1, 2), 1, '.011'), (Fr(1, 2), 2, '.011'), (2, 2, '.002')):
        tot = arb(0)
        for j in range(1, 7):
            tot += Ms[j] * ((pi * ts[j]) ** nu / 6 + Lj(j, R(sst)) * (pi * as_[j - 1]) ** nu)
        tot = delta * tot / FACT[nu]
        check('perturbation cap s*=%s nu=%d < %s' % (sst, nu, bound), tot < R(bound), str(tot))
    tv = R('2.29') * sum(Ms[1:], arb(0)) / 6 + R('.017') * sum(ms, arb(0))
    check('direct TV < 80', tv < 80, str(tv))
    a0 = as_[0]
    for (mm, Y, nu) in ((0, Fr(1, 2), 0), (1, Fr(1), 1), (1, Fr(1), 2), (3, Fr(1), 2), (4, Fr(3, 2), 2)):
        q = 29 + nu; Yb = R(Y)
        ratio = pi * a0 * Yb / (q + 1)
        assert ratio < 1
        inner = 80 * (pi * Yb) ** q
        inner += R('2.3') * sum((Ms[j] * Lj(j, arb(mm)) * (pi * as_[j - 1] * Yb) ** q for j in range(1, 7)), arb(0))
        inner += R('.017') * sum((LP(j, arb(mm)) * (pi * as_[j] * Yb) ** q for j in range(7)), arb(0))
        E = Yb ** (-nu) / FACT[q] / (1 - ratio) * inner
        check('Taylor truncation (m,Y,nu)=(%d,%s,%d) < .001' % (mm, Y, nu), E < R('.001'), str(E))
    s0 = R('41.5')
    E2 = R('.00005') * 77 + R('2.3') * sum((Ms[j] * Lj(j, s0) * (pi * as_[j - 1]) ** 2 for j in range(1, 7)), arb(0)) \
         + R('.017') * sum((LP(j, s0) * (pi * as_[j]) ** 2 for j in range(7)), arb(0))
    d77 = sum((Ms[j] * (pi * ts[j]) ** 2 / 6 for j in range(1, 7)), arb(0))
    check('sum M_j (pi t_j)^2/6 < 77', d77 < 77, str(d77))
    check('tail E\'\' bound < .006', E2 < R('.006'), str(E2))
    # rational lower sums (exact fractions; only n<=40)
    s0f = Fr(83, 2)
    for i, need in ((1, Fr(12, 1000)), (2, Fr(16, 1000))):
        sg = SIG[i]
        u = {n: sg * c for n, (c, d) in LIST[i].items() if n <= 40}
        v = {n: sg * d for n, (c, d) in LIST[i].items() if n <= 40}
        tot = sg * CONST[i] + min(sum(v.values()), 0) / s0f
        for n in u:
            tot += min(u[n], 0) / (s0f - n) ** 2 + min(v[n], 0) * n / (s0f * (s0f - n))
        check('rational lower sum i=%d > %s' % (i, need), tot > need, '%.10f' % float(tot))
    # P lemma: midpoint quotients
    for lo_, hi_ in ((0, 1), (1, 3), (3, 4), (4, 7), (7, 9), (9, 12)):
        mid = arb(lo_ + hi_) / 2
        Pm = arb(1)
        for a in A: Pm *= (2 * (pi * (mid - a) / 12).sin()) ** 2
        qv = Pm / (mid - lo_) ** 2
        check('P midpoint quotient gap[%d,%d] > .68' % (lo_, hi_), qv > R('.68'), str(qv))
    # Schur / inverse arithmetic (exact rationals)
    L_ = Fr(9, 100) + Fr(15, 1000); sch = Fr(28, 100) + Fr(11, 100) * Fr(37, 10)
    Vf = max(32 + (1 + Fr(37, 10)) * Fr(11, 100) * 32 / (1 - sch), (1 + Fr(37, 10)) / (1 - sch))
    dl = Fr(5, 10 ** 10)
    check('Schur arithmetic', L_ < Fr(11, 100) and sch < 1 and Vf < 90 and dl * (1 + 90 * 28) < 1, 'sch=%.4f Vf=%.3f' % (float(sch), float(Vf)))
    c91 = 90 + 2521 * 90 * dl / (1 - 2521 * dl); c2540 = 2521 / (1 - 2521 * dl)
    check('inverse constants 91, 2540', c91 < 91 and c2540 < 2540, '%.6f %.4f' % (float(c91), float(c2540)))
    Rf = Fr(11, 10 ** 8); Rt = Fr(229, 100) * dl + Fr(17, 1000) * Fr(32, 10 ** 9)
    check('coefficient distance < 3e-5', 2 * (91 * Rf + 2540 * Rt) < Fr(3, 10 ** 5), '%.10g' % float(2 * (91 * Rf + 2540 * Rt)))
    norms = {i: sum(abs(c) + abs(d) for c, d in LIST[i].values()) for i in (1, 2)}
    tails = {i: sum(abs(c) + abs(d) for n, (c, d) in LIST[i].items() if n > 40) for i in (1, 2)}
    check('table norms <2.29, tails n>40 <2e-5', all(norms[i] < Fr(229, 100) and tails[i] < Fr(2, 10 ** 5) for i in (1, 2)),
          '%s %s %s %s' % (float(norms[1]), float(norms[2]), float(tails[1]), float(tails[2])))

def part_residual():
    print('== residuals at all 51 finite nodes ==')
    t0 = time.time(); worst = arb(0)
    cf = {i: coefs_all(LIST[i]) for i in (1, 2)}
    for m in nodes(100):
        for i in (1, 2):
            X0 = K_moment(cf[3 - i], CONST[3 - i], m, 0); X1 = K_moment(cf[3 - i], CONST[3 - i], m, 1)
            c, d = LIST[i][m]
            rc = R(c) + X0 / Qn(m)
            rd = R(d) + (1 if (i == 1 and m == 1) else 0) + (X1 - Dn(m) * X0) / Qn(m)
            worst = worst.max(abs(rc)).max(abs(rd))
            if not (abs(rc) < R('3e-10') and abs(rd) < R('3e-10')):
                print('  residual too big m=%d i=%d rc=%s rd=%s' % (m, i, rc, rd))
    check('all 204 residuals < 3e-10', worst < R('3e-10'), 'max ' + str(worst) + ' time %.0fs evals %d' % (time.time() - t0, NEVAL[0]))

def part_matrix():
    print('== finite matrix blocks (basis columns), own inverse ==')
    t0 = time.time()
    J = [1, 3, 4]; E = [n for n in nodes(100) if n >= 7]; E16 = [n for n in E if n <= 16]
    def column(n, typ, rows):
        L = {n: (Fr(1), Fr(0)) if typ == 'c' else (Fr(0), Fr(1))}
        cf = coefs_all(L)
        out_c, out_d = [], []
        for m in rows:
            X0 = K_moment(cf, 0, m, 0); X1 = K_moment(cf, 0, m, 1)
            out_c.append(-X0 / Qn(m)); out_d.append((Dn(m) * X0 - X1) / Qn(m))
        return out_c + out_d
    colsJ = [column(n, 'c', J + E16) for n in J] + [column(n, 'd', J + E16) for n in J]
    nr = len(J + E16)
    # D = S_JJ, ordered (c1,c3,c4,d1,d3,d4)
    def rowsel(col, nodeset):
        allr = J + E16
        idx = [allr.index(m) for m in nodeset]
        return [col[q] for q in idx] + [col[nr + q] for q in idx]
    D = arb_mat([[rowsel(colsJ[cidx], J)[r] for cidx in range(6)] for r in range(6)])
    SE16J = [[rowsel(colsJ[cidx], E16)[r] for cidx in range(6)] for r in range(2 * len(E16))]
    nE16 = max(sum((abs(SE16J[r][c]) for r in range(len(SE16J))), arb(0)) for c in range(6))
    check('||S_{E cap[7,16],J}|| < .088', nE16 < R('.088'), str(nE16))
    colsE = [column(n, 'c', J) for n in E] + [column(n, 'd', J) for n in E]   # 6 x 96 block U=S_{J,E}
    U = arb_mat([[colsE[cidx][r] for cidx in range(len(colsE))] for r in range(6)])
    for eta in (1, -1):
        Ae = arb_mat([[(1 if r == c else 0) - eta * D[r, c] for c in range(6)] for r in range(6)])
        V = Ae.inv()
        nV = max(sum((abs(V[r, c]) for r in range(6)), arb(0)) for c in range(6))
        VU = V * U
        nVU = max(sum((abs(VU[r, c]) for r in range(6)), arb(0)) for c in range(VU.ncols()))
        check('eta=%+d ||(I-eta D)^-1|| < 32' % eta, nV < 32, str(nV))
        check('eta=%+d ||(I-eta D)^-1 S_JE|| < 3.7' % eta, nVU < R('3.7'), str(nVU))
    print('  matrix time %.0fs evals %d' % (time.time() - t0, NEVAL[0]))

def part_bernstein(lo=0, hi=40):
    print('== Bernstein coefficients on half-gaps (centers %d..%d) ==' % (lo, hi))
    t0 = time.time()
    N0 = [0] + nodes(43)
    cf = {i: coefs_all(LIST[i]) for i in (1, 2)}
    groups = {}
    thr = lambda m: R('.74') if m == 0 else R('.18') if m == 1 else R('.02') if m in (3, 4) else R('.009')
    worst_margin = None
    for idx, m in enumerate(N0):
        if m > 40 or m < lo or m > hi: continue
        ys = [Fr(N0[idx + 1] - m, 2)] + ([Fr(N0[idx - 1] - m, 2)] if m > 0 else [])
        for i in (1, 2):
            cases = []
            for y in ys:
                if i == 1 and (m + y < 1 or m == 0): continue  # interval must lie in [1, inf)
                nu = 0 if m == 0 else (1 if (i == 1 and m == 1) else 2)
                cases.append((y, nu))
            if not cases: continue
            ordmax = 28 + max(nu for _, nu in cases)
            U = p_taylor(LIST[i], CONST[i], m, ordmax)
            X = [K_moment(cf[3 - i], CONST[3 - i], m, l) for l in range(ordmax + 1)]
            for y, nu in cases:
                T = [SIG[i] * R(y) ** r * (U[r + nu] + X[r + nu]) for r in range(29)]
                Bmin = None
                for k in range(29):
                    Bk = sum((arb(comb(k, r)) / comb(28, r) * T[r] for r in range(k + 1)), arb(0))
                    Bmin = Bk if Bmin is None else Bmin.min(Bk)
                ok = Bmin > thr(m)
                key = (i, m)
                groups[key] = Bmin if key not in groups else groups[key].min(Bmin)
                check('Bernstein i=%d m=%d y=%s nu=%d: min_k B_k > %s' % (i, m, y, nu, thr(m)), ok, str(Bmin))
    print('  bernstein time %.0fs evals %d' % (time.time() - t0, NEVAL[0]))
    json.dump({'%d,%d' % k: str(v) for k, v in groups.items()}, open(SCR + '/bern_%d_%d.json' % (lo, hi), 'w'))

if __name__ == '__main__':
    what = sys.argv[1]
    if what == 'scalar': part_scalar()
    elif what == 'residual': part_residual()
    elif what == 'matrix': part_matrix()
    elif what == 'bernstein': part_bernstein(int(sys.argv[2]), int(sys.argv[3]))
    print('SUMMARY %s: %d pass, %d fail %s' % (what, len(OK), len(BAD), BAD))

