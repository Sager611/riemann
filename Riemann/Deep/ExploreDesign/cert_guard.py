"""cert_guard.py c m eps prec K2 [maxit]  -- even sector of Weil's form on [-L/2, L/2], L = ln c.
Theorem C (DERIVATION.md sec. 4) on the P = mL lattice with guard-band patches (sec. 5), python-flint balls.
Steps: float pre-scan of w_k = mu(t_k) + kappa^(t_k) up to kscan; patch LP (float) at forced nodes above the low block;
coefficients rounded to dyadics (exact data); Arb re-evaluation of every flagged node (float w < tau + FLAG); analytic tail beyond kscan;
Woodbury (exact arb solve) and an Arb Cholesky of Id - W - nu Id.  Prints one JSON line."""
import math, sys, json, time, resource
import numpy as np
from flint import arb, acb, arb_mat, ctx
sys.path.insert(0, '.')
from wsym import mu as mu_f, vonmangoldt_table, sigma
from scipy.optimize import linprog

c = int(sys.argv[1]); m = int(sys.argv[2]); epss = sys.argv[3]; ctx.prec = int(sys.argv[4]); K2 = int(sys.argv[5])
MAXIT = int(sys.argv[6]) if len(sys.argv) > 6 else 3
T0 = time.time()
def LOG(msg):
    sys.stderr.write(f'[{time.time()-T0:7.1f}s {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:6.0f}MB] {msg}\n'); sys.stderr.flush()
Lf = math.log(c); Pf = m * Lf
L = arb(c).log(); P = m * L; PI = arb.pi()
eps = arb(2) ** (-int(epss[3:])) if epss.startswith('2^-') else arb(epss)
deff = arb(1) / 2; tau = m * deff / (m - 1); tauf = float(tau.mid())
FLAG = 0.01
VM = vonmangoldt_table(c)

def mu_arb(t):
    t = arb(t)
    v = -PI.log() + acb(arb(1) / 4, t / 2).digamma().real
    cl, sl = (t * L).cos(), (t * L).sin()
    N = 0
    while True:
        b = 2 * N + arb(1) / 2
        w = arb(c) ** (-b)
        v += 2 * w * (b * cl - t * sl) / (b * b + t * t)
        N += 1
        if float(w.mid()) < 2.0 ** (-ctx.prec - 20):
            break
    v += arb(0, 4 * float(w.mid()) * (1 + abs(float(t.mid()))))   # tail of the b_n series (geometric, ratio c^-2)
    for n, lam in VM:
        p = int(round(math.exp(lam)))
        v -= 2 * arb(p).log() / arb(n).sqrt() * (t * arb(n).log()).cos()
    return v

# ---- patches: list of (t_d, D, xs, a, b) with dyadic coefficients ----
Gf = Pf / 2 - Lf; MP = 12
Dp = Gf / (MP + 1); xsf = Lf + np.arange(1, MP + 1) * Dp
patches = []

def kap_f(t):
    t = np.asarray(t, dtype=float); out = np.zeros_like(t)
    for td, a, b in patches:
        for sgn in (1, -1):
            u = t - sgn * td
            v = u * Dp / 2
            env = Dp * np.where(np.abs(v) < 1e-12, 1.0, (np.sin(v) / np.where(v == 0, 1, v)) ** 2)
            out += env * ((np.cos(np.outer(u, xsf)) @ a) + sgn * (np.sin(np.outer(u, xsf)) @ b))
    return out

def kap_arb(t):
    t = arb(t); out = arb(0)
    D = (P / 2 - L) / (MP + 1)
    for td, a, b in patches:
        for sgn in (1, -1):
            u = t - sgn * arb(td)
            v = u * D / 2
            env = D * ((v.sin() / v) ** 2 if abs(float(v.mid())) > 1e-30 else arb(1))
            s = arb(0)
            for i in range(MP):
                x = L + (i + 1) * D
                s += arb(a[i]) * (u * x).cos() + sgn * arb(b[i]) * (u * x).sin()
            out += env * s
    return out

def wf(k):
    t = 2 * math.pi * np.asarray(k, dtype=float) / Pf
    return mu_f(t, c) + kap_f(t)

def thr(k):
    return np.where(np.asarray(k) % m == 0, 0.0, tauf)

# ---- scan end: for t >= T, w >= log(T/2pi) - log(pi)?  use Re psi(1/4+it/2) increasing; mu >= psi_low - log pi - 2 sigma - bnd - |kap| ----
sig = sigma(c)
def lower_bound_f(t):
    return float((acb(arb(1)/4, arb(t)/2).digamma().real - arb.pi().log()).lower()) - 2 * sig - 4 * c ** -0.5 * (1 + 1.0 / t)
tlow = math.pi * K2 / Lf            # low block height (Chekmachev's odd K2)
T = tlow
while lower_bound_f(T) < tauf + 1.0:
    T *= 1.05
kscan = int(math.ceil(T * Pf / (2 * math.pi)))
klow = int(math.floor(tlow * Pf / (2 * math.pi)))
J = K2   # sublattice-0 exact block j <= J (Chekmachev's J = K2)

def forced_scan():
    out = []
    CH = 200000
    for a0 in range(1, kscan + 1, CH):
        k = np.arange(a0, min(kscan + 1, a0 + CH))
        w = wf(k); th = thr(k)
        sel = (w < th + FLAG) & ~((k <= klow) & (k % m != 0)) & ~((k % m == 0) & (k // m <= J))
        out.extend(k[sel].tolist())
    return out

def design_patch(kd):
    td = 2 * math.pi * kd / Pf
    W = 25.0
    ka = int(math.ceil((td - W) * Pf / (2 * math.pi))); kb = int(math.floor((td + W) * Pf / (2 * math.pi)))
    k = np.arange(max(ka, 1), kb + 1); t = 2 * math.pi * k / Pf
    base = wf(k); th = thr(k) + 2 * FLAG
    u = t - td; v = u * Dp / 2
    env = Dp * np.where(np.abs(v) < 1e-12, 1.0, (np.sin(v) / np.where(v == 0, 1, v)) ** 2)
    Bm = np.hstack([env[:, None] * np.cos(np.outer(u, xsf)), env[:, None] * np.sin(np.outer(u, xsf))])
    n = Bm.shape[1]
    outer = np.abs(u) > W - 6
    A = np.vstack([np.hstack([-Bm, np.ones((len(t), 1))]), np.hstack([Bm[outer], np.zeros((outer.sum(), 1))]), np.hstack([-Bm[outer], np.zeros((outer.sum(), 1))])])
    bb = np.concatenate([base - th, np.full(outer.sum(), 0.2), np.full(outer.sum(), 0.2)])
    cobj = np.zeros(n + 1); cobj[-1] = -1
    r = linprog(cobj, A_ub=A, b_ub=bb, bounds=[(-200, 200)] * n + [(None, 0.5)], method='highs')
    if r.status != 0 or r.x[-1] < 0:
        return None
    q = 2.0 ** 20
    a = np.round(r.x[:MP] * q) / q; b = np.round(r.x[MP:n] * q) / q
    return (td, a, b)

info = {'c': c, 'm': m, 'eps': epss, 'prec': ctx.prec, 'K2': K2, 'kscan_P': kscan, 'tscan': T, 'klow': klow}
LOG(f'kscan {kscan} klow {klow}')
forced0 = forced_scan()
LOG(f'forced before {len(forced0)}')
info['forced_high_before'] = len(forced0)
def local_forced(td, R=80.0):
    ka = int(math.ceil((td - R) * Pf / (2 * math.pi))); kb = int(math.floor((td + R) * Pf / (2 * math.pi)))
    k = np.arange(max(ka, 1), kb + 1)
    k = k[~((k <= klow) & (k % m != 0)) & ~((k % m == 0) & (k // m <= J))]
    w = wf(k); return set(k[w < thr(k) + FLAG].tolist())
nacc = nrej = 0
for kd in sorted(forced0, reverse=True):
    td = 2 * math.pi * kd / Pf
    if any(abs(td - x[0]) < 3.0 for x in patches):
        continue
    before = local_forced(td)
    if kd not in before:
        continue
    pch = design_patch(kd)
    if pch is None:
        nrej += 1; continue
    patches.append(pch)
    after = local_forced(td)
    if len(after) < len(before) and after <= before:
        nacc += 1
    else:
        patches.pop(); nrej += 1
info['patch_accepted'] = nacc; info['patch_rejected'] = nrej
LOG(f'patch loop done: {len(patches)} patches')
forced = forced_scan()
info['patches'] = len(patches); info['forced_high_after_float'] = len(forced)
# Arb re-check of flagged nodes: those still below threshold go to K / E
Kh, Eh = [], []
for k in forced:
    w = mu_arb(2 * PI * k / P) + kap_arb(2 * PI * k / P)
    th = 0 if k % m == 0 else tau
    if not (w > th):
        (Eh if k % m == 0 else Kh).append(k)
info['forced_high_after_arb'] = len(Kh) + len(Eh)
LOG(f'arb recheck: {len(forced)} flagged, {len(Kh)+len(Eh)} kept')
import os
if os.environ.get('CG_SCAN_ONLY'):
    info['Kh'] = len(Kh); info['Eh'] = len(Eh); info['time_s'] = round(time.time() - T0, 1); print(json.dumps(info)); sys.exit(0)
# tail of patches beyond kscan: |kap(t)| <= sum_d A_d Dp (2/(Dp |t - t_d|))^2 for t > tscan; included in the +1.0 slack of lower_bound_f
if patches:
    tail = 2 * sum(Dp * (np.abs(a).sum() + np.abs(b).sum()) * (2 / (Dp * max(1.0, T - td))) ** 2 for td, a, b in patches)
else:
    tail = 0.0
info['patch_tail_at_tscan'] = tail
assert tail < 0.5, 'patch tail too large at tscan'
# ---------------- certificate (Theorem C, even sector) ----------------
# neighbour rule (Chekmachev Def. 6.1, transplanted): E gets j within +-NBJ spacings of every forced K node; every forced E node
# gets the r >= 1 nodes within +-NBK spacings into K and j +- 3 into E.  Spacing 2 pi/L.
NBJ, NBK = 32, 1.01
Kset = set(k for k in range(1, klow + 1) if k % m != 0) | set(Kh)
Eset = set(range(0, J + 1)) | set(k // m for k in Eh)
for k in Kh:
    jc = (2 * math.pi * k / Pf) * Lf / (2 * math.pi)
    Eset |= set(range(max(0, int(jc) - NBJ), int(jc) + NBJ + 2))
jlow = np.arange(0, J + 1); wlow = wf(m * jlow)
Eneg = [m * int(j) for j in jlow[wlow < FLAG]] + list(Eh)   # every sublattice-0 node with w < 0 (low or high) gets neighbours
for k in Eneg:
    j = k // m
    Eset |= set(range(max(0, j - 3), j + 4))
    for kk in range(m * j - int(NBK * m) - 1, m * j + int(NBK * m) + 2):
        if kk > 0 and kk % m != 0:
            Kset.add(kk)
Klist = sorted(Kset); Elist = sorted(Eset)
tK = [2 * PI * k / P for k in Klist]
wK = [mu_arb(t) + kap_arb(t) for t in tK]
sE = [2 * PI * j / L for j in Elist]
wE = [mu_arb(2 * PI * (m * j) / P) + kap_arb(2 * PI * (m * j) / P) for j in Elist]
de = deff - eps
dl = [w / m + de for w in wE]
NE = [(1 / L).sqrt() if j == 0 else (2 / L).sqrt() for j in Elist]
Ep = [i for i in range(len(Elist)) if dl[i] > 0]
bad = [i for i in range(len(Elist)) if not (dl[i] > 0)]
Kp = [i for i in range(len(Klist)) if wK[i] > tau]
Km = [i for i in range(len(Klist)) if not (wK[i] > tau)]
info.update({'K': len(Klist), 'E': len(Elist), 'Kplus': len(Kp), 'Kminus': len(Km), 'bad_even': len(bad)})
sh, ch = (L / 4).sinh(), (L / 4).cosh()
def ip_phi_c(i, t):
    s = sE[i]
    return NE[i] * (((t - s) * L / 2).sin() / (t - s) + ((t + s) * L / 2).sin() / (t + s))
def ip_c_c(t1, t2, same):
    if same:
        return L / 2 + (t1 * L).sin() / (2 * t1)
    return ((t1 - t2) * L / 2).sin() / (t1 - t2) + ((t1 + t2) * L / 2).sin() / (t1 + t2)
def ip_phi_v(i):
    s = sE[i]
    return NE[i] * 2 * ((sh / 2) * (s * L / 2).cos() + s * ch * (s * L / 2).sin()) / (arb(1) / 4 + s * s)
def ip_c_v(t):
    return 2 * ((sh / 2) * (t * L / 2).cos() + t * ch * (t * L / 2).sin()) / (arb(1) / 4 + t * t)
vv = L / 2 + (L / 2).sinh()
gam = {i: 1 / dl[i] - 1 / de for i in Ep}
# column descriptors: ('c', idx_in_K) or ('phi', idx_in_E) or ('v',)
Ucols = [('c', i) for i in Kp] + [('v',)]
Dvals = [2 * (wK[i] - tau) / P for i in Kp] + [arb(2)]
Bcols = [('c', i) for i in Km] + [('phi', i) for i in bad]
Bsc = [(2 * (tau - wK[i]) / P).sqrt() for i in Km] + [(de - dl[i]).sqrt() for i in bad]
# coefficient vectors on E' (for the gamma correction) and Gram with identity part
Epl = list(Ep)
def coefE(col):
    if col[0] == 'c':
        t = tK[col[1]]; return [ip_phi_c(i, t) for i in Epl]
    if col[0] == 'v':
        return [ip_phi_v(i) for i in Epl]
    return [arb(0)] * len(Epl)     # phi_j with j bad is orthogonal to E'
def plain(c1, c2):
    if c1[0] == 'c' and c2[0] == 'c':
        return ip_c_c(tK[c1[1]], tK[c2[1]], c1[1] == c2[1])
    if c1[0] == 'v' and c2[0] == 'v':
        return vv
    if c1[0] == 'phi' and c2[0] == 'phi':
        return arb(1) if c1[1] == c2[1] else arb(0)
    if c1[0] == 'v' or c2[0] == 'v':
        o = c2 if c1[0] == 'v' else c1
        return ip_c_v(tK[o[1]]) if o[0] == 'c' else ip_phi_v(o[1])
    a, b = (c1, c2) if c1[0] == 'phi' else (c2, c1)
    return ip_phi_c(a[1], tK[b[1]])
def gram(colsA, colsB):
    gv = [gam[i] for i in Epl]
    CA = arb_mat(len(colsA), len(Epl), [x for col in colsA for x in coefE(col)]) if Epl else None
    CB = arb_mat(len(colsB), len(Epl), [x * gv[r] for col in colsB for r, x in enumerate(coefE(col))]) if Epl else None
    M = arb_mat(len(colsA), len(colsB), [plain(x, y) / de for x in colsA for y in colsB])
    if Epl:
        M = M + CA * CB.transpose()
    return M
LOG(f'sets K {len(Klist)} E {len(Elist)} U {len(Ucols)} B {len(Bcols)}')
t1 = time.time()
X = gram(Bcols, Bcols); Y = gram(Ucols, Bcols); S = gram(Ucols, Ucols)
for i in range(len(Ucols)):
    S[i, i] += 1 / Dvals[i]
nb = len(Bcols)
for p in range(nb):
    for q in range(nb):
        X[p, q] *= Bsc[p] * Bsc[q]
for p in range(len(Ucols)):
    for q in range(nb):
        Y[p, q] *= Bsc[q]
LOG('grams done')
Z = S.solve(Y)
LOG('solve done')
Wm = X - Y.transpose() * Z
A = arb_mat(nb, nb)
for p in range(nb):
    for q in range(nb):
        A[p, q] = (1 if p == q else 0) - Wm[p, q]
import mpmath
mpmath.mp.prec = ctx.prec
if nb:
    Mm = mpmath.matrix(nb, nb)
    for p_ in range(nb):
        for q_ in range(nb):
            Mm[p_, q_] = (mpmath.mpf(A[p_, q_].mid().str(int(ctx.prec / 3.3), radius=False, more=True)) + mpmath.mpf(A[q_, p_].mid().str(int(ctx.prec / 3.3), radius=False, more=True))) / 2
    ev = mpmath.eigsy(Mm, eigvals_only=True)
    lam = float(min(ev)); info['eig_low3'] = [mpmath.nstr(x, 6) for x in sorted(ev)[:3]]
    info['maxabsA'] = float(max(abs(Mm[p_, q_]) for p_ in range(nb) for q_ in range(nb)))
else:
    lam = 1.0
info['lam_min_est'] = lam
def ball_cho_pd(Mx, n):
    """interval Cholesky; success => every symmetric matrix with lower triangle in the balls is positive definite"""
    Lc = [[arb(0)] * n for _ in range(n)]
    for j in range(n):
        d = Mx[j, j] - sum((Lc[j][k] * Lc[j][k] for k in range(j)), arb(0))
        if not (d > 0):
            return False
        Lc[j][j] = d.sqrt()
        for i in range(j + 1, n):
            Lc[i][j] = (Mx[i, j] - sum((Lc[i][k] * Lc[j][k] for k in range(j)), arb(0))) / Lc[j][j]
    return True
ok = False; nu = None
if nb == 0:
    ok = True
else:
    for frac in (0.5, 0.125, 0.01):
        if lam <= 0:
            break
        nu = arb(lam * frac)
        Bn = arb_mat(nb, nb, [A[p, q] - (nu if p == q else 0) for p in range(nb) for q in range(nb)])
        if ball_cho_pd(Bn, nb):
            ok = True; break
info.update({'order': nb, 'U': len(Ucols), 'rigorous': ok, 'nu': None if nu is None else float(nu.mid()),
             'time_s': round(time.time() - T0, 1), 'time_matrix_s': round(time.time() - t1, 1)})
print(json.dumps(info))
