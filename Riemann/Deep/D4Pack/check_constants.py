#!/usr/bin/env python3
"""check_constants.py -- recheck the quantitative constants of arXiv 2609.25120v2 with mpmath.
Every number is recomputed from its definition (not from the paper's closed forms),
then compared with the value printed in the TeX (line numbers of D4.tex v2)."""
import re, sys, itertools
import numpy as np
import mpmath as mp

mp.mp.dps = 30
pi = mp.pi
TEX = sys.argv[1] if len(sys.argv) > 1 else None
out = []
def rep(name, got, claimed, tol, line, trunc=False):
    ok = abs(mp.mpf(got) - mp.mpf(claimed)) <= tol
    if trunc:  # paper prints the leading digits followed by \ldots
        nd = len(claimed.split('.')[1])
        ok = mp.floor(mp.mpf(got) * 10**nd) == mp.floor(mp.mpf(claimed) * 10**nd + mp.mpf('0.5'))
    out.append(ok)
    print(("OK  " if ok else "FAIL") + f"  {name}: computed {mp.nstr(got, 12)}  paper {claimed}  (TeX l.{line})")

# ---------------------------------------------------------------- caps and lenses on S^3
def C_closed(r):
    return pi * (2 * r - mp.sin(2 * r))
def C_int(r):  # 2 pi int_{cos r}^1 t * 2 arccos(cos r / t) dt  (Hopf coordinates)
    c = mp.cos(r)
    return 2 * pi * mp.quad(lambda t: 2 * t * mp.acos(c / t), [c, 1])
def Lam_int(r, g):  # two caps of radius r, centres g apart
    c = mp.cos(r)
    if 2 * r <= g:
        return mp.mpf(0)
    t0 = c / mp.cos(g / 2)
    return 2 * pi * mp.quad(lambda t: t * (2 * mp.acos(c / t) - g), [t0, 1])
def Lam_paper(r, g):  # closed form printed in the proof of thm:certificate
    return 4 * pi * ((r - g / 2) / 2 - (mp.sin(2 * r) - mp.sin(g)) / 4
                     - mp.tan(g / 2) / 2 * (mp.sin(r) ** 2 - mp.sin(g / 2) ** 2))

for r in [mp.mpf('0.3'), mp.mpf('0.6'), mp.mpf('1.2')]:
    assert abs(C_int(r) - C_closed(r)) < mp.mpf('1e-20')
print("cap measure: integral == pi(2r - sin 2r) at 3 radii")
mx = 0
for r, g in [(0.55, 1.0), (0.6, 1.05), (0.62, 0.9), (0.66, 1.2)]:
    r = mp.mpf(r); g = mp.mpf(g)
    mx = max(mx, abs(Lam_int(r, g) - Lam_paper(r, g)))
print("lens closed form vs own integral, max diff:", mp.nstr(mx, 3))
out.append(mx < 1e-15)

# ---------------------------------------------------------------- covering table
def r_m(m):
    return mp.findroot(lambda r: 2 * r - mp.sin(2 * r) - 2 * pi / m, (mp.mpf('0.1'), mp.mpf('1.56')), solver='bisect')
def cov(m):
    return pi * m / 3 * mp.tan(r_m(m)) ** 3
tab = {}
if TEX:
    s = open(TEX).read().split('\n')
    for ln in range(2110, 2135):
        for a, b, c in re.findall(r'(\d+)\s*&\s*\$([0-9.]+)\$\s*&\s*\$(?:\\mathbf\{)?([0-9.]+)\}?\$', s[ln]):
            tab[int(a)] = (b, c, ln + 1)
bad = 0
for m in range(5, 25):
    rm = r_m(m)
    deg = rm * 180 / pi
    b = cov(m)
    if m in tab:
        d_c, b_c, ln = tab[m]
        okd = abs(deg - mp.mpf(d_c)) <= mp.mpf('6e-7')
        okb = abs(b - mp.mpf(b_c)) <= mp.mpf('6e-7')
        if not (okd and okb):
            bad += 1
            print(f"FAIL covering m={m}: r_m={mp.nstr(deg,10)} (paper {d_c}), bound={mp.nstr(b,10)} (paper {b_c})")
print(f"covering table: {len(tab)} rows parsed, {bad} mismatches;  m=22 -> {mp.nstr(cov(22),10)}, m=23 -> {mp.nstr(cov(23),10)}, m=24 -> {mp.nstr(cov(24),10)}")
out.append(bad == 0 and len(tab) == 20)
vals = [cov(m) for m in range(5, 25)]
out.append(all(vals[i] > vals[i + 1] for i in range(len(vals) - 1)))
print("covering bound strictly decreasing on 5..24:", out[-1])
# direct check of the covering bound against an independent formula: the layer-cake integral
m = 22; R = r_m(m)
lc = (2 * pi**2 + mp.quad(lambda r: (2 * pi**2 - m * C_closed(r)) * 4 * mp.sec(r)**4 * mp.tan(r), [0, R])) / 4
rep("covering bound m=22 via layer-cake integral", lc, '8.046376', mp.mpf('6e-7'), 2125)

# ---------------------------------------------------------------- radii
r_star = mp.asin(1 / mp.sqrt(3)); r_4 = mp.asin(mp.sqrt(mp.mpf(3) / 8)); r23 = r_m(23)
rep("r_* (deg)", r_star * 180 / pi, '35.2644', mp.mpf('5e-5'), 4750)
rep("r_4 (deg)", r_4 * 180 / pi, '37.7612', mp.mpf('5e-5'), 4751)
rep("2 r_23 (deg)", 2 * r23 * 180 / pi, '69.2213', mp.mpf('5e-5'), 3489)
rep("2 r_* (deg)", 2 * r_star * 180 / pi, '70.5288', mp.mpf('5e-5'), 3490)

# ---------------------------------------------------------------- A_*, A_4, omega
def A_of(R, m=23):
    return (2 * pi**2 + mp.quad(lambda r: (2 * pi**2 - m * C_closed(r)) * 4 * mp.sec(r)**4 * mp.tan(r), [0, R])) / 4
def omega(g, R):
    if 2 * R <= g:
        return mp.mpf(0)
    return mp.quad(lambda r: Lam_int(r, g) * mp.tan(r) * mp.sec(r)**4, [g / 2, R])
Astar = A_of(r_star)
rep("A_*", Astar, '7.907144430', mp.mpf('2e-9'), 4879)
closed = 9 * pi**2 / 8 - 207 * pi / 8 * mp.atan(1 / mp.sqrt(2)) + 253 * pi / (12 * mp.sqrt(2))
rep("A_* closed form (eq:A-star-closed)", closed, str(Astar), mp.mpf('1e-20'), 5167)
rep("8 - A_*", 8 - Astar, '0.0928555702', mp.mpf('2e-10'), 5088)
A4 = A_of(r_4)
rep("A_4", A4, '7.647558', mp.mpf('2e-6'), 4801)
rep("8 - A_4", 8 - A4, '0.352442', mp.mpf('2e-6'), 4887)
g60 = pi / 3
w_star = omega(g60, r_star); w_23 = omega(g60, r23); w_4 = omega(g60, r_4)
rep("omega(60deg) at r_*", w_star, '0.00144540', mp.mpf('6e-9'), 4973, trunc=True)
rep("omega(60deg) at r_23", w_23, '0.000922', mp.mpf('6e-7'), 3488)
rep("omega_4(60deg)", w_4, '0.00565947', mp.mpf('6e-9'), 4854)
rep("88 omega_4(60deg)", 88 * w_4, '0.498033', mp.mpf('2e-6'), 4854)
rep("deletion pair-only bound A_*+88 omega - 8 (spare)", Astar + 88 * w_star - 8, '0.034', mp.mpf('6e-4'), 4866)

# ---------------------------------------------------------------- omega_3 (three caps, pairwise 60 deg) by quadrature
kap = 2 / np.sqrt(6.0)
# orthonormal frame of the 3-space of c1,c2,c3 with g = (c1+c2+c3)/|.| along axis 0
c = np.array([[1, 1, 0], [1, 0, 1], [0, 1, 1]], float) / np.sqrt(2)  # pairwise inner products 1/2
g = c.sum(0); g /= np.linalg.norm(g)
e1 = c[0] - (c[0] @ g) * g; e1 /= np.linalg.norm(e1); e2 = np.cross(g, e1)
phi = np.linspace(0, 2 * np.pi, 4001)[:-1]
W = np.outer(np.cos(phi), e1) + np.outer(np.sin(phi), e2)
S = W @ c.T  # (nphi, 3)
Rr = np.sqrt(kap**2 + S**2); dl = np.arctan2(S, kap)
def A3(ca):  # area on S^2 of three caps with cos(radius) = ca, ca <= kap
    bmax = (dl + np.arccos(np.clip(ca / Rr, -1, 1))).min(axis=1)
    return (1 - np.cos(bmax)).mean() * 2 * np.pi
def Lam3(r):
    cr = np.cos(r)
    psimax = np.arccos(min(1.0, cr / kap))
    x, w = np.polynomial.legendre.leggauss(80)
    psi = psimax * x; wts = psimax * w
    return sum(wi * np.cos(p)**2 * A3(cr / np.cos(p)) for p, wi in zip(psi, wts))
xr, wr = np.polynomial.legendre.leggauss(60)
a, b = float(r_star), float(r_4)
rr = (b - a) / 2 * xr + (a + b) / 2
w3 = sum(wi * (b - a) / 2 * Lam3(r) * np.tan(r) / np.cos(r)**4 for r, wi in zip(rr, wr))
print(f"omega_3(60,60,60) = {w3:.6e}  (paper 5.65e-5, TeX l.4856)")
out.append(abs(w3 - 5.65e-5) < 1e-6)
del_bound = A4 + 88 * w_4 - 84 * mp.mpf(w3)
rep("deletion three-point bound A_4+88w4-84w3 (uses own omega_3)", del_bound, '8.140848', mp.mpf('3e-6'), 4853)

# ---------------------------------------------------------------- D4 pair / triangle counts
R24 = [np.array(v, float) / np.sqrt(2) for v in
       [tuple((s1 if k == i else s2 if k == j else 0) for k in range(4))
        for i in range(4) for j in range(i + 1, 4) for s1 in (1, -1) for s2 in (1, -1)]]
G = np.array([[x @ y for y in R24] for x in R24])
from collections import Counter
cnt = Counter(round(G[i, j], 6) for i in range(24) for j in range(i + 1, 24))
print("D4 pair inner products:", dict(cnt), "(paper 96,72,96,12, TeX l.2856)")
out.append(cnt[0.5] == 96 and cnt[0.0] == 72 and cnt[-0.5] == 96 and cnt[-1.0] == 12)
keep = list(range(1, 24))
tp = sum(1 for i, j in itertools.combinations(keep, 2) if abs(G[i, j] - .5) < 1e-9)
tt = sum(1 for i, j, k in itertools.combinations(keep, 3) if abs(G[i, j] - .5) < 1e-9 and abs(G[i, k] - .5) < 1e-9 and abs(G[j, k] - .5) < 1e-9)
print(f"deletion: {tp} tight pairs (paper 88), {tt} tight triangles (paper 84)")
out.append(tp == 88 and tt == 84)
print("\nALL CHECKS PASS" if all(out) else f"\n{out.count(False)} CHECK(S) FAILED")

