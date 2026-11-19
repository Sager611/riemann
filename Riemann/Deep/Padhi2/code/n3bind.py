"""padhi2 step A: can Bellotti's N3 (Table 2 rows, uniform bound) ever be the minimum in Padhi's N+?
Cellwise and rigorous for the published constants: on every Lemma D cell (k, lk) (delta in ((k-1)/1e4, k/1e4], log T in
(lk-0.5, lk] for lk <= 31, (lk-1, lk] above), each N3 candidate is >= its value at the lower-left corner (it increases in
delta and in log T), and the cell's N4 value is a constant (certified logN_up). If N3_lower > N4_up on every cell where both
apply, N3 never attains the min there, for the printed constants and a fortiori for any larger (corrected) C.
Also: the largest log T1 over the published scan (N4 exists only for log T <= 40), and missing cells."""
import math, json, glob, os, sys
sys.path.insert(0, os.path.dirname(__file__))
import common as CM
from flint import arb
S = CM.S
rows = [(math.log(a) if j == 0 else a, b, float(al), C, B) for j, (a, b, al, C, B) in enumerate(CM.ORIG['ROWS'])]
rows.append((math.log(3e12), 1e9, 0.9927, arb('1.62e11'), arb('1.448')))      # uniform (Thm 1.1)
rows.append((math.log(3e12), 29.0, 0.9927, arb('37341.72'), arb('14.160')))   # T <= e^29 (same as row 1)
LKS = [29, 29.5, 30, 30.5, 31] + list(range(32, 41))
worst = None; nboth = 0; missing = []
for lk in LKS:
    lo_c = lk - 0.5 if lk <= 31 else lk - 1.0
    for k in range(1, 301):
        key = '%d,%s' % (k, lk if isinstance(lk, int) else (int(lk) if float(lk).is_integer() else float(lk)))
        dlo, dhi = (k - 1)/1e4, k/1e4
        for j, (lo, hi, al, C, B) in enumerate(rows):
            l0, l1 = max(lo_c, lo), min(float(lk), hi)
            if l1 <= l0 or dlo >= 1 - al + 1e-12: continue        # no overlap in log T or in delta (s >= alpha needs delta <= 1-alpha)
            if key not in S.CELL: missing.append((key, j)); continue
            n3lo = C.log() + B*arb(dlo)*arb(l0)
            gap = float(n3lo.lower()) - float(S.CELL[key].upper()); nboth += 1
            if worst is None or gap < worst[0]: worst = (gap, key, j, float(n3lo.lower()), float(S.CELL[key].upper()))
print('cell/row pairs compared:', nboth, ' missing N4 cells where N3 applies:', sorted(set(m[0] for m in missing))[:20], len(missing))
print('min over pairs of (N3 lower - N4 upper) in log N:', worst)
# largest log T1 used by the published scan: T1 = c/eps with eps as in q_arb (same float recipe)
mx = (0, None)
for f in glob.glob(os.path.join(S.HERE, 'scan_arb_*.txt')):
    for line in open(f):
        p = line.split(); L1, L2 = float(p[0]), float(p[1]); cc = int(p[3][2:]); fr = float(p[4][2:])
        tol_f = S.RU.eps_max(L2 + 1)/(1 + 2/L1); eps_f = fr*tol_f*math.sqrt(cc)/(0.7946*1.02)
        for _ in range(3):
            lam = math.log(10*math.sqrt(cc)/eps_f); kap = (L2 + lam)/(L1 - lam); eps_f = fr*tol_f*math.sqrt(cc)/(0.7946*1.02*kap)
        lT1 = math.log(cc/eps_f)
        if lT1 > mx[0]: mx = (lT1, (L1, L2))
print('max log T1 over the published scan: %.4f at %s' % mx)
# zero-free edge vs the rows' alpha0: N3 is evaluated only where s < edge(log T)
for lt in [28.73, 29, 29.5, 30, 31, 32, 34, 36, 38, 40, 41]:
    print('log T %.2f: edge delta = lambda/l = %.6f' % (lt, float((1 - S.edge(lt)).mid())))

# missing N4 cells where an N3 row applies: are they inside the certified zero-free region (never evaluated)?
bad = []
for key, j in missing:
    k, lk = key.split(','); k = int(k); lk = float(lk)
    zf = 1 - S.edge(lk)                      # zero-free delta at the cell's top height (smallest over the cell)
    if not (arb(k)/10**4 < zf): bad.append((key, j, float(zf.mid())))
print('missing cells (with an N3 row) NOT inside the zero-free region:', len(bad), bad[:10])

# per N3 candidate: where it is evaluated at all (non-zero-free, log T <= max log T1) and its smallest margin over N4 there
LT1MAX = mx[0]
names = ['row%d (e^%s, e^%s]' % (j + 1, ('log3e12' if j == 0 else '%g' % rows[j][0]), '%g' % rows[j][1]) for j in range(13)] + ['uniform 1.62e11', 'T<=e^29 (=row1)']
for j, (lo, hi, al, C, B) in enumerate(rows):
    top = min(hi, LT1MAX)
    if top <= lo: print('%-28s never evaluated (log T > %.4f)' % (names[j], LT1MAX)); continue
    zf_top = float((1 - S.edge(top)).lower())
    if 1 - al <= zf_top: print('%-28s inside the zero-free region (delta <= %.4f <= edge delta %.6f at log T %.4f)' % (names[j], 1 - al, zf_top, top)); continue
    g = None
    for lk in LKS:
        lo_c = lk - 0.5 if lk <= 31 else lk - 1.0
        l0, l1 = max(lo_c, lo), min(float(lk), top)
        if l1 <= l0: continue
        for k in range(1, 301):
            dlo = (k - 1)/1e4
            if dlo >= 1 - al + 1e-12: continue
            key = '%d,%s' % (k, lk if isinstance(lk, int) else (int(lk) if float(lk).is_integer() else float(lk)))
            if key not in S.CELL: continue
            gap = float((C.log() + B*arb(dlo)*arb(l0)).lower()) - float(S.CELL[key].upper())
            if g is None or gap < g[0]: g = (gap, key)
    print('%-28s evaluated for log T in (%.4f, %.4f], delta in (edge, %.4f]; min log(N3/N4) = %.4f at cell %s' % (names[j], lo, top, 1 - al, g[0], g[1]))
