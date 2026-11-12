# an_det.py: tables and growth exponents from det_*.out
import re, math, glob, collections
rows = []
for fn in sorted(glob.glob('det_*.out'), key=lambda s: int(re.findall(r'\d+', s)[0])):
    X = int(re.findall(r'\d+', fn)[0]); fam = None
    for line in open(fn):
        if line.startswith('# X='):
            fam = int(re.search(r'famsize=(\d+)', line).group(1)); continue
        if line.startswith('#') or not line.strip(): continue
        p = line.split(); N, sig, V, tp, kill, der = int(p[0]), float(p[1]), float(p[2]), float(p[3]), int(p[4]), int(p[5])
        d = dict(x.split('=') for x in p[6:]); d = {k: float(v) for k, v in d.items()}
        rows.append(dict(X=X, fam=fam, N=N, sig=sig, V=V, tp=tp, kill=kill, der=der, **d))
def key(r):
    return (r['sig'], 'V=2' if r['V'] == 2 else 'V=sqrtX/2', r['tp'], r['kill'], r['der'])
names = {(0.8333, 'V=sqrtX/2', 0.0, 1, 0): 'main s=5/6 killed', (0.8333, 'V=sqrtX/2', 0.0, 0, 0): 's=5/6 POLE',
         (0.55, 'V=sqrtX/2', 0.0, 1, 0): 's=0.55 killed', (0.55, 'V=2', 0.0, 1, 0): 's=0.55 V=2 killed',
         (0.8333, 'V=sqrtX/2', 5.0, 1, 0): "s=5/6 t'=5 killed", (0.8333, 'V=sqrtX/2', 0.0, 1, 1): 's=5/6 d/dt killed',
         (0.9, 'V=2', 0.0, 0, 0): 's=0.9 V=2 POLE'}
by = collections.defaultdict(list)
for r in rows: by[key(r)].append(r)
for k, rs in by.items():
    print(f"\n== {names.get(k, k)}   R1 = LHS/(famsize*sum|b|^2), R2 = LHS/(X^1.5 N^(1-2s)), triv/gen shares, R1ng = (LHS-GEN)/(famsize*sum|b|^2)")
    print(" X      N      N/X     R1      R2      triv%   gen%    R1ng")
    for r in sorted(rs, key=lambda r: (r['X'], r['N'])):
        R1 = r['LHS'] / (r['fam'] * r['NB']); R2 = r['LHS'] / (r['X'] ** 1.5 * r['N'] ** (1 - 2 * r['sig']))
        R1ng = (r['LHS'] - r['GEN']) / (r['fam'] * r['NB'])
        print(f"{r['X']:6d} {r['N']:6d} {r['N']/r['X']:7.4f} {R1:7.3f} {R2:8.4f} {100*r['TRIV']/r['LHS']:6.1f} {100*r['GEN']/r['LHS']:6.1f} {R1ng:7.3f}")
    # exponent in N at fixed X over N in [X^0.5, X^0.95]: slope of log R1
    for X in sorted(set(r['X'] for r in rs)):
        pts = [(math.log(r['N']), math.log(r['LHS'] / (r['fam'] * r['NB']))) for r in rs if r['X'] == X and X ** 0.5 <= r['N'] <= X ** 0.95]
        if len(pts) >= 3:
            mx = sum(a for a, _ in pts) / len(pts); my = sum(b for _, b in pts) / len(pts)
            sl = sum((a - mx) * (b - my) for a, b in pts) / sum((a - mx) ** 2 for a, _ in pts)
            print(f"   X={X}: slope d log R1 / d log N on [X^0.5, X^0.95] = {sl:+.3f} ({len(pts)} pts)")
    # exponent in X at fixed N/X (nearest power of 2 to X/4 and X/16)
    for frac in (1/64, 1/16, 1/4, 1.0):
        pts = []
        for X in sorted(set(r['X'] for r in rs)):
            cand = [r for r in rs if r['X'] == X]
            r = min(cand, key=lambda r: abs(math.log(r['N'] / (frac * X))))
            if abs(math.log(r['N'] / (frac * X))) < 0.75: pts.append((math.log(X), math.log(r['LHS'] / (r['fam'] * r['NB']))))
        if len(pts) >= 3:
            mx = sum(a for a, _ in pts) / len(pts); my = sum(b for _, b in pts) / len(pts)
            sl = sum((a - mx) * (b - my) for a, b in pts) / sum((a - mx) ** 2 for a, _ in pts)
            print(f"   N ~ X*{frac:.4f}: slope d log R1 / d log X = {sl:+.3f} ({len(pts)} pts)")
