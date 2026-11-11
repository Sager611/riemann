# Parse gram_win logs: window table, trend of min mu vs log t, tail fractions, candidate zero-structure classes.
import sys, re, math, glob
files = sys.argv[1:]
W = {}; H = {}; C = []
num = r'(-?[0-9.]+(?:e[-+]?[0-9]+)?)'
for fn in files:
    for line in open(fn):
        if line.startswith('WIN '):
            id_ = line.split()[1]
            kv = dict(re.findall(r'(\w+)=([^\s]+)', line))
            t0, t1 = map(float, re.search(r't=\[([^,]+),([^\]]+)\]', line).groups())
            mm = re.search(r'minmu=(\S+) n=(\S+) t=(\S+) Z=(\S+) Zp=(\S+) Zl=(\S+) Zr=(\S+) minnorm=(\S+) n=(\S+) \(mu=(\S+)\) maxabsZ=(\S+) n=(\S+) iso_bigZ=(\S+) secs=(\S+)', line).groups()
            W[id_] = dict(file=fn, n0=int(kv['n0']), cnt=int(kv['cnt']), t0=t0, t1=t1, thp=float(kv['thetap']), bad=int(kv['bad']), iso=int(kv['iso']),
                          minmu=float(mm[0]), argn=int(mm[1]), argt=float(mm[2]), Z=float(mm[3]), Zp=float(mm[4]), Zl=float(mm[5]), Zr=float(mm[6]),
                          minnorm=float(mm[7]), nnorm=int(mm[8]), maxz=float(mm[10]), bigz=int(mm[12]), secs=float(mm[13]))
        elif line.startswith('HIST '):
            id_ = line.split()[1]; parts = line.split(':')
            a = list(map(int, parts[1].split()[:9])); b = list(map(int, parts[2].split()[:9])); H[id_] = (a, b)
        elif line.startswith('CAND '):
            d = dict(re.findall(r'(\w+)=([^\s]+)', line)); d['id'] = line.split()[1]; d['ctx'] = re.search(r'Z\[-2\.\.2\]=(\S+)', line).group(1); C.append(d)
ids = sorted(W, key=lambda i: W[i]['t0'])
print(f'{len(ids)} windows; Gram points {sum(W[i]["cnt"] for i in ids):,}; isolated bad {sum(W[i]["iso"] for i in ids):,}; cpu {sum(W[i]["secs"] for i in ids):.0f} s')
print('| id | t (start) | theta\'(t) | Gram points | isolated bad | min mu | at n | min mu/theta\' | 4/theta\' | margin | max abs Z |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for i in ids:
    w = W[i]
    print(f'| {i} | {w["t0"]:.4e} | {w["thp"]:.3f} | {w["cnt"]:,} | {w["iso"]:,} | {w["minmu"]:.4f} | {w["argn"]} | {w["minnorm"]:.4f} | {4/w["thp"]:.4f} | {w["minnorm"]-4/w["thp"]:.4f} | {w["maxz"]:.1f} |')
# trend: least squares of min mu and min mu/theta' against log10 t, plus per-iso-sample correction
def lsq(xs, ys):
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    sxx = sum((x-mx)**2 for x in xs); sxy = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
    b = sxy/sxx; a = my - b*mx; r = sxy/math.sqrt(sxx*sum((y-my)**2 for y in ys)); return a, b, r
xs = [math.log10(W[i]['t0']) for i in ids]
for key in ('minmu', 'minnorm'):
    a, b, r = lsq(xs, [W[i][key] for i in ids]); print(f'TREND {key} = {a:.4f} + {b:.4f} log10 t   (r = {r:.3f})')
a, b, r = lsq(xs, [W[i]['minmu']/math.log(W[i]['iso']) for i in ids]); print(f'TREND minmu/log(iso) = {a:.4f} + {b:.4f} log10 t (r = {r:.3f})')
# tail fractions by height decade (pooled)
AT = [4.0, 4.25, 4.5, 4.75, 5.0, 5.5, 6.0, 7.0, 8.0]; NT = [0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.80]
print('POOLED tail fractions per height bin: iso, then P(mu<x) for x in', AT, ' | P(mu/theta\'<x) for x in', NT)
bins = {}
for i in ids:
    k = math.floor(2*math.log10(W[i]['t0']))/2; bins.setdefault(k, []).append(i)
for k in sorted(bins):
    iso = sum(W[i]['iso'] for i in bins[k]); A = [sum(H[i][0][j] for i in bins[k]) for j in range(9)]; Bn = [sum(H[i][1][j] for i in bins[k]) for j in range(9)]
    thp = sum(W[i]['thp'] for i in bins[k])/len(bins[k])
    print(f'10^{k:.1f}: n_win={len(bins[k])} iso={iso:,} thp={thp:.2f} abs:', ' '.join(f'{c/iso:.1e}' for c in A), '| norm:', ' '.join(f'{c/iso:.1e}' for c in Bn))
# candidates: zero structure
cls = {}
for d in C:
    zl, zr = int(d['zeros_left']), int(d['zeros_right']); key = 'zero-free block' if zl+zr == 0 else ('2 zeros one side' if {zl, zr} == {0, 2} else f'{zl}+{zr} zeros')
    cls.setdefault(key, []).append(float(d['norm']))
print('CAND classes (norm < NTHR or mu < THR): count, min norm:', {k: (len(v), round(min(v), 4)) for k, v in cls.items()})
clsm = {}
for d in C:
    zl, zr = int(d['zeros_left']), int(d['zeros_right']); key = 'zero-free' if zl+zr == 0 else ('confined pair' if {zl, zr} == {0, 2} else f'{zl}+{zr}')
    clsm.setdefault(key, []).append((float(d['mu']), d['n'], d['id']))
print('CAND classes min abs mu:', {k: (len(v), min(v)) for k, v in clsm.items()})
C.sort(key=lambda d: float(d['mu']))
print('LOWEST mu candidates:')
for d in C[:15]:
    print(d['id'], d['n'], d['t'], 'mu', d['mu'], 'norm', d['norm'], 'Z', d['Z'], 'Z[-2..2]', d['ctx'], 'zl/zr', d['zeros_left'], d['zeros_right'])
C.sort(key=lambda d: float(d['norm']))
print('LOWEST normalised candidates:')
for d in C[:10]:
    print(d['id'], d['n'], d['t'], 'mu', d['mu'], 'norm', d['norm'], 'zl/zr', d['zeros_left'], d['zeros_right'])

# ---- tail model: log P(mu < x | isolated bad) = a + b x (+ c log10 t), Poisson-weighted LSQ on HIST counts x in XS ----
XS = [4.75, 5.0, 5.5, 6.0]
rows = []
for i in ids:
    A = H[i][0]
    for x in XS:
        j = AT.index(x); rows.append((x, math.log10(W[i]['t0']), A[j], W[i]['iso']))
def fit(rows, use_t):
    # weighted least squares on y = log(c/I) with weight c (bins pooled per (x, half-decade) to avoid zeros)
    pool = {}
    for x, lt, c, I in rows:
        k = (x, math.floor(2*lt)/2); p = pool.setdefault(k, [0, 0, 0.0]); p[0] += c; p[1] += I; p[2] += lt*I
    pts = [(x, p[2]/p[1], math.log(p[0]/p[1]), p[0]) for (x, _), p in pool.items() if p[0] > 0]
    import itertools
    # normal equations
    cols = (lambda x, lt: [1.0, x, lt]) if use_t else (lambda x, lt: [1.0, x])
    m = len(cols(0, 0)); M = [[0.0]*m for _ in range(m)]; v = [0.0]*m
    for x, lt, y, w in pts:
        f = cols(x, lt)
        for a in range(m):
            v[a] += w*f[a]*y
            for b in range(m): M[a][b] += w*f[a]*f[b]
    # solve
    for c in range(m):
        piv = M[c][c]
        for r in range(c+1, m):
            q = M[r][c]/piv
            for k in range(c, m): M[r][k] -= q*M[c][k]
            v[r] -= q*v[c]
    sol = [0.0]*m
    for c in range(m-1, -1, -1):
        sol[c] = (v[c] - sum(M[c][k]*sol[k] for k in range(c+1, m)))/M[c][c]
    return sol, pts
sol0, pts = fit(rows, False); sol1, _ = fit(rows, True)
print(f'TAIL fit (pooled heights): ln P(mu<x) = {sol0[0]:.3f} + {sol0[1]:.3f} x  ->  P(mu<4) = {math.exp(sol0[0]+4*sol0[1]):.2e}, P(mu<4.2155) = {math.exp(sol0[0]+4.2155*sol0[1]):.2e}')
print(f'TAIL fit with height: ln P(mu<x) = {sol1[0]:.3f} + {sol1[1]:.3f} x + {sol1[2]:.3f} log10 t')
for lt in (8, 9, 10, 11, 11.7):
    print(f'   at t=10^{lt}: P(mu<4) = {math.exp(sol1[0]+4*sol1[1]+sol1[2]*lt):.2e}')

# ---- min mu against log10 t controlling for sample size: min mu = a + b log10 t + c ln(iso) ----
def lsqn(X, y):
    m = len(X[0]); M = [[sum(r[a]*r[b] for r in X) for b in range(m)] for a in range(m)]; v = [sum(r[a]*yy for r, yy in zip(X, y)) for a in range(m)]
    for c in range(m):
        for r in range(c+1, m):
            q = M[r][c]/M[c][c]
            for k in range(c, m): M[r][k] -= q*M[c][k]
            v[r] -= q*v[c]
    sol = [0.0]*m
    for c in range(m-1, -1, -1): sol[c] = (v[c] - sum(M[c][k]*sol[k] for k in range(c+1, m)))/M[c][c]
    res = [yy - sum(s*xx for s, xx in zip(sol, r)) for r, yy in zip(X, y)]
    return sol, math.sqrt(sum(e*e for e in res)/(len(y)-m))
X = [[1.0, math.log10(W[i]['t0']), math.log(W[i]['iso'])] for i in ids]
for key in ('minmu', 'minnorm'):
    sol, sd = lsqn(X, [W[i][key] for i in ids])
    print(f'TREND2 {key} = {sol[0]:.3f} + {sol[1]:.4f} log10 t + {sol[2]:.4f} ln(iso)   (resid sd {sd:.3f}); at iso=1e5: t=1e8 -> {sol[0]+8*sol[1]+sol[2]*math.log(1e5):.3f}, t=1e11 -> {sol[0]+11*sol[1]+sol[2]*math.log(1e5):.3f}, t=5e11 -> {sol[0]+math.log10(5e11)*sol[1]+sol[2]*math.log(1e5):.3f}')
