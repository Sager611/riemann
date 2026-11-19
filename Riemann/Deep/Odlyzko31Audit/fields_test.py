#!/usr/bin/env python3
"""odlaudit falsification (ii): fields with small rd against every bound of the certified odlyzko31 table, plus table sanity
(ct <= Odlyzko 1976 GRH Table 1, ct >= classical, ct < Prop D limits). Own code. Usage:
  uv run --no-project python3 fields_test.py > logs/fields_test.log
"""
import json, math, csv, bisect
SRC = '/Users/adriansager/playground/riemann/Riemann/Deep/Odlyzko31/'
cert = json.load(open(SRC + 'certified.json'))
R = [e for e in cert['results'] if e['kind'] == 'table']
tab = {w: {(e['n'], e['r1']): e['cert_rd'] for e in R if e['w'] == w} for w in ('classical', 'ct', 'ab')}
print('table sizes', {w: len(v) for w, v in tab.items()})

def best(w, N, r1):
    """best certified bound valid for a degree-N field with r1 real places: exact signature entry (n <= 12), else
    max over tabulated n' <= N of the r1 = 0 entry (every signature), and for totally real fields also of the r1 = n' entry."""
    t = tab[w]; b = t.get((N, r1), 0.0)
    for (n, s), v in t.items():
        if n <= N and (s == 0 or (r1 == N and s == n)):
            b = max(b, v)
    return b

ACT, ART = 23.3878942309, 65.6439656548   # consts_mp.log (this audit)
viol = []
# --- sanity of the table itself
odl = {}
for row in csv.reader(l for l in open(SRC + 'data/odlyzko1976_tables12.csv') if not l.startswith('#')):
    if row[0] == 'n': continue
    odl[int(row[0])] = tuple(map(float, row[1:]))   # grh_tr, grh_tc, unc_tr, unc_tc
worst_grh, worst_cls, worst_lim, cnt = 9, 9, 9, 0
for (n, r1), v in tab['ct'].items():
    lim = ART ** (r1 / n) * ACT ** ((n - r1) / n)
    worst_lim = min(worst_lim, lim / v - 1)
    if v >= lim: viol.append(('above Prop D limit', n, r1, v, lim))
    c = tab['classical'].get((n, r1))
    if c is not None:
        worst_cls = min(worst_cls, v / c - 1)
    if n in odl and r1 in (0, n):
        g = odl[n][0] if r1 == n else odl[n][1]
        cnt += 1; worst_grh = min(worst_grh, g / v - 1)
        if v > g: viol.append(('above GRH table', n, r1, v, g))
print(f'ct entries vs Prop D limit A_r^(r1/n) A_c^(2r2/n): min (limit/ct - 1) = {worst_lim:.2e}')
print(f'ct vs own classical recomputation: min (ct/classical - 1) = {worst_cls:.2e}')
print(f'ct vs Odlyzko 1976 GRH Table 1 ({cnt} common rows): min (GRH/ct - 1) = {worst_grh:.3e}')
gain = sorted(((n, r1), tab['ct'][(n, r1)] / tab['classical'][(n, r1)] - 1) for (n, r1) in tab['ct'] if r1 in (0, n) and n in (6, 7, 8, 9, 10, 12, 20, 50, 100, 1000))
print('ct/classical - 1 at selected (n, r1):', ', '.join(f'{k}:{v:.4%}' for k, v in gain))
for n in (6, 7, 8, 9, 10, 12, 20, 100, 1000):
    if n in odl:
        print(f'  n={n}: ct tc {tab["ct"].get((n, 0))} vs Odlyzko unc tc {odl[n][3]} (+{tab["ct"].get((n, 0)) / odl[n][3] - 1:.3%}); ct tr {tab["ct"].get((n, n))} vs unc tr {odl[n][2]} (+{tab["ct"].get((n, n)) / odl[n][2] - 1:.3%})')

def report(name, fields, ws=('ct',)):
    """fields: iterable of (label, N, r1, log rd)"""
    for w in ws:
        m, arg, k = 9.0, None, 0
        for lab, N, r1, lrd in fields:
            if N > 1000000: continue
            b = best(w, N, r1)
            if b <= 0: continue
            k += 1; mar = math.exp(lrd) / b - 1
            if mar < m: m, arg = mar, (lab, N, r1, round(math.exp(lrd), 5), b)
            if mar <= 0: viol.append((name, w, lab, N, r1, math.exp(lrd), b))
        print(f'{name} [{w}]: {k} fields, tightest margin rd/bound - 1 = {m:.3%} at {arg}')

# --- minimal discriminants n <= 8, every signature (from memory of the Voight / LMFDB tables), and the agent's data file
mins = {(2, 0): 3, (2, 2): 5, (3, 1): 23, (3, 3): 49, (4, 0): 117, (4, 2): 275, (4, 4): 725, (5, 1): 1609, (5, 3): 4511, (5, 5): 14641,
        (6, 0): 9747, (6, 2): 28037, (6, 4): 92779, (6, 6): 300125, (7, 1): 184607, (7, 3): 612233, (7, 5): 2306599, (7, 7): 20134393,
        (8, 0): 1257728, (8, 2): 4286875, (8, 4): 15243125, (8, 6): 65106259, (8, 8): 282300416}
lm = json.load(open(SRC + 'data/lmfdb_min.json'))
agree = [(n, n - 2 * r2, d) for n, r2, gal, lab, d in lm if not gal and (n, n - 2 * r2) in mins and mins[(n, n - 2 * r2)] != d]
print('memory list vs data/lmfdb_min.json (non-Galois minima): disagreements', agree)
report('minimal fields n<=8', [(f'|d|={d}', n, r1, math.log(d) / n) for (n, r1), d in mins.items()])
report('lmfdb_min.json', [(lab, n, n - 2 * r2, math.log(d) / n) for n, r2, gal, lab, d in lm])
# --- quadratic fields and their Hilbert class fields (degree 2h, rd = sqrt|D|)
Q = [tuple(map(int, l.split(','))) for l in open('/Users/adriansager/playground/riemann/.scratch/odlaudit/quadh.csv')]
report('imag. quadratic Hilbert class fields', [(f'D={D},h={h}', 2 * h, 0, 0.5 * math.log(-D)) for D, h in Q if D < 0])
report('real quadratic Hilbert class fields', [(f'D={D},h={h}', 2 * h, 2 * h, 0.5 * math.log(D)) for D, h in Q if D > 0])
report('quadratic fields', [(f'D={D}', 2, 0 if D < 0 else 2, 0.5 * math.log(abs(D))) for D, h in Q], ws=('ct', 'ab'))
# --- cyclotomic fields Q(zeta_m) and maximal real subfields, m <= 200000, via local conductor counts of Dirichlet characters
M = 200000
spf = list(range(M + 1))
for p in range(2, int(M ** 0.5) + 1):
    if spf[p] == p:
        for q in range(p * p, M + 1, p):
            if spf[q] == q: spf[q] = p
def local(p, e):
    """[(count, sum log f) for even, odd] over characters mod p^e"""
    ev, od = [1, 0.0], [0, 0.0]
    if p == 2:
        if e >= 2: od[0] += 1; od[1] += 2 * math.log(2)
        for k in range(3, e + 1):
            c = 2 ** (k - 3); ev[0] += c; ev[1] += c * k * math.log(2); od[0] += c; od[1] += c * k * math.log(2)
    else:
        if e >= 1:
            ev[0] += (p - 3) // 2; ev[1] += (p - 3) // 2 * math.log(p); od[0] += (p - 1) // 2; od[1] += (p - 1) // 2 * math.log(p)
        for k in range(2, e + 1):
            c = (p ** k - p ** (k - 1) - (p ** (k - 1) - p ** (k - 2))) // 2
            ev[0] += c; ev[1] += c * k * math.log(p); od[0] += c; od[1] += c * k * math.log(p)
    return ev, od
cyc, cycp = [], []
bad_check = 0
for m in range(3, M + 1):
    if m % 4 == 2: continue
    x, fac = m, {}
    while x > 1:
        p = spf[x]; fac[p] = fac.get(p, 0) + 1; x //= p
    E, O = [1, 0.0], [0, 0.0]
    for p, e in fac.items():
        ev, od = local(p, e)
        nE = [E[0] * ev[0] + O[0] * od[0], E[1] * ev[0] + E[0] * ev[1] + O[1] * od[0] + O[0] * od[1]]
        nO = [E[0] * od[0] + O[0] * ev[0], E[1] * od[0] + E[0] * od[1] + O[1] * ev[0] + O[0] * ev[1]]
        E, O = nE, nO
    phi = E[0] + O[0]
    lrd = math.log(m) - sum(math.log(p) / (p - 1) for p in fac)
    if abs((E[1] + O[1]) / phi - lrd) > 1e-9: bad_check += 1
    cyc.append((f'Q(zeta_{m})', phi, 0, lrd))
    if E[0] >= 2: cycp.append((f'Q(zeta_{m})+', E[0], E[0], E[1] / E[0]))
print('cyclotomic conductor-discriminant self-check failures:', bad_check)
report('cyclotomic Q(zeta_m), m<=2e5', cyc, ws=('ct', 'ab'))
report('real cyclotomic Q(zeta_m)+, m<=2e5', cycp, ws=('ct', 'ab'))
print('VIOLATIONS:', len(viol), viol[:10])
