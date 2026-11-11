# compare.py -- comparison table (floating point; the claims are in logs/cert_ch_results*.json, certified in Arb).
# record row : smallest published Mertens78 all-x constant C' among rows with theta_i <= theta
# baseline   : sup_x min(x, CH(x), C'_i x^theta_i, C_i x^theta_i [x >= x_a,i]) / x^theta  (published bounds combined)
import json, math, sys
LN10 = math.log(10)
rows = json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/Mertens78/results.json'))
R0 = [(r['theta'], r['allx']['Cprime_log10'], r['best']['log10C'], r['best']['log10xa']) for r in rows]
import os
RB = [(r['theta'], r['allx']['Cprime_log10'], r['best']['log10C'], r['best']['log10xa'])
      for r in (json.load(open('logs/bprime_results.json')) if os.path.exists('logs/bprime_results.json') else [])]
R = R0 + RB if '--with-bp-rows' in sys.argv else R0
def prevL(L):     # log10 of the best published bound at x = 10^L
    v = [L, math.log10(3/(math.pi*1e10)*10**min(L, 300)) if L < 300 else math.log10(3/(math.pi*1e10)) + L]
    ch = math.log10(3/(math.pi*1e10)) + L
    ch = ch + math.log10(1 + 11.39*math.pi*1e10/3*10**(-L/2)) if L < 200 else ch
    v = [L, ch]
    for th, cp, ca, xa in R:
        v.append(cp + th*L)
        if L >= xa: v.append(ca + th*L)
    return min(v)
Ls = [k*0.01 for k in range(0, 500000)] + [5000*1.0005**k for k in range(1, 15000)]
P = [(L, prevL(L)) for L in Ls]
def base(th): return max(p - th*L for L, p in P)
def rowenv(th):
    c = [cp for t, cp, _, _ in R0 if t <= th + 1e-12]; return min(c) if c else None
files = [f for f in sys.argv[1:] if not f.startswith('--')] or ['logs/cert_ch_results.json']
# crossover check: CH vs the record row 1.84e6 x^0.9382 (TREE N10 says 10^263.5)
L = 200.0
while math.log10(3/(math.pi*1e10)) + L < math.log10(1.84e6) + 0.9382*L: L += 0.01
print("crossover CH vs 1.84e6 x^0.9382: x = 10^%.1f" % L)
cols = []
for f in files:
    d = {}
    for r in json.load(open(f)):
        d[round(r['theta'], 4)] = (r['log10c'], r['c'])
    cols.append(d)
ths = sorted(set().union(*[set(d) for d in cols]), reverse=True)
def fmt(c):
    m, e = c.split('e'); m = int(m); e = int(e)
    while m >= 1000: m = (m + 9)//10; e += 1
    k = len(str(m)) - 1; return "%s.%s·10^%d" % (str(m)[0], str(m)[1:], e + k)
hdr = "| θ | record row C′ | published bounds combined | " + " | ".join("new (%s)" % f.split('results')[-1].strip('_.json') or "7/8" for f in files) + " | gain vs row | gain vs combined |"
print(hdr); print("|" + "---|"*(5 + len(files)))
for th in ths:
    re = rowenv(th); b = base(th)
    vals = [d.get(th) for d in cols]
    best = min(v[0] for v in vals if v)
    print("| %.4f | %s | 10^%.2f | %s | %s | %s |" % (th, ("10^%.2f" % re) if re is not None else "-", b,
          " | ".join(fmt(v[1]) if v else "-" for v in vals),
          ("%.3gx" % 10**(re - best)) if re is not None else "-", ("%.3gx" % 10**(b - best))))
