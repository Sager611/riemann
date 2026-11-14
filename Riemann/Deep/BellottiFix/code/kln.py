
import json, re, math
L = open('/Users/adriansager/playground/riemann/.scratch/bellotti2/tex/final.tex').read().split('\n')
rows = []
for ln in L[140:170]:
    if '&' in ln and '%' in ln:
        c = [x.strip().strip('$').replace('\\','') for x in ln.replace('\\\\','').split('&')]
        rows.append(c)
ind = json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/BellottiFix/logs/indep_all.json'))
tab = json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/BellottiFix/code/tables.json'))
def lT(s):
    s = s.strip(); return math.log(3e12) if s.startswith('3') else float(s[5:-1])
bounds = [(lT(r[0]), lT(r[1])) for r in tab['CB']]
P = lambda n: float(tab['CB'][n-1][3]) if 'cdot' not in tab['CB'][n-1][3] else float(tab['CB'][n-1][3].split('\\cdot')[0])*10**int(re.search(r'\{(\d+)\}|\^(\d+)', tab['CB'][n-1][3]).group(0).strip('{}^'))
for r in rows:
    s, T, imp = r[0], r[1], float(r[2].rstrip('%').replace('%',''))/100
    lt = math.log(1e13) if '10^{13}' in T or '10^13' in T else float(re.search(r'exp\(([\d.]+)\)', T).group(1))
    n = next(i+1 for i,(a,b) in enumerate(bounds) if a < lt <= b + 1e-12)
    pr = P(n); f2 = float(ind['%d_flaws2' % n]['C'])/pr; fa = float(ind['%d_all' % n]['C'])/pr
    print('sigma %s T %s row %d alpha0 %s printed impr %.1f%% -> both flaws %.1f%%, all fixes %.1f%%' % (s, T, n, tab['CB'][n-1][2], 100*imp, 100*(1-(1-imp)*f2), 100*(1-(1-imp)*fa)))
