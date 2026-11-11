# exploratory (floating point, untrusted): where can RH-to-T and CH help?
import json, math
R = json.load(open('/Users/adriansager/playground/riemann/Riemann/Deep/Mertens78/results.json'))
LN10 = math.log(10)
def lse(*v):
    m = max(v); return m + math.log(sum(math.exp(x - m) for x in v))
certs = []
for d in R:
    certs.append(dict(name=d['name'], s=float(d['sigma1']), a=d['a'], lC=d['logCinf'], J0=d['J0'],
        th=d['theta'], lCa=d['best']['log10C']*LN10, xa=d['best']['log10xa']*LN10, lCp=d['allx']['Cprime_log10']*LN10))
Tp = 3.0e12
import sys
BETA = 1/4 if 'nosqf' in sys.argv else 3/(2*math.pi**2)  # T' <= Platt-Trudgian height (refine later)
def tailint(la, tau, T):
    # int_T^inf t^a min(1/t, tau/t^2) dt
    if tau <= T: return math.log(tau) + (la-1)*math.log(T) - math.log(1-la)
    return math.log((tau**la - T**la)/la + tau**la/(1-la)) if tau < 1e300 else la*math.log(tau) + math.log(1/la + 1/(1-la))
def m78(c, L, lh):
    lX = L + math.log1p(math.exp(lh - L)/2); ltau = math.log(2) + lX - lh; a = c['a']
    return lse(lh - math.log(4), 0.0, c['s']*lX + math.log(c['J0']/(2*math.pi)),
               c['s']*lX + c['lC'] + math.log(2*(1/a + 1/(1-a))/(2*math.pi)) + a*ltau)
def hyb(c, L, lh, T=Tp):
    lX = L + math.log1p(math.exp(lh - L)/2); ltau = math.log(2) + lX - lh; a = c['a']
    tau = math.exp(min(ltau, 690))
    ti = tailint(a, tau, T) if ltau < 690 else a*ltau + math.log(1/a + 1/(1-a))
    return lse(lh + math.log(BETA), math.log(0.0268) + lX/2, c['s']*lX + c['lC'] + math.log(2/(2*math.pi)) + ti)
def opt(f, c, L):
    lo, hi = 0.5*L, L - 0.7
    best = None
    for k in range(200):
        lh = lo + (hi - lo)*k/199
        v = f(c, L, lh)
        if best is None or v < best[0]: best = (v, lh)
    # refine
    lo, hi = best[1] - (hi-lo)/199, best[1] + (hi-lo)/199
    for k in range(60):
        lh = lo + (hi-lo)*k/59; v = f(c, L, lh)
        if v < best[0]: best = (v, lh)
    return best
def ch(L): return lse(math.log(3/(math.pi*1e10)) + L, math.log(11.39) + L/2)
def prev(L):
    v = [L, ch(L)]
    for c in certs:
        v.append(c['lCp'] + c['th']*L)
        if L >= c['xa']: v.append(c['lCa'] + c['th']*L)
    return min(v)
import sys
thetas = [0.9409, 0.94, 0.9382, 0.937, 0.935, 0.9325, 0.93, 0.9296, 0.9277, 0.925, 0.92, 0.9192, 0.915, 0.9121, 0.91, 0.9092, 0.905, 0.9024, 0.9, 0.8992]
Ls = [LN10*k for k in range(1, 4500, 3)] + [LN10*10**(3.66 + 0.01*k) for k in range(300)]
Ls2 = [LN10*k for k in range(1, 4500, 1)]
# precompute per-L bests
rows = []
for L in Ls:
    m = min(opt(m78, c, L)[0] for c in certs)
    hy = min(opt(hyb, c, L)[0] for c in certs)
    rows.append((L, prev(L), ch(L), m, hy))
def sup(th, idx):
    return max(min(r[i] for i in idx) - th*r[0] for r in rows)/LN10
fam = lambda th: min(c['lCp']/LN10 for c in certs if c['th'] <= th + 1e-12) if any(c['th'] <= th + 1e-12 for c in certs) else None
def argsup(th, idx):
    v = [(min(r[i] for i in idx) - th*r[0], r[0]/LN10) for r in rows]; return max(v)[1]
print('theta  family_row  prev(rows+CH)  CH+M78opt  CH+M78opt+hybrid   (log10 c)')
for th in thetas:
    print(th, fam(th), round(sup(th, [1]), 3), round(sup(th, [1, 3]), 3), round(sup(th, [1, 3, 4]), 3), 'at x=1e%.0f' % argsup(th, [1,3,4]))
