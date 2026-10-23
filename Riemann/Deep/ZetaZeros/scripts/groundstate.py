
# groundstate.py -- try to build periodic zero configurations that beat the paper's certificates
# (task 3). For a periodic configuration (period: p points, labels d, gaps g) we compute the
# asymptotic per-period defect
#   delta = [tr phi_c(G) - tr D^2]_per_period - (x1*n1 + x2*n2) + b*len      (distinct, c = 67/20)
#   delta = [tr phi_2(U) - n]_per_period - x*n + b*len                       (line, c = 2)
# A negative delta (beyond numerical error) would contradict the certificate inequality
# (eq:certificate / eq:certificate2) and hence the proof. Per-period values are obtained from two
# chunk lengths (boundary terms cancel).
import json, sys, time, itertools
import numpy as np
which = sys.argv[1]; seed = int(sys.argv[2]); budget = float(sys.argv[3]) if len(sys.argv)>3 else 200
rng = np.random.default_rng(seed)
A='src/anc/'
d = json.load(open(A+{'first':'certificate/inputs.json','second':'certificate_line/inputs.json'}[which]))
a = np.array(d['window_coefficients'],float)/d['denominator']
if which=='first':
    x1, x2 = 0.054350, 0.069433; b = 6*0.006219; c = 3.35; labels_allowed = (1,2)
else:
    x1 = x2 = 0.008722; b = 7*0.000627; c = 2.0; labels_allowed = (1,)
def K(t):
    s = np.sinc(t)
    for j,aj in enumerate(a,1): s = s + aj/2*(np.sinc(t-j)+np.sinc(t+j))
    return s
def phi(lam): return np.where(lam<=c, lam**2, 2*c*lam-c*c)
def chunk_energy(g, lab, M):
    p = len(g); L = g.sum()
    y = (np.concatenate([[0.0], np.cumsum(g)[:-1]])[None,:] + L*np.arange(M)[:,None]).ravel()
    dl = np.tile(np.array(lab,float), M)
    U = K(y[:,None]-y[None,:])
    G = np.sqrt(dl)[:,None]*U*np.sqrt(dl)[None,:]
    lam = np.linalg.eigvalsh(G)
    val = phi(lam).sum() - (dl**2).sum()
    n1 = (dl==1).sum(); n2 = (dl==2).sum()
    return val - (x1*n1 + x2*n2) + b*(y[-1]-y[0])
def defect(g, lab):
    p = len(g); M1 = max(4, int(60/p)); M2 = M1 + max(4, int(60/p))
    return (chunk_energy(g,lab,M2) - chunk_energy(g,lab,M1))/(M2-M1)/p   # per point
def refine(g, lab, lo, step=0.08, minstep=2e-4):
    best = defect(g, lab)
    while step > minstep:
        improved = False
        for k in rng.permutation(len(g)):
            for s in (step, -step):
                h = g.copy(); h[k] = max(lo, h[k]+s)
                v = defect(h, lab)
                if v < best - 1e-12: best, g, improved = v, h, True
        if not improved: step /= 2
    return best, g
t0 = time.time(); res = []
while time.time()-t0 < budget:
    p = int(rng.integers(1, 7))
    lab = tuple(int(rng.choice(labels_allowed)) for _ in range(p))
    lo = 0.85 if rng.random() < 0.7 else 0.25          # also allow close pairs (clipping computed)
    g = lo + rng.exponential(0.5, size=p) + (0.15 if lo>0.5 else 0.6)*rng.random(p)
    v, g = refine(g, lab, lo)
    res.append((v, lab, np.round(g,4), lo))
res.sort(key=lambda z: z[0])
print(which, 'configs tried', len(res), 'time %.0fs' % (time.time()-t0))
for v, lab, g, lo in res[:10]:
    dens = sum(lab)/g.sum(); nd = len(lab)/sum(lab)
    print('defect/point %+.4e  labels %s gaps %s  lo %.2f  (density %.3f, Nd/N %.3f)' % (v, ''.join(map(str,lab)), g, lo, dens, nd))
print('NEGATIVE DEFECT FOUND' if res and res[0][0] < -1e-6 else 'no negative defect')

