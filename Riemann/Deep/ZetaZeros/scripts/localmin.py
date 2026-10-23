
# localmin.py -- numerical search for frames violating the local inequalities (eq:local, eq:local8)
import json, itertools, sys, time
import numpy as np
A = 'src/anc/'
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv)>2 else 1)
which = sys.argv[1] if len(sys.argv)>1 else 'first'
path = {'first':'certificate/inputs.json','second':'certificate_line/inputs.json','third':'certificate_third/inputs.json'}[which]
d = json.load(open(A+path)); D = d['denominator']
a = np.array(d['window_coefficients'],float)/D
npts = d.get('points',7); ng = npts-1
alpha = d['alpha']/D
W = np.array(d['weights'],float)/D          # [lab][pair]
kap = np.array(d['h'],float)/D              # [label-1][site]
cob = np.array(d['cob_knots'],float)/D      # [r][lab][knot]
pairs = [(i,j) for i in range(npts) for j in range(i+1,npts)]
assert W.shape[1]==len(pairs), (W.shape, len(pairs))
knots = 0.75 + np.arange(cob.shape[2])/8.0
theta = 0.85
def K(t):
    s = np.sinc(t)
    for j,aj in enumerate(a,1): s = s + aj/2*(np.sinc(t-j)+np.sinc(t+j))
    return s
def B(r, lab, g):
    return np.interp(g, knots, cob[r][lab], left=0.0, right=0.0)
def lab2(di,dj): return (di-1)+(dj-1)
labelled = which=='first'
labelings = list(itertools.product([1,2],repeat=npts)) if labelled else [tuple([1]*npts)]
def slack(G, dl):
    # G: (m, ng) gaps
    y = np.concatenate([np.zeros((G.shape[0],1)), np.cumsum(G,axis=1)],axis=1)
    val = alpha*G.sum(axis=1)
    for p,(i,j) in enumerate(pairs):
        val = val + W[lab2(dl[i],dl[j])][p]*K(y[:,j]-y[:,i])**2
    for r in range(ng):
        val = val + B(r, lab2(dl[r],dl[r+1]), G[:,r])
    return val - sum(kap[dl[i]-1][i] for i in range(npts))
def refine(g, dl, step=0.05, minstep=1e-6):
    best = slack(g[None,:], dl)[0]
    while step > minstep:
        # all +-step moves on each coordinate, plus random directions
        cand = [g + s*e for e in np.eye(ng) for s in (step,-step)]
        dirs = rng.normal(size=(2*ng, ng)); dirs /= np.linalg.norm(dirs,axis=1)[:,None]
        cand += [g + step*v for v in dirs]
        C = np.maximum(np.array(cand), theta)
        v = slack(C, dl); k = np.argmin(v)
        if v[k] < best - 1e-15: best, g = v[k], C[k]
        else: step /= 2
    return best, g
t0=time.time(); results=[]
for dl in labelings:
    m = 40000 if labelled else 400000
    G = theta + rng.exponential(0.35, size=(m, ng))
    G[: m//4] = theta + rng.uniform(0, 1.6, size=(m//4, ng))
    v = slack(G, dl)
    idx = np.argsort(v)[:12 if labelled else 60]
    best = (np.inf, None)
    for k in idx:
        b, g = refine(G[k].copy(), dl)
        if b < best[0]: best = (b, g)
    results.append((best[0], dl, best[1]))
results.sort(key=lambda z: z[0])
print(which, 'labelings', len(labelings), 'time %.1fs' % (time.time()-t0))
for b, dl, g in results[:8]:
    print('min slack %.3e  labels %s  gaps %s' % (b, ''.join(map(str,dl)), np.round(g,4)))
print('NEGATIVE FOUND' if results[0][0] < 0 else 'no negative slack found')

