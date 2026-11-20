#!/usr/bin/env python3
# dbnhigh lp.py: per-row source design LP in the float model of model.py.
# maximize F = c0 v(0) - sum_k beta_k |v(u_k)| - O(q*)   over signed weights of a dictionary of probe designs,
# subject to the boundary kernel R(d, sqrt q') <= 4/(d^2+4q') on a (d, q') grid, q' in [0, cap], and the far field d^2 R -> <= 4.
# usage: lp.py ROWS CHAIN_JSON row dict_spec [tail]
#   dict_spec: comma list of  base  |  ray  |  jets:h0:kmax[:shifts]  (k-th height-derivative designs (1/k!)(-d/dh)^k at h0, k=1..kmax)
import sys, json, math, time, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.optimize import linprog
from scipy import sparse
from model import dt_coeffs, load_rows, v_row, R_row, c0, L0, NMAX
T0 = time.time()
rows = load_rows(sys.argv[1]); ch = json.load(open(sys.argv[2]))
i = int(sys.argv[3]); spec = sys.argv[4].split(','); tail = len(sys.argv) > 5 and sys.argv[5] == 'tail'
row = rows[i]; t = row['t0']; cap = row['cap']
g = ch['rows'][str(i)]; qs = 0.5*(g['q0'] + g['q1']); ys = math.sqrt(qs)
# ---- band points
d, ln = dt_coeffs(row['t1'])
n = np.arange(2, NMAX + 1); beta = 0.5*d[n]*np.exp(-(0.5 + t*c0)*ln[n]); un = ln[n]
keep = beta > 0; n, beta, un = n[keep], beta[keep], un[keep]
NEX = int(os.environ.get('NEX', '150'))
ex = n <= NEX
U = list(un[ex]); W = list(beta[ex]); HEAD = list(n[ex] <= 4096)
def binit(lo, hi, w):
    sel = (un > lo) & (un <= hi)
    edges = np.arange(lo, hi + w, w)
    idx = np.digitize(un[sel], edges)
    for k in np.unique(idx):
        m = idx == k; bw = beta[sel][m]; U.append(float(np.sum(bw*un[sel][m])/np.sum(bw))); W.append(float(np.sum(bw))); HEAD.append(bool(np.all(n[sel][m] <= 4096)))
binit(math.log(NEX), math.log(4096), float(os.environ.get('BW1', '0.02'))); binit(math.log(4096), math.log(NMAX) + 1e-9, float(os.environ.get('BW2', '0.05')))
if tail:  # PNT continuation beyond the native cutoff (t=0: full zeta'/zeta series), weight (1/2) e^{u/2} e^{-t u (L0-u)/4} du
    for a in np.arange(math.log(NMAX), float(os.environ.get('TMAX', '40')), 0.1):
        um = a + 0.05; U.append(float(um)); W.append(float(0.5*math.exp(um/2)*math.exp(-t*um*max(L0 - um, 0)/4)*0.1)); HEAD.append(False)
U = np.array(U); W = np.array(W); HEAD = np.array(HEAD)
# ---- kernel grid
D = np.concatenate([np.arange(0, 5, 0.2), np.arange(5, 30, 1.0), np.arange(30, 301, 15.0)])
QP = [cap*f for f in (0.0, 0.25, 0.5, 0.75, 1.0)]
KD, KB, KR = [], [], []
for qp in QP:
    dd = D[D > 0] if qp == 0 else D
    KD.append(dd); KB.append(np.full(dd.shape, math.sqrt(qp))); KR.append(4/(dd*dd + 4*qp))
KD = np.concatenate(KD); KB = np.concatenate(KB); KR = np.concatenate(KR)
# ---- dictionary: each column -> (v(0), v(U), R(KD,KB), O(qs), far-field coefficient, label)
cols = []
def Kh(h, dd, b): return (h - b)/(dd*dd + (h - b)**2) + (h + b)/(dd*dd + (h + b)**2)
def add_node(h, a, lab):
    v0 = 1.0; vU = np.exp(-h*U/2)*np.cos(a*U/2)
    R = 0.5*(Kh(h, KD + a, KB) + Kh(h, KD - a, KB))
    O = 0.5*(Kh(h, a, ys) + Kh(h, -a, ys))
    cols.append((v0, vU, R, O, 2*h, lab))
def add_jet(h0, k, a, lab):
    lf = math.lgamma(k + 1)
    vU = np.exp(k*np.log(U/2) - h0*U/2 - lf)*np.cos(a*U/2)
    def Rz(dd, b):
        return 0.5*sum(np.real((h0 - b + 1j*(dd + s))**(-k - 1) + (h0 + b + 1j*(dd + s))**(-k - 1)) for s in (a, -a))
    R = Rz(KD, KB); O = float(Rz(np.zeros(1), ys)[0])
    cols.append((1.0 if k == 0 else 0.0, vU, R, O, 2*h0 if k == 0 else (-2.0 if k == 1 else 0.0), lab))
for sp in spec:
    p = sp.split(':')
    if p[0] == 'base':
        hmin = max(float(os.environ.get('HMIN', '1.2')), 1.5 - 2*t*c0 + 1e-6, math.sqrt(cap) + 1e-3)
        hstep = float(os.environ.get('HSTEP', '0.2')); astep = float(os.environ.get('ASTEP', '0.195'))
        HS = sorted(set(h for h, a, w in row['nodes'])) if os.environ.get('HSET') == 'row' else list(np.arange(math.ceil(hmin/hstep - 1e-9)*hstep, 6.01, hstep))
        for h in HS:
            if h < hmin: continue
            for a in np.arange(0, 1.951, astep): add_node(float(h), float(a), 'n%.4f/%.4f' % (h, a))
    elif p[0] == 'ray':
        v0 = float(v_row(row, np.zeros(1))[0]); vU = v_row(row, U); R = R_row(row, KD, KB); O = float(R_row(row, np.zeros(1), ys)[0])
        FF = float(1e8*R_row(row, np.array([1e4]), 0.0)[0])
        cols.append((v0, vU, R, O, FF, 'RAY'))
    elif p[0] == 'jets':
        h0 = float(p[1]); kmax = int(p[2]); shifts = [float(s) for s in p[3].split('/')] if len(p) > 3 else [0.0]
        kstep = int(p[4]) if len(p) > 4 else 1
        for a in shifts:
            for k in sorted(set([1] + list(range(kstep, kmax + 1, kstep)))): add_jet(h0, k, a, 'j%.2f/%d/%.3f' % (h0, k, a))
ENV = os.environ.get('ENV') == '1'
if ENV:  # Raymond-style native: per-column envelope charges on labels n > 4096 (and the tail); head keeps signed sums
    NAT = ~HEAD
    def envcharge(c):
        if c[5] == 'RAY':
            e = np.zeros(int(NAT.sum()))
            for h, a, w in row['nodes']: e += abs(w)*np.abs(np.exp(-h*U[NAT]/2)*np.cos(a*U[NAT]/2))
            if row['ov'] is not None:
                e *= row['ov'][0]
                for p, bs in row['ov'][1]:
                    for mm, bm in enumerate(bs): e += abs(bm*2*p**(mm - 1))*(U[NAT]/2)**mm*np.exp(-p*U[NAT]/2)/math.factorial(mm)
            return float(np.dot(W[NAT], e))
        return float(np.dot(W[NAT], np.abs(c[1][NAT])))
    NENV = np.array([envcharge(c) for c in cols])
    cols = [(c[0], c[1][HEAD], c[2], c[3], c[4], c[5]) for c in cols]
    U, W, HEAD = U[HEAD], W[HEAD], HEAD[HEAD]
nc = len(cols); nb = len(U); nk = len(KD)
scale = np.array([1/max(abs(c[0]), np.max(np.abs(c[1])), np.max(np.abs(c[2])), 1e-300) for c in cols])
V0 = np.array([c[0] for c in cols])*scale; VU = np.array([c[1] for c in cols]).T*scale; RK = np.array([c[2] for c in cols]).T*scale
OO = np.array([c[3] for c in cols])*scale; FF = np.array([c[4] for c in cols])*scale
# LP: x (nc, free), s (nb, >= 0); minimize -(c0 V0 - OO).x + W.s ; kernel rows added by cutting planes from a fine grid
def colR(j, dd, b):
    lab = cols[j][5]
    if lab == 'RAY': return R_row(row, dd, b)
    if lab[0] == 'n':
        h, a = [float(z) for z in lab[1:].split('/')]; return 0.5*(Kh(h, dd + a, b) + Kh(h, dd - a, b))
    h0, k, a = lab[1:].split('/'); h0 = float(h0); k = int(k); a = float(a)
    return 0.5*sum(np.real((h0 - b + 1j*(dd + sg))**(-k - 1) + (h0 + b + 1j*(dd + sg))**(-k - 1)) for sg in (a, -a))
Df = np.concatenate([np.arange(0, 5, 0.01), np.arange(5, 40, 0.05), np.arange(40, 400, 1.0), np.arange(400, 3000, 20.0)])
FD, FB, FR = [], [], []
for qp in np.linspace(0, cap, 9):
    dd = Df[Df > 0] if qp == 0 else Df
    FD.append(dd); FB.append(np.full(dd.shape, math.sqrt(qp))); FR.append(4/(dd*dd + 4*qp))
FD = np.concatenate(FD); FB = np.concatenate(FB); FR = np.concatenate(FR)
def fineR(xs):
    Rv = np.zeros(FD.shape)
    for j in range(nc):
        if abs(xs[j]) > 0: Rv += xs[j]*colR(j, FD, FB)
    return Rv
cvec = np.concatenate([-(c0*V0 - OO), np.ones(nb)])
Ib = sparse.identity(nb, format='csr')
VUW = VU*W[:, None]
Aband = sparse.vstack([sparse.hstack([sparse.csr_matrix(VUW), -Ib]), sparse.hstack([sparse.csr_matrix(-VUW), -Ib]),
                   sparse.hstack([sparse.csr_matrix(FF[None, :]), sparse.csr_matrix((1, nb))])]).tocsr()
bband = np.concatenate([np.zeros(2*nb), [4.0]])
if ENV:  # z_l >= |x_l| (scaled units), cost NENV_l * scale_l * z_l
    Ic = sparse.identity(nc, format='csr')
    Aband = sparse.vstack([sparse.hstack([Aband, sparse.csr_matrix((Aband.shape[0], nc))]),
                           sparse.hstack([Ic, sparse.csr_matrix((nc, nb)), -Ic]), sparse.hstack([-Ic, sparse.csr_matrix((nc, nb)), -Ic])]).tocsr()
    bband = np.concatenate([bband, np.zeros(2*nc)])
    cvec = np.concatenate([cvec, NENV*scale])
del VU
XB = float(os.environ.get('XB', '200'))
bounds = [(-XB, XB)]*nc + [(0, None)]*nb + ([(0, None)]*nc if ENV else [])
KRows = RK; KRhs = KR; good = None
for it in range(int(os.environ.get('NIT', '6'))):
    A = sparse.vstack([Aband, sparse.hstack([sparse.csr_matrix(KRows), sparse.csr_matrix((KRows.shape[0], nb + (nc if ENV else 0)))])]).tocsr()
    bvec = np.concatenate([bband, KRhs])
    print('row %d t %.5f q* %.5f cap %.4f | cols %d band %d kernel %d nnz %d | it %d %.1fs' % (i, t, qs, cap, nc, nb, KRows.shape[0], A.nnz, it, time.time() - T0), flush=True)
    res = linprog(cvec, A_ub=A, b_ub=bvec, bounds=bounds, method=os.environ.get('LPM', 'highs-ds'))
    del A; import gc; gc.collect()
    if res.status != 0:
        print('status', res.status, res.message[:80], '(keeping the last solved iteration)', flush=True)
        if good is not None: res = good
        break
    good = res
    Rv = fineR(res.x[:nc]*scale); viol = (Rv - FR)/FR
    bad = np.argsort(-viol)[:int(os.environ.get('NCUT', '200'))]; bad = bad[viol[bad] > 1e-6]
    print('  F_lp %.6f  max fine violation %.2e (%d pts)' % (-res.fun, float(np.max(viol)), len(bad)), flush=True)
    if len(bad) == 0: break
    newrows = np.array([[colR(j, FD[bad], FB[bad])[m]*scale[j] for j in range(nc)] for m in range(len(bad))]) if False else np.stack([colR(j, FD[bad], FB[bad])*scale[j] for j in range(nc)], axis=1)
    KRows = np.vstack([KRows, newrows]); KRhs = np.concatenate([KRhs, FR[bad]])
print('status', res.status, res.message[:60], '%.1fs' % (time.time() - T0), flush=True)
x = res.x[:nc]*scale; s = res.x[nc:nc + nb]/W
v0 = sum(x[j]*cols[j][0] for j in range(nc)); O = sum(x[j]*cols[j][3] for j in range(nc))
vU = sum(x[j]*cols[j][1] for j in range(nc)); band = np.abs(vU)*W
head = float(np.sum(band[HEAD])); nat = float(np.sum(band[~HEAD]))
if ENV: nat = float(np.dot(NENV, np.abs(x)))
F = c0*v0 - head - nat - O
bandgap = float(np.sum(band) - np.dot(W, s)); bandviol = float(np.max(np.abs(vU) - s))
# Raymond's design in the same model
vR = v_row(row, U); bR = np.abs(vR)*W; v0R = float(v_row(row, np.zeros(1))[0]); OR = float(R_row(row, np.zeros(1), ys)[0])
FRm = c0*v0R - float(np.sum(bR)) - OR - (float(NENV[[c[5] for c in cols].index('RAY')]) if ENV and 'RAY' in [c[5] for c in cols] else 0.0)
# fine kernel re-check of the LP design (relative to the bound)
worst = float(np.min((FR - fineR(x))/FR))
# arc-lemma check: own pair O(q) = R(0, sqrt q) of the LP design on [0, cap] (must be nondecreasing in q)
qg = np.linspace(0, cap, 201); Oq = np.array([sum(x[j]*float(np.atleast_1d(colR(j, np.zeros(1), math.sqrt(qq)))[0]) for j in range(nc) if abs(x[j]) > 0) for qq in qg])
Omono = float(np.min(np.diff(Oq)))
top = sorted(range(nc), key=lambda j: -abs(x[j])*max(abs(cols[j][0]), np.max(np.abs(cols[j][1]))))[:8]
out = dict(O_min_step=Omono, O_at_0=float(Oq[0]), O_at_cap=float(Oq[-1]), ENV=ENV, F_lp=float(-res.fun), bandgap=bandgap, bandviol=bandviol, row=i, t=t, qs=qs, spec=spec, tail=tail, status=int(res.status), F=F, main=c0*v0, head=head, native=nat, O=O, v0=v0,
           FR_model=FRm, FR_file=g['S0'] - g['O0'], kernel_worst_rel=worst, ncols=nc,
           top=[(cols[j][5], float(x[j])) for j in top], seconds=time.time() - T0)
print(json.dumps(out), flush=True)
