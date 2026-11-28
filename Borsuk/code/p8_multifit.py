# Fit a multi-cell power diagram (k affine functions per colour class, cell = argmax, class = union of its cells) to a
# sampled 5-colouring; evaluate class diameters on a dense sample.  LP: each point's sub-cluster function beats every
# function of the other classes by a margin (same-class functions compete freely).
import numpy as np, sys, json
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
sys.path.insert(0, 'Borsuk/code')
from p6_pd import body, interior
from p6_survey import boundary, d2
spec = json.loads(open(sys.argv[1]).read()); colfile = sys.argv[2]; kper = int(sys.argv[3]); seed = int(sys.argv[4])
rng = np.random.default_rng(seed); rho, P, cuts = body(spec)
Z = np.load(colfile)
if len(sys.argv) > 5: Z = Z[np.random.default_rng(0).choice(len(Z), min(len(Z), int(sys.argv[5])), replace=False)]
S = Z[:, :4]; col = Z[:, 4].astype(int); n = len(S)
# sub-clusters by k-means within each class
sub = np.zeros(n, int); lab_of = []
for j in range(5):
    idx = np.flatnonzero(col == j); Y = S[idx]
    cen = Y[rng.choice(len(Y), kper, replace=False)]
    for it in range(30):
        a = np.argmin(d2(Y, cen), 1)
        for q in range(kper):
            if (a == q).any(): cen[q] = Y[a == q].mean(0)
    sub[idx] = len(lab_of) + a; lab_of += [j] * kper
K = len(lab_of); lab_of = np.array(lab_of)
def fit(margin=2e-3):
    rows, cols, vals = [], [], []; r = 0
    for i in range(n):
        f = sub[i]
        for g in np.flatnonzero(lab_of != col[i]):
            for q in range(4): rows += [r, r]; cols += [5 * g + q, 5 * f + q]; vals += [S[i, q], -S[i, q]]
            rows += [r, r, r]; cols += [5 * g + 4, 5 * f + 4, 5 * K + i]; vals += [1.0, -1.0, -1.0]; r += 1
    A = coo_matrix((vals, (rows, cols)), shape=(r, 5 * K + n)).tocsr()
    res = linprog(np.r_[np.zeros(5 * K), np.ones(n)], A_ub=A, b_ub=np.full(r, -margin), bounds=[(-50, 50)] * (5 * K) + [(0, None)] * n, method='highs')
    z = res.x[:5 * K].reshape(K, 5); return z[:, :4], z[:, 4], res.x[5 * K:]
G, c, sl = fit()
def classes(X): return lab_of[np.argmax(X @ G.T + c[None], 1)]
Xd = np.concatenate([boundary(rho, P, cuts, 30000, rng), interior(rho, P, cuts, 6000, rng)])
cl = classes(Xd); diam = []
for j in range(5):
    Y = Xd[cl == j]; m = 0.0
    for i0 in range(0, len(Y), 400): m = max(m, float(d2(Y[i0:i0 + 400], Y).max()))
    diam.append(round(np.sqrt(m), 4))
print(json.dumps({'kper': kper, 'K': K, 'misfit': int((sl > 1e-7).sum()), 'n': n, 'class_diam_dense': diam}), flush=True)
np.save('Borsuk/.scratch/multifit_%d.npy' % kper, np.c_[G, c, lab_of])
