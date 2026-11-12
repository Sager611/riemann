# t6_cheeger.py (explore39): Cheeger and nearest-neighbour (Muckenhoupt) Markov minorants in the smallest Markov case,
# piecewise-constant cell basis on [-a, a], c <= rho (full form Markov). Z-matrix H; killing k_i = row sums (signed).
import sys, numpy as np, math, json
from scipy.integrate import quad
from scipy.special import digamma
c = float(sys.argv[1]); a = math.log(c)/2
D0 = float(digamma(0.25)) - math.log(math.pi)
kap = lambda z: math.exp(z/2)/(math.exp(z) - math.exp(-z))
K = lambda r: math.atan(math.exp(-r/2)) + 0.5*math.log(1/math.tanh(r/4))
G2 = lambda r: quad(K, r, 60, limit=200)[0]
for N in [int(x) for x in sys.argv[2:]]:
    h = 2*a/N; xs = -a + h*np.arange(N + 1)
    g2 = np.array([G2(m*h) for m in range(N + 2)])
    iu = quad(lambda u: u*kap(u), 0, h)[0]
    H = np.zeros((N, N))
    for i in range(N):
        H[i, i] = D0 + 2*K(h) + 2*iu/h
        for j in range(i + 1, N):
            m = j - i; v = -(g2[m + 1] - 2*g2[m] + g2[m - 1])/h
            H[i, j] = H[j, i] = v
    Ep = np.array([2*(math.exp(xs[i+1]/2) - math.exp(xs[i]/2)) for i in range(N)]); Em = np.array([-2*(math.exp(-xs[i+1]/2) - math.exp(-xs[i]/2)) for i in range(N)])
    H += (np.outer(Ep, Em) + np.outer(Em, Ep))/h
    off = H - np.diag(np.diag(H)); assert off.max() <= 1e-12, off.max()
    e, V = np.linalg.eigh(H); lam0 = e[0]; v = np.abs(V[:, 0])
    W = -off; k = H.sum(1); t = max(0.0, -k.min()); kt = k + t; d = float(np.max(kt + W.sum(1)))
    # sweep cuts: level sets of the Perron vector
    order = np.argsort(-v); best = np.inf; inA = np.zeros(N, bool)
    for r in range(N):
        inA[order[r]] = True
        A = np.where(inA)[0]; nA = len(A)
        cut = kt[A].sum() + W[np.ix_(A, np.where(~inA)[0])].sum() if nA < N else kt.sum()
        best = min(best, cut/nA)
    ch_lo = best**2/(2*d) - t
    # nearest-neighbour minorant (Muckenhoupt/Chen birth-death setting): keep killing, drop jumps beyond neighbours
    Hnn = np.diag(k.copy()); 
    for i in range(N - 1):
        w = W[i, i + 1]; Hnn[i, i] += w; Hnn[i + 1, i + 1] += w; Hnn[i, i + 1] = Hnn[i + 1, i] = -w
    lnn = np.linalg.eigvalsh(Hnn)[0]
    print(json.dumps({'c': c, 'N': N, 'lam0_cell': round(float(lam0), 6), 'min_killing': round(float(k.min()), 4), 'shift_t': round(t, 4),
                      'max_degree_d': round(d, 4), 'h_sweep(>=h_Cheeger)': round(float(best), 4), 'cheeger_lower<=': round(float(ch_lo), 4),
                      'buser-type_upper h-t': round(float(best - t), 4), 'lam0_nn_minorant': round(float(lnn), 4)}), flush=True)
