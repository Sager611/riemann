# explore55: CUE / function-field model of the CGG statistic. L(u) = prod(1 - a_j u), zeros u_j = conj(a_j),
# "zeta'(rho)" <-> u L'(u) at u_j, mollifier M(u) = sum_{k<=K} ((K-k)/K)^i h_k(a) u^k (1/L = sum h_k u^k), K = theta*N.
import numpy as np, sys, json
from moll import bound_complex, bound_real
def haar(N, rng):
    Z = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2)
    Q, R = np.linalg.qr(Z); d = np.diag(R); return Q*(d/abs(d))
def sample(N, K, J, rng):
    a = np.linalg.eigvals(haar(N, rng)); u = np.conj(a)
    X = 1 - np.outer(u, a); np.fill_diagonal(X, 1)
    Lp = -a*np.exp(np.sum(np.log(X), 1))                    # L'(u_j)
    p = np.array([np.sum(a**m) for m in range(1, K+1)]); h = np.zeros(K+1, complex); h[0] = 1
    for k in range(1, K+1):
        h[k] = np.dot(p[:k], h[k-1::-1])/k
    k = np.arange(K+1); xx = (K-k)/K
    if J < 0:   # optimal polynomials: BHB 1D and 2D (theta = K/N)
        th = K/N; Wt = np.stack([(-th*xx**2 + (1+th)*xx)*h, (((0.5+th)*xx - th*xx**2/2)/(0.5+th/2))*h], 1)
    else:
        Wt = np.stack([xx**i*h for i in range(1, J+1)], 1)   # (K+1, J)
    M = (u[:, None]**k[None, :]) @ Wt                       # (N, J)
    return (u*Lp)[:, None]*M
def run(N, theta, J, nzero, seed=1):
    rng = np.random.default_rng(seed); K = max(1, int(round(theta*N))); Ys = []
    while sum(len(y) for y in Ys) < nzero:
        Ys.append(sample(N, K, J, rng))
    Y = np.concatenate(Ys); r = dict(N=N, theta=theta, K=K, nz=len(Y))
    if J < 0:
        r['CGG_PBHB'] = bound_complex(Y[:, [0]]); r['2D_P2D'] = bound_real(Y[:, [1]]); r['2D_PBHB'] = bound_real(Y[:, [0]])
    for JJ in sorted({1, max(J, 1)}):
        r[f'CGG_J{JJ}'] = bound_complex(Y[:, :JJ]); r[f'2D_J{JJ}'] = bound_real(Y[:, :JJ])
    y1 = Y[:, 0]; z = y1/(y1.mean()/abs(y1.mean()))
    r['mu/N'] = float(abs(y1.mean())/N); r['A/N2'] = float(np.mean(abs(z)**2)/N**2)
    r['ReB/N2'] = float(np.mean(z**2).real/N**2); r['ImB/N2'] = float(np.mean(z**2).imag/N**2)
    r['par'] = float(np.var(z.real)/(np.var(z.real)+np.var(z.imag)))
    return r
if __name__ == '__main__':
    for N in [int(x) for x in sys.argv[1].split(',')]:
        for th in [float(x) for x in sys.argv[2].split(',')]:
            print(json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in run(N, th, int(sys.argv[3]), int(sys.argv[4])).items()}), flush=True)
