# T5 (spacing record -> band-limited statement): Landau-density-constrained Wirtinger.
# R(kappa; I, config) = inf <omega^2> over F with spectrum in I vanishing on a periodic zero set of density 1 whose max gap is kappa.
# Bloch: F = sum_k c_k e^{i w_k tau}, w_k = q + 2 pi k / P in I; constraints F(x_j) = 0 in one period; minimize over q.
import numpy as np
from scipy.linalg import null_space
PI = np.pi
def R_config(x, P, I, nq=241):
    best = np.inf
    for q in np.linspace(0, 2*PI/P, nq, endpoint=False):
        k = np.arange(int(np.floor((I[0] - q)*P/(2*PI))) - 1, int(np.ceil((I[1] - q)*P/(2*PI))) + 2)
        w = q + 2*PI*k/P; w = w[(w >= I[0] - 1e-12) & (w <= I[1] + 1e-12)]
        C = np.exp(1j*np.outer(x, w))
        Nn = null_space(C, rcond=1e-10)
        if Nn.shape[1] == 0: continue
        M = Nn.conj().T @ np.diag(w**2) @ Nn
        best = min(best, np.linalg.eigvalsh(M)[0])
    return best
print('kappa  n=P  small gap  I            R        Wirtinger (pi/kappa)^2   ratio')
for kap in [2.0, 3.0, 3.324, 4.0]:
    for n in [4, 6, 10, 20]:
        if n <= kap + 0.5: continue
        x = np.concatenate([[0.0], kap + (n - kap)*np.arange(0, n - 1)/(n - 1)])   # zeros 0, kappa, then n-2 more, all small gaps equal
        sg = (n - kap)/(n - 1)
        for name, I in [('[-pi,pi]', (-PI, PI)), ('[-1.5pi,pi]', (-1.5*PI, PI)), ('[-2pi,2pi]', (-2*PI, 2*PI))]:
            r = R_config(x, n, I)
            print(f'{kap:5.3f} {n:4d} {sg:8.3f}  {name:12s} {r:.5f}  {(PI/kap)**2:.5f}  {r/(PI/kap)**2:.4f}', flush=True)
