# e1phase dls_count.py — the 2D spacing count B_a that the Bombieri–Iwaniec double large sieve needs for the phase-kept dual sum
#   kernel e(n qbar/p) exp(4 pi i sqrt(xi n/(pq))) = e(x_r . y_n),  x_r = (qbar/p, 2 (pq)^{-1/2}),  y_n = (n, sqrt(xi n))
#   B_a = #{(r,r'): |alpha-alpha'| <= 1/Y1, |x2-x2'| <= 1/Y2},  Y1 = N/2, Y2 = (sqrt(2)-1) sqrt(xi N)/2  (half-ranges)
# Block scaling as at the binding block (theta, nu): P=Q=y, N = y^{n0/nu}, X0 = y^{(1-theta)/nu}; xi from X0 = 4 pi sqrt(xi N/(PQ)).
# Control: same (pq) values, alpha replaced by uniform random reals.  dls_count2: argv[4] = 'primes' restricts p, q to primes.
import sys, math
import numpy as np
th = float(sys.argv[1]); nu = float(sys.argv[2])
rng = np.random.default_rng(2)
for y in [int(v) for v in sys.argv[3].split(',')]:
    P = Q = y; N = y**((2*nu+1-2*th)/nu); X0 = y**((1-th)/nu); xi = (X0/(4*math.pi))**2*P*Q/N
    pr, qr, al = [], [], []
    isp = lambda m: m > 1 and all(m % k for k in range(2, int(m**0.5)+1))
    rgp = [m for m in range(P, 2*P) if len(sys.argv) < 5 or isp(m)]
    for p in rgp:
        for q in rgp:
            if math.gcd(p, q) == 1: pr.append(p); qr.append(q); al.append(pow(q, -1, p)/p)
    pq = np.array(pr, float)*np.array(qr, float); al = np.array(al); x2 = 2/np.sqrt(pq); R = len(pq)
    Y1 = N/2; Y2 = (math.sqrt(2)-1)*math.sqrt(xi*N)/2
    o = np.argsort(x2); x2s = x2[o]; als = al[o]; alr = rng.random(R)[o]
    lo = np.searchsorted(x2s, x2s - 1/Y2, 'left'); hi = np.searchsorted(x2s, x2s + 1/Y2, 'right')
    Ba = Br = 0
    for i in range(R):
        Ba += int(np.count_nonzero(np.abs(als[lo[i]:hi[i]] - als[i]) <= 1/Y1))
        Br += int(np.count_nonzero(np.abs(alr[lo[i]:hi[i]] - alr[i]) <= 1/Y1))
    nb = float(np.mean(hi - lo))
    print(('primes ' if len(sys.argv) > 4 else 'all    ') + 'theta=%.3f nu=%.4f y=%d R=%d N=%.0f X0=%.0f  mean x2-neighbours=%.1f  B_a/R: arithmetic %.3f  random-alpha %.3f  (heuristic 1+nb*2/Y1 = %.3f)  X2Y2=%.2f' % (
        th, nu, y, R, N, X0, nb, Ba/R, Br/R, 1+nb*2/Y1, (1/math.sqrt(P*Q)/2)*Y2), flush=True)
