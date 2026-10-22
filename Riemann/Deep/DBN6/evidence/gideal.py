# "Ideal" physical derivative bound (non-rigorous): G0(t,h) = (1/4) sum_{2<=n<=N} d_t(n) log n n^{-sigma_B},
# where d_t = D log_* b_t (formal convolution log-derivative of the heat-weighted series; d_t >= 0 by Raymond's lemma),
# sigma_B = (1+h)/2 + t L/4, N = 690988 (native cutoff at X0 = 6000000185827).  Phase-aligned worst case, first branch only.
import numpy as np, mpmath as mp, json, sys
N = 690988
L = float(mp.log(mp.mpf(6000000185827)/(4*mp.pi)))
n = np.arange(N+1, dtype=np.float64); ln = np.zeros(N+1); ln[1:] = np.log(n[1:])
def dcoef(t):
    b = np.exp(t*ln*ln/4)
    d = np.zeros(N+1); acc = np.zeros(N+1)
    for k in range(2, N+1):
        dk = ln[k]*b[k] - acc[k]
        d[k] = dk
        if dk != 0.0 and 2*k <= N:
            m = N//k
            acc[2*k::k][:m-1] += dk*b[2:m+1]
    return d
out = {}
ts = [0.0, 0.01, 0.02, 0.03, 0.045, 0.06, 0.075, 0.09, 0.105, 0.12]
hs = [0.65, 0.75, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9]
for t in ts:
    d = dcoef(t)
    row = []
    for h in hs:
        sB = (1+h)/2 + t*L/4
        row.append(0.25*float(np.sum(d[2:]*ln[2:]*np.exp(-sB*ln[2:]))))
    out[t] = row
    print('t=%.3f ' % t + ' '.join('h=%.2f:%.3f' % (h, g) for h, g in zip(hs, row)), flush=True)
json.dump({'ts': ts, 'hs': hs, 'G': [out[t] for t in ts]}, open('gideal.json', 'w'))
