# Audit numerics on lambda_n, n <= 1e5: FFT table (rerun of LiRate fft_lambda.py) vs independent zeros-based values,
# and the quantitative claims c1, 0.289, kappa3, 7/4 (Thm 4 windows), Thm 2 window sums, Cor 4.
# usage: python3 audit_num.py lam_r999.txt lamZ.npy trend_r999.txt
import sys, numpy as np
from math import log, sqrt, pi, exp
EG = 0.57721566490153286061
lam = np.zeros(100001)
for l in open(sys.argv[1]):
    if l[0] == '#': continue
    p = l.split(); lam[int(p[0])] = float(p[1])
Z = np.load(sys.argv[2]); lamZ, errZ = Z[0], Z[1]
d = np.abs(lam - lamZ)
print('A. FFT vs zeros: max|diff| = %.2e (n=%d), max |diff|/bound = %.3f, max|diff| n<=1000 = %.2e' % (d.max(), d.argmax(), (d[1:]/errZ[1:]).max(), d[:1001].max()))
n = np.arange(100001, dtype=float); n[0] = 1
kap = 1 + log(2*pi) - EG
t = 0.5*n*(np.log(n) - kap); t[0] = 0
R = lam - t; R[0] = 0
c1 = sqrt(3/5)*(4/5)**2/(2*pi*sqrt(2))
cs = np.cumsum(R**2)
print('B1. c1=%.5f; sum_{n<=N} R^2/(c1 N^1.5 sqrt(log N)):' % c1, ' '.join('N=%d:%.3f' % (N, cs[N]/(c1*N**1.5*sqrt(log(N)))) for N in (3, 10, 30, 100, 300, 1000, 3000, 10**4, 3*10**4, 10**5)))
r4 = np.abs(R[2:])/(n[2:]*np.log(n[2:]))**0.25
print('B2. max |R|/(n log n)^(1/4) = %.4f at n=%d; per decade max:' % (r4.max(), r4.argmax()+2), ' '.join('%d-%d:%.3f' % (a, 10*a, r4[a-2:10*a-1].max()) for a in (10, 100, 1000, 10**4)), '; frac n in [5e4,1e5] with ratio>0.289: %.3f' % (r4[5*10**4-2:].__gt__(0.289).mean()))
k3 = (2/(3*pi))*(pi**2/12)**0.25
m = np.arange(100, 100001)
q = np.abs(R[m])/(k3*np.sqrt(m)*np.log(m)/np.sqrt(np.log(np.log(m))))
print('B3. kappa3=%.4f; max_{100<=n<=1e5} |R|/(kappa3 sqrt n log n/sqrt(loglog n)) = %.3f at n=%d; max over n in [16,100): %.3f' % (k3, q.max(), m[q.argmax()], max(abs(R[j])/(k3*sqrt(j)*log(j)/sqrt(log(log(j)))) for j in range(16, 100))))
def bump(x, a, b):
    y = np.zeros_like(x); s = (x > a) & (x < b); u = (2*x[s] - a - b)/(b - a); y[s] = np.exp(-1/(1 - u*u)); return y
print('B4. Thm 4: smoothed mean of R minus 7/4, and (dev)*N/log N:')
for (a, b, wt) in ((0.5, 1, 0), (0.9, 1, 0), (0.98, 1, 0), (0.1, 1, 0), (0.3, 1, 5), (0.001, 0.002, 0), (0.01, 0.011, 0)):
    out = []
    for N in (10**3, 10**4, 10**5):
        x = np.arange(1, 100001)/N; w = bump(x, a, b)*(x**wt)
        if w.sum() == 0: continue
        dev = (w*R[1:]).sum()/w.sum() - 1.75; out.append('N=%d dev=%+.2e dev*N/logN=%+.2f' % (N, dev, dev*N/log(N)))
    print('  (%g,%g) x^%d: ' % (a, b, wt) + '; '.join(out))
for N in (10**3, 10**4, 10**5):
    print('  sharp [1,N] mean R - 7/4 = %+.4f (N=%d)' % (R[1:N+1].mean() - 1.75, N))
print('B5. Thm 2: S2(N) = sum_{N<n<=2N} R^2 / (N^1.5 sqrt(log N)) [and / loglog N]:', ' '.join('N=%d:%.3f[%.3f]' % (N, (cs[2*N]-cs[N])/(N**1.5*sqrt(log(N))), (cs[2*N]-cs[N])/(N**1.5*sqrt(log(N))*log(log(N)))) for N in (100, 500, 2000, 5000, 20000, 50000)))
trend = {}
for l in open(sys.argv[3]):
    if l[0] == '#': continue
    p = l.split(); trend[int(p[0])] = float(p[1])
lb = np.array([trend.get(j, 0.0) for j in range(100001)])
lt = lam - lb
bad = [j for j in range(2, 100001) if abs(lt[j]) > 2*log(j)]
print('B6. Cor 4: lambda_bar - t - 3/4 + 1/(24n) at n=10,1e5: %.2e %.2e; first n>=2 with |lt|>2log n: %d (lt=%.6f, 2log n=%.6f); count n in [2,1e5]: %d; max_{30<=n} |lt|/log n = %.4f at n=%d' % (lb[10]-t[10]-0.75+1/240, lb[10**5]-t[10**5]-0.75+1/(24e5), bad[0], lt[bad[0]], 2*log(bad[0]), len(bad), max(abs(lt[j])/log(j) for j in range(30, 100001)), max(range(30, 100001), key=lambda j: abs(lt[j])/log(j))))
for a in (3, 4):
    print('  a=%d first fails at n=%d' % (a, next(j for j in range(2, 100001) if abs(lt[j]) > a*log(j))))
np.save(sys.argv[2].replace('lamZ', 'R_fft'), R)
