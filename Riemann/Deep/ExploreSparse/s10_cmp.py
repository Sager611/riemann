# explore61 s10: measured kappa(theta) = L (pi_N^2 - om^2) (s3_rate.npz) vs zero-formula prediction (logs/s5_zero_3000_40k.log).
import re, math, numpy as np
z = np.load('logs/s3_rate.npz'); om2 = z['om2']
TH = ['2/3','3/4','3/5','4/5','5/6','4/7','5/7','5/8','2/5','3/7','3/8','2/7','2/9','3/10']
pred = {m.group(1): float(m.group(2)) for m in re.finditer(r'^(\d+/\d+):.*-> ([0-9.]+)', open('logs/s5_zero_3000_40k.log').read(), re.M)}
Ns = np.arange(100, 2001); L = np.log(Ns)
print('theta  kappa(N=160)  kappa(640)  kappa(2000)  fit k_inf + b/L (N>=100)   predicted   ratio k(2000)/pred  k_inf/pred')
for i, t in enumerate(TH):
    pi2 = z[t.replace('/', '_')]
    k = np.log(np.arange(1, len(pi2) + 1)) * (pi2 - om2[i])
    A = np.vstack([np.ones_like(L), 1 / L]).T; c = np.linalg.lstsq(A, k[99:2000], rcond=None)[0]
    print(f'{t:5s} {k[159]:.5f} {k[639]:.5f} {k[1999]:.5f}   {c[0]:.5f} {c[1]:+.5f}   {pred[t]:.5f}   {k[1999]/pred[t]:.4f}  {c[0]/pred[t]:.4f}')
