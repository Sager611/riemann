# explore38 s5: (i) out-of-sample test: predict C2(N) = L^2 d_N^2 - C L at N = 2560, 5120 from B(N) + flat deficit, vs 2611.03924's V2;
# (ii) regression c2o ~ a + beta*c2B on [100,1280]; (iii) asymptotics of c2B from s3 runs at large N.
import numpy as np, math
OUT = '/Users/adriansager/playground/riemann/.scratch/explore38/'
EG = 0.57721566490153286061; C = 2 + EG - math.log(4 * math.pi)
Z = np.load(OUT + 'compare.npz'); N = Z['N']; L = Z['L']; c2B = Z['c2B']; c2o = Z['c2o']
m = N >= 100
dl = c2o - c2B
print(f'deficit on [100,1280]: mean {dl[m].mean():.5f} sd {dl[m].std():.5f}; fit a + b/L: ', np.polyfit(1 / L[m], dl[m], 1)[::-1].round(5))
X = np.vstack([np.ones(m.sum()), c2B[m]]).T
co = np.linalg.lstsq(X, c2o[m], rcond=None)[0]; res = c2o[m] - X @ co
print(f'regression c2o = {co[0]:+.5f} + {co[1]:.4f} c2B, residual sd {res.std():.5f} (raw sd of c2o {c2o[m].std():.5f})')
X2 = np.vstack([np.ones(m.sum()), c2B[m], 1 / L[m]]).T
co2 = np.linalg.lstsq(X2, c2o[m], rcond=None)[0]; res2 = c2o[m] - X2 @ co2
print(f'regression c2o = {co2[0]:+.5f} + {co2[1]:.4f} c2B + {co2[2]:+.4f}/L, residual sd {res2.std():.5f}')
big = np.load(OUT + 'scan_big600.npy')
for row in big:
    n_, L_, BK, Bd, td = row
    B = BK + td
    cb = L_ * L_ * B - C * L_
    pred_flat = cb + dl[m].mean()
    pred_reg = co[0] + co[1] * cb
    pred_reg2 = co2[0] + co2[1] * cb + co2[2] / L_
    tag = ''
    if int(round(n_)) == 2560: tag = '  2611.03924 V2: +0.0028 (V2-V1 ~ +0.0003)'
    if int(round(n_)) == 5120: tag = '  2611.03924 V2: -0.0091'
    print(f'N={n_:.4g} L={L_:.3f} B L={B*L_:.6f} c2B={cb:+.5f} pred C2 (flat deficit) {pred_flat:+.5f} (2-par regression) {pred_reg:+.5f} (3-par) {pred_reg2:+.5f}{tag}')
