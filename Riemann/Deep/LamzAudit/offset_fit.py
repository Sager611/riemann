# Reads logs/num_pairs_zeros{1,6}.log; tabulates J - S_g0(Z)/N against 1/log T and the GUE+local-density model.
import re, math
import numpy as np
rows = []
for f in ('logs/num_pairs_zeros1.log', 'logs/num_pairs_zeros6.log'):
    for blk in re.split(r'\nT=', open(f).read())[1:]:
        T = float(blk.split()[0])
        obs = float(re.search(r'S_g0\(Z\)/N=([0-9.]+)', blk).group(1))
        mod = float(re.search(r'model S_g0\(Z\)/N=([0-9.]+)', blk).group(1))
        rows.append((T, math.log(T), obs, mod))
rows = sorted(set(rows)); J = 1.000100147643
print(' T          logT    obs      deficit  deficit*logT  obs-model  (obs-model)*log^2T')
for T, LT, o, m in rows:
    print(f'{T:10.1f} {LT:7.3f} {o:.6f} {J-o:.5f}  {(J-o)*LT:.4f}       {o-m:+.5f}   {(o-m)*LT*LT:+.3f}')
L = np.array([r[1] for r in rows]); r = np.array([r[2] - r[3] for r in rows]); dfc = J - np.array([x[2] for x in rows])
for k in (6, 4):
    (b, c), *_ = np.linalg.lstsq(np.vstack([1/L[-k:]**2, np.ones(k)]).T, r[-k:], rcond=None)
    print(f'fit (obs-model) = b/log^2T + c over last {k} points: b={b:+.3f}, c={c:+.5f}')
    (a, b, c), *_ = np.linalg.lstsq(np.vstack([1/L[-k:], 1/L[-k:]**2, np.ones(k)]).T, dfc[-k:], rcond=None)
    print(f'fit deficit = a/logT + b/log^2T + c over last {k} points: a={a:+.3f}, b={b:+.3f}, c={c:+.5f}  (predicted a=(1+log 2pi)(2f(0)-J)={(1+math.log(2*math.pi))*(2*0.7892-J):.3f})')
