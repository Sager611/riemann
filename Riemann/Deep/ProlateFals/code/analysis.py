# analysis of the even negative spectrum, lambda = sqrt 2 (Dirac 2D: t_n = 2 E_n)
import math, sys
import mpmath as mp
mp.mp.dps = 30
D = sys.argv[1] if len(sys.argv) > 1 else '.'
def load_eig(fn):
    out = []
    for l in open(fn):
        if l.startswith('#'): continue
        p = l.split(); out.append((int(p[0]), float(p[1])))
    return out
ev = load_eig(f'{D}/e0_1000.txt')
lam = math.sqrt(2); kap = 2*math.pi*lam**2
def smooth(t): return t/(2*math.pi)*(math.log(t/(2*math.pi)) - 1)
print('n0 check: first even index j =', ev[0][0], ' (n = j, so n0 = 0)')
# 1. location test: theta(t_n)/pi - (n - 3/8)
print('\n# B26 Thm B: e_n = vartheta(t_n)/pi - (n - 3/8)')
rows = []
for j, E in ev:
    t = 2*E
    th = float(mp.siegeltheta(t))
    rows.append((j, t, th/math.pi - (j - 0.375)))
for j, t, e in rows:
    if j in (1,2,3,5,10,20,50,100,200,500,1000,1500) or j == rows[-1][0]:
        print(f'n={j:5d} t_n={t:14.8f} e_n={e: .6e}  e_n*t^(2/3)={e*t**(2/3): .4f}  e_n*t/log t={e*t/math.log(t): .4f}')
# local exponent fits of |e_n| ~ C t^-p on windows
print('\n# local exponent p in |e_n| ~ C t^-p (least squares on log-log, windows of n)')
def fit(xs, ys):
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    sxx = sum((x-mx)**2 for x in xs); sxy = sum((x-mx)*(y-my) for x, y in zip(xs, ys))
    b = sxy/sxx; a = my - b*mx
    res = [y-(a+b*x) for x, y in zip(xs, ys)]
    se = math.sqrt(sum(r*r for r in res)/(n-2)/sxx) if n > 2 else float('nan')
    return a, b, se
for lo, hi in ((10,50),(50,200),(200,600),(600,1516)):
    xs = [math.log(t) for j,t,e in rows if lo <= j <= hi]
    ys = [math.log(abs(e)) for j,t,e in rows if lo <= j <= hi]
    a, b, se = fit(xs, ys)
    print(f'n in [{lo},{hi}]: p = {-b:.4f} +- {se:.4f}')
# 2. counting law: R_n = (n - 1/2) - smooth(t_n); fit R = c + b log t
print('\n# counting law residual R_n = N(t_n^-)+1/2 - (t/2pi)(log(t/2pi)-1) = n - 1/2 - smooth(t_n)')
R = [(t, j - 0.5 - smooth(t)) for j, t, e in rows]
for lo, hi in ((1,1516),(1,100),(100,1516),(500,1516),(1000,1516)):
    xs = [math.log(t) for (t, r), (j, _, _) in zip(R, rows) if lo <= j <= hi]
    ys = [r for (t, r), (j, _, _) in zip(R, rows) if lo <= j <= hi]
    a, b, se = fit(xs, ys)
    print(f'n in [{lo},{hi}] (t in [{2*ev[lo-1][1]:.1f},{2*ev[hi-1][1]:.1f}]): b = {b: .5f} +- {se:.5f}  (CM22: -0.159155, B26: 0), c = {a: .4f}')
for j in (1,10,100,1000,1516):
    t, r = R[j-1]; print(f'  n={j}: t={t:.3f} R={r: .6f}  CM22 shape -(log t)/2pi = {-math.log(t)/(2*math.pi): .4f}')
# 3. zeta comparison, first 100
print('\n# first 100: t_n vs zeta zeros gamma_n')
zs = [float(mp.zetazero(n).imag) for n in range(1, 101)]
for n in (1,2,3,4,5,10,20,50,100):
    print(f'n={n:3d} t_n={rows[n-1][1]:10.4f} gamma_n={zs[n-1]:10.4f} diff={rows[n-1][1]-zs[n-1]: .4f}')
md = max(abs(rows[n-1][1]-zs[n-1]) for n in range(1,101))
print('max |t_n - gamma_n| n<=100:', round(md, 4))
open(f'{D}/tn_even_sqrt2.txt','w').write(''.join(f'{j} {t:.12f} {e:.6e}\n' for j,t,e in rows))
