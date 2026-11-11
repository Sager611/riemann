# analysis 2: large-t windows, split of Theta-Phi, Heywood-type phase integral
import math, sys, glob
import mpmath as mp
mp.mp.dps = 30
D = sys.argv[1] if len(sys.argv) > 1 else '.'
lam = math.sqrt(2); kap = 2*math.pi*lam**2
def smooth(t): return t/(2*math.pi)*(math.log(t/(2*math.pi)) - 1)
def load(fn):
    return [(int(l.split()[0]), float(l.split()[1])) for l in open(fn) if not l.startswith('#')]
ev = load(f'{D}/e0_1000.txt')
win = {}
for tag in ('1e4', '5e4', '1e5'):
    win[tag] = load(f'{D}/ew_{tag}.txt')
print('# windows: location error e_n = vartheta(t_n)/pi - (n - 3/8) and count residual R_n = n - 1/2 - smooth(t_n)')
allr = []
for j, E in ev + sum(win.values(), []):
    t = 2*E
    e = float(mp.siegeltheta(t)/mp.pi) - (j - 0.375)
    allr.append((j, t, e, j - 0.5 - smooth(t)))
for tag, w in win.items():
    for j, E in w:
        t = 2*E; r = [x for x in allr if x[0] == j and abs(x[1]-t) < 1e-9][0]
        print(f'  n={j:7d} t_n={t:.9f} e_n={r[2]: .6e} e_n*t/log t={r[2]*t/math.log(t): .3f} R_n={r[3]: .6f}')
# exponent between windows using per-window mean e
def mean_e(lst): return sum(x[2] for x in lst)/len(lst), sum(x[1] for x in lst)/len(lst)
groups = [[x for x in allr if 1400 <= x[1] <= 2000]] + [[x for x in allr if abs(x[1]/2 - c) < 3] for c in (10000, 50000, 100000)]
ms = [mean_e(g) for g in groups]
print('\n# local decay exponent p of |e| between scales (t ~ 1.7e3, 2e4, 1e5, 2e5)')
for (e1, t1), (e2, t2) in zip(ms, ms[1:]):
    print(f'  t {t1:9.0f} -> {t2:9.0f}: p = {math.log(abs(e1)/abs(e2))/math.log(t2/t1):.4f}  (B26 bound needs p >= 2/3; log t/t gives p -> 1)')
# count law: local log-slope b between scales, CM22 predicts -1/(2 pi) = -0.159155
print('\n# count residual, mean per scale, and local slope b = dR/dlog t (CM22: -0.159155; B26: -> 0)')
rs = [(sum(x[3] for x in g)/len(g), sum(x[1] for x in g)/len(g)) for g in groups]
for (r1, t1), (r2, t2) in zip(rs, rs[1:]):
    print(f'  t {t1:9.0f} -> {t2:9.0f}: R {r1: .5f} -> {r2: .5f}, b = {(r2-r1)/math.log(t2/t1): .5f}')
print(f'  CM22 would need R(2e5) - R(1.7e3) = {-math.log(2e5/1.7e3)/(2*math.pi):.4f}; observed {rs[-1][0]-rs[0][0]:.4f}')
# split Theta - Phi = [I(E) - (E log(4E/kap) - E)] + [eta(inf) + pi/4]
print('\n# split of Theta - Phi (lambda = sqrt2) on the log grid')
def I_exact(E):
    E = mp.mpf(E); k = mp.mpf(kap)
    f = lambda t: E**2/(mp.sqrt(E**2 + k**2*mp.cosh(t)**2) + k*mp.cosh(t))
    tc = mp.acosh(max(E/k, 1))
    return mp.quad(f, [0, tc/2, tc, tc+2, tc+6, tc+20, mp.inf])
def J_heywood(E):
    # (1/pi) J = Heywood-type count: int_0^inf sqrt(E^2+V) - sqrt(V), V = kap^2 cosh^2 + (coth^2 - 2)/4 = -q (CM22 Eq. 34)
    E = mp.mpf(E); k = mp.mpf(kap)
    V = lambda y: k**2*mp.cosh(y)**2 + (mp.coth(y)**2 - 2)/4
    f = lambda y: E**2/(mp.sqrt(E**2 + V(y)) + mp.sqrt(V(y)))
    tc = mp.acosh(max(E/k, 1))
    return mp.quad(f, [0, 1/E, 10/E, 0.5, tc/2, tc, tc+2, tc+6, tc+20, mp.inf])
rows = [tuple(map(float, l.split()[:2])) for l in open(f'{D}/lg1.txt')]
print('  E           Theta-Phi     I-Iasym       eta+pi/4      (J-I)+0.5logE   Heywood count - true (phase units)')
for E, T in rows:
    if E < 3: continue
    Iex = float(I_exact(E)); Ias = E*math.log(4*E/kap) - E
    Phi = Ias + math.pi/4
    eta = T - math.pi/2 - Iex          # eta(inf)
    J = float(J_heywood(E))
    print(f'  {E:10.2f} {T-Phi: .6e} {Iex-Ias: .6e} {eta+math.pi/4: .6e} {J-Iex+0.5*math.log(E): .6f}  {J - (T - math.pi/4): .5f}')
