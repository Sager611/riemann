# t2_count.py: on-line zeros of DH below the off-line pairs (sign changes of the real function Lambda(1/2+it)), vs N(T) formula.
import sys, json
sys.path.insert(0, '.')
import mpmath as mp
import ecommon as E
mp.mp.dps = 15
T1 = float(sys.argv[1]) if len(sys.argv) > 1 else 120.0
dt = 0.02
prev = None; zs = []
n = int(T1 / dt)
for k in range(1, n + 1):
    t = k * dt
    v = mp.re(E.dh(mp.mpc(0.5, t)))
    if prev is not None and v * prev < 0: zs.append(round(t - dt / 2, 3))
    prev = v
def NT(T):
    return T / (2 * mp.pi) * mp.log(5 * T / (2 * mp.pi * mp.e)) + mp.mpf(1) / 8
out = {'T': T1, 'n_signchanges': len(zs), 'first5': zs[:5]}
for g in (85.699348, 114.163343):
    out['below_%g' % g] = sum(1 for z in zs if z < g)
    out['NT_%g' % g] = float(NT(g))
print(json.dumps(out))
