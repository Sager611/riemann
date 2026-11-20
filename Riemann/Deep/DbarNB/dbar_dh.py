# explore59 dbar_dh.py : Davenport-Heilbronn zeros (on-line by sign changes of Lambda, off-line refined) and U_DH(L)
import numpy as np, mpmath as mp, time, sys
from dbar_U import U
mp.mp.dps = 20
t0 = time.time()
kap = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
a = [1, kap, -kap, -1]
def f(s):
    return mp.power(5, -s) * sum(a[r - 1] * mp.zeta(s, mp.mpf(r) / 5) for r in range(1, 5))
def Lam(t):
    s = mp.mpf('0.5') + 1j * t
    return mp.exp((s + 1) / 2 * mp.log(5 / mp.pi) + mp.loggamma((s + 1) / 2)) * f(s)
T = float(sys.argv[1]) if len(sys.argv) > 1 else 180.0
ts = np.arange(0.05, T, 0.04)
vals = []
for t in ts:
    v = Lam(t); vals.append(complex(v))
vals = np.array(vals)
print('max |Im Lambda|/|Lambda| =', float(np.max(np.abs(vals.imag) / (np.abs(vals) + 1e-300))))
re = vals.real
on = []
for i in range(len(ts) - 1):
    if re[i] == 0 or re[i] * re[i + 1] < 0:
        r = mp.findroot(lambda x: mp.re(Lam(x)), (ts[i], ts[i + 1]), solver='anderson')
        on.append(float(r))
on = np.array(on)
guess = [0.808517 + 85.699348j, 0.650830 + 114.163343j, 0.574356 + 166.479306j, 0.724258 + 176.702461j]
off = []
for g0 in guess:
    if g0.imag < T:
        r = mp.findroot(f, mp.mpc(g0)); off.append(complex(r))
print('on-line zeros found:', len(on), ' first:', on[:4], ' off-line (refined):', [f'{z.real:.6f}+{z.imag:.6f}i' for z in off])
NT = lambda x: x / (2 * np.pi) * np.log(5 * x / (2 * np.pi * np.e)) + 7.0 / 8   # smooth zero count, conductor 5 (approx)
print(f'smooth N({T}) ~ {NT(T):.1f}; on-line {len(on)} + 2*off-line {2*len(off)} = {len(on)+2*len(off)}')
np.save('logs/dh_zeros_on.npy', on); np.save('logs/dh_zeros_off.npy', np.array(off))
Zon = np.concatenate([on, -on])
Zoff = []
for z in off:
    b = z.real - 0.5; g = z.imag
    Zoff += [g - 1j * b, -g - 1j * b, g + 1j * b, -g + 1j * b]
Zoff = np.array(Zoff)
Con = float(np.sum(1 / (0.25 + Zon ** 2)))
print(f'DH on-line sum 1/|rho|^2 (|gamma|<{T}) = {Con:.6f}')
for L in (7.155, 10.0, 20.0, 40.0, 80.0, 160.0):
    u_on = U(Zon, L); u_all = U(np.concatenate([Zon, Zoff]), L)
    print(f'L={L:7.3f}  L*U_on={L*u_on:.6f}  L*U_all={L*u_all:.6f}  U_all={u_all:.3e}')
fl = sum(1 / (((0.25 + (z.imag ** 2 - (z.real - .5) ** 2)) / (2 * (z.real - .5)) + 0.5 + (z.real - .5))) for z in off) * 2
print(f'floor (sum over right zeros of 1/K_inf) ~ {fl:.3e}')
print('time', round(time.time() - t0, 1), 's')
