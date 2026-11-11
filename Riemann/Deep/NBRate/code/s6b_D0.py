# nbrate: diagonal pair constant D0 = sum_rho K(rho, 1-rho)/|zeta'(rho)|^2 for the first M zeros (BCF ramp second-order term)
import mpmath as mp, numpy as np, math, sys, time
t0 = time.time(); mp.mp.dps = 15
M = int(sys.argv[1]) if len(sys.argv) > 1 else 40
O = '/Users/adriansager/playground/riemann/.scratch/nbrate/'
c = 0.25; T = 300.0; h = 0.01
ts = np.arange(0.0, T + h / 2, h)
z1 = np.array([complex(mp.zeta(mp.mpc(c, t))) for t in ts]); z2 = np.array([complex(mp.zeta(mp.mpc(1 - c, t))) for t in ts])
np.savez(O + 'zgrid.npz', ts=ts, z1=z1, z2=z2)
print(f'grid {time.time()-t0:.0f}s', flush=True)
tt = np.concatenate([-ts[:0:-1], ts]); Zs = np.concatenate([np.conj(z1[:0:-1]), z1]); Z1s = np.concatenate([z2[:0:-1], np.conj(z2)])
s = c + 1j * tt; base = Zs * Z1s / (s * (1 - s)); wts = np.ones_like(tt) * h; wts[0] = wts[-1] = h / 2
tot = 0.0; S = 0.0
print(' m  gamma  term=2Re K/|z\'|^2  term*g^2/log(g/2pi)  partial D0   partial sum log(g/2pi)/g^2')
for m in range(1, M + 1):
    r = complex(mp.zetazero(m)); dz = complex(mp.zeta(mp.mpc(r.real, r.imag), derivative=1))
    K = complex(np.sum(base / ((r - s) ** 4) * wts) / (2 * math.pi))   # K(rho, 1-rho): (rho-s)^2 (1-rho-1+s)^2 = (s-rho)^4
    term = 2 * (K / (dz * dz.conjugate())).real
    g = r.imag; tot += term; S += math.log(g / (2 * math.pi)) / g ** 2
    if m <= 10 or m % 5 == 0:
        print(f'{m:3d} {g:8.3f} {term:+.4e} {term*g*g/math.log(g/(2*math.pi)):+.4f} {tot:+.6f} {S:.6f}', flush=True)
gM = r.imag
tail = (1 / (2 * math.pi)) * (math.log(gM / (2 * math.pi)) ** 2 + 2 * math.log(gM / (2 * math.pi)) + 2) / gM
print(f'tail of sum log(g/2pi)/g^2 beyond gamma_M (density log(t/2pi)/2pi): {tail:.6f}; full S ~ {S+tail:.6f}')
print(f'done {time.time()-t0:.0f}s')
