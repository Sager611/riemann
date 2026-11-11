# nbrate: pair term D(N) of the BCF ramp, truncated to the first M zero pairs.
# D(N) = sum_{rho1,rho2} K(rho1,rho2) N^{rho1+rho2-1} / (zeta'(rho1) zeta'(rho2)),
# K(r1,r2) = (1/2pi i) int_{(c)} zeta(s) zeta(1-s) ds / ((r1-s)^2 (r2-1+s)^2 s(1-s)),  0 < c < 1/2 (here c = 1/4).
import mpmath as mp, numpy as np, math, sys, time
t0 = time.time()
mp.mp.dps = 15
M = int(sys.argv[1]) if len(sys.argv) > 1 else 6
c = 0.25; T = 260.0; h = 0.01
ts = np.arange(0.0, T + h / 2, h)
z1 = np.array([complex(mp.zeta(mp.mpc(c, t))) for t in ts])        # zeta(c+it), t>=0
z2 = np.array([complex(mp.zeta(mp.mpc(1 - c, t))) for t in ts])    # zeta(1-c+it)
print(f'grid done {time.time()-t0:.1f}s', flush=True)
# full grid t in [-T,T]: zeta(c+it) for t<0 is conj(zeta(c-it)); zeta(1-s) = zeta(1-c-it) = conj(zeta(1-c+it))
tt = np.concatenate([-ts[:0:-1], ts])
Zs = np.concatenate([np.conj(z1[:0:-1]), z1])                # zeta(s), s = c+it
Z1s = np.concatenate([z2[:0:-1], np.conj(z2)])               # zeta(1-s) = zeta(1-c-it): t<0 -> zeta(1-c+i|t|)
s = c + 1j * tt
base = Zs * Z1s / (s * (1 - s))
wts = np.ones_like(tt) * h; wts[0] = wts[-1] = h / 2     # trapezoid (analytic, fast decay)
def K(r1, r2):
    return complex(np.sum(base / ((r1 - s) ** 2 * (r2 - 1 + s) ** 2) * wts) / (2 * math.pi))
zs = [complex(mp.zetazero(n)) for n in range(1, M + 1)]
zp = [complex(mp.zeta(mp.mpc(r.real, r.imag), derivative=1)) for r in zs]
rho = zs + [z.conjugate() for z in zs]; dz = zp + [z.conjugate() for z in zp]
Kd = {}
for i, r1 in enumerate(rho):
    for j, r2 in enumerate(rho):
        Kd[(i, j)] = K(r1, r2)
print(f'K done {time.time()-t0:.1f}s', flush=True)
# diagonal (r2 = conj r1 = 1 - r1): contributes to the constant D0
D0 = sum(Kd[(i, (i + M) % (2 * M))] / (dz[i] * dz[(i + M) % (2 * M)]) for i in range(2 * M))
print('D0 partial sums (first m zeros):')
for m in range(1, M + 1):
    idx = list(range(m)) + list(range(M, M + m))
    s0 = sum(Kd[(i, (i + M) % (2 * M))] / (dz[i] * dz[(i + M) % (2 * M)]) for i in idx)
    print(f'  m={m}: gamma_m={zs[m-1].imag:.4f}  D0_m = {s0.real:+.6f} {s0.imag:+.1e}i   term_m = {(2*(Kd[(m-1,m-1+M)]/(dz[m-1]*dz[m-1+M]))).real:+.3e}')
def Dtr(N):
    L = math.log(N); tot = 0
    for i in range(2 * M):
        for j in range(2 * M):
            ph = np.exp(L * (rho[i] + rho[j] - 1))
            tot += Kd[(i, j)] * ph / (dz[i] * dz[j])
    return tot
np.save('/Users/adriansager/playground/riemann/.scratch/nbrate/Ktab.npy', np.array([[Kd[(i, j)] for j in range(2 * M)] for i in range(2 * M)]))
np.save('/Users/adriansager/playground/riemann/.scratch/nbrate/rho_dz.npy', np.array([rho, dz]))
# compare with Gram numerics R(N) = L^2 Q_BCF - C L - J
O = '/Users/adriansager/playground/riemann/.scratch/nbrate/'
bcf = np.load(O + 'bcf_rig_mid.npy')
EG = 0.57721566490153286061; C = 2 + EG - math.log(4 * math.pi); J = 3.562663456536367655
print('N, R(N)=L^2Q-CL-J (Gram, rigorous), D_trunc(N) (first M pairs), difference')
for N in [50, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100, 1200, 1280]:
    L = math.log(N); R = L * L * bcf[N] - C * L - J; Dt = Dtr(N)
    print(f'{N:5d} {R:+.5f} {Dt.real:+.5f} ({Dt.imag:+.0e}i) diff {R - Dt.real:+.5f}', flush=True)
# largest off-diagonal amplitudes
amp = []
for i in range(M):
    for j in range(2 * M):
        if j == (i + M) % (2 * M): continue
        a = Kd[(i, j)] / (dz[i] * dz[j]); fr = (rho[i] + rho[j] - 1).imag
        amp.append((abs(a), fr, i, j))
amp.sort(reverse=True)
print('largest off-diagonal |K/(zeta zeta)| (each also has its conjugate partner):', [(round(a, 6), round(f, 3)) for a, f, i, j in amp[:8]])
print(f'done {time.time()-t0:.1f}s')
