# explore59 dbar_pair.py : (a) close pair -> no 1/|zeta'| blow-up; (b) trivial-zero penalty nu_k; (c) off-line moved pair floor
import numpy as np, math, time
from dbar_U import Kblock, U, zz
t0 = time.time()
K = 400
g = zz[:K].copy()
Cfull = 2 + 0.5772156649015329 - math.log(4 * math.pi)
r1 = 0.25 + g[0] ** 2
print('== (a) zero gamma_1 split into gamma_1 +- eta/2 (both signs of t); predicted limits: m=1+1 ->', f'{np.sum(2/(0.25+g**2)) + 2/r1:.6f}', ' m=2 (double) ->', f'{np.sum(2/(0.25+g**2)) + 6/r1:.6f}')
for L in (7.155, 10.0, 20.0, 40.0, 80.0):
    row = []
    for eta in (2.0, 1.0, 0.3, 0.1, 0.03, 0.01, 1e-3, 1e-4):
        gg = np.concatenate([[g[0] - eta / 2, g[0] + eta / 2], g[1:]])
        Zs = np.concatenate([gg, -gg])
        u, c, cond = U(Zs, L, True)
        row.append(f'{eta:g}:{L*u:.5f}')
    print(f'L={L:6.2f} L*U by eta  ' + '  '.join(row) + f'   (cond at 1e-4 {cond:.1e})')
print('   J_-1 proxy for the split model: |F\'(rho)|^-2 ~ eta^-2 -> 1e8 at eta=1e-4 (BCF needs sum |zeta\'(rho)|^-2 << T^{3/2-d})')
print('== (b) trivial-zero penalty nu_k = |1-E_low(-2k)| e^{-L(2k+1/2)} / |zeta\'(-2k)| for E_low of band L\'=(1-eps)L')
zp = lambda k: math.factorial(2 * k) * (1 + 2 ** (-2 * k - 1) + 3 ** (-2 * k - 1)) / (2 * (2 * math.pi) ** (2 * k))
Zs = np.concatenate([g, -g])
for L in (10.0, 20.0, 40.0):
    for eps in (0.0, 0.1, 0.25):
        Lp = (1 - eps) * L
        u, c, cond = U(Zs, Lp, True)
        out = []
        for k in range(1, 7):
            zk = 1j * (2 * k + 0.5)
            Ek = (Kblock([zk], Zs, Lp) @ c)[0]
            # log of nu_k to avoid overflow
            lnu = math.log(abs(1 - Ek)) - L * (2 * k + 0.5) - math.log(zp(k))
            out.append(f'k{k}:{lnu/math.log(10):+.1f}')
        print(f'L={L:5.1f} eps={eps:4.2f} L\'*U={Lp*u:.5f}  log10 nu_k: ' + ' '.join(out))
print('== (c) off-line quadruple: gamma_1 moved to beta=1/2+b (z = +-g1 +- i b); L*U should grow ~L (floor), DH-like b=0.3085 at g=85.699')
for b, gm in ((0.3085, 85.699348), (0.1, 85.699348), (0.01, 85.699348), (0.3085, 14.134725)):
    zs_off = np.array([gm - 1j * b, gm + 1j * b, -gm - 1j * b, -gm + 1j * b])
    gg = g[np.abs(g - gm) > 0.5]
    Zs = np.concatenate([gg, -gg, zs_off])
    row = []
    for L in (7.155, 10.0, 20.0, 40.0, 80.0):
        row.append(f'L={L:g}:{L*U(Zs, L):.5f}')
    # evaluator floor for one right zero: 1/K(z,z) as L->inf
    zr = gm - 1j * b
    print(f'b={b} gamma={gm}  L*U: ' + '  '.join(row) + f'   floor 1/K_inf(z,z)~{1/ (( (0.25 + (gm**2 - b**2))/(2*b) + (0.5 + b)) ):.2e} per right zero')
print('time', round(time.time() - t0, 1), 's')
