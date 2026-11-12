# explore38 s6: asymptotic second-order coefficient of the Christoffel bound: c2B(inf) = -(S1 + P0),
# S1 = sum_rho |c_rho|^2 (Phi'(gamma) - L) = sum (2 theta'(gamma) + 3/|rho|^2)/|rho|^2,
# P0 = Re sum_{a != b} conj(c_a) c_b i/(gamma_b - gamma_a)  (non-oscillating part of the off-diagonal kernel), c = (rho-1)/rho^2.
import numpy as np, math
OUT = '/Users/adriansager/playground/riemann/.scratch/explore38/'
zc = np.load(OUT + 'zeros_cache.npz'); g = zc['g']; th1 = zc['th1']
for K in [200, 500, 1000, 2000, 3000]:
    gp = g[:K]; gam = np.concatenate([gp, -gp]); rho = 0.5 + 1j * gam
    c = (rho - 1) / rho ** 2
    A = np.concatenate([2 * th1[:K], 2 * th1[:K]]) + 3 / np.abs(rho) ** 2
    S1 = float(np.sum(np.abs(c) ** 2 * A))
    # tail of S1: zeros K..22491 (file) + integral
    gt = g[K:]; S1t = 2 * float(np.sum((np.log(gt / (2 * np.pi)) + 3 / (0.25 + gt ** 2)) / (0.25 + gt ** 2)))
    T = g[-1]; u = np.linspace(0, 60, 60001); t = T * np.exp(u); lg = np.log(t / (2 * np.pi))
    S1t += float(np.trapezoid(2 * (lg / (2 * np.pi)) * lg / t, u))
    P0 = 0.0
    for r0 in range(0, len(gam), 100):   # block 100 (block 500 peaked at 255 MB)
        dg = gam[None, :] - gam[r0:r0 + 100, None]
        with np.errstate(divide='ignore', invalid='ignore'):
            Ksub = 1j / dg
        Ksub[~np.isfinite(Ksub)] = 0
        P0 += float(np.real(np.conj(c[r0:r0 + 100]) @ (Ksub @ c)))
    print(f'K={K:5d} gamma_K={gp[-1]:9.2f}  S1_K={S1:.6f} S1_tail={S1t:.6f} S1={S1+S1t:.6f}  P0_K={P0:+.6f}  c2B(inf) ~ {-(S1+S1t+P0):+.6f}', flush=True)
