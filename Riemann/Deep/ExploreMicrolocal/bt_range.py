# explore47: optimized triangle (Bérard/CCM) bound B_tri(t) = min_Δ[log(t/2π)/(2Δ) + TB(Δ)] vs unconditional explicit S-bounds
import sys, os, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from btlib import prime_powers, TB
n, lam, _ = prime_powers(int(math.exp(17)) + 1); logn = np.log(n.astype(float))
Df = np.arange(1.0, 17.001, 0.05); TBf = np.array([TB(d, logn, lam, n) for d in Df])
print('Δ  TB(Δ):', ' '.join(f'{d:.0f}:{TBf[np.argmin(abs(Df-d))]:.3f}' for d in range(2, 18)))
for e in [3, 4, 6, 8, 10, 12, 20, 21, 30, 50, 100]:
    t = 10.0 ** e; L = math.log(t / (2 * math.pi)); tot = L / (2 * Df) + TBf; j = int(np.argmin(tot))
    lt = math.log(t); llt = math.log(lt)
    tr = 0.112 * lt + 0.278 * llt + 2.510; hsw = 0.1038 * lt + 0.2573 * llt + 9.3675
    print(f't=1e{e}: Δ_tri={Df[j]:.2f} (2loglog t={2*llt:.2f}) B_tri={tot[j]:.3f} CCM¼={0.25*lt/llt:.3f} B_tri/(log t/loglog t)={tot[j]*llt/lt:.3f} | Trudgian14={tr:.2f} HSW22={hsw:.2f}')
