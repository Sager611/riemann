# explore62: Riemann-Siegel (zdil.zrs) vs Euler-Maclaurin (oscz_copy) vs mpmath siegelz at a few heights
import numpy as np, mpmath as mp, zdil, oscz_copy as O
for t0 in (1e4+0.37, 1e5+0.61, 1e6+0.123):
    t = np.array([t0]); a = zdil.zrs(t)[0]; b = float(mp.siegelz(t0)); c = O.Z('zeta', t)[0][0] if t0 < 2e5 else float('nan')
    print(f"t={t0}: RS={a:.10f} mp={b:.10f} EM={c:.10f} |RS-mp|={abs(a-b):.1e}")
