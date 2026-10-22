import hbmodel as H, numpy as np
cases = [('HB X, b=1/25, dt=1/400', 60000083951.5, 1/25, 1/400),
         ('X=6e12, b=1/25', 6e12, 1/25, 1/200),
         ('PM Table1 X=5e12+194858, b=0.16733^2', 5e12+194858, 0.16733**2, 1/200),
         ('Raymond line wall X0=6000000185827, b=1e-8', 6000000185827, 1e-8, 1/200),
         ('X=6e10 with line wall (hypothetical), b=1e-8', 60000083951.5, 1e-8, 1/200)]
for name, X, b, dt in cases:
    val, rows = H.chain(X, b, dt=dt)
    r0 = rows[0]; rl = rows[-1]
    print('%-50s Lambda<=%.5f  t_end=%.5f  k0=%.3f kend=%.3f C0=%.4f G0=%.4f steps=%d' % (name, val, val-b/2, r0[7], rl[7], r0[3], r0[4], len(rows)))
