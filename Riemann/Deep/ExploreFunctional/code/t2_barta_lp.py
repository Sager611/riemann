# t2_barta_lp.py (explore39): best Barta bound in the Legendre span (LP), Markov ranges, and the flip minorant.
import sys, numpy as np, math, json, time
sys.path.insert(0, '.scratch/explore39'); from xform import Form, barta_lp
rho = 1.324717957244746
t0 = time.time()
cases = [(1.2, 'even', False), (1.2, 'odd', False), (rho, 'even', False), (1.5, 'even', False), (rho**2, 'even', False),
         (1.5, 'odd', False), (2.0, 'odd', False)]
for c, sec, fl in cases:
    F = Form(c)
    for n in [12, 20]:
        Qm, p = F.Q(n, sec); e, V = np.linalg.eigh(Qm)
        b, v, xs = barta_lp(F, n, sec, flip=fl, v0=V[:, 0])
        print(json.dumps({'c': round(c, 6), 'sec': sec, 'n': n, 'galerkin': round(float(e[0]), 7), 'barta_lp': (round(b, 7) if b is not None else None), 't': round(time.time() - t0, 1)}), flush=True)
