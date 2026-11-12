# t3_flip.py (explore39): term-wise Beurling-Deny minorant Q_flip = Q - F (anti-Markov parts flipped), per parity.
# lambda_0(Q) >= lambda_0(Q_flip) and Q_flip is Markov, so Barta evaluates it; where lambda_0(Q_flip) < 0 no Markov comparison of the
# explicit-formula terms certifies the window.
import sys, numpy as np, math, json, time
sys.path.insert(0, '.scratch/explore39'); from xform import Form, barta_lp
kind = sys.argv[1] if len(sys.argv) > 1 else 'zeta'
rho = 1.324717957244746
t0 = time.time()
ev = [1.8, 1.9, 2.0, 2.1, 2.2, 2.3, 2.5, 3.0, 5.0, 10.0, 20.0, 47.0]
od = [2.02, 2.05, 2.1, 2.2, 2.5, 3.0, 4.0, 5.0, 10.0, 20.0, 47.0]
if kind == 'dh':
    ev = [2.5, 3.0, 3.2, 3.5, 4.0, 5.0, 10.0, 20.0, 31.0]; od = [2.02, 2.05, 2.1, 2.2, 2.5, 3.0, 5.0, 10.0, 20.0, 31.0]
for sec, cs in (('even', ev), ('odd', od)):
    for c in cs:
        F = Form(c, kind); n = 20 if c < 10 else 28
        Qm, p = F.Q(n, sec); Fm, info = F.flip(n, sec)
        e = np.linalg.eigvalsh(Qm); ef, Vf = np.linalg.eigh(Qm - Fm)
        out = {'kind': kind, 'sec': sec, 'c': c, 'n': n, 'lam0_Q': float('%.6g' % e[0]), 'lam0_Qflip': float('%.6g' % ef[0]),
               'normF': float('%.4g' % np.linalg.norm(Fm, 2))}
        out.update({k: float('%.5g' % v) for k, v in info.items()})
        if ef[0] > -0.05 and c < 4:
            b, v, xs = barta_lp(F, n, sec, flip=True, v0=Vf[:, 0]); out['barta_Qflip'] = (float('%.6g' % b) if b is not None else None)
        out['t'] = round(time.time() - t0, 1)
        print(json.dumps(out), flush=True)
