"""explore53: scan h_J (J-th order local criterion, joint with the amplifier) over kappa for several J at once.
usage: scanJ.py r theta d Jmax kappa1 kappa2 ...   -> logs/scanJ_r{r}_d{d}.json (appends)"""
import sys, os, json
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import joint as Jt
from sdpcut import sdp_cut
from scipy.optimize import minimize_scalar
NIT = int(os.environ.get('NIT', '2'))
ROUNDS = int(os.environ.get('ROUNDS', '6'))
MARGIN = float(os.environ.get('MARGIN', '0'))
JMIN = int(os.environ.get('JMIN', '3'))

r = int(sys.argv[1]); th = float(sys.argv[2]); d = int(sys.argv[3]); Jmax = int(sys.argv[4])
kaps = [float(x) for x in sys.argv[5:]]
v0 = 1.3617
fn = os.path.join(HERE, 'logs', 'scanJ_r%d_d%d.json' % (r, d))
res = json.load(open(fn)) if os.path.exists(fn) else []
for kap in kaps:
    C0 = Jt.central_forms(r, th, d, kap, 2 * Jmax, v0)
    fw = lambda dv: -Jt.lmax(*[Jt.shift(C0[:3], dv, kap)[i] for i in (0, 2)])[0]
    rw = minimize_scalar(fw, bounds=(-0.01, 0.01), method='bounded', options={'xatol': 1e-9})
    h1 = -rw.fun
    line = {'kappa': kap, 'h1': h1}
    for J in range(JMIN, Jmax + 1):
        CJ = C0[:2 * J + 1]
        dv = rw.x
        _, y = Jt.lmax(*[Jt.shift(CJ, dv, kap)[i] for i in (0, 2)])
        for it in range(NIT):
            C = Jt.shift(CJ, dv, kap)
            n0 = y @ C[0] @ y
            m = [y @ C[2 * j] @ y / n0 for j in range(J + 1)]
            val, a, ph, mnv, rds = sdp_cut(m, J, rounds=ROUNDS, margin=MARGIN)
            line['minvalid%d' % J] = mnv
            def fa(dvv):
                Cs = Jt.shift(CJ, dvv, kap)
                return a[0] * Jt.lmax(Cs[0], sum(a[j] * Cs[2 * j] for j in range(1, J + 1)))[0]
            rr = minimize_scalar(fa, bounds=(dv - 0.004, dv + 0.004), method='bounded', options={'xatol': 1e-9})
            dv = rr.x; hJ = -rr.fun
            Cs = Jt.shift(CJ, dv, kap)
            _, y = Jt.lmax(Cs[0], sum(a[j] * Cs[2 * j] for j in range(1, J + 1)))
        line['h%d' % J] = hJ; line['a%d' % J] = list(map(float, a)); line['phi%d' % J] = {str(k): v for k, v in ph.items()}
        line['v%d' % J] = v0 + dv
        print('kappa=%.4f J=%d: h1=%.6f hJ=%.6f (SDP %.3g)' % (kap, J, h1, hJ, val), flush=True)
    res.append(line)
    json.dump(res, open(fn, 'w'), indent=1)
