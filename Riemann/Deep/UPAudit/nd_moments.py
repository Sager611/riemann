# N-D: the M3 threshold q < kappa/theta in the toy model. One planted pair at depth a/L over the rescaled zeta background;
# J_q(a) := int_{|t-gamma_p| <= 6/L} (|M_toy(1/2+it)|^q - |M_bg(1/2+it)|^q) dt  (excess q-moment caused by the zero).
# Bin a of the density profile carries ~ c0 N e^{-kappa a} zeros, so the q-moment ~ N sum_a e^{-kappa a} J_q(a):
# finite uniformly in L iff the growth rate g_q := lim d log J_q / da satisfies g_q < kappa. Lemma Z predicts g_q <= q*theta.
import numpy as np, sys, json, math
sys.path.insert(0, '/Users/adriansager/playground/riemann/Riemann/Deep/UPAudit')
from nc_toy import Toy, background, planted
qs = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0]
res = {}
for L in ([20.0] if len(sys.argv) > 1 else [44.0]):  # toy is scale invariant in units of 1/L; L only bounds a <= L/2
    m = Toy(L, du=0.02); t0 = 1000.0; bg = background(L, t0, 12.0)
    ts = t0 + np.arange(-60, 61)/(10*L); dt = 1/(10*L)
    Mbg = np.array([m.M(bg, 0.5+1j*t) for t in ts])
    rows = []
    for a in range(1, int(L/2)):
        zs = np.concatenate([bg, planted('single', L, t0, a)])
        Mz = np.array([m.M(zs, 0.5+1j*t) for t in ts])
        J = {q: float(np.sum(Mz**q - Mbg**q)*dt) for q in qs}
        rows.append(dict(a=a, logMmax=float(np.log(Mz.max())), **{f'J{q}': J[q] for q in qs}))
        print(json.dumps({k: (round(v, 5) if isinstance(v, float) else v) for k, v in rows[-1].items()}), flush=True)
    res[str(L)] = rows
    # growth rates over the last third of a-range, and partial sums for kappa = 8/7, 1.6
    A = np.array([r['a'] for r in rows]); sel = A >= A.max()*2/3
    for q in qs:
        J = np.array([max(r[f'J{q}'], 1e-300) for r in rows]); g = np.polyfit(A[sel], np.log(J[sel]), 1)[0]
        out = dict(L=L, q=q, growth=round(float(g), 3), q_theta=round(q*m.th, 3))
        for kap in (8/7, 1.6):
            ps = np.cumsum(np.exp(-kap*A)*J); out[f'partial_k{kap:.3f}'] = [float(f'{x:.3g}') for x in ps[[len(ps)//3, 2*len(ps)//3, -1]]]
        print('SUMMARY', json.dumps(out), flush=True); res[f'summary_{L}_{q}'] = out
if len(sys.argv) == 1: json.dump(res, open('/Users/adriansager/playground/riemann/Riemann/Deep/UPAudit/logs/nd_moments.json', 'w'), indent=1)
