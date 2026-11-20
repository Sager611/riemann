"""gapjinf: collect validated h_J (form JSONs), crossings kappa_J by linear interpolation in kappa (3.358, 3.40),
and J -> infinity extrapolation of the increments (geometric and power-law fits).  usage: extrap.py"""
import json, os, glob
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
H = {}
for f in glob.glob(os.path.join(HERE, 'logs', 'form_k*_J*.json')):
    F = json.load(open(f)); H[(round(F['kappa_eff'], 4), F['J'])] = F['h']
ks = sorted(set(k for k, _ in H)); Js = sorted(set(j for _, j in H))
for k in ks:
    print('kappa %.4f: ' % k + '  '.join('J%d %.7f' % (j, H[(k, j)]) for j in Js if (k, j) in H))
def extrap(hs, Js):
    d = np.diff(hs); r = d[1:] / d[:-1]
    out = {'increments': list(d), 'ratios': list(r)}
    rr = r[-1] if len(r) else np.nan
    out['geom_tail'] = d[-1] * rr / (1 - rr) if 0 < rr < 1 else np.nan
    # power law d_J = C J^-alpha from the last two increments
    if len(d) >= 2 and d[-1] > 0 and d[-2] > 0:
        al = np.log(d[-2] / d[-1]) / np.log(Js[-1] / Js[-2])
        C = d[-1] * Js[-1] ** al
        out['alpha'] = al
        out['power_tail'] = sum(C * J ** -al for J in range(Js[-1] + 1, 4000)) if al > 1 else np.inf
    return out
for k in ks:
    Jk = [j for j in Js if (k, j) in H and j >= 4]
    if len(Jk) >= 3:
        hs = np.array([H[(k, j)] for j in Jk])
        e = extrap(hs, Jk)
        print('kappa %.4f, J=%s: increments %s ratios %s; tail geom %.5f, power (alpha %.2f) %.5f -> h_inf in [%.5f, %.5f]' % (
            k, Jk, ' '.join('%.5f' % x for x in e['increments']), ' '.join('%.3f' % x for x in e['ratios']), e['geom_tail'],
            e.get('alpha', np.nan), e.get('power_tail', np.nan), hs[-1] + e['geom_tail'], hs[-1] + e.get('power_tail', np.nan)))
if len(ks) >= 2:
    k1, k2 = ks[0], ks[-1]
    for j in Js:
        if (k1, j) in H and (k2, j) in H:
            h1, h2 = H[(k1, j)], H[(k2, j)]
            sl = (h2 - h1) / (k2 - k1)
            print('J=%d: dh/dkappa = %.4f, crossing kappa_J = %.4f' % (j, sl, k1 + (1 - h1) / sl))
