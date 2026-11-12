# t1_sign_barta.py (explore39): (a) sign structure of the reduced jump kernels; (b) Barta (Allegretto-Piepenbrink) lower bounds
# lambda_0 >= inf_x (H phi)(x)/phi(x) in the Markov ranges (even: c <= rho; odd: c <= 2), phi = Galerkin ground state.
import sys, numpy as np, math, json, time
sys.path.insert(0, '.scratch/explore39'); from xform import Form
rho = 1.324717957244746
print('ln rho = %.6f' % math.log(rho))
# (a) kernel signs on a grid of [0,a]^2 (off-diagonal), continuous parts
for kind in ['zeta', 'dh']:
    for c in [1.2, rho, 1.4, 1.755, 2.0, 2.5, 3.0, 5.0, 10.0]:
        F = Form(c, kind); a = F.a; xs = np.linspace(0, a, 301)[1:]; X, Y = np.meshgrid(xs, xs, indexing='ij'); off = np.abs(X - Y) > 1e-9
        je = F.jev(X, Y)[off]; jo = F.jodd(X, Y)[off]
        neg_e = float(np.mean(je < 0)); neg_o = float(np.mean(jo < 0))
        at = [(round(math.exp(h)), round(w, 4)) for h, w in F.atoms]
        anti_e = [n for n, w in at if w < 0]; anti_o_dir = [n for n, w in at if w < 0]; anti_o_ref = [n for n, w in at if w > 0]
        print(json.dumps({'kind': kind, 'c': round(c, 6), 'min_jev': round(float(je.min()), 5), 'frac_jev<0': round(neg_e, 4),
                          'min_jodd': round(float(jo.min()), 5), 'atoms': at[:8], 'even_antiMarkov_atoms': anti_e,
                          'odd_antiMarkov_direct': anti_o_dir, 'odd_antiMarkov_reflected': anti_o_ref}), flush=True)
# (b) Barta
def barta(F, sec, n=24, flip=False, npts=500):
    Qm, p = F.Q(n, sec)
    if flip:
        Fm, info = F.flip(n, sec); Qm = Qm - Fm
    e, V = np.linalg.eigh(Qm); v = V[:, 0]
    a = F.a; u = (np.arange(npts) + 0.5)/npts; xs = a*np.sin(np.pi*u/2)
    if sec == 'odd': xs = xs[xs > 1e-4*a]
    Hp, phi = F.Hphi(v, sec, xs, flip=flip)
    ph = phi(xs)
    if np.median(ph) < 0: ph = -ph; Hp = -Hp
    ok = bool(np.all(ph > 0))
    r = Hp/ph; k = int(np.argmin(r))
    return {'gal': float(e[0]), 'gal1': float(e[1]), 'barta': float(r.min()), 'argmin_x/a': round(float(xs[k]/a), 4), 'phi>0': ok, 'phi_min/max': round(float(ph.min()/ph.max()), 4), 'r_median': float(np.median(r))}
t0 = time.time()
for c in [1.1, 1.2, 1.3, rho]:
    F = Form(c)
    for sec in ['even', 'odd']:
        for n in [16, 24]:
            out = barta(F, sec, n); out.update({'c': round(c, 6), 'sec': sec, 'n': n}); print(json.dumps(out), flush=True)
for c in [1.5, 1.8, 2.0]:
    F = Form(c)
    for n in [16, 24]:
        out = barta(F, 'odd', n); out.update({'c': c, 'sec': 'odd', 'n': n}); print(json.dumps(out), flush=True)
print('sec', round(time.time() - t0, 1))
