
import mpmath as mp, itertools
mp.mp.dps = 150; w = None
for eta, f, L in itertools.product(('1e-12','1e-8','1e-3','0.1','0.3','0.5'), ('1','-1','0.5','-0.5','0.99','-0.99'), ('20','23.03','60','111.1')):
    eta = mp.mpf(eta); u = eta*mp.mpf(f); L = mp.mpf(L)
    D = u/L - mp.log(1 + mp.log(1 + u)/L)
    b = eta**2/(2*L)*(1 + 2*eta/(3*(1 - eta)) + (1 + eta)**2/L + 2*eta*(1 + eta)**3/(3*L**2*(1 - eta*(1 + eta)/L)))
    r = 1 - D/b
    if w is None or r < w[0]: w = (r, float(eta), float(f), float(L))
print('Lemma 13: 144 cases, min (1 - D/bound) =', mp.nstr(w[0], 6), 'at (eta, u/eta, L) =', w[1:])
