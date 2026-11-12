# explore44 c1/c2: smallest case (one prime) exact MGF vs Herbst; Lipschitz constants L(sigma); Bohr-lift Monte Carlo for zeta and DH.
import numpy as np, mpmath as mp, json, time, sys
t0 = time.time()
def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
out = {}
# ---- c1: one prime p = 2 (smallest case)
print('c1: one prime, f = -log|1 - a e^{i th}|, a = 2^-sigma')
th = np.linspace(0, 2*np.pi, 200001)[:-1]
for sig in [0.6, 0.75]:
    a = 2.0**(-sig)
    f = -np.log(np.abs(1 - a*np.exp(1j*th)))
    fp = -a*np.sin(th)/(1 - 2*a*np.cos(th) + a*a)
    Lre = a/(1 - a*a); var_exact = sum(a**(2*n)/(2*n*n) for n in range(1, 400))
    print(' sigma', sig, 'mean f', f.mean(), 'max|f\'| num', np.abs(fp).max(), 'a/(1-a^2)', Lre, 'var num', f.var(), 'exact', var_exact)
    rows = []
    for lam in [0.5, 1, 2, 4, 8]:
        exact = float(mp.hyp2f1(lam/2, lam/2, 1, a*a)); num = np.exp(lam*f).mean()
        herbst = np.exp(lam*lam*Lre*Lre/2); gauss = np.exp(lam*lam*var_exact/2)
        rows.append((lam, exact, num, herbst, gauss))
        print('  lam %4.1f  E e^{lam f}: hyp2f1 %.6f num %.6f | Herbst bound %.6f | Gaussian(var) %.6f | Herbst ok %s' % (lam, exact, num, herbst, gauss, exact <= herbst))
    # entropy check of LSI for g = e^{lam f/2}: Ent(g^2) <= 2 E|g'|^2
    for lam in [1, 4, 8]:
        g2 = np.exp(lam*f); Eg2 = g2.mean(); ent = (g2*np.log(g2)).mean() - Eg2*np.log(Eg2)
        dg2 = ((lam/2)*fp)**2*g2
        print('  LSI lam', lam, 'Ent(g^2) %.6f <= 2E|g\'|^2 %.6f : %s' % (ent, 2*dg2.mean(), ent <= 2*dg2.mean()))
    out['c1_%g' % sig] = rows
# ---- L(sigma)^2 with primes to 1e7 + integral tail
P = primes_upto(10**7).astype(float)
print('primes', len(P), round(time.time() - t0, 1))
Ls = {}
for sig in [0.55, 0.6, 0.7, 0.75, 0.8, 0.9, 1.0]:
    x = P**(-sig)
    Lre2 = np.sum((x/(1 - x*x))**2); Lim2 = np.sum((x/(1 - x))**2)
    # tail p > 1e7: sum p^{-2 sig} ~ int_{1e7}^inf u^{-2sig}/log u du
    tail = float(mp.quad(lambda u: u**(-2*sig)/mp.log(u), [1e7, 1e9, 1e12, mp.inf]))
    var = 0.5*np.sum(x*x) + 0.25*np.sum(x**4)/2 + 0.5*tail   # Var log|zeta_lift| = (1/2) sum_p sum_n p^{-2n sig}/n^2
    Ls[sig] = dict(Lre2=Lre2 + tail, Lim2=Lim2 + tail, var=var)
    print(' sigma %.2f  L_re^2 %.5f  L_arg^2 %.5f  Var log|zeta_lift| %.5f  ratio L^2/var %.3f  tail(p>1e7) %.5f' % (sig, Lre2 + tail, Lim2 + tail, var, (Lre2 + tail)/var, tail))
out['L'] = {str(k): v for k, v in Ls.items()}
json.dump(out, open('/Users/adriansager/playground/riemann/.scratch/explore44/c1.json', 'w'), default=float)
print('done', round(time.time() - t0, 1))
