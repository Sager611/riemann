# explore44 c2: Bohr-lift Monte Carlo at sigma for zeta and Davenport-Heilbronn (DH = c1 L(s,chi) + conj(c1) L(s,chibar), chi mod 5, chi(2) = i).
import numpy as np, json, time, sys
t0 = time.time(); rng = np.random.default_rng(44)
def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
PMAX = 30000
P = primes_upto(PMAX); Pf = P.astype(float)
chi = {1: 1, 2: 1j, 3: -1j, 4: -1, 0: 0}
cp = np.array([chi[int(p) % 5] for p in P]); 
kk = (np.sqrt(10 - 2*np.sqrt(5)) - 2)/(np.sqrt(5) - 1); c1 = (1 - 1j*kk)/2
Pbig = primes_upto(10**7).astype(float); Pbig = Pbig[Pbig > PMAX]
res = {}
for sig in [0.6, 0.7, 0.8]:
    x = Pf**(-sig)
    s2 = float(np.sum(Pbig**(-2*sig)))           # tail variance E|sum_{p>PMAX} p^-sig z_p|^2 (to 1e7; beyond negligible here)
    NZ, ND, CH = 20000, 40000, 100
    lz = []; ld = []; ldarg = []
    for i in range(0, NZ, CH):
        z = np.exp(2j*np.pi*rng.random((CH, len(P))))
        g = -np.log(np.abs(1 - x*z)).sum(1) + rng.normal(0, np.sqrt(s2/2), CH)
        lz.append(g)
    for i in range(0, ND, CH):
        z = np.exp(2j*np.pi*rng.random((CH, len(P))))
        G1 = -np.log(1 - cp*x*z).sum(1); G2 = -np.log(1 - np.conj(cp)*x*z).sum(1)
        m = (P % 5 != 0)
        # tails: independent complex Gaussians, E|.|^2 = s2 each (sum chi^2(p) p^-2sig ~ 0)
        G1 = G1 + (rng.normal(0, np.sqrt(s2/2), CH) + 1j*rng.normal(0, np.sqrt(s2/2), CH))
        G2 = G2 + (rng.normal(0, np.sqrt(s2/2), CH) + 1j*rng.normal(0, np.sqrt(s2/2), CH))
        F = c1*np.exp(G1) + np.conj(c1)*np.exp(G2)
        ld.append(np.log(np.abs(F)))
    lz = np.concatenate(lz); ld = np.concatenate(ld)
    rs = [0.5, 1, 1.5, 2, 3, 4, 5, 6, 7, 8]
    r = dict(tail_var=s2, zeta_mean=float(lz.mean()), zeta_mean_se=float(lz.std()/np.sqrt(NZ)), zeta_var=float(lz.var()),
             dh_mean=float(ld.mean()), dh_mean_se=float(ld.std()/np.sqrt(ND)), dh_var=float(ld.var()),
             zeta_low={str(q): float(np.mean(lz < -q)) for q in rs}, zeta_up={str(q): float(np.mean(lz > q)) for q in rs},
             dh_low={str(q): float(np.mean(ld - ld.mean() < -q)) for q in rs}, dh_low_abs={str(q): float(np.mean(ld < -q)) for q in rs},
             zeta_q=[float(v) for v in np.quantile(lz, [1e-3, 1e-2, 0.5, 0.99, 0.999])], dh_q=[float(v) for v in np.quantile(ld, [1e-4, 1e-3, 1e-2, 0.5, 0.99])])
    res[str(sig)] = r
    np.save('/Users/adriansager/playground/riemann/.scratch/explore44/lift_zeta_%g.npy' % sig, lz.astype(np.float32))
    np.save('/Users/adriansager/playground/riemann/.scratch/explore44/lift_dh_%g.npy' % sig, ld.astype(np.float32))
    print('sigma', sig, 'tail var', round(s2, 5), '| zeta: mean %.4f +- %.4f var %.4f | DH: mean (Jessen phi) %.4f +- %.4f var %.4f' % (r['zeta_mean'], r['zeta_mean_se'], r['zeta_var'], r['dh_mean'], r['dh_mean_se'], r['dh_var']), round(time.time() - t0, 1), flush=True)
    print('   P(log|zeta| < -r):', {k: round(v, 5) for k, v in r['zeta_low'].items()})
    print('   P(log|DH| < -r):  ', {k: round(v, 5) for k, v in r['dh_low_abs'].items()})
json.dump(res, open('/Users/adriansager/playground/riemann/.scratch/explore44/c2.json', 'w'))
print('done', round(time.time() - t0, 1))
