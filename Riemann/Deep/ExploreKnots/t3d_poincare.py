# t3d_poincare.py — Lawrence–Zagier: the Poincare sphere Sigma(2,3,5) series chi_60 (odd, 8 residues) = (psi + conj psi)/2-type
# combination of two odd characters of conductor 60; its zeros right of the line, and those of the scalar-FE eigencombinations.
import numpy as np, mpmath as mp, math, json
exec(open('t5_family.py').read().split('single = []')[0])
src = open('t3_zeros.py').read()
exec(src[src.index('# Bernoulli numbers'):src.index('# sanity')])
mp.mp.dps = 30
N, f60 = brieskorn((2, 3, 5))
units, chars = char_table(N)
comps = [(sum(f60[n]*np.conj(c[n]) for n in units)/len(units), c) for c in chars]
comps = [(a, c) for a, c in comps if abs(a) > 1e-9]
(a1, c1), (a2, c2) = comps
print('coefficients:', np.round(a1, 6), np.round(a2, 6), ' conjugate pair:', np.allclose(c1, np.conj(c2)), ' odd:', np.isclose(c1[N-1], -1))
psi = c1
tau = sum(psi[n]*np.exp(2j*np.pi*n/N) for n in range(N)); eps = tau/(1j*np.sqrt(N))
print('root number eps =', np.round(eps, 12), '|eps| =', abs(eps))
G = 2*np.real(eps**-0.5*psi); H = 2*np.imag(eps**-0.5*psi)
fam = {'Poincare_chi60': f60, 'Poincare_G_plus': G, 'Poincare_H_minus': H}
# check FE numerically: Lambda(s) = (N/pi)^((s+1)/2) Gamma((s+1)/2) L(s) ; G(1-s) = G(s), H(1-s) = -H(s)
def Lam(s, a): return (mp.mpf(N)/mp.pi)**((s+1)/2)*mp.gamma((s+1)/2)*L_mp(s, N, a)
s0 = mp.mpc(0.3, 7.1)
print('FE check G:', mp.nstr(abs(Lam(s0, G)/Lam(1-s0, G)), 12), ' H:', mp.nstr(Lam(s0, H)/Lam(1-s0, H), 12))
res = {}
sig = np.arange(0.53, 1.80, 0.02); tt = np.arange(0.5, 60, 0.04)
for name, a in fam.items():
    Z = np.empty((len(sig), len(tt)))
    for i, sg in enumerate(sig):
        for j0 in range(0, len(tt), 4000): Z[i, j0:j0+4000] = np.abs(L_np(sg + 1j*tt[j0:j0+4000], N, a))
    zs = []
    for i in range(1, len(sig)-1):
        row = Z[i]
        m = (row[1:-1] < row[:-2]) & (row[1:-1] < row[2:]) & (row[1:-1] < Z[i-1, 1:-1]) & (row[1:-1] < Z[i+1, 1:-1])
        for j in np.nonzero(m)[0]+1:
            if Z[i, j] > 0.2: continue
            try: z = mp.findroot(lambda s: L_mp(s, N, a), mp.mpc(sig[i], tt[j]), tol=1e-20, maxsteps=50)
            except Exception: continue
            if abs(L_mp(z, N, a)) < 1e-15 and z.real > 0.5+1e-8 and 0 < z.imag < 60:
                zz = (round(float(z.real), 10), round(float(z.imag), 8))
                if all(abs(zz[0]-w[0])+abs(zz[1]-w[1]) > 1e-6 for w in zs): zs.append(zz)
    zs.sort(key=lambda w: w[1]); res[name] = zs
    print(f'{name}: {len(zs)} zeros with Re s > 1/2, t < 60; first:', zs[:5], flush=True)
json.dump(res, open('t3d_out.json', 'w'), indent=1)
