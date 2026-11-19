# Numerical test of Lemmas 2-7 for COMPLEX characters (Littlewood78's sanity.py did real ones only): for every nonprincipal
# chi mod q in {5, 7, 11, 13} compute exactly L(1,chi) = -(1/q) sum chi(a) psi(a/q), L'(1,chi) = (1/q) sum chi(a)(log q psi(a/q)
# - gamma_1(a/q)), Z by Lemma 4, the smoothed prime sums P_0 (L=10, h=1/2) and P_1 (L_1=8, h_1=1/2), and check
#   (W1) residual r1 := Re L'/L(1) + Re P_1 = -sum Re V_1(rho-1) - sum_triv V_1  lies in [-kappa_1 Z, kappa_1 Z + T_1];
#   (W0) residual r0 := log|L(1)| - Re P_0 = sum Re G(rho) + sum_triv G  lies in [-kappa_0 Z - T_0, kappa_0 Z];  Z > 0.
# Run: uv run --no-project --with numpy --with mpmath python3 cplxsanity.py
import math, numpy as np, mpmath as mp
mp.mp.dps = 25
N = int(math.exp(15)) + 1
s = np.ones(N + 1, dtype=bool); s[:2] = False
for p in range(2, int(N**0.5) + 1):
    if s[p]: s[p*p::p] = False
ps = np.nonzero(s)[0]; del s
nn, lg = [], []
for p in ps.tolist():
    pk = p
    while pk <= N: nn.append(pk); lg.append(math.log(p)); pk *= p
nn = np.array(nn, dtype=np.int64); lam = np.array(lg); ln = np.log(nn.astype(float))
def v(u, h): return np.where(u <= 1, 1.0, np.where(u < 1 + h, (1 + h - u)/h, 0.0))
B = 0.87338196661320752; b = 1 - B
L, h, L1, h1 = 10.0, 0.5, 8.0, 0.5
g0 = lambda be: (math.exp(-(1+h)*L*(1-be))/(1+h) + math.exp(-L*(1-be)))/(h*L*L)
g1 = lambda be: (math.exp(-(1+h1)*L1*(1-be)) + math.exp(-L1*(1-be)))/(h1*L1)
k0 = g0(B)*B/b + g0(0.5); k1 = g1(B)*B/b + g1(0.5)
T0 = math.exp(-L)/(L*(1 - math.exp(-2*L))); T1 = math.exp(-L1)/(1 - math.exp(-2*L1))
w0 = lam*v(ln/L, h)/(nn*ln); w1 = lam*v(ln/L1, h1)/nn
print(f'kappa_0 = {k0:.5f}, kappa_1 = {k1:.5f}, T_0 = {T0:.2e}, T_1 = {T1:.2e}; {len(nn)} prime powers <= e^15')
worst0 = worst1 = 0.0; cnt = 0; bad = 0
for q in (5, 7, 11, 13):
    g = next(g for g in range(2, q) if all(pow(g, (q-1)//f, q) != 1 for f in range(2, q) if (q-1) % f == 0 and all(f % d for d in range(2, f))))
    ind = {pow(g, j, q): j for j in range(q - 1)}
    r = (nn % q).astype(np.int64); indarr = np.array([ind.get(int(t), -1) for t in range(q)])
    for k in range(1, q - 1):
        chi_n = np.where(r == 0, 0, np.exp(2j*np.pi*k*indarr[r]/(q - 1)))
        P0 = complex(np.sum(chi_n*w0)); P1 = complex(np.sum(chi_n*w1))
        chi = lambda a: mp.e**(2j*mp.pi*k*ind[a]/(q - 1))
        L1v = -sum(chi(a)*mp.digamma(mp.mpf(a)/q) for a in range(1, q))/q
        dL1 = sum(chi(a)*(mp.log(q)*mp.digamma(mp.mpf(a)/q) - mp.stieltjes(1, mp.mpf(a)/q)) for a in range(1, q))/q
        a_par = k % 2                                     # chi(-1) = chi(g^{(q-1)/2}) = (-1)^k
        Z = float(mp.re(dL1/L1v)) + 0.5*math.log(q/math.pi) + 0.5*float(mp.digamma((1 + a_par)/2))
        r1 = float(mp.re(dL1/L1v)) + P1.real; r0 = math.log(float(abs(L1v))) - P0.real
        ok = Z > 0 and -k1*Z <= r1 <= k1*Z + T1 and -k0*Z - T0 <= r0 <= k0*Z
        bad += not ok; cnt += 1
        worst1 = max(worst1, abs(r1)/(k1*Z)); worst0 = max(worst0, abs(r0)/(k0*Z))
        if k in (1, (q - 1)//2): print(f'q={q} k={k} ({"real" if 2*k == q-1 else "complex"}, a={a_par}): |L(1)|={float(abs(L1v)):.6f} Z={Z:.5f} r1={r1:+.2e} (<= k1 Z = {k1*Z:.2e}) r0={r0:+.2e} (<= k0 Z = {k0*Z:.2e}) {"OK" if ok else "FAIL"}')
print(f'{cnt} characters, {bad} failures; max |r1|/(kappa_1 Z) = {worst1:.3f}, max |r0|/(kappa_0 Z) = {worst0:.3f} (GRH zeros: small)')
