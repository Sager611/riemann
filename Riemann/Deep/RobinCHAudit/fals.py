#!/usr/bin/env python3
# fals.py (robinchaudit) -- falsification tests of RobinCH Lemma 11_k ((B2_k), (B1'), m_c = E|V - c|, MAD_k) and Lemma 13
# at corner parameters, by a method different from RobinCH/test_b2k.py: F_phi(rho) = E F(rho; L' + V) expanded in Taylor
# series of G(y) = e^{-sy} F(rho; y) around L' (exact E_n closed forms) against the moments E[V^n e^{sV}] = d^n/ds^n M(s),
# M(s) = ((e^{sa} - 1)/(sa))^k. Non-rigorous (mpmath, 60 digits); a ratio true/bound > 1 would falsify.
#   uv run --no-project --with mpmath python3 fals.py
import itertools, math, time
import mpmath as mp
mp.mp.dps = 60
t0 = time.time()

def Fy(s, y):       # F(rho; y) = int_y^inf e^{su} (1/u + 1/u^2) du = E_1(-sy) + E_2(-sy)/y
    z = -s*y
    return mp.expint(1, z) + mp.expint(2, z)/y

def Ij(j, z):        # m^{(j)}(z) = int_0^1 t^j e^{zt} dt
    if abs(z) < 1:
        return sum(z**i/(mp.factorial(i)*(i + j + 1)) for i in range(90))
    v = (mp.exp(z) - 1)/z
    for i in range(1, j + 1): v = (mp.exp(z) - i*v)/z
    return v

def Mders(s, a, k, N):   # [M^{(n)}(s) for n <= N], M(s) = m(sa)^k: Taylor series of m around z0 = sa, power k, exact
    z0 = s*a
    c = [Ij(j, z0)/mp.factorial(j) for j in range(N + 1)]
    p = [mp.mpc(1)] + [mp.mpc(0)]*N
    for _ in range(k):
        p = [sum(p[i]*c[n - i] for i in range(n + 1)) for n in range(N + 1)]
    return [p[n]*mp.factorial(n)*a**n for n in range(N + 1)]

def Fphi(s, Lp, lam, k, N=10):   # E F(rho; L'+V) = sum_n (-1)^n [L'^-n E_{n+1}(-sL') + (n+1) L'^{-n-1} E_{n+2}(-sL')] M^{(n)}(s)
    a = mp.mpf(lam)/k
    Md = Mders(s, a, k, N)
    z = -s*Lp; tot = mp.mpc(0); last = 0
    for n in range(N + 1):
        t = (-1)**n*(Lp**(-n)*mp.expint(n + 1, z) + (n + 1)*Lp**(-n - 1)*mp.expint(n + 2, z))*Md[n]/mp.factorial(n)
        tot += t; last = t
    return tot, last

def MADk(k):
    return 2*sum((-1)**j*math.comb(k, j)*(mp.mpf(k)/2 - j)**(k + 1) for j in range(k//2 + 1))/mp.factorial(k + 1)

worst2 = (0, None); worst1 = (0, None); n = 0; maxtrunc = (0, None); close = []
for k, gam, L, beta, fac, cen in itertools.product((1, 2, 40, 57, 120), ('14.134725', '3000000000001', '1e16'),
                                                    ('16.11809565', '111.1'), ('0.5', '0.87338196661320752', '0.999'),
                                                    (math.pi, 2.0, 4.0, 0.0), (False, True)):
    gam = mp.mpf(gam); L = mp.mpf(L); beta = mp.mpf(beta)
    lam = min(fac*k/gam, mp.mpf('1e-3')) if fac > 0 else mp.mpf('1e-3')
    if fac > 0 and gam < 100: lam = mp.mpf(fac)*k/gam       # deliberately outside lam <= 1e-3 (the proof does not use it)
    s = beta - 1 + 1j*gam; c = lam/2 if cen else 0; Lp = L - c
    if Lp <= 0 or lam/Lp > 0.05: continue                                  # Taylor series in lam/L' would converge too slowly
    Fp, last = Fphi(s, Lp, lam, k)
    xp = mp.exp((beta - 1)*Lp); H = 1/Lp + 1/Lp**2; Hp = 1/Lp**2 + 2/Lp**3; S = abs(s)
    b2 = xp*(H*(2*k/(lam*S))**k/S + Hp*lam/(2*S) + 2*Hp/S**2)
    mc = lam*MADk(k)/k if cen else lam/2
    b1 = xp*H*mc
    r2 = abs(Fp)/b2; r1 = abs(Fp - Fy(s, L))/b1; n += 1
    case = (k, float(gam), float(L), float(beta), round(float(lam*gam/k), 3), cen)
    if abs(last)/abs(Fp) > maxtrunc[0]: maxtrunc = (abs(last)/abs(Fp), case)
    if r2 > 1 - mp.mpf('1e-6'): close.append((r2, case, abs(last)/abs(Fp)))
    if r2 > worst2[0]: worst2 = (r2, case)
    if r1 > worst1[0]: worst1 = (r1, case)
print(f'(B2_k): {n} cases, worst |F_phi|/bound = {mp.nstr(worst2[0], 8)} (1 - ratio = {mp.nstr(1 - worst2[0], 6)}) at (k, gamma, L, beta, lam*gamma/k, centred) = {worst2[1]}')
close.sort(key=lambda t: -t[0]); print(f'   {len(close)} cases with 1 - ratio < 1e-6 (all k <= 40, lam*gamma/k = pi); tightest 3:')
for r2, case, tr in close[:3]: print(f'   near-tight: 1 - ratio = {mp.nstr(1 - r2, 6)} at {case}, truncation {mp.nstr(tr, 3)}')
print(f"(B1'): worst |F_phi - F|/bound = {mp.nstr(worst1[0], 8)} at {worst1[1]};  max last Taylor term / |F_phi| = {mp.nstr(maxtrunc[0], 3)} at {maxtrunc[1]}; {time.time()-t0:.0f}s")
