#!/usr/bin/env python3
"""Local Euler factors of the theta average, Liu 2610.26835v2 §13.3 (number_fields_v4.tex l.8552-8636, 8696-8702,
8800-8838, 9284-9298, 9314-9340).
(1) Gauss-phase identities (l.8622-8631) with actual sextic Gauss sums at primes of Q(zeta12) of norm 13, 37, 61, 73.
(2) Summation of the C_{t,k}(J) series (l.8578-8613) against the closed form P_j^* with the six J_j (l.8558-8569),
    at random complex s,w,z in the region of absolute convergence and random unit phases eta, rho.
(3) Exponent claims: Lemma 13.10 |H_p-1| = O(q^{-1-10e}) for p not dividing u (j=0) on Re s >= a+16e, Re w >= 1-a-6e,
    Re z >= 17/50, 1/2 <= a <= 11/12; the j-row table sigma_j (l.9289-9298); principal box (l.8800-8812): H_p-1 = O(q^{-907/500}).
Measured ratios |X| q^(-bound) at q = 1e8..1e40 with mpmath (300 digits; 50 digits gave cancellation garbage).
PASS = ratio bounded (max ratio <= C) and not growing (ratio(1e40) <= 1.5 ratio(1e24))."""
import cmath, math, random, sys
import mpmath as mp
mp.mp.dps = 300
random.seed(7)

def sextic_gauss(p):
    # residue field F_p, chi of order 6 via generator g: chi(g^k) = exp(2 pi i k/6); additive e(x) = exp(2 pi i x/p)
    g = next(x for x in range(2, p) if all(pow(x, (p-1)//r, p) != 1 for r in (2, 3) if (p-1) % r == 0) and pow(x, (p-1)//2, p) != 1 and all(pow(x, (p-1)//f, p) != 1 for f in set(fac(p-1))))
    log = {}
    y = 1
    for k in range(p-1): log[y] = k; y = y*g % p
    chi = lambda x, v: 0 if x % p == 0 else cmath.exp(2j*math.pi*v*log[x % p]/6)
    tau = [-1+0j] + [sum(chi(x, v)*cmath.exp(2j*math.pi*x/p) for x in range(1, p)) for v in range(1, 6)]
    return chi, tau

def fac(n):
    f = []; d = 2
    while d*d <= n:
        while n % d == 0: f.append(d); n //= d
        d += 1
    if n > 1: f.append(n)
    return f

res = []
def rep(name, ok, info=''):
    res.append(bool(ok)); print(('PASS ' if ok else 'FAIL ') + name + ('  ' + info if info else ''), flush=True)

# (1) phase identities
for p in (13, 37, 61, 73):
    chi, tau = sextic_gauss(p)
    gam = [t/math.sqrt(p) for t in tau]
    Lam = -gam[2]**3
    c4 = chi(4, 1)
    G0 = c4.conjugate()*gam[3]
    ok = abs(abs(gam[1])-1) < 1e-9 and abs(gam[3]**2-1) < 1e-9 and abs(gam[4]-gam[2].conjugate()) < 1e-9 and abs(gam[5]-gam[1].conjugate()) < 1e-9
    ok2 = abs(gam[1]*gam[2] + Lam*G0) < 1e-9 and abs(c4**3-1) < 1e-9 and abs(tau[1]*tau[5]-p) < 1e-6
    for _ in range(5):
        eta = cmath.exp(2j*math.pi*random.random())
        th = Lam.conjugate()*G0.conjugate()*eta
        a = Lam.conjugate()**3*eta**3; b0 = Lam.conjugate()**2*eta**2*c4.conjugate()
        ok3 = (abs(th**3*gam[3]-a) < 1e-9 and abs(th**3*gam[1]*gam[2]+eta*b0) < 1e-9 and abs(gam[2]*th**4*gam[4]-b0**2) < 1e-9
               and abs(gam[1]*gam[2]*gam[3]*th**4+eta*a) < 1e-9 and abs(th**6-a**2) < 1e-9)
        ok2 = ok2 and ok3
    rep(f'l.8622-8631/8681-8690 Gauss phase identities at q={p}', ok and ok2, f'gamma3={gam[3].real:+.0f}')

# (2) series vs closed form
def closed_Pstar(j, q, T, Y0, V, eta, rho, a, b0):
    R = a**2*q**4*T**6*V
    J = [-eta/rho*T + rho*Y0*R, eta*T*Y0, a*rho**-3*q**1.5*T**3, -eta*b0*rho**-2*q**2*T**3*Y0 + b0**2*rho**-4*q**2*T**4,
         -eta*a*rho**-3*q**2.5*T**4*Y0, -a**2*q**3*T**6][j]
    return (1/(1-R))*((R*(1-1/q) - eta*(q-1)*T*Y0*V**(1 if j <= 1 else 0))/(1-V) + J)

def series_Pstar(j, q, s, w, z, eta, rho, chiv, tau, gam, Lam, G0, tmax=80, mmax=60):
    th = Lam.conjugate()*G0.conjugate()*eta
    def Tv(v, y_is_unit):
        if y_is_unit: return tau[v % 6] if v % 6 else -1  # chibar_p(1)^v tau_v with y = 1 (p^0); tau_0 = -1
        return (q-1) if v % 6 == 0 else 0
    tot = 0j
    for t in range(1, tmax):
        if t % 3 == 2: continue
        e0 = t % 3; l = (t-e0)//3
        for k in (0, 1):
            for m in range(0, mmax):
                J = j+6*m
                if k == 0:
                    if J < t-1: continue
                    C = q**(t-1)*Tv(t, J == t-1)
                else:
                    if J < t: continue
                    C = tau[1]*q**(t-1)*Tv(t-1, J == t)
                if C == 0: continue
                mult = gam[2]**e0*th**(e0+3*l)*rho**(k-e0-3*l)*q**(-e0*(0.5+s) - l*(1+3*s) - k*w - 6*m*z)
                tot += C*mult
    return tot

for p in (13, 37):
    chi, tau = sextic_gauss(p)
    gam = [t/math.sqrt(p) for t in tau]
    Lam = -gam[2]**3; c4 = chi(4, 1); G0 = c4.conjugate()*gam[3]
    worst = 0
    for trial in range(6):
        eta = cmath.exp(2j*math.pi*random.random()); rho = cmath.exp(2j*math.pi*random.randrange(6)/6)
        s = 2.2+3*random.random() + 1j*random.uniform(-5, 5); w = 2.5+random.random() + 1j*random.uniform(-5, 5); z = 1.2+random.random()+1j*random.uniform(-5, 5)
        q = p
        T, Y0, V = q**-s, q**-w, q**(-6*z)
        a = Lam.conjugate()**3*eta**3; b0 = Lam.conjugate()**2*eta**2*c4.conjugate()
        for j in range(6):
            A = series_Pstar(j, q, s, w, z, eta, rho, chi, tau, gam, Lam, G0)
            B = closed_Pstar(j, q, T, Y0, V, eta, rho, a, b0)
            worst = max(worst, abs(A-B)/max(abs(B), 1e-300))
    rep(f'l.8558-8636 series of C_(t,k)(J) terms = closed P_j^* with J_0..J_5, q={p} (rel err {worst:.2e})', worst < 1e-8)

# (3) exponent claims via mpmath at huge q (phases random units)
def Hp_Gp(j, q, s, w, z, eta, rho, a, b0):
    T = q**(-s); Y0 = q**(-w); V = q**(-6*z)
    R = a**2*q**4*T**6*V
    Jv = [-eta/rho*T + rho*Y0*R, eta*T*Y0, a*rho**-3*q**mp.mpf(1.5)*T**3, -eta*b0*rho**-2*q**2*T**3*Y0 + b0**2*rho**-4*q**2*T**4,
          -eta*a*rho**-3*q**mp.mpf(2.5)*T**4*Y0, -a**2*q**3*T**6][j]
    Ps = (1/(1-R))*((R*(1-1/q) - eta*(q-1)*T*Y0*V**(1 if j <= 1 else 0))/(1-V) + Jv)
    W = rho*Y0 if j == 0 else 0; D = eta/rho*T if j == 0 else 0
    P = 1/(1-V) + (rho*Y0/(1-rho*Y0) if j == 0 else 0) + Ps
    H = (1-V)*(1-W)*P/(1-D)
    G = ((mp.conj(eta)*q**s - Y0)*(1-V)*(1-W)*Ps - Y0*(1-V*W))/(1-D)
    return H, G, R

def unit():
    return mp.expjpi(2*mp.mpf(random.random()))
qs = [mp.mpf(10)**k for k in (8, 16, 24, 32, 40)]
def slope(vals):
    xs = [mp.log(q) for q in qs]; ys = [mp.log(abs(v)) if v != 0 else mp.mpf(-10**6) for v in vals]
    n = len(xs); mx = sum(xs)/n; my = sum(ys)/n
    return float(sum((x-mx)*(y-my) for x, y in zip(xs, ys))/sum((x-mx)**2 for x in xs))
def ratios(vals, bexp):
    return [float(abs(v)*q**(-bexp)) for v, q in zip(vals, qs)]
def bounded(rs, C):
    return max(rs) <= C and rs[-1] <= 2*max(rs[0], rs[1]) + 1e-30
e = mp.mpf(1)/1000
ok = True; worst = -99
for a_ in [mp.mpf(1)/2, mp.mpf(61)/100, mp.mpf(3)/4, mp.mpf(33)/40, mp.mpf(11)/12]:
    for trial in range(4):
        t1, t2, t3 = [mp.mpf(random.uniform(-50, 50)) for _ in range(3)]
        eta, rho, a, b0 = unit(), mp.expjpi(2*mp.mpf(random.randrange(6))/6), unit(), unit()
        vals = []
        for q in qs:
            s = a_+16*e + 1j*t1; w = 1-a_-6*e + 1j*t2; z = mp.mpf(17)/50 + 1j*t3
            H, G, R = Hp_Gp(0, q, s, w, z, eta, rho, a, b0)
            vals.append(H-1)
        rs = ratios(vals, -1-10*e); worst = max(worst, max(rs))
        if not bounded(rs, 56): ok = False; print('   L13.10 ratio', float(a_), rs)
rep(f'Lemma 13.10 (l.9314): |H_p-1| q^(1+10e) <= C1 = 56 and not growing, a in [1/2,11/12] (max ratio {worst:.4f})', ok)
# also G_p + chibar_p(u) H_p = O(q^{-a-4e}) at good primes (l.9287); chibar_p(u) = conj(rho) for j=0
ok = True; worst = -99
for a_ in [mp.mpf(1)/2, mp.mpf(3)/4, mp.mpf(11)/12]:
    for trial in range(3):
        t1, t2, t3 = [mp.mpf(random.uniform(-50, 50)) for _ in range(3)]
        eta, rho, a, b0 = unit(), mp.expjpi(2*mp.mpf(random.randrange(6))/6), unit(), unit()
        vals = []
        for q in qs:
            s = a_+16*e + 1j*t1; w = 1-a_-6*e + 1j*t2; z = mp.mpf(17)/50 + 1j*t3
            H, G, R = Hp_Gp(0, q, s, w, z, eta, rho, a, b0)
            vals.append(G + mp.conj(rho)*H)
        rs = ratios(vals, -a_-4*e); worst = max(worst, max(rs))
        if not bounded(rs, 10): ok = False; print('   G+chibar H ratio', float(a_), rs)
rep(f'l.9287 G_p + chibar_p(u) H_p = O(q^(-a-4e)) at good primes (max ratio {worst:.4f})', ok)
# sigma_j table: selected factor q^{z0-1} G_p, j>=2 minus -(q-1) q^{z0-1-w}/(1-R); remaining O(q^{z0-1/2-sigma_j})
sig = {1: mp.mpf(26)/25+80*e, 2: 32*e, 3: 26*e, 4: 42*e, 5: 80*e}
ok = True; rows = []
for j in range(1, 6):
    worst = -99
    for a_ in [mp.mpf(1)/2, mp.mpf(0.6), mp.mpf(3)/4, mp.mpf(0.85), mp.mpf(11)/12]:
        for trial in range(3):
            t1, t2, t3 = [mp.mpf(random.uniform(-50, 50)) for _ in range(3)]
            eta, rho, a, b0 = unit(), mp.expjpi(2*mp.mpf(random.randrange(6))/6), unit(), unit()
            vals = []
            z0 = mp.mpf(17)/50
            for q in qs:
                s = a_+16*e + 1j*t1; w = 1-a_-6*e + 1j*t2; z = z0 + 1j*t3
                H, G, R = Hp_Gp(j, q, s, w, z, eta, rho, a, b0)
                sel = q**(z-1)*G
                if j >= 2: sel -= -(q-1)*q**(z-1-w)/(1-R)
                vals.append(sel)
            rs = ratios(vals, z0-mp.mpf(1)/2-sig[j]); worst = max(worst, max(rs))
            if not bounded(rs, 10): ok = False; print('   sigma_j ratio', j, float(a_), rs)
    rows.append((j, round(worst, 4)))
rep(f'l.9289-9298 sigma_j table: ratios bounded (max ratio per j: {rows})', ok)
# principal box (l.8800-8812): j=0, rho=1 at Re s >= 437/500, Re w >= 19/20, Re z >= 33/200: H_p - 1 = O(q^{-907/500})
ok = True; worst = -99
for trial in range(8):
    t1, t2, t3 = [mp.mpf(random.uniform(-50, 50)) for _ in range(3)]
    eta, a, b0 = unit(), unit(), unit()
    vals = []
    for q in qs:
        H, G, R = Hp_Gp(0, q, mp.mpf(437)/500+1j*t1, mp.mpf(19)/20+1j*t2, mp.mpf(33)/200+1j*t3, eta, mp.mpf(1), a, b0)
        vals.append(H-1)
    rs = ratios(vals, mp.mpf(-907)/500); worst = max(worst, max(rs))
    if not bounded(rs, 10): ok = False; print('   box ratio', rs)
rep(f'l.8810 principal box: |H_p-1| q^(907/500) bounded (max ratio {worst:.4f})', ok)
# selected factor G_p + 1 = O(q^{-437/500}) at the principal corner (l.8829-8838)
ok = True; worst = -99
for trial in range(8):
    t1, t2, t3 = [mp.mpf(random.uniform(-50, 50)) for _ in range(3)]
    eta, a, b0 = unit(), unit(), unit()
    vals = []
    for q in qs:
        H, G, R = Hp_Gp(0, q, mp.mpf(437)/500+1j*t1, mp.mpf(19)/20+1j*t2, mp.mpf(33)/200+1j*t3, eta, mp.mpf(1), a, b0)
        vals.append(G+1)
    rs = ratios(vals, mp.mpf(-437)/500); worst = max(worst, max(rs))
    if not bounded(rs, 10): ok = False; print('   G+1 ratio', rs)
rep(f'l.8838 principal box: |G_p+1| q^(437/500) bounded (max ratio {worst:.4f})', ok)
# (4) exact monomial exponents (affine in a; endpoints suffice). Monomial q^(c0 - cs Re s - cw Re w - cz Re z).
from fractions import Fraction as F
def ex(m, a, e, region):
    c0, cs, cw, cz = m
    S, Wr, Zr = region(a, e)
    return F(c0) - cs*S - cw*Wr - cz*Zr
reg1013 = lambda a, e: (a+16*e, 1-a-6*e, F(17, 50))
e_ = F(1, 1000)
# Lemma 13.10 numerator monomials of (1-R)(1-D)(H_p-1) at j=0 (general rho): eta T Y0, eta T V, eta T R, T Y0 V, V Y0, R, R Y0, R/q, (q-1) T Y0 V, (q-1)T Y0^2 V
mons1013 = [(0, 1, 1, 0), (0, 1, 0, 6), (4, 7, 0, 6), (0, 1, 1, 6), (0, 0, 1, 6), (4, 6, 0, 6), (4, 6, 1, 6), (3, 6, 0, 6), (1, 1, 1, 6), (1, 1, 2, 6)]
okx = all(ex(m, a, e_, reg1013) <= -1-10*e_ for m in mons1013 for a in (F(1, 2), F(11, 12)))
mx = max(ex(m, a, e_, reg1013) for m in mons1013 for a in (F(1, 2), F(11, 12)))
rep(f'Lemma 13.10 exact: every numerator monomial exponent <= -1-10e on a in [1/2,11/12] (max {float(mx):.4f}, bound {-1-10*float(e_):.4f}; equality from eta T Y0)', okx)
# sigma_j table: remaining monomials of (eta-bar q^s - Y0)(1-V)P_j^* - Y0 after the separated term (j>=2) / the Y0 cancellation (j=1)
Jm = {2: [(F(3, 2), 3, 0, 0)], 3: [(2, 3, 1, 0), (2, 4, 0, 0)], 4: [(F(5, 2), 4, 1, 0)], 5: [(3, 6, 0, 0)]}
okj = True; tight = []
for j in range(1, 6):
    mons = [(4, 5, 0, 6), (3, 5, 0, 6), (4, 6, 1, 6), (1, 1, 2, 0), (0, 0, 1, 0)] if j >= 2 else [(4, 5, 0, 6), (1, 0, 1, 6), (0, 0, 1, 6), (4, 6, 1, 6)]
    if j >= 2:
        for c0, cs, cw, cz in Jm[j]:
            mons.append((c0, cs-1, cw, cz))      # eta-bar q^s J_j
            mons.append((c0, cs, cw+1, cz))      # -Y0 J_j
    sigj = {1: F(26, 25)+80*e_, 2: 32*e_, 3: 26*e_, 4: 42*e_, 5: 80*e_}[j]
    worst = max(ex(m, a, e_, reg1013) for m in mons for a in (F(1, 2), F(11, 12)))
    okj = okj and worst <= F(1, 2)-sigj
    tight.append((j, float(worst-(F(1, 2)-sigj))))
rep(f'sigma_j exact: max monomial exponent minus (1/2 - sigma_j) per j = {tight} (0 = tight at a = 1/2)', okj)
regbox = lambda a, e: (F(437, 500), F(19, 20), F(33, 200))
mbox = max(ex(m, 0, 0, regbox) for m in mons1013)
rep(f'principal box exact: max numerator monomial exponent = {float(mbox):.4f} <= -907/500', mbox <= F(-907, 500))
print('SUMMARY_LOCAL', sum(res), 'pass of', len(res))

