#!/usr/bin/env python3
# test_b2k.py -- numerical test of Lemma 11_k (per-zero bounds for the k-box ramp, one-sided c = 0 and centred c = lam/2).
#   uv run --no-project --with mpmath python3 test_b2k.py
# F_phi(rho) = E[F(L' + V)], F(y) = int_y^inf e^{s u} H(u) du = (1+s) E1(-s y) + e^{s y}/y, s = rho - 1, L' = L - c,
# V = (lam/k) * (sum of k uniforms on [0,1]); compared with
# (B2_k) |F_phi| <= e^{(beta-1)L'} [H(L')(2k/lam)^k/|s|^(k+1) + |H'(L')| lam/(2|s|) + 2|H'(L')|/|s|^2]
# (B1')  |F_phi - F(L)| <= e^{(beta-1)L'} H(L') m_c,  m_c = lam/2 (c = 0) or (lam/k) MAD_k (c = lam/2).
# Also checks MAD_k against direct numerical integration of E|S - k/2|.
import mpmath as mp, time
mp.mp.dps = 30
def H(u): return 1/u + 1/u**2
def Hp(u): return -1/u**2 - 2/u**3
def F(s, y): return (1 + s)*mp.e1(-s*y) + mp.exp(s*y)/y
def ih(k, z):                                   # Irwin-Hall density on [0, k]
    if z <= 0 or z >= k: return mp.mpf(0)
    return sum((-1)**j*mp.binomial(k, j)*(z - j)**(k - 1) for j in range(int(mp.floor(z)) + 1))/mp.factorial(k - 1)
def Fphi(s, Lp, lam, k):
    g = lambda z: ih(k, z)*F(s, Lp + lam*z/k)  # V = lam z/k, density in z
    return mp.quad(g, [mp.mpf(j) for j in range(k + 1)], maxdegree=10)
def MAD(k):
    return 2*sum((-1)**j*mp.binomial(k, j)*(mp.mpf(k)/2 - j)**(k + 1) for j in range(k//2 + 1))/mp.factorial(k + 1)
for k in (1, 2, 4, 8):
    num = mp.quad(lambda z: abs(z - mp.mpf(k)/2)*ih(k, z), [mp.mpf(j) for j in range(k + 1)])
    print('MAD', k, mp.nstr(MAD(k), 12), 'numerical', mp.nstr(num, 12))
t0 = time.time(); worst = {}; n = 0
# direct check of the closed form for F
s = mp.mpc(-0.5, 40); y = mp.mpf(23)
direct = mp.quad(lambda u: mp.exp(s*u)*H(u), mp.linspace(y, y + 40, 400) + [mp.inf])
print('closed form check:', mp.nstr(abs(direct - F(s, y))/abs(F(s, y)), 5))
for beta in (0.5, 0.875, 0.97):
    for gam in (30, 300, 3000):
        for L in (23, 70):
            for k in (1, 2, 4, 8):
                for r in (0.2, 0.6, 2.0):          # r = 2k/(lam gamma)
                    lam = mp.mpf(2*k)/(r*gam)
                    for cen in (0, 1):
                        c = lam/2 if cen else mp.mpf(0)
                        s = mp.mpc(beta - 1, gam); Lp = L - c
                        fp = Fphi(s, Lp, lam, k)
                        a = abs(s); ex = mp.exp((beta - 1)*Lp)
                        b2 = ex*(H(Lp)*(2*k/lam)**k/a**(k + 1) + abs(Hp(Lp))*lam/(2*a) + 2*abs(Hp(Lp))/a**2)
                        mc = lam/2 if cen == 0 else lam*MAD(k)/k
                        b1 = ex*H(Lp)*mc
                        q2 = abs(fp)/b2; q1 = abs(fp - F(s, mp.mpf(L)))/b1
                        n += 1
                        for key, q in (('B2k', q2), ('B1', q1)):
                            if key not in worst or q > worst[key][0]: worst[key] = (q, beta, gam, L, k, r, cen)
print(f'{n} cases, {time.time()-t0:.0f}s')
for key, v in worst.items():
    print(key, 'worst ratio true/bound =', mp.nstr(v[0], 6), 'at beta,gamma,L,k,r,centred =', v[1:])
# large k (the certificates use k up to about 60): higher precision for the Irwin-Hall density
mp.mp.dps = 80
for beta in (0.5, 0.97):
    for gam in (300, 3000):
        for k in (16, 32):
            for r in (0.3, 0.6):
                lam = mp.mpf(2*k)/(r*gam); c = lam/2; L = 70
                s = mp.mpc(beta - 1, gam); Lp = L - c
                fp = Fphi(s, Lp, lam, k)
                a = abs(s); ex = mp.exp((beta - 1)*Lp)
                b2 = ex*(H(Lp)*(2*k/lam)**k/a**(k + 1) + abs(Hp(Lp))*lam/(2*a) + 2*abs(Hp(Lp))/a**2)
                b1 = ex*H(Lp)*lam*MAD(k)/k
                q2 = abs(fp)/b2; q1 = abs(fp - F(s, mp.mpf(L)))/b1; n += 1
                for key, q in (('B2k', q2), ('B1', q1)):
                    if q > worst[key][0]: worst[key] = (q, beta, gam, L, k, r, 1)
                print(f'  large k: beta={beta} gamma={gam} k={k} r={r}: B2k ratio {mp.nstr(q2, 4)}, B1 ratio {mp.nstr(q1, 4)}', flush=True)
print(f'{n} cases in total, {time.time()-t0:.0f}s')
print('PASS' if all(v[0] < 1 for v in worst.values()) else 'FAIL')
