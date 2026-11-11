# theta_vs_aclg.py: Chekmachev's empirical scale Theta_c (Zenodo 23848597; "found empirically", eps_c >= 0.97 Theta_c proved for
# 3 <= c <= 31, 45, 47) against ACLG's min-max upper bounds eps^N(1) (2611.13129 Table logbound l. 684-697, exact rungs D_1 l. 640-656).
import mpmath as mp
mp.mp.dps = 50
def Phi(u): return 2*mp.fsum((2*mp.pi**2*n**4*mp.e**(9*u/2) - 3*mp.pi*n**2*mp.e**(5*u/2))*mp.e**(-mp.pi*n**2*mp.e**(2*u)) for n in range(1, 30))
def Theta(c):
    c = mp.mpf(c); a = mp.log(c)/2; return (2*mp.pi*Phi(a)/(2*mp.pi*c - mp.mpf(9)/2))**4
# (c, D_1 = -log delta_4 exact, ACLG ratio eps^N/delta_4 at largest N)
rows = [(5, 41.185, 7.06), (9, 88.680, 8.11), (22, 247.933, 9.17), (24, 272.669, 9.31), (28, 322.233, 9.76), (32, 371.892, 11.16),
        (36, 421.623, 10.46), (40, 471.411, 11.12), (50, 596.064, 9.66), (100, 1221.251, 12.27), (150, 1847.741, 13.83),
        (200, 2474.762, 18.07), (250, 3102.075, 16.08), (300, 3729.573, 19.02)]
print('| lambda^2 | Theta_c/delta_4 | 0.97 Theta_c/delta_4 | ACLG eps^N(1)/delta_4 (upper bound) | 0.97 Theta_c <= eps^N ? |')
print('|---|---|---|---|---|')
for c, D1, r in rows:
    t = Theta(c)*mp.e**mp.mpf(D1)
    print(f'| {c} | {float(t):.3f} | {0.97*float(t):.3f} | {r} | {"yes" if 0.97*float(t) <= r else "NO"} |')
