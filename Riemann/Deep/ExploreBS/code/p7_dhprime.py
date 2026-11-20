# p7_dhprime.py (explore57): the zero of DH' (derivative of the Davenport-Heilbronn function) behind phi'_DH(85.70) = -9.50:
# Newton on f'(s) from starts left and right of the line near t = 85.7; reports the zero, its distance delta' = beta' - 1/2, and
# the Poisson term (beta'-1/2)/((beta'-1/2)^2 + (t - gamma')^2) at t = 85.70 (negative if the zero is left of the line).
import sys, mpmath as mp
sys.path.insert(0, __file__.rsplit('/', 1)[0])
mp.mp.dps = 20
kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
AR = [(mp.mpf(1)/5, 1), (mp.mpf(2)/5, kap), (mp.mpf(3)/5, -kap), (mp.mpf(4)/5, -1)]
def fd(s, k):
    L5 = mp.log(5); tot = 0
    for j in range(k + 1):
        tot += mp.binomial(k, j)*(-L5)**(k - j)*sum(a*mp.zeta(s, r, derivative=j) for r, a in AR)
    return mp.power(5, -s)*tot
found = []
for s0 in [mp.mpc(x, y) for x in (0.2, 0.35, 0.45, 0.6, 0.9) for y in (84.8, 85.4, 85.7, 86.0, 86.6)]:
    try:
        z = mp.findroot(lambda s: fd(s, 1), s0, tol=1e-25, maxsteps=60)
    except Exception:
        continue
    if 80 < z.imag < 92 and all(abs(z - w) > 1e-8 for w in found): found.append(z)
found.sort(key=lambda z: float(z.imag))
t = mp.mpf('85.70')
for z in found:
    d = z.real - 0.5
    print(f"DH' zero {mp.nstr(z, 10)}: beta'-1/2 = {mp.nstr(d, 6)}, Poisson term at t = 85.70: {mp.nstr(d/(d**2 + (t - z.imag)**2), 6)}")
for z0 in [mp.mpc('0.8085', '85.6993')]:
    z = mp.findroot(lambda s: fd(s, 0), z0, tol=1e-25)
    print('DH zero near 0.8085+85.6993i:', mp.nstr(z, 12), ' its reflection 1 - conj:', mp.nstr(1 - mp.conj(z), 12))
