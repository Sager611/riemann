# Independent float check of the c = 0.001 hinge numbers (DERIVATION 3.7). Run: uv run --no-project --with mpmath python3 sanity_mp.py
import mpmath as mp
mp.mp.dps = 30
S2 = mp.sqrt(2); cot = mp.cot(1/S2)
def K(y): return (mp.cosh(mp.pi*y) + S2*mp.pi*y*cot*mp.sinh(mp.pi*y))/(1+2*mp.pi**2*y**2)
def Kint(y):  # direct quadrature of 2 int_0^{1/2} f0 cosh(2 pi x y)
    f0 = lambda x: mp.cos(S2*x)/(S2*mp.sin(1/S2))
    return 2*mp.quad(lambda x: f0(x)*mp.cosh(2*mp.pi*x*y), [0, 0.5])
c = mp.mpf('0.001')
def K0r(u):
    f0 = lambda x: mp.cos(S2*x)/(S2*mp.sin(1/S2))
    return 2*mp.quad(lambda x: f0(x)*mp.cos(2*mp.pi*u*x), [0, 0.5])
lam = 2*(1-K0r(2*c)); thb = mp.pi**2*c**2/2; s = 3*lam
print('lam', lam, ' K(iy) closed vs quad at y=2.3:', K(2.3), Kint(2.3))
H = lambda y: 1 - (1-thb)**2*4*K(y)**2/(1+K(2*y))
# root of H - s
r = mp.findroot(lambda y: H(y)-s, 1.2)
print('rho_max (float root of H - 3lam):', r)
rho = mp.mpf('1.2036')
best = max((H(y)-s)/(y-rho) for y in [rho + mp.mpf(k)/2000 for k in range(1, 40000)])
print('sup h/(y-rho) on grid:', best, ' (cert: <= 0.670124)')
R = 2*mp.pi*rho; th = mp.mpf(4)/7
B = mp.log(1 + (1-mp.exp(-2*R))/(mp.exp(2*th*R)-1))/(4*mp.pi)
print('B(rho):', B, ' (cert: <= 1.404e-05)')
p = (2 - lam - (mp.mpf(1)/2 + cot/S2) - best*B)/(1-lam)
print('p(0.001) float:', p, ' (cert >= 0.672487)')
