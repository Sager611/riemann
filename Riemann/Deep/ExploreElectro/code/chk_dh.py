import mpmath as mp, sys
sys.path.insert(0, '.')
import ecommon as E
mp.mp.dps = 30
kap = E.kappa()
A = {1: mp.mpf(1), 2: kap, 3: -kap, 4: mp.mpf(-1), 0: mp.mpf(0)}
def theta(x):
    s = mp.mpf(0); n = 1
    while True:
        ex = mp.exp(-mp.pi*n*n*x/5); s += A[n%5]*n*ex
        if ex < mp.mpf(10)**-40 and n > 3: break
        n += 1
    return s
Psi = lambda v: 2*mp.exp(1.5*v)*theta(mp.exp(2*v))
for v in (0.1, 0.3, 0.7):
    print('Psi(v), Psi(-v):', mp.nstr(Psi(v), 20), mp.nstr(Psi(-v), 20))
I0 = 2*mp.quad(Psi, [0, 1, 2, 4])
print('2 int Psi =', mp.nstr(I0, 20), ' dh(1/2) =', mp.nstr(E.dh(mp.mpf(0.5)), 20))
for z in (3, 7.5):
    Iz = 2*mp.quad(lambda v: Psi(v)*mp.cos(z*v), [0, 0.5, 1, 2, 4])
    print('z', z, mp.nstr(Iz, 20), mp.nstr(E.dh(mp.mpc(0.5, z)), 20))
print('min Psi on grid [0,3]:', mp.nstr(min(Psi(mp.mpf(k)/100) for k in range(0, 301)), 10))
