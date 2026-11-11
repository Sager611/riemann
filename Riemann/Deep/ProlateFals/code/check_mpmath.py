# independent check: mpmath (arbitrary precision Taylor ODE solver) vs prolate.c raw Pruefer phase at finite X
import mpmath as mp, sys, math, subprocess
mp.mp.dps = 25
lam = mp.sqrt(2); om = 2*mp.pi*lam
def raw_phase(E, X):
    E = mp.mpf(E); X = mp.mpf(X)
    # Frobenius start
    h1 = mp.mpf('1e-3')
    V0 = om**2*lam**2 + E**2
    U = [mp.mpf(1)]; Um1 = 0; Um2 = 0
    u = mp.mpf(1); up = mp.mpf(0)
    Un, Unm1, Unm2 = mp.mpf(1), mp.mpf(0), mp.mpf(0)
    for n in range(0, 200):
        Unp1 = -((n*(n+1) + V0)*Un + 2*om**2*lam*Unm1 + om**2*Unm2)/(2*lam*(n+1)**2)
        u += Unp1*h1**(n+1); up += (n+1)*Unp1*h1**n
        Unm2, Unm1, Un = Unm1, Un, Unp1
    x0 = lam + h1
    # ODE: p u'' + 2x u' + V u = 0
    f = mp.odefun(lambda x, y: [y[1], -(2*x*y[1] + (om**2*x**2 + E**2)*y[0])/(x**2 - lam**2)], x0, [u, up])
    # unwrap phase on a grid
    k = lambda x: mp.sqrt(E**2 + om**2*x**2)
    al = mp.atan2(k(x0)*u, mp.sqrt(x0**2-lam**2)*up); W = 0; phi = al
    N = int(mp.ceil((X - x0)*(om + E)*2)) + 200
    for i in range(1, N+1):
        x = x0 + (X - x0)*i/N
        y = f(x)
        a2 = mp.atan2(k(x)*y[0], mp.sqrt(x**2-lam**2)*y[1])
        d = a2 - al
        while d > mp.pi: d -= 2*mp.pi
        while d <= -mp.pi: d += 2*mp.pi
        phi += d; al = a2
    return phi - om*X
for E, X in ((3, 12), (10, 12)):
    r = raw_phase(E, X)
    out = subprocess.check_output(['/Users/adriansager/playground/riemann/.scratch/prolatefals/prolate', '1.4142135623730951', 'theta', str(E), repr(float(om*X))]).decode()
    rawc = float(out.split('raw=')[1].split()[0])
    print(f'E={E} X={X}: mpmath raw = {mp.nstr(r, 15)}  prolate.c raw = {rawc:.15f}  diff = {float(r) - rawc:.2e}')
