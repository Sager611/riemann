import sys, time; sys.path.insert(0, __file__.rsplit('/', 1)[0])
# adversarial (c): H_s(z) = 4 int_0^inf Phi e^{s u^2} cos(zu) du (H_0 = Xi; Polymath t = 4 s, z_P = 2 z) at s < 0: the gamma_9, gamma_10 pair leaves the line
from kern import ev, TZ, xi
import mpmath as mp
mp.mp.dps = 25
t0 = time.time()
NODES = mp.linspace(0, 4, 41)
def H(s, z): return 4*mp.quad(lambda u: ev(TZ, u)*mp.e**(s*u*u)*mp.cos(z*u), NODES, method='gauss-legendre')
print('check H_0(48) =', mp.nstr(H(0, mp.mpf(48)), 10), ' Xi(48) =', mp.nstr(mp.re(xi(mp.mpf(0.5) + 48j)), 10), flush=True)
for s, guess in [(mp.mpf('-0.4'), None), (mp.mpf('-0.5'), mp.mpc(49.6, 0.25)), (mp.mpf('-0.6'), mp.mpc(49.8, 0.5))]:
    xs = [mp.mpf('48.6') + mp.mpf(i)/10 for i in range(0, 21)]
    vals = [H(s, x) for x in xs]
    ch = [mp.nstr((xs[i] + xs[i+1])/2, 5) for i in range(20) if vals[i]*vals[i+1] < 0]
    print('s=%s (t=%s): sign changes of H_s on [48.6, 50.6]:' % (s, 4*s), ch, flush=True)
    if guess is not None:
        z = mp.findroot(lambda z: H(s, z), guess, tol=1e-18)
        print('    non-real zero z = %s (z_P = %s), |H_s(z)| = %s, |H_s(Re z)| = %s' % (mp.nstr(z, 12), mp.nstr(2*z, 10), mp.nstr(abs(H(s, z)), 3), mp.nstr(abs(H(s, mp.re(z))), 3)), flush=True)
print('time %.1f' % (time.time() - t0))
