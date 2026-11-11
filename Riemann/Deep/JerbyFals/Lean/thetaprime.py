# Evidence (not proof) for the named hypothesis ThetaDerivPos: theta'(t) = Re psi(1/4+it/2)/2 - log(pi)/2 > 0 for t >= 7.
import mpmath as mp
mp.mp.dps = 30
def tp(t): return mp.re(mp.digamma(mp.mpf(1)/4 + 1j*mp.mpf(t)/2))/2 - mp.log(mp.pi)/2
ts = mp.findroot(tp, 6.3)
print('theta\' zero t* =', mp.nstr(ts, 15), ' theta(t*) =', mp.nstr(mp.siegeltheta(ts), 12))
grid = [7 + k*0.01 for k in range(0, 1001)] + [17 + k for k in range(0, 2000)] + [10**e for e in range(4, 13)]
mn = min((tp(t), t) for t in grid)
print('min theta\' on grid (t>=7):', mp.nstr(mn[0], 10), 'at t =', mn[1])
for t in [7, 568009639.698978187581685, 906089969.752693821200281, 1155454823.107585085442780]:
    print('t =', t, " theta'(t) =", mp.nstr(tp(t), 15), ' vs siegeltheta deriv', mp.nstr(mp.siegeltheta(t, derivative=1), 15))
# monotonicity of Re psi(1/4+iy) in y: derivative = -Im psi'(1/4+iy) > 0 on a grid
md = min(-mp.im(mp.psi(1, mp.mpf(1)/4 + 1j*mp.mpf(y))) for y in [0.01*k for k in range(1, 2001)])
print('min d/dy Re psi(1/4+iy) on y in (0,20]:', mp.nstr(md, 10))
# check Lean theta against mpmath siegeltheta: integral of theta' from 0 to t
for t in [7, 50]:
    print('t =', t, ' int_0^t theta\' =', mp.nstr(mp.quad(tp, [0, t]), 20), ' siegeltheta =', mp.nstr(mp.siegeltheta(t), 20))
