# p1_prufer.py (explore57): Prufer angle of Hardy's Z w.r.t. the theta-oscillator (y'/theta')' + theta' y = 0.
# Z = r sin(phi), Z'/theta' = r cos(phi);  phi' = (theta'(Z'^2 - Z Z'') + theta'' Z Z')/(Z'^2 + theta'^2 Z^2).
# Checks: (1) crossings of k*pi by phi = sign changes of Z (count), (2) phi + theta + arg zeta'(1/2+it) = const mod 2pi,
# (3) floor: min phi' vs (1/2) log 2 (RH + Speiser prediction), (4) same for Davenport-Heilbronn (DH) on [18, TDH].
# usage: python3 p1_prufer.py zeta|dh|ident a b step dps
import sys, time, math
import mpmath as mp
WHICH = sys.argv[1]; A_ = float(sys.argv[2]); B_ = float(sys.argv[3]); h = float(sys.argv[4]); mp.mp.dps = int(sys.argv[5])
LOG2H = 0.5*math.log(2)
t_start = time.time()

def zeta_data(t):
    Z = mp.siegelz(t); Z1 = mp.siegelz(t, derivative=1); Z2 = mp.siegelz(t, derivative=2)
    th1 = mp.siegeltheta(t, derivative=1); th2 = mp.siegeltheta(t, derivative=2)
    return float(Z), float(Z1), float(Z2), float(th1), float(th2)

kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
AR = [(mp.mpf(1)/5, 1), (mp.mpf(2)/5, kap), (mp.mpf(3)/5, -kap), (mp.mpf(4)/5, -1)]
def fdh(s, k):
    # k-th s-derivative of f(s) = 5^{-s} sum_r a_r zeta(s, r/5)
    L5 = mp.log(5); tot = 0
    for j in range(k + 1):
        hj = sum(a*mp.zeta(s, r, derivative=j) for r, a in AR)
        tot += mp.binomial(k, j)*(-L5)**(k - j)*hj
    return mp.power(5, -s)*tot
def dh_data(t):
    s = mp.mpc(0.5, t)
    vt = mp.im(mp.loggamma(mp.mpf(0.75) + 0.5j*t)) + 0.5*t*mp.log(5/mp.pi)
    v1 = 0.5*mp.re(mp.digamma(mp.mpf(0.75) + 0.5j*t)) + 0.5*mp.log(5/mp.pi)
    v2 = -0.25*mp.im(mp.psi(1, mp.mpf(0.75) + 0.5j*t))
    F = fdh(s, 0); F1 = 1j*fdh(s, 1); F2 = -fdh(s, 2)
    e = mp.expj(vt)
    g = e*F; g1 = e*(1j*v1*F + F1); g2 = e*(1j*v2*F - v1**2*F + 2j*v1*F1 + F2)
    return float(mp.re(g)), float(mp.re(g1)), float(mp.re(g2)), float(v1), float(v2), float(abs(mp.im(g)) + abs(mp.im(g1)))

def scan(name, data, a, b, h):
    n = int((b - a)/h) + 1
    phi_prev = None; phi0 = None; kcross = 0; signch = 0; zprev = None
    minp = (1e9, None); negs = 0; worst_im = 0.0; acc = 0.0; pprev = None
    for i in range(n):
        t = a + i*h
        d = data(t)
        Z, Z1, Z2, w, w1 = d[:5]
        if len(d) > 5: worst_im = max(worst_im, d[5])
        ang = math.atan2(Z, Z1/w)
        if phi_prev is None: phi = ang; phi0 = ang
        else:
            dlt = ang - (phi_prev % (2*math.pi) if False else math.atan2(math.sin(phi_prev), math.cos(phi_prev)))
            dlt = (dlt + math.pi) % (2*math.pi) - math.pi
            phi = phi_prev + dlt
            kcross += math.floor(phi/math.pi) - math.floor(phi_prev/math.pi)
        p = (w*(Z1*Z1 - Z*Z2) + w1*Z*Z1)/(Z1*Z1 + w*w*Z*Z)
        if pprev is not None: acc += 0.5*h*(p + pprev)
        pprev = p
        if p < minp[0]: minp = (p, t)
        if p < 0: negs += 1
        if zprev is not None and Z*zprev < 0: signch += 1
        zprev = Z; phi_prev = phi
    print(f'{name}: [{a},{b}] step {h}: sign changes of Z = {signch}, net k*pi crossings of phi = {kcross}, '
          f'(phi_end - phi_0)/pi = {(phi - phi0)/math.pi:.4f}, int phi\'/pi = {acc/math.pi:.4f}', flush=True)
    print(f'{name}: min phi\' = {minp[0]:.6f} at t = {minp[1]:.3f}; grid points with phi\' < 0: {negs}; (1/2)log2 = {LOG2H:.6f}'
          + (f'; max |Im| check = {worst_im:.2e}' if worst_im else ''), flush=True)
    return minp

# identity (2): phi + theta + arg zeta'(1/2+it) = const mod 2pi at sample points
if WHICH == 'ident':
    vals = []
    for t in [20.5, 37.3, 61.1, 99.9, 143.2, 201.7, 280.4, 1000.3, 5000.7]:
        Z, Z1, Z2, w, w1 = zeta_data(t)
        ang = math.atan2(Z, Z1/w)
        zp = mp.zeta(mp.mpc(0.5, t), derivative=1)
        c = (ang + float(mp.siegeltheta(t)) + float(mp.arg(zp))) % (2*math.pi)
        vals.append(round(c, 12))
    print('identity phi + theta + arg zeta\' mod 2pi:', vals, '(3pi/2 =', round(1.5*math.pi, 12), ')', flush=True)
elif WHICH == 'zeta':
    scan('zeta', zeta_data, A_, B_, h)
    print('zeta N(b) - N(a) via mp.nzeros:', mp.nzeros(B_) - mp.nzeros(A_), ' elapsed', round(time.time() - t_start, 1), flush=True)
else:
    scan('DH', dh_data, A_, B_, h)
    print('elapsed', round(time.time() - t_start, 1), flush=True)
