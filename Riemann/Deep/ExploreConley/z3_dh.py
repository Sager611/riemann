# z3_dh.py — explore66 adversarial test: Davenport–Heilbronn's off-line pair 0.808517 + 85.699348 i under DH's own heat flow.
# H^DH_0(z) = Lambda_f(1/2 + i z/2)/8, Lambda_f(s) = (5/pi)^{(s+1)/2} Gamma((s+1)/2) [c1 L(s,chi) + c2 L(s,chibar)] (as explore7).
# Same exact Taylor-coefficient heat flow as z1/z2: landing time t* (double real zero), du/dt at t*, and the Conley data of a
# symmetric block D = {|w| <= r0} for t in [0, 0.2]: zero counts of H_t, H_t', exit set of grad-descent |H_t|^2, 1D signs.
import mpmath as mp, json, time
mp.mp.dps = 40
T0 = time.time()
LOG = '/Users/adriansager/playground/riemann/Riemann/Deep/ExploreConley/logs/'
kap = (mp.sqrt(10 - 2*mp.sqrt(5)) - 2)/(mp.sqrt(5) - 1)
c1 = (1 - 1j*kap)/2; c2 = mp.conj(c1)
chi = [0, 1, 1j, -1j, -1]; chib = [0, 1, -1j, 1j, -1]
def H0(z):
    s = mp.mpf(1)/2 + 1j*z/2
    return (5/mp.pi)**((s+1)/2)*mp.gamma((s+1)/2)*(c1*mp.dirichlet(s, chi) + c2*mp.dirichlet(s, chib))/8
x0 = mp.mpf(2)*mp.mpf('85.699348485378')
print('reality check H0(x0+0.1)/|.| :', mp.nstr(H0(x0 + mp.mpf('0.1'))/abs(H0(x0 + mp.mpf('0.1'))), 6))
zoff = mp.findroot(H0, mp.mpc(x0, 2*(mp.mpf('0.808517182457') - mp.mpf(1)/2)))
c = mp.re(zoff)
print('off-line zero of H0 (z units):', mp.nstr(zoff, 16), ' s =', mp.nstr(mp.mpf(1)/2 + 1j*zoff/2, 16), flush=True)
R = mp.mpf(3); M = 128; K = M//2
vals = [H0(c + R*mp.expjpi(2*mp.mpf(j)/M)) for j in range(M)]
a = [sum(vals[j]*mp.expjpi(-2*mp.mpf(j*k)/M) for j in range(M))/(M*R**k) for k in range(K)]
sc = a[0]; a = [x/sc for x in a]
print('max |Im a_k| / max|a_k|:', mp.nstr(max(abs(mp.im(x)) for x in a)/max(abs(x) for x in a), 3), ' tail |a_k|R^k:', [mp.nstr(abs(a[k])*R**k, 3) for k in (K-3, K-1)], ' time', round(time.time()-T0, 1), flush=True)
a = [mp.re(x) for x in a]
json.dump({'c': str(c), 'R': str(R), 'M': M, 'a': [str(x) for x in a]}, open(LOG + 'taylor_dh.json', 'w'))
fact = [mp.factorial(n) for n in range(K+1)]
def coeffs(t):
    t = mp.mpf(t)
    return [sum(a[m+2*j]*(-t)**j*fact[m+2*j]/(fact[j]*fact[m]) for j in range((K-1-m)//2+1)) for m in range(K)]
def ev(b, w):
    s = mp.mpc(0)
    for x in reversed(b): s = s*w + x
    return s
def der(b): return [b[m+1]*(m+1) for m in range(len(b)-1)]
b0 = coeffs(0)
rts = mp.polyroots(b0[::-1], maxsteps=400, extraprec=200)
inside = sorted([r for r in rts if abs(r) < 2.6], key=lambda r: (mp.re(r), mp.im(r)))
print('zeros of H0^DH with |w| < 2.6 (w = z - c):', [mp.nstr(r, 8) for r in inside])
crs = mp.polyroots(der(b0)[::-1], maxsteps=400, extraprec=200)
print("zeros of H0^DH' with |w| < 2.6:", [mp.nstr(r, 8) for r in sorted([r for r in crs if abs(r) < 2.6], key=lambda r: mp.re(r))], flush=True)
def crit(t):
    b = coeffs(t); bd = der(b)
    ws = mp.findroot(lambda w: ev(bd, w), mp.mpf(0))
    return mp.re(ev(b, ws)), ws
ts = mp.findroot(lambda t: crit(t)[0], mp.mpf('0.176'), tol=mp.mpf(10)**-30)
vs, wss = crit(ts)
print('EXACT landing (double real zero) t* =', mp.nstr(ts, 15), ' at w* =', mp.nstr(wss, 12), ' de Bruijn bound y0^2/2 =', mp.nstr(mp.im(zoff)**2/2, 10), flush=True)
def pairu(t):
    b = coeffs(t); ws_ = crit(t)[1]
    if t > ts: d = mp.sqrt(2*(t-ts)); g1, g2 = ws_-d, ws_+d
    else: d = mp.sqrt(2*(ts-t)); g1, g2 = ws_-1j*d, ws_+1j*d
    z1 = mp.findroot(lambda w: ev(b, w), g1); z2 = mp.findroot(lambda w: ev(b, w), g2)
    return ((z2-z1)/2)**2
for dl in ['1e-6', '1e-4']:
    dl = mp.mpf(dl); print('du/dt at t* (delta=%s): %s' % (mp.nstr(dl, 2), mp.nstr((pairu(ts+dl)-pairu(ts-dl))/(2*dl), 12)))
json.dump({'tstar': str(ts), 'zoff': str(zoff), 'c': str(c)}, open(LOG + 'dh_tstar.json', 'w'))
print('time', round(time.time()-T0, 1))
