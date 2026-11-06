# duvall.py (explore23, T1): Duvall/Abel "zeta star".
# (a) 1-D acoustic-cutoff profile from the smooth count (Wu-Sprung in acoustic dress), closed form + forward check.
# (b) Duvall law for the family of L-functions: pi(n+alpha)/omega = F(q*omega), F(w) = (1/2)log(w/2pi e); Abel inversion -> c(r) = 2 r ln(R/r).
# (c) Duvall diagram from zeros of real primitive Dirichlet L-functions (q = 1,3,4,5,7,8,8',11,12,13) and DH (full vs on-line only).
import sys, json, numpy as np, mpmath as mp
sys.path.insert(0, __file__.rsplit('/', 1)[0])
from dh import hurwitz
PI = np.pi
out = {}
# (a) tau(omega_c) = (1/pi)[(pi/2 - phi0) log(omega_c/2pi) + int_{phi0}^{pi/2} log sin], sin phi0 = 2pi/omega_c ; E = omega^2, V = omega_c^2
def tau_of_wc(wc):
    phi0 = np.arcsin(2*PI/wc)
    ph = np.linspace(phi0, PI/2, 20001)
    return ((PI/2 - phi0)*np.log(wc/(2*PI)) + np.trapezoid(np.log(np.sin(ph)), ph))/PI
wcs = np.unique(np.concatenate([np.linspace(2*PI*1.0000001, 30, 400), np.geomspace(30, 1e6, 600)]))
taus = np.array([tau_of_wc(x) for x in wcs])
out['a_monotone_tau'] = bool(np.all(np.diff(taus) > 0))
out['a_asymptote_wc_over_4pi_e2tau_at_1e6'] = float(wcs[-1]/(4*PI*np.exp(2*taus[-1])))
# forward WKB: m(E) = (1/pi) int_0^{tau_t} sqrt(E - wc(tau)^2) dtau, should equal Nbar(omega) - Nbar(2pi), Nbar = (w/2pi) log(w/2pi e)
def Nbar(w): return w/(2*PI)*np.log(w/(2*PI*np.e))
chk = []
for om in [10.0, 50.0, 200.0, 1000.0, 5000.0]:
    sel = wcs <= om
    tt = taus[sel]; ww = wcs[sel]
    # integrate in tau: sqrt(om^2 - wc^2) d tau, plus the last segment to the turning point
    val = np.trapezoid(np.sqrt(np.maximum(om**2 - ww**2, 0)), tt)
    tau_t = tau_of_wc(om)
    # tail from last grid point to turning point: sqrt(om^2-wc^2) ~ sqrt(2 om (om - wc)), integrate analytically via linearization
    dw = om - ww[-1]; dtau = tau_t - tt[-1]
    val += (2/3)*np.sqrt(2*om*dw)*dtau
    chk.append([om, float(val/PI), float(Nbar(om) - Nbar(2*PI))])
out['a_forward_check_[omega, WKB m, Nbar(omega)-Nbar(2pi)]'] = chk
# (b) Duvall forward check for c(r) = 2 r ln(1/r) (R = 1), a = c/r = 2 zeta, zeta = ln(1/r), surface cut a_s
a_s = 0.01
def Fduv(wv):
    z = np.linspace(a_s/2, wv/2, 200001)
    return np.trapezoid(np.sqrt(np.maximum(1/(2*z)**2 - 1/wv**2, 0)), z)
out['b_duvall_[w, F numeric, 0.5 ln(2w/(e a_s))]'] = [[wv, float(Fduv(wv)), float(0.5*np.log(2*wv/(np.e*a_s)))] for wv in [1.0, 10.0, 100.0, 1000.0]]
# (c) real primitive characters
def leg(r, p):
    v = pow(r % p, (p - 1)//2, p)
    return 0 if r % p == 0 else (1 if v == 1 else -1)
chars = {'1': (1, {1: 1}), '-3': (3, {r: leg(r, 3) for r in range(1, 3)}), '-4': (4, {1: 1, 3: -1}), '5': (5, {r: leg(r, 5) for r in range(1, 5)}),
         '-7': (7, {r: leg(r, 7) for r in range(1, 7)}), '8': (8, {1: 1, 3: -1, 5: -1, 7: 1}), '-8': (8, {1: 1, 3: 1, 5: -1, 7: -1}),
         '-11': (11, {r: leg(r, 11) for r in range(1, 11)}), '12': (12, {1: 1, 5: -1, 7: -1, 11: 1}), '13': (13, {r: leg(r, 13) for r in range(1, 13)})}
TM = 300.0
ts = np.arange(2.0, TM, 0.02)
def varth(t, q, a):
    return np.array([float(mp.im(mp.loggamma((mp.mpf(0.5) + a + 1j*mp.mpf(x))/2))) for x in t]) + 0.5*t*np.log(q/PI)
res = {}
allpts = []
for name, (q, ch) in chars.items():
    a = 0 if ch.get(q - 1, 1) == 1 else 1
    s = 0.5 + 1j*ts
    L = q**(-s)*sum(c*hurwitz(s, r/q, int(TM/2) + 60) for r, c in ch.items() if c != 0)
    if q == 1:
        Z = (np.exp(1j*varth(ts, 1, 0))*L).real   # theta for zeta (no s(s-1) factor needed for sign changes)
    else:
        Z = (np.exp(1j*varth(ts, q, a))*L).real
    idx = np.nonzero(np.sign(Z[:-1])*np.sign(Z[1:]) < 0)[0]
    g = ts[idx] - Z[idx]*(ts[idx + 1] - ts[idx])/(Z[idx + 1] - Z[idx])
    n = np.arange(1, len(g) + 1)
    vt = varth(g, q, a)
    c = np.median(n - 0.5 - vt/PI)
    y = PI*(n - 0.5 - c)/g                    # Duvall ordinate pi(n + alpha)/omega
    Fw = 0.5*np.log(q*g/(2*PI*np.e))          # F(w), w = q*omega
    r_ = y - Fw
    res[name] = {'q': q, 'parity': a, 'zeros': int(len(g)), 'alpha(c)': round(float(c), 4), 'max|S| est': round(float(np.max(np.abs(n - 0.5 - c - vt/PI))), 3),
                 'resid_rms_g>50': float(np.sqrt(np.mean(r_[g > 50]**2))), 'resid_rms_g>150': float(np.sqrt(np.mean(r_[g > 150]**2)))}
out['c_duvall_family'] = res
# DH (conductor 5, odd): full list vs on-line only, on the same curve
D = json.load(open('dh_zeros.json'))
on = np.array(D['online']); off = np.array(D['offline'])
for lab, g in [('dh_full', np.sort(np.concatenate([on, off[:, 1], off[:, 1]]))), ('dh_online_only', on)]:
    n = np.arange(1, len(g) + 1); vt = varth(g, 5, 1)
    c = np.median((n - 0.5 - vt/PI)[g < 80])
    dev = n - 0.5 - c - vt/PI
    out['c_' + lab] = {'zeros': int(len(g)), 'alpha(c)': round(float(c), 4), 'dev_at_[100,200,300,419]': [round(float(dev[np.searchsorted(g, x) - 1]), 3) for x in [100, 200, 300, 419]],
                       'resid_rms_g>150': float(np.sqrt(np.mean((PI*(n - 0.5 - c)/g - 0.5*np.log(5*g/(2*PI*np.e)))[g > 150]**2)))}
json.dump(out, open('duvall.json', 'w'), indent=1)
print(json.dumps(out)[:6000])
