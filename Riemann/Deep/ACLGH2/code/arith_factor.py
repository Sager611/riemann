# arith_factor.py: F_S = prod_{p in S} E_theta[exp(delta_p(theta))], delta_p = 2 log p Re(a e^{i theta}/(1 - a e^{i theta})), a = p^{-1/2}:
# the semiclassical count of coupled directions of the semilocal pair at scale t is N_sc = (1/2pi) int D_t(tau)_+ dtau = 2 t^2 F_S
# under equidistribution of (tau log p mod 2pi)_p (Kronecker-Weyl); D_t = 2 log t - 2 vartheta' + sum_p delta_p (total delay).
# Prints log F_S for S = {inf} u {p <= x}, compared with sum (log p)^2/p and (1/2) log^2 x.
import math, json, numpy as np
X = 10**6
sieve = np.ones(X + 1, dtype=bool); sieve[:2] = False
for i in range(2, int(X**0.5) + 1):
    if sieve[i]: sieve[i*i::i] = False
ps = np.nonzero(sieve)[0]
th = np.linspace(0, 2*math.pi, 2048, endpoint=False); c = np.cos(th)
logF = 0.0; s2 = 0.0; out = []; marks = {3, 5, 7, 11, 13, 23, 47, 97, 199, 499, 997, 1999, 4999, 9973, 19997, 49999, 99991, 199999, 499979, 999983}
for p in ps:
    a = p**-0.5; L = math.log(p)
    e = float(np.mean(np.exp(2*L*(a*c - a*a)/(1 - 2*a*c + a*a))))   # trapezoid on the circle: spectrally accurate
    logF += math.log(e); s2 += L*L/p
    if int(p) in marks or p == ps[-1]:
        out.append({'x': int(p), 'S': int(np.searchsorted(ps, p)) + 2, 'logF': round(logF, 3), 'F': float('%.4g' % math.exp(logF)),
                    'sum_log2p_over_p': round(s2, 3), 'half_log2x': round(0.5*math.log(p)**2, 3), 'logF_over_logx': round(logF/math.log(p), 3)})
for o in out: print(json.dumps(o))
p = 9973; a = p**-0.5; L = math.log(p)
print('check p=9973: quadrature', float(np.mean(np.exp(2*L*(a*c - a*a)/(1 - 2*a*c + a*a)))), 'second-order series 1 + (log p)^2/(p-1):', 1 + L*L/(p - 1))
