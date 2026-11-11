# semicl.py: semiclassical count of coupled (near-self-dual) directions of the semilocal pair at scale t:
# N_sc(T) = (1/2pi) int_{|tau|<=T} D(tau)_+ dtau, D = psi' = total delay of the symbol phi_t = t^{2 i tau} u_inf prod_p u_p:
# D = 2 log t - 2 vartheta'(tau) + sum_p delta_p(tau), delta_p = 2 log p (a cos(tau L) - a^2)/(1 - 2 a cos(tau L) + a^2), a = p^-1/2, L = log p.
# Also the Cramer/large-deviation estimate N_ld = prod_p E_theta[exp(delta_p(theta))] (theta uniform; Kronecker-Weyl).
# Usage: semicl.py primes Tmax [t]
import sys, math, json, numpy as np
from scipy.special import digamma
primes = [] if sys.argv[1] == '-' else [int(q) for q in sys.argv[1].split(',')]
Tmax = float(sys.argv[2]); tsc = float(sys.argv[3]) if len(sys.argv) > 3 else 1.0
def D(tau):
    out = 2*math.log(tsc) - (np.real(digamma(0.25 + 0.5j*tau)) - math.log(math.pi))     # -2 vartheta'(tau)
    for p in primes:
        a = p**-0.5; L = math.log(p); c = np.cos(tau*L)
        out = out + 2*L*(a*c - a*a)/(1 - 2*a*c + a*a)
    return out
h = 0.01; edges = [10, 30, 100, 300, 588.75, 600, 885, 1000, 1200, 3000, 1e4, 3e4, 1e5, 3e5, 1e6]
acc = 0.0; t0 = 0.0; out = {}
for T in [e for e in edges if e <= Tmax] + [Tmax]:
    if T <= t0: continue
    for c0 in np.arange(t0, T, 2e5*h):
        tau = np.arange(c0, min(T, c0 + 2e5*h), h) + h/2
        acc += np.clip(D(tau), 0, None).sum()*h
    t0 = T; out[str(T)] = round(acc/math.pi, 2)     # (1/2pi) * 2 (both signs of tau)
th = np.linspace(0, 2*math.pi, 4096, endpoint=False)
ld = 1.0
for p in primes:
    a = p**-0.5; L = math.log(p); c = np.cos(th)
    ld *= float(np.mean(np.exp(2*L*(a*c - a*a)/(1 - 2*a*c + a*a))))
M = sum(2*math.log(p)/(math.sqrt(p) - 1) for p in primes)
print(json.dumps({'primes': primes, 't': tsc, 'M_S': round(M, 3), 'horizon_2pi_e^M': round(2*math.pi*math.exp(M + 2*math.log(tsc)), 1), 'N_sc_cum': out, 'N_ld_prod': round(ld, 3)}))
