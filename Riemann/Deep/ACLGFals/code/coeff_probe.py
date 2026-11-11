# coeff_probe.py: decay of the Hankel coefficients phi_m (m > 0) of the semilocal symbol, S = {inf} vs {inf, 2}; diagnostics only.
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS'): os.environ[k] = '1'
import sys, math, json, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from hankel import symbol, chi_cut, circle_coeffs, tail_coeffs
mmax = int(sys.argv[1]); T = float(sys.argv[2]); T1 = float(sys.argv[3])
out = {}
for primes in ([], [2]):
    rex = sum(math.log(p) for p in primes)
    sym = lambda tau: symbol(tau, 1.0, primes)
    c1, Nth = circle_coeffs(lambda tau: chi_cut(tau, T1/2)*sym(tau), T1, mmax, rex)
    c2, npt = tail_coeffs(lambda tau: (1 - chi_cut(tau, T1/2))*sym(tau), T1/2, T, mmax, rex)
    ph = (c1 + c2).real
    # also negative-index coefficients: phi_{-m} = coefficient of conj side; compute by fun = conj(sym) (since (conj phi)_m = phi_{-m})
    d1, _ = circle_coeffs(lambda tau: chi_cut(tau, T1/2)*np.conj(sym(tau)), T1, mmax, rex)
    d2, _ = tail_coeffs(lambda tau: (1 - chi_cut(tau, T1/2))*np.conj(sym(tau)), T1/2, T, mmax, rex)
    phm = (d1 + d2).real
    key = 'S=' + ','.join(map(str, primes))
    out[key] = {'pos': {m: float(ph[m]) for m in [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 1500, 1999] if m < mmax},
                'neg': {m: float(phm[m]) for m in [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 1500, 1999] if m < mmax},
                'tail_only_pos_max_m>100': float(np.max(np.abs(c2.real[100:]))), 'circle_pos_max_m>100': float(np.max(np.abs(c1.real[100:])))}
print(json.dumps(out, indent=0))
