# make_cert.py: optimal SDP + eta * strictly interior SDP (margin) -> certificate JSON (decimal Hermite coefficients)
# usage: OMP_NUM_THREADS=1 uv run --no-project --with cvxpy --with numpy --with scipy python3 make_cert.py d s lam d_int margin eta out.json
#   s as fraction string (e.g. 11/10), lam decimal string (e.g. 0.7852)
import sys, json, math, gc, numpy as np
from fractions import Fraction
from sdp import solve, check
from scipy.special import gammaln
d = int(sys.argv[1]); sstr = sys.argv[2]; s = float(Fraction(sstr)); lamstr = sys.argv[3]; lam = float(lamstr)
d_int = int(sys.argv[4]); margin = float(sys.argv[5]); eta = float(sys.argv[6]); out = sys.argv[7]
v, a, st = solve(d, s, lam); gc.collect()
print('optimal SDP', d, sstr, lam, '->', v, st, check(d, s, lam, a), flush=True)
vi, ai, sti = solve(d_int, s, lam, margin=margin); gc.collect()
print('interior SDP', d_int, 'margin', margin, '->', vi, sti, check(d_int, s, lam, ai), flush=True)
ai = np.concatenate([ai, np.zeros(d - d_int)])
ac = (1-eta)*a + eta*ai
print('combined (eta=%g)' % eta, check(d, s, lam, ac), flush=True)
B = []
for k in range(d+1):
    logc = 0.25*math.log(2*math.pi) - 0.5*(2*k*math.log(2) + gammaln(2*k+1) + 0.5*math.log(math.pi))
    B.append('%.20e' % (ac[k]*math.exp(logc)))
json.dump(dict(d=d, s=sstr, lam=lamstr, theta='399993/700000', d_int=d_int, margin=margin, eta=eta, sdp_value=v, interior_value=vi,
               note='f(x) = P(sqrt(2 pi) x / s) exp(-pi x^2/s^2), P(u) = sum_k B[k] H_{2k}(u) (physicists Hermite)', B=B), open(out, 'w'), indent=1)
print('wrote', out)
