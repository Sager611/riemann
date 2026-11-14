# Sharper rigorous bound c = sup_{y>=1} S(y)/sqrt(y), S(y) = sum_{m<=y} mu^2(m) log(y/m)/sqrt(m),
# giving |G(1/2+it)| <= sum |theta_j theta_k| [j,k]^(-1/2) <= c^2 P W/log^2 W  (c = 4 is the crude integral bound).
# y <= Y0: exact sup on every interval between consecutive squarefree m (S = A log y - B there; critical point y* = e^(2+B/A)).
# y >= Y0: S(y) = int_1^y Q(u) (-f'(u)) du, Q(u) <= 6u/pi^2 + 3 sqrt(u) (elementary, u >= 1) gives
#          S(y) <= (6/pi^2)(4 sqrt y - 4 - log y) + 3 (log y + log^2 y / 4); (log y + log^2 y/4)/sqrt y decreases for y > 17.
import math, numpy as np
Y0 = 2 * 10**6
q = np.ones(Y0 + 2, dtype=bool); q[0] = False
for p in range(2, int((Y0 + 1) ** 0.5) + 1):
    q[p*p::p*p] = False
m = np.nonzero(q)[0].astype(np.float64)          # squarefree 1 <= m <= Y0+1
A = np.cumsum(1 / np.sqrt(m)); B = np.cumsum(np.log(m) / np.sqrt(m))
lo = m[:-1]; hi = m[1:]; A = A[:-1]; B = B[:-1]  # S(y) = A log y - B on [lo, hi)
r_lo = (A * np.log(lo) - B) / np.sqrt(lo)
r_hi = (A * np.log(hi) - B) / np.sqrt(hi)        # left limit at hi (upper end of the interval)
ys = np.exp(2 + B / A); inside = (ys >= lo) & (ys < hi)
r_st = np.where(inside, (A * np.log(np.where(inside, ys, 1.0)) - B) / np.sqrt(np.where(inside, ys, 1.0)), -np.inf)
best = np.maximum(np.maximum(r_lo, r_hi), r_st)
k = int(np.argmax(best)); c_small = float(best[k])
L = math.log(Y0)
c_tail = 24 / math.pi**2 + 3 * (L + L * L / 4) / math.sqrt(Y0)
c = max(c_small * (1 + 1e-9) + 1e-9, c_tail)
print('squarefree <= Y0:', len(m) - 1, ' sup_{y<=Y0} S/sqrt(y) = %.9f near y = %.0f' % (c_small, lo[k]))
print('ratio at Y0 = %.6f (asymptote 24/pi^2 = %.6f); tail bound for y >= Y0: %.6f' % (float(r_lo[-1]), 24 / math.pi**2, c_tail))
P = 1.47252045
print('c = %.6f, c^2 = %.6f, c^2 P <= %.4f (crude 16P = 23.5604)' % (c, c * c, c * c * P))
