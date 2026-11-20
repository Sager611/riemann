# ceilaudit6: copy of ACLGH2/code/pockets.py that dumps every pocket [tau_peak, area, width] on [T0, T1].
# pockets.py: coupling pockets {tau : D(tau) > 0} of the semilocal symbol at t = 1 (D = total delay, see semicl.py) on [T0, T1]:
# per pocket its centre, width, max D and area A = int D_+ (semiclassical mode count A/2pi). Prints summary and the largest pockets.
# Usage: pockets.py primes T0 T1 [h]
import sys, math, json, numpy as np
from scipy.special import digamma
primes = [int(q) for q in sys.argv[1].split(',')]; T0 = float(sys.argv[2]); T1 = float(sys.argv[3])
h = float(sys.argv[4]) if len(sys.argv) > 4 else 0.02
def D(tau):
    out = -(np.real(digamma(0.25 + 0.5j*tau)) - math.log(math.pi))
    for p in primes:
        a = p**-0.5; L = math.log(p); c = np.cos(tau*L)
        out = out + 2*L*(a*c - a*a)/(1 - 2*a*c + a*a)
    return out
pk = []; cur = None
for c0 in np.arange(T0, T1, 1e6*h):
    tau = np.arange(c0, min(T1, c0 + 1e6*h), h)
    d = D(tau); pos = d > 0
    if not pos.any() and cur is None: continue
    idx = np.nonzero(np.diff(np.concatenate([[0], pos.astype(np.int8), [0]])))[0]
    for s, e in zip(idx[0::2], idx[1::2]):
        seg = d[s:e]
        piece = (float(tau[s]), float(tau[e - 1]), float(seg.sum()*h), float(seg.max()), float(tau[s + np.argmax(seg)]))
        if cur is not None and s == 0 and abs(piece[0] - cur[1] - h) < 1.5*h:
            cur = (cur[0], piece[1], cur[2] + piece[2], max(cur[3], piece[3]), cur[4] if cur[3] >= piece[3] else piece[4])
        else:
            if cur is not None: pk.append(cur)
            cur = piece
        if e < len(tau): pk.append(cur); cur = None
if cur is not None: pk.append(cur)
print(json.dumps([[round(p[4], 3), round(p[2], 5), round(p[1]-p[0], 3)] for p in pk]))
