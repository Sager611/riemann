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
A = np.array([p[2] for p in pk])
print(json.dumps({'primes': primes, 'T0': T0, 'T1': T1, 'n_pockets': len(pk), 'total_area_over_2pi': round(float(A.sum())/(2*math.pi), 3),
                  'n_area_gt_2pi': int((A > 2*math.pi).sum()), 'area_quantiles': [round(float(x), 4) for x in np.quantile(A, [0.5, 0.9, 0.99])] if len(A) else [],
                  'examples_q50_q90': [{'tau_peak': round(p[4], 2), 'width': round(p[1] - p[0], 3), 'Dmax': round(p[3], 3), 'area': round(p[2], 3)}
                        for q in (0.5, 0.9) for p in [pp for pp in pk if abs(pp[2] - float(np.quantile(A, q))) < 0.1*float(np.quantile(A, q))][:3]],
                  'largest': [{'tau_peak': round(p[4], 2), 'width': round(p[1] - p[0], 3), 'Dmax': round(p[3], 3), 'area': round(p[2], 3)} for p in sorted(pk, key=lambda p: -p[2])[:12]]}))
