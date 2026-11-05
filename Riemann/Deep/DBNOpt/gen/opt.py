import time, math, sys
from dbn import *
th = 5458637291332547 / 6250000000000000; q0n = (2 * th - 1) ** 2
def run(q0, h, mU=16, mV=1, verbose=False):
    b = 1 / 25; t = 0.0; q = q0; rows = []
    while True:
        s = scalars_f(t, t + h); C, G, ep0 = CG(s)
        lam, mu, k = lp_best(q, C, G, mU, mV)
        qn = (q + 2 / k) * math.exp(-k * h) - 2 / k
        rows.append((t, h, q, C, G, lam, mu, k, qn, ep0))
        if qn <= b:
            dd = math.log((2 + k * q) / (2 + k * b)) / k
            # recompute with the shorter final step
            s = scalars_f(t, t + dd); C, G, ep0 = CG(s); lam, mu, k = lp_best(q, C, G, mU, mV)
            dd = math.log((2 + k * q) / (2 + k * b)) / k
            rows[-1] = (t, dd, q, C, G, lam, mu, k, b, ep0)
            return t + dd, rows
        t += h; q = qn
t0 = time.time()
for h in [float(x) for x in sys.argv[1:]]:
    tf, rows = run(q0n, h)
    tf9, rows9 = run(9 / 16, h)
    print("h=%.6f rows=%d Lambda(new)=%.6f Lambda(9/16)=%.6f  maxeps0=%.4f  %.1fs" % (h, len(rows), tf + 0.02, tf9 + 0.02, max(r[9] for r in rows), time.time() - t0), flush=True)
