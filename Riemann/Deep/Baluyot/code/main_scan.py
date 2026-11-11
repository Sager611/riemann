# Part 2a: exact diagonal M1 (thesis (3.3.9), I=1) vs the thesis's main term 1 + O(y^{1-2sigma})
# and vs its residue formula M11 = 1/W(y)^2 + A y^{1-2sigma}/W(y)^2 (3.5.9 with F1 = 1).
import sys, time, json
import numpy as np
from mollifier import coeffs, main_terms

out = []
ys = [float(v) for v in sys.argv[1].split(',')] if len(sys.argv) > 1 else [1e2, 1e3, 1e4, 1e5]
print("%8s %7s %6s %9s %10s %10s %9s %10s %10s %9s" % ("y", "sigma", "c", "W(y)", "M1", "M11form", "y^(1-2s)",
                                                        "(M1-1)/yp", "(M1-M11)/yp", "S2"))
for y in ys:
    L = np.log(y)
    sigmas = [0.5 + c / L for c in (1, 2, 4)] + [s for s in (0.6, 0.75) if s >= 0.5 + 1 / L]
    N = int(20 * y)
    for sigma in sigmas:
        t0 = time.time()
        a, mu, primes, W = coeffs(y, sigma, N)
        M1, S2 = main_terms(a, mu, primes, sigma)
        Wy = W(np.array([y]))[0]
        yp = y ** (1 - 2 * sigma)
        M11 = (1 + W.A * yp) / Wy ** 2
        r = dict(y=y, sigma=sigma, c=(sigma - 0.5) * L, W=Wy, M1=M1, M11=M11, yp=yp,
                 r1=(M1 - 1) / yp, r2=(M1 - M11) / yp, S2=S2, N=N, sec=time.time() - t0)
        out.append(r)
        print("%8.0e %7.4f %6.2f %9.5f %10.6f %10.6f %9.2e %10.4f %10.4f %9.5f" % (
            y, sigma, r['c'], Wy, M1, M11, yp, r['r1'], r['r2'], S2), flush=True)
json.dump(out, open('../logs/main_scan.json', 'w'), indent=1)
