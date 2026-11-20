# probe (a2): the archive's mean-square regime Z^2 < d (QH.meanSq_abs needs floor(Z)^2 < q):
# mean over primes d in (1e6, 2e6] of |S(Z,y)|^2, S = sum_{n<=Z} Lambda(n) min(1,y/n) xi(n) (n/d), Z = 1000,
# versus the exact 'square-pairing' diagonal Diag(y) = sum_{nm = square} c_n c_m (c_n = Lambda(n) min(1,y/n) |xi(n)|^2-type)
import numpy as np, math, time
t0 = time.time()
Z, D0, D1 = 1000, 10**6, 2 * 10**6
sv = np.ones(D1 + 1, dtype=bool); sv[:2] = False
for i in range(2, int(D1 ** 0.5) + 1):
    if sv[i]: sv[i*i::i] = False
allP = np.nonzero(sv)[0].astype(np.int64)
P = allP[allP <= Z]
D = allP[(allP > D0)]
pp = []   # (n, log p, idx, k)
for i, p in enumerate(P):
    n, k = int(p), 1
    while n <= Z:
        pp.append((n, math.log(p), i, k)); n *= int(p); k += 1
pp.sort()
N = np.array([x[0] for x in pp], dtype=np.float64); LP = np.array([x[1] for x in pp]); IB = np.array([x[2] for x in pp]); K = np.array([x[3] for x in pp])
Ni = N.astype(np.int64)
ys = np.array([10, 20, 50, 100, 200, 500, 1000], dtype=np.float64)
W = np.minimum(1.0, ys[:, None] / N[None, :])          # (ny, npp)
xi4 = np.where(Ni % 2 == 0, 0.0, np.where(Ni % 4 == 1, 1.0, -1.0))
def diag(xi):
    out = []
    for r in range(len(ys)):
        c = LP * W[r] * (xi if xi is not None else 1.0)
        # same prime, a+b even
        same = 0.0; evsum_by_p = {}
        byp = {}
        for j in range(len(c)):
            byp.setdefault(IB[j], []).append((K[j], c[j]))
        tot_even = 0.0; sq_even = 0.0
        for i, lst in byp.items():
            for (a, ca) in lst:
                for (b, cb) in lst:
                    if (a + b) % 2 == 0: same += ca * cb
            ev = sum(cv for (a, cv) in lst if a % 2 == 0)
            tot_even += ev; sq_even += ev * ev
        out.append(same + tot_even ** 2 - sq_even)
    return np.array(out)
print('#d', len(D), 'd in', int(D[0]), int(D[-1]), '#prime powers <= Z', len(N), flush=True)
for name, xi in (('M=1', None), ('M=4', xi4)):
    acc = np.zeros(len(ys)); mx = []
    for d in D:
        d = int(d); e = (d - 1) // 2
        b = P % d; res = np.ones_like(b)
        while e:
            if e & 1: res = (res * b) % d
            b = (b * b) % d; e >>= 1
        chi = np.where(res == 1, 1.0, -1.0)
        a = LP * chi[IB] ** K
        if xi is not None: a = a * xi
        S = W @ a
        acc += S * S
        mx.append(np.max(np.abs(S) / np.sqrt(ys * np.log(ys))))
    emp = acc / len(D); dg = diag(xi)
    print(name, 'y:', ' '.join('%d' % v for v in ys))
    print(name, 'mean |S|^2/(y log y):', ' '.join('%.3f' % v for v in emp / (ys * np.log(ys))))
    print(name, 'square-pairing diagonal/(y log y):', ' '.join('%.3f' % v for v in dg / (ys * np.log(ys))))
    print(name, 'ratio empirical/diagonal:', ' '.join('%.4f' % v for v in emp / dg))
    print(name, 'max_y |S|/sqrt(y log y) quantiles 0/.5/.9/.99/1:', ' '.join('%.3f' % v for v in np.quantile(mx, [0, .5, .9, .99, 1])), ' time %.1fs' % (time.time() - t0), flush=True)
