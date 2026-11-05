# t2b_scaling.py — explore6 T2 follow-up: the rigidity of M sits at shifts/windows H ~ X (scales with X), and the
# random-sign control averaged over 32 seeds. No FFT: Q(H), R(H), V(H) from cumulative sums.
#   Q(H) = 2 sum_{1<=h<=H} C_X(h) / C(0),  C_X(h) = sum_{n<=X-h} a_n a_{n+h}
#   R(H) = (D(0) + 2 sum_{1<=h<=H} D(h)) / D(0),  D(h) = sum_n a_n a_{n+h} (X+1-n-h)   [R(X-1) = sum_{x<=X} S(x)^2 / D(0)]
#   V(H) = mean_{x in [X/4, X/2]} (S(x+H)-S(x))^2 / (6H/pi^2)
import numpy as np, math, sys, time
t0 = time.time()
LOGX = 22; XM = 1 << LOGX
def mobius(N):
    mu = np.ones(N + 1, np.int8); acc = np.zeros(N + 1, np.float32)
    R = int(math.isqrt(N)) + 1
    sv = np.ones(R + 1, bool); sv[:2] = False
    for i in range(2, int(R**0.5) + 1):
        if sv[i]: sv[i*i::i] = False
    for p in np.nonzero(sv)[0].tolist():
        mu[::p] *= -1; acc[::p] += math.log(p); mu[::p*p] = 0
    n = np.arange(N + 1); big = acc < (np.log(np.maximum(n, 1)) - 1.0); mu[big] *= -1
    mu[0] = 0
    return mu
mu = mobius(XM)
assert int(mu[1:10**6 + 1].astype(int).sum()) == 212
def qrv(a, X, Hs):
    a = a[:X + 1].astype(np.float64)
    S = np.cumsum(a)                                  # S[x] = sum_{n<=x}
    w = a * (X + 1 - np.arange(X + 1)); w[0] = 0
    W = np.cumsum(w)
    C0 = float((a * a).sum()); D0 = float((a * w).sum())
    n = np.arange(1, X + 1)
    out = []
    for H in Hs:
        top = np.minimum(n + H, X)
        q = float((a[1:] * (S[top] - S[n])).sum())
        r = float((a[1:] * (W[top] - W[n])).sum())
        xs = np.arange(X // 4, X // 2 + 1)
        v = float(((S[xs + H] - S[xs])**2).mean()) / (6 * H / math.pi**2) if H <= X // 2 else float('nan')
        out.append((2 * q / C0, (D0 + 2 * r) / D0, v))
    la = float((S[1000:]**2 / np.arange(1000, X + 1, dtype=np.float64)**2).sum()) / math.log(X / 1000)
    return out, la
print("# mu: rigidity location scales with X (columns: H = X/2^j)")
js = [1, 2, 3, 4, 6, 8, 10, 12]
for LX in (18, 20, 22):
    X = 1 << LX
    res, la = qrv(mu, X, [X >> j for j in js])
    print(f"X = 2^{LX}: log-mean M^2/x = {la:.4f}")
    print("   j (H=X/2^j): " + "  ".join(f"{j:>6d}" for j in js))
    print("   Q(H)       : " + "  ".join(f"{q:6.3f}" for q, r, v in res))
    print("   R(H)       : " + "  ".join(f"{r:6.3f}" for q, r, v in res))
    print("   V(H)/GC    : " + "  ".join(f"{v:6.3f}" for q, r, v in res))
print()
NS = 32
print(f"# random signs on squarefrees, X = 2^22, {NS} seeds: mean (sd)")
X = XM; Hs = [X >> j for j in js]
allr = []; las = []
for seed in range(NS):
    rng = np.random.default_rng(100 + seed)
    eps = (mu != 0) * (2 * rng.integers(0, 2, X + 1) - 1); eps[0] = 0
    res, la = qrv(eps.astype(np.int8), X, Hs); allr.append(res); las.append(la)
A = np.array(allr)   # seeds x H x 3
print(f"log-mean S^2/x: {np.mean(las):.3f} ({np.std(las):.3f})   [6/pi^2 = {6/math.pi**2:.4f}]")
for k, name in enumerate(("Q(H)", "R(H)", "V(H)/GC")):
    print(f"   {name:<10s}: " + "  ".join(f"{A[:, i, k].mean():6.3f}({A[:, i, k].std():.2f})" for i in range(len(js))))
print(f"seconds {time.time()-t0:.1f}")
