# t2_chowla_fft.py — explore6 T2: where in the shift range h does the size of M(X) live?
# mu(n), n <= X = 2^22; C(h) = sum_{n<=X-h} mu(n) mu(n+h) by FFT; identities
#   M(X)^2 = C(0) + 2 sum_{h>=1} C(h)        and   sum_{x<=X} M(x)^2 = D(0) + 2 sum_{h>=1} D(h),
#   D(h) = sum_n mu(n) mu(n+h) (X+1-n-h);  short-interval variance V(H) vs Good-Churchhouse 6H/pi^2;
#   log-average of M(x)^2/x vs 6/pi^2 and vs Ng's constant 2 sum_rho 1/|rho zeta'(rho)|^2 (300 zeros, mpmath).
# Same statistics for a random-sign control on the squarefrees.
# Run: OMP_NUM_THREADS=1 perl -e 'alarm 1800; exec @ARGV' uv run --no-project --with numpy --with mpmath python3 t2_chowla_fft.py
import numpy as np, math, time, sys
t0 = time.time()
LOGX = int(sys.argv[1]) if len(sys.argv) > 1 else 22
X = 1 << LOGX
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
mu8 = mobius(X)
assert int(mu8[1:1001].astype(int).sum()) == 2 and (X < 10**6 or int(mu8[1:10**6+1].astype(int).sum()) == 212)
def stats(a, label):
    # a[1..X] (a[0] = 0)
    a = a.astype(np.float64)
    L2 = 1 << (LOGX + 1)
    A = np.fft.rfft(a[1:], L2)
    C = np.rint(np.fft.irfft(A * np.conj(A), L2)[:X])          # C[h] = sum_n a_n a_{n+h}
    del A
    Mx = np.cumsum(a)                                             # Mx[x] = sum_{n<=x}
    MX = Mx[X]
    ident = C[0] + 2 * C[1:].sum() - MX**2
    w = a[1:] * (X + 1 - np.arange(1, X + 1))                     # w(m) = (X+1-m) a_m
    Aa = np.fft.rfft(a[1:], L2); Aw = np.fft.rfft(w, L2)
    D = np.fft.irfft(np.conj(Aa) * Aw, L2)[:X]                    # D[h] = sum_n a_n w_{n+h}
    del Aa, Aw
    msq = float((Mx[1:]**2).sum())
    ident2 = D[0] + 2 * D[1:].sum() - msq
    out = [f"## {label}: X = 2^{LOGX} = {X}, M(X) = {int(MX)}, C(0) = {int(C[0])}"]
    out.append(f"identity M(X)^2 = C(0)+2 sum C(h): residual {ident:.3e};  sum_x M(x)^2 = D(0)+2 sum D(h): residual/msq {ident2/msq:.2e}")
    Sabs = C[0] + 2 * np.abs(C[1:]).sum()
    out.append(f"triangle route: sqrt(C(0)+2 sum|C(h)|) = {math.sqrt(Sabs):.4g} = X^{math.log(math.sqrt(Sabs))/math.log(X):.4f}  vs |M(X)| = {abs(MX):.0f} = X^{math.log(max(abs(MX),1))/math.log(X):.4f}")
    rmsC = math.sqrt(float((C[1:X//2]**2).mean())); out.append(f"rms C(h), 1<=h<X/2: {rmsC:.4g} = X^{math.log(rmsC)/math.log(X):.4f}")
    out.append("H        Q(H)/C(0)=2sum_{h<=H}C/C(0)   R(H)/D(0)=(D(0)+2sum_{h<=H}D)/D(0)   V(H)/(6H/pi^2)  [x in [X/4,X/2]]")
    cC = np.cumsum(C[1:]); cD = np.cumsum(D[1:])
    for j in list(range(0, LOGX, 2)) + [LOGX - 1]:
        H = 1 << j
        H = min(H, X - 1)
        Q = 2 * cC[H - 1] / C[0]; Rr = (D[0] + 2 * cD[H - 1]) / D[0]
        if H <= X // 2:
            xs = np.arange(X // 4, X // 2 + 1); inc = Mx[xs + H] - Mx[xs]; V = float((inc**2).mean()) / (6 * H / math.pi**2)
        else:
            V = float('nan')
        out.append(f"2^{j:<5d}  {Q: .4f}                        {Rr: .4f}                              {V:.4f}")
    out.append(f"H = X-1: Q = {2*cC[-1]/C[0]: .4f} (= M^2/C(0) - 1),  R = {(D[0]+2*cD[-1])/D[0]: .4f} (= sum M^2 / D(0))")
    xs = np.arange(1000, X + 1)
    la = float((Mx[xs]**2 / xs.astype(np.float64)**2).sum()) / math.log(X / 1000)
    out.append(f"log-average of M(x)^2/x over [1e3, X]: {la:.4f}   (6/pi^2 = {6/math.pi**2:.4f})")
    return "\n".join(out)
res = [stats(mu8, "mu")]
rng = np.random.default_rng(2026)
eps = (mu8 != 0) * (2 * rng.integers(0, 2, X + 1) - 1); eps[0] = 0
res.append(stats(eps.astype(np.int8), "random signs on squarefrees (control)"))
try:
    import mpmath as mp
    mp.mp.dps = 20
    s = mp.mpf(0)
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 300
    for k in range(1, K + 1):
        r = mp.zetazero(k); d = mp.zeta(r, derivative=1)
        s += 2 / abs(r * d)**2               # rho and its conjugate
    res.append(f"Ng constant sum over all zeros rho of 1/|rho zeta'(rho)|^2 (= log-mean of M(x)^2/x under RH + simple zeros), first {K} conjugate pairs: {float(s):.5f}")
except Exception as e:
    res.append("mpmath part failed: " + str(e))
res.append(f"seconds {time.time()-t0:.1f}")
print("\n\n".join(res))
