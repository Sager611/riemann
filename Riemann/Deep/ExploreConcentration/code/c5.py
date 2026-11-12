# explore44 c5: Weissler-Bayart shift needed for H^q <= H^2, with small primes removed; and the [T,2T] failure of hypercontractivity
# for a resonant Dirichlet polynomial of length N = T^theta (4th moment over [T,2T] vs the T^inf (Haar) 4th moment).
import numpy as np, time
t0 = time.time()
def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
P = primes_upto(2*10**6).astype(float)
print('Bayart: ||F(.+s)||_q <= ||F||_2 needs p^-s <= sqrt(2/q) for every prime in the support; s(q) = log(q/2)/(2 log p_min)')
for q in [3, 4, 6, 8, 12]:
    row = ['q=%d s(q; p>=2)=%.3f' % (q, np.log(q/2)/(2*np.log(2)))]
    for pmin in [11, 101, 1009, 10007, 1000003]:
        s = np.log(q/2)/(2*np.log(pmin)); sm = P[P < pmin]
        cost = np.exp(-np.sum(np.log(1 - sm**(-(0.5 + s)))))   # bound for the removed Euler factors at real part 1/2 + s
        row.append('p>=%d: s=%.3f cost=%.2e' % (pmin, s, cost))
    print('  ' + ' | '.join(row))
print('resonance: D(t) = sum_{n<=N} n^{-1/2 - i(t - t0)}, t in [T, 2T], t0 = 1.5 T')
T = 2000.0; t = np.arange(T, 2*T, 0.01); t0r = 1.5*T
for N in [20, 45, 100, 200, 400, 1000]:
    n = np.arange(1, N + 1, dtype=float); ln = np.log(n); w = n**-0.5
    m4 = 0.0; m2 = 0.0
    for i in range(0, len(t), 2000):
        D = (w[None, :]*np.exp(-1j*np.outer(t[i:i+2000] - t0r, ln))).sum(1)
        a = np.abs(D)**2; m2 += a.sum(); m4 += (a*a).sum()
    m2 /= len(t); m4 /= len(t)
    # Haar (T^inf) moments: E|D|^2 = sum 1/n ; E|D|^4 = sum_m r_N(m)^2/m, r_N(m) = #{(a,b) <= N : ab = m}
    r = {}
    for a_ in range(1, N + 1):
        for b_ in range(1, N + 1):
            r[a_*b_] = r.get(a_*b_, 0) + 1
    h2 = float(np.sum(1/n)); h4 = float(sum(v*v/m for m, v in r.items()))
    print('  N=%5d theta=%.2f  [T,2T]: E|D|^2=%.3f E|D|^4=%.2f | Haar: E|D|^2=%.3f E|D|^4=%.2f | ratio4=%.2f  hypercontr. bound (q/2)^? fails if ratio4 >> 1' % (N, np.log(N)/np.log(T), m2, m4, h2, h4, m4/h4), flush=True)
print('done', round(time.time() - t0, 1))
