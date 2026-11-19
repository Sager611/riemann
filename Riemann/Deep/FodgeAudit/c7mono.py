# FodgeAudit c7mono.py: the key step of Prop. C.7's proof: with N = s y, L = s l and the integers (M, K, J) fixed, B/s^2 is non-increasing in s
# (B = squared block bound of C.2+C.3, C.6 or C.8). Floating-point test on random (y, l, M, K, J) and s-grids from s_k to 50 s_k (mpmath, 30 digits).
import random, mpmath as mp
mp.mp.dps = 30
pi = mp.pi
def B1(s, y, l, M):            # C.2 (flat w) + C.3: (L+M-1)/M (L + 2 sum (1-m/M) min(L, rho_m))
    N, L = s*y, s*l; Np = N - mp.mpf(1)/2; W = pi*(y+l)**3; mu = 2*(N+L)**1.5/L*(1/mp.sqrt(Np) - 1/mp.sqrt(Np+L))
    rho = lambda m: mu*L*mp.sqrt(m/W) + mp.mpf('4.524')*mp.sqrt(W/m) + mp.mpf('1.0004') + mp.mpf('2.24')*m/W
    return (L+M-1)/M*(L + 2*mp.fsum((1-mp.mpf(m)/M)*min(L, rho(m)) for m in range(1, M)))
def B2(s, y, l, M, K):         # C.6 with flat outer weights
    N, L = s*y, s*l; Np = N - mp.mpf(1)/2; V = pi*s*(y+l)**4/3; mu = (N+L)**2/(Np*(Np+L))
    Q = lambda m: L/K + mp.mpf(8)/15*mu*L*mp.sqrt(m/V)*mp.sqrt(K) + mp.mpf(8)/3*mp.mpf('4.527')*mp.sqrt(V/m)/mp.sqrt(K) + mp.mpf('1.0006') + mp.mpf('2.272')/3*(m/V)*K
    return (L+M-1)/M*(L + 2*mp.fsum((1-mp.mpf(m)/M)*min(L, mp.sqrt((L+K-1)*Q(m))) for m in range(1, M)))
def B3(s, y, l, M, K, J):      # C.8 with flat weights at all levels
    N, L = s*y, s*l; Np = N - mp.mpf(1)/2; V3 = pi*s**2*(y+l)**5/12; mu3 = 2*(N+L)**2.5/(3*L)*(Np**-1.5 - (Np+L)**-1.5)
    def tau(m, k):
        Q = L/J + mp.mpf(8)/15*mu3*L*mp.sqrt(m*k/V3)*mp.sqrt(J) + mp.mpf(8)/3*mp.mpf('4.531')*mp.sqrt(V3/(m*k))/mp.sqrt(J) + mp.mpf('1.001') + mp.mpf('2.31')/3*(m*k/V3)*J
        return min(L, mp.sqrt((L+J-1)*Q))
    def sig(m): return min(L, mp.sqrt((L+K-1)*(L/K + 2/mp.mpf(K)*mp.fsum((1-mp.mpf(k)/K)*tau(m, k) for k in range(1, K)))))
    return (L+M-1)/M*(L + 2*mp.fsum((1-mp.mpf(m)/M)*sig(m) for m in range(1, M)))
random.seed(13); worst = {1: -1e9, 2: -1e9, 3: -1e9}; n = 0
for it in range(120):
    nu_r = random.uniform(0.004, 3); y = mp.mpf(nu_r); l = y*(mp.mpf(1)/80 + random.uniform(0, 0.0008))
    sk = max(mp.mpf(2*10**4)/y, 2*pi*y*y)*mp.mpf(random.uniform(1, 3))
    ss = [sk*f for f in (1, 1.001, 1.01, 1.1, 1.5, 3, 10, 50)]
    for kind in (1, 2, 3):
        M = random.randint(2, 12); K = random.randint(1, 8); J = random.randint(1, 6)
        f = (lambda s: B1(s, y, l, M)) if kind == 1 else ((lambda s: B2(s, y, l, M, K)) if kind == 2 else (lambda s: B3(s, y, l, M, K, J)))
        vals = [f(s)/s**2 for s in ss]; n += 1
        for a, b in zip(vals, vals[1:]): worst[kind] = max(worst[kind], float((b - a)/a))   # must be <= 0
print(f'C.7 monotonicity: {n} (y, l, integers) samples x 8 values of s in [s_k, 50 s_k]; max relative increase of B/s^2 between consecutive s: '
      f'C.2+C.3 {worst[1]:.3e}, C.6 {worst[2]:.3e}, C.8 {worst[3]:.3e} (all must be <= 0)')
