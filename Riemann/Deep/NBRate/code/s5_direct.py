# nbrate: independent check of BCF ramp values by direct integration on the critical line
# I(N) = (1/2pi) int |1 - zeta V_N(1/2+it)|^2 dt/(1/4+t^2),  V_N = sum_{n<=N} mu(n)(1 - log n/log N) n^-s
import mpmath as mp, math, sys
mp.mp.dps = 15
def mobius(n):
    r = 1; p = 2; m = n
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    if m > 1: r = -r
    return r
def run(N, T, h):
    L = math.log(N)
    co = [(n, mobius(n) * (1 - math.log(n) / L)) for n in range(1, N + 1)]
    co = [(n, w) for (n, w) in co if w != 0]
    lg = [math.log(n) for (n, w) in co]
    tot = 0.0; t = 0.0
    # Simpson on [0,T], symmetric integrand (real coefficients) -> factor 2
    import cmath
    def f(t):
        s = complex(0.5, t)
        z = complex(mp.zeta(s))
        V = sum(w * cmath.exp(-s * l) for (n, w), l in zip(co, lg))
        return abs(1 - z * V) ** 2 / (0.25 + t * t)
    n = int(T / h); n += n % 2
    hh = T / n
    acc = f(0.0) + f(T)
    for i in range(1, n):
        acc += (4 if i % 2 else 2) * f(i * hh)
    I = 2 * acc * hh / 3 / (2 * math.pi)
    # tail t > T: mean value |zeta V|^2 ~ (log(t/2pi)+2 gamma) sum w^2/n  (diagonal), so tail ~ (1/pi) int_T^inf (1 + that) dt/t^2
    S2 = sum(w * w / n for (n, w) in co)
    tail = (1 / math.pi) * ((1 + S2 * (math.log(T / (2 * math.pi)) + 2 * 0.5772156649 + 1)) / T)
    return I, tail
for N, T, h in [(2, 400, 0.02), (100, 400, 0.02)]:
    I, tail = run(N, T, h)
    print(f'N={N}: int_0^T part = {I:.6f}, tail estimate ~ {tail:.6f}, total ~ {I+tail:.6f}', flush=True)
