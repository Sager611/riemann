# nflp: tails for the L^p table of nbk.c (K = Q(i), Q(sqrt5)); reads run1e8.txt (nbk.c output, M = 1e8).
# (1) rigorous tail bound: |v(m)| <= V0 + V1 sqrt(m) from |Delta_K(u)| <= C1 sqrt(u) + C0 (Q(i): lattice squares; Q(sqrt5): hyperbola method,
#     |S_chi| <= 1), |lambda_k(I_m)| <= c_K + |Delta(m)|/k + |Delta(m/k)|, and sum_{m>M} (V0+V1 sqrt m)^p/m^2 <= 2^{p-1}(V0^p/M + V1^p M^{p/2-1}/(1-p/2)).
# (2) heuristic tail: geometric continuation of the last decade increment (ratio = last/previous increment).
# Run: python3 nbk_tail.py   (stdlib only)
import math, re
FIELD = {1: ("Q(i)", lambda n: (0, 1, 0, -1)[n % 4], math.pi / 4, math.sqrt(2) * math.pi / 4, (math.pi / 2 + 1) / 4),
         2: ("Q(sqrt5)", lambda n: (0, 1, -1, -1, 1)[n % 5], 2 * math.log((1 + 5 ** 0.5) / 2) / 5 ** 0.5, 5.0, 0.5)}
def bK(k, chi):
    v, n, p = 1, k, 2
    while n > 1:
        if p * p > n: p = n
        e = 0
        while n % p == 0: n //= p; e += 1
        if e:
            c = chi(p); v *= (-(1 + c) if e == 1 else (c if e == 2 else 0))
            if v == 0: return 0
        p += 1
    return v
runs = {}
for line in open("run1e8.txt"):
    m = re.match(r"type=(\d) N=(\d+) m=(\d+) S32/31=([\d.]+) S1.1=([\d.]+) S1.5=([\d.]+) S2=([\d.]+)", line)
    if m:
        t, N, mm = int(m.group(1)), int(m.group(2)), int(m.group(3))
        runs.setdefault((t, N), {})[mm] = [float(m.group(i)) for i in (4, 5, 6, 7)]
ps = [32 / 31, 1.1, 1.5, 2.0]
M = 10 ** 8
print("field  N   p      norm(M=1e8)   heuristic-tail norm   rigorous upper (if finite)")
for (t, N), d in sorted(runs.items()):
    name, chi, cK, C1, C0 = FIELD[t]
    ck = [(k, bK(k, chi) * (1 - k / N)) for k in range(1, N + 1)]
    V0 = 1 + sum(abs(c) * (cK + C0 / k + C0) for k, c in ck if k > 1)
    V1 = C1 * sum(abs(c) * (1 / k + 1 / math.sqrt(k)) for k, c in ck if k > 1)
    for q, p in enumerate(ps):
        S6, S7, S8 = d[10 ** 6][q], d[10 ** 7][q], d[10 ** 8][q]
        i1, i2 = S7 - S6, S8 - S7
        rho = i2 / i1 if i1 > 0 else 0
        tail_h = i2 * rho / (1 - rho) if 0 < rho < 1 else float("nan")
        line = f"{name:8s} {N:4d} {p:6.4f}  {S8 ** (1 / p):.7f}   {(S8 + tail_h) ** (1 / p):.7f} (ratio {rho:.3f})"
        if p < 2:
            Tr = 2 ** (p - 1) * (V0 ** p / M + V1 ** p * M ** (p / 2 - 1) / (1 - p / 2))
            line += f"   [{S8 ** (1 / p):.5f}, {(S8 + Tr) ** (1 / p):.5f}]"
        print(line)
