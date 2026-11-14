# E1Audit: remaining rows of perterm_check.py (its 1200 s cap hit during T = 1e7 under background QoS): d = 0 exactly and the N0 < 1 regime.
import math
from perterm_check import row, Tj
print("d = 0 exactly (alpha=0.01+0.02i, beta=-0.01-0.02i):", flush=True)
for T, th, pq, j in [(1e5, 0.505, 30, 0), (1e5, 0.505, 30, 1), (1e6, 0.55, 30, 0)]:
    row(T, th, pq, j, complex(0.01, 0.02), complex(-0.01, -0.02))
print("N0 < 1 regime (p=q=1), complex shifts alpha=0.05+0.1i, beta=-0.03+0.07i:", flush=True)
al, be = complex(0.05, 0.1), complex(-0.03, 0.07)
for T, th in [(1e6, 0.55), (1e6, 0.6), (1e6, 0.75)]:
    H = T ** th; N0 = T / H ** 2; v = Tj(T, H, 1, [1, 2, 5], 0, al, be)
    bd = H * (2 * math.sqrt(2 * math.pi) * math.sqrt(T)) ** -0.5
    print(f"  T={T:.0e} th={th} N0={N0:.3g}: |T(n)|/(H X0(1)^-1/2) at n=1,2,5: {abs(v[0])/bd:.2e}, {abs(v[1])/bd:.2e}, {abs(v[2])/bd:.2e}", flush=True)
