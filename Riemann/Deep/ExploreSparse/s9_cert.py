# explore61 s9: certified sampling defect.  For theta = a/b not of the form 1/k and every N:
#   dist(rho_theta, V_N)^2 >= dist(rho_theta, W1)^2 >= min_beta [ sum_{m<=M} f_m(beta) + beta^2 ]   (drop m > M: f_m >= 0)
# f_m(beta) = |I_m| (q - q^2) [jump of floor(theta/x) at x = theta/n in I_m, fraction q]  - 2 beta * B_m + beta^2 * c_m,
#   B_m = -(log(x/l) - q log(r/l)) (one jump; 0 if none),  c_m = 1 - m(m+1) log^2(1+1/m) = |1/x - avg|^2_{I_m}.
# Minimum over beta of A - 2 beta B + beta^2 (1 + sum c_m) is A - B^2/(1 + sum c_m): exact in Arb ball arithmetic (python-flint).
# usage: uv run --no-project --with python-flint python3 s9_cert.py M
import sys
from flint import arb, ctx
ctx.prec = 120
M = int(sys.argv[1])
def cert(a, b):
    A = arb(0); B = arb(0); Cc = arb(1)
    for m in range(1, M + 1):
        l, r = arb(1) / (m + 1), arb(1) / m
        Cc += 1 - m * (m + 1) * (arb(1) + arb(1) / m).log() ** 2
        # jumps x = (a/b)/n in the open interval (l, r): 1/x = n b / a in (m, m+1)
        for n in range(max(1, (m * a) // b), ((m + 1) * a) // b + 2):
            if m * a < n * b < (m + 1) * a:  # strictly inside
                x = arb(a) / (b * n)
                q = (x - l) / (r - l)
                A += (q - q * q) / (m * (m + 1))
                B += -((x / l).log() - q * (r / l).log())
    return A - B * B / Cc, A, B, Cc
for a, b in [(2, 3), (3, 4), (2, 5), (3, 5), (2, 9)]:
    lo, A, B, Cc = cert(a, b)
    print(f'theta={a}/{b}: certified dist(rho_theta, V_N)^2 >= {lo.lower().str(10, radius=False)} for every N (M={M} intervals; ball {lo.str(12)})', flush=True)
