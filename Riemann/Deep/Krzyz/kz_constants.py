#!/usr/bin/env python3
"""Independent Arb checks of the scalar constants used by the proof:
 - T_lo = -log(1-2/e)/2 versus the binary64 literal used in the kernels,
 - the binary64 threshold literal versus 2/e,
 - the Poisson/Cauchy mass cutoffs of Table thi (upper Riemann sums, P_r decreasing on [0,pi]),
 - the q caps from T_n(q) <= n*T_hi + n - 1."""
from flint import arb, ctx
ctx.prec = 96
two_e = 2 / arb(1).exp()
tlo = -(1 - two_e).log() / 2
lit_tlo = arb(0.66544663410202)          # binary64 value of the literal in fixed.hpp
lit_thr = arb(0.73575888234288456)       # binary64 value of the threshold literal
print('T_lo =', tlo, ' literal <= T_lo:', bool(lit_tlo < tlo), ' gap', (tlo - lit_tlo).mid())
print('2/e =', two_e, ' threshold literal < 2/e:', bool(lit_thr < two_e), ' gap', (two_e - lit_thr).mid())

def poisson_bound(n, t, r, M):
    """r^-n (1/2pi) int exp(-t P_r) <= r^-n (1/M) sum_{i=1..M} exp(-t P_r(i pi/M)) (integrand increasing on [0,pi])."""
    r = arb(r); t = arb(t); s = arb(0); pi = arb.pi()
    one_m_r2 = 1 - r * r; one_p_r2 = 1 + r * r
    for i in range(1, M + 1):
        th = pi * i / M
        P = one_m_r2 / (one_p_r2 - 2 * r * th.cos())
        s += (-t * P).exp()
    return s / M / r ** n

rows = [(3, 3.539, 0.8456, 0.734882), (4, 4.689, 0.8790, 0.734844), (5, 5.839, 0.9009, 0.735144),
        (6, 6.989, 0.9161, 0.735487), (7, 8.149, 0.9280, 0.735332)]
for n, t, r, claim in rows:
    b = poisson_bound(n, t, r, 20000)
    print(f'n={n} t={t} r={r}: bound <= {b.upper()}  claimed {claim}  claim_ok={bool(b < arb(claim))}  below_2/e={bool(b < two_e)}')
# The binary64 cutoffs used in code are T_hi = 3.54,4.69,5.84,6.99,8.15; certified point is T_hi - 0.001.
# monotone in t, so it suffices that the certified point is <= the binary64 literal (it is, by 0.001).
for n, T, cap in [(5, 5.84, 1.4), (6, 6.99, 1.31), (7, 8.15, 1.3)]:
    q = arb(cap); a, b = arb(1), q
    for _ in range(n - 1):
        a, b = b, 2 * q * b - a
    print(f'n={n}: T_n({cap}) = {b.mid()} > n*T_hi+n-1 = {n*T+n-1}: {bool(b > arb(n*T+n-1))}')
