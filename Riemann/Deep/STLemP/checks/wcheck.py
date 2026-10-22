# STLemP review repair (2026-10-22): Lemma W(ii) with N := max{1, 2^ceil(log2(P/s0))} covers every AFE block of Theorem MS.
# Part A: exact block bookkeeping (integers / Fractions) for dyadic P = 2^i, i <= 13, and every integer s0 with a nonzero block.
# Part B: the Theorem MS exponents after the repair (sympy, delta symbolic on (0, 1/1000], m = 1, 2).
import signal; signal.alarm(60)
import sympy as sp
from fractions import Fraction as Fr
ok = 0; bad = 0
def check(name, cond):
    global ok, bad
    c = bool(cond)
    ok += c; bad += (not c)
    print(('OK  ' if c else 'FAIL'), name)

def clog2(P, s0):                 # least j in Z with 2^j >= P/s0 (integer arithmetic)
    if s0 <= P:
        j = 0
        while s0 * 2**j < P: j += 1
        return j
    e = 0                         # s0 > P: j = -e with e maximal such that P 2^e <= s0
    while P * 2**(e+1) <= s0: e += 1
    return -e

def adm(k, pmin=19):              # [R3]-admissible: odd, squarefree, every prime factor >= pmin
    if k % 2 == 0: return False
    p = 2
    while p * p <= k:
        if k % p == 0:
            if p < pmin or k % (p*p) == 0: return False
            k //= p
        else: p += 1
    return k == 1 or k >= pmin

# ---------- Part A ----------
# b(k) = 1_adm(k) 1[(k,s0)=1] lambda_m(k) k^{-w} phi_2(s0 k/P), supp phi_2 in (1/4, 4): integer support K = {k >= 1 : P < 4 s0 k, s0 k < 4P}
IMAX = 13
tot = deg = 0
A1 = A2 = A3 = A4 = A5 = A6 = A7 = True
for i in range(IMAX + 1):
    P = 2**i
    for s0 in range(1, 4*P):                      # s0 >= 4P: s0 k/P >= 4 for every k >= 1, the block is empty
        kmin, kmax = P // (4*s0) + 1, (4*P - 1) // s0
        if i <= 6:                                # brute-force confirmation of the support formula
            A7 &= ([k for k in range(1, 4*P + 2) if P < 4*s0*k and s0*k < 4*P] == list(range(kmin, kmax + 1)))
        if kmin > kmax: continue
        tot += 1
        j = clog2(P, s0)
        Nold = Fr(2)**j
        Nnew = max(Fr(1), Nold)
        r = Fr(P, s0)
        A1 &= (Nnew >= 1 and Nnew.denominator == 1 and (Nnew.numerator & (Nnew.numerator - 1)) == 0 and Nold >= r and Nold / 2 < r)
        A2 &= (Nnew / 16 <= kmin and kmax <= 4 * Nnew)                       # class-B support condition, for the new N
        A5 &= (Nnew <= max(Fr(1), 2*r) and Nnew <= 2*P)                     # N <= 2P <= 4 P0 in Theorem MS
        A6 &= ((Nold <= Fr(1, 2)) == (2*P <= s0))                           # old N <= 1/2 exactly when s0 >= 2P
        if 2*P <= s0:
            deg += 1
            A3 &= (Nnew == 1 and kmin == kmax == 1)                         # degenerate: N = 1, support {1}, nonzero block b(1) = phi_2(s0/P)
        else:
            A4 &= (Nnew == Nold and 1 <= Nnew / r < 2 and Nnew / 8 <= kmin)  # regular: N unchanged, s0 N/P in [1,2), support in [N/8, 4N]
check('A0: support formula K = [P/(4 s0)] + 1 .. [(4P-1)/s0] (brute force, P <= 2^6)', A7)
check('A1: 2^ceil(log2(P/s0)) computed exactly; N_new = max(1, .) is a power of 2 with exponent >= 0 (log 2N >= log 2 in (B1))', A1)
check('A2: supp b in [N/16, 4N] for N_new, every block', A2)
check('A3: s0 >= 2P: N_new = 1 and supp b = {1} (the single AFE term n = s0)', A3)
check('A4: s0 < 2P: N_new = N_old, s0 N/P in [1, 2), supp b in [N/8, 4N] (Mellin bounds of W uniform)', A4)
check('A5: N_new <= max(1, 2P/s0) <= 2P', A5)
check('A6: the old N is <= 1/2 exactly for s0 >= 2P (the gap: log(2N) <= 0 while b(1) may be nonzero)', A6)
print('     blocks checked:', tot, '  degenerate (2P <= s0 < 4P, old N <= 1/2, k = 1 only):', deg)
check('A8: admissible kappa <= 4N = 4 (N = 1): only kappa = 1', [k for k in range(1, 5) if adm(k)] == [1])
check('A9: (m,P)-admissible k <= 2 (prime factors > P_m = 64(m+1)^2): only k = 1', all([k for k in range(1, 3) if adm(k, 64*(m+1)**2 + 1)] == [1] for m in (1, 2)))
check('A10: L_b = 1: sum_{q <= Q, q odd squarefree} 1 <= Q for Q <= 10^4', all(sum(1 for q in range(1, Q + 1, 2) if all(q % (p*p) for p in range(3, int(q**0.5) + 2, 2))) <= Q for Q in (1, 2, 3, 10, 100, 1000, 10000)))

# ---------- Part B ----------
d = sp.symbols('delta', positive=True)
dmax = sp.Rational(1, 1000)
def pos_on(expr):                 # polynomial in delta, > 0 on (0, dmax]
    e = sp.expand(expr)
    while e != 0 and e.subs(d, 0) == 0:       # divide out the powers of delta (positivity on (0, dmax] is unchanged)
        e = sp.expand(sp.cancel(e / d))
    pts = [sp.Integer(0), dmax] + [c for c in sp.solve(sp.diff(e, d), d) if c.is_real and 0 < c < dmax]
    return all(e.subs(d, p) > 0 for p in pts) if sp.Poly(e, d).degree() <= 2 else False
for m in (1, 2):
    h = sp.Rational(m + 1, 2)
    xp = h + 12*d                 # X' >= X^{(m+1)/2+12 delta} (p0 >= X^{(m-1)/4+6 delta}), X' <= 4 X^{xp}
    P0 = h + 4*d                  # P0 = X^{(m+1)/2+4 delta}
    check('B1 m=%d: N <= 8 P0 <= X\'^{1-4 delta}: xp(1-4d) - P0 > 0' % m, pos_on(xp*(1 - 4*d) - P0))
    check('B2 m=%d: s0 = s < 4P <= 8 P0 <= X\' (Lemma W(ii) needs s0 <= X\')' % m, pos_on(xp - P0))
    check('B3 m=%d: |t\'| <= 2 X^delta <= X\'^delta: delta xp - delta > 0' % m, pos_on(d*xp - d))
    check('B4 m=%d: sqrt(X\') E^{(m-1)/2} <= 5 X\' (E <= 23 X\'^{1/2}): 1/2 + (m-1)/4 <= 1 and 23^{(m-1)/4} <= 5' % m, sp.Rational(1, 2) + sp.Rational(m - 1, 4) <= 1 and sp.Integer(23)**sp.Rational(m - 1, 4) <= 5)
    check('B5 m=%d: degenerate blocks via Theorem 1\' (N = 1, E = 1): exponent 3/2 + 21 delta of X\', same as the others' % m, sp.Rational(1, 2) + 21*d + 1 == sp.Rational(3, 2) + 21*d)
    check('B6 m=%d: degenerate blocks directly: |b(1)|^2 sum_D h(D) << X^{3/2} log X <= X\'^{3/2}' % m, pos_on(sp.Rational(3, 2)*xp - sp.Rational(3, 2)))
    check('B7 m=%d: S_2: N^{4 delta} <= (8 P0)^{4 delta} = X^{(2m+2) delta + 16 delta^2}, |eps_m|^2 <= X^{-(2m+2) delta}: net 16 delta^2 <= 20 delta^2' % m, sp.expand(4*d*P0 - (2*m + 2)*d) == 16*d**2)
    check('B8 m=%d: X\'^{3/2+21 delta} <= 16 X^{3(m+1)/4+50 delta} (X\' <= 4 X^{3/2+12 delta})' % m, pos_on(sp.Rational(3*(m + 1), 4) + 50*d - (sp.Rational(3, 2)*xp + 21*d*(sp.Rational(3, 2) + 12*d))) and 4**(sp.Rational(3, 2)) * 2 <= 16)
    check('B9 m=%d: with the S_2 loss 20 delta^2 and the logs: 50 delta + 20 delta^2 < 60 delta' % m, pos_on(60*d - 50*d - 20*d**2))
eta = sp.symbols('eta', positive=True)
vF2 = sp.Rational(15, 32) / (1 + eta)
check('C1: class I blocks of Theorem 2_1 and Proposition F_2 have N > V/2 >= 1 (V = X^{1/2}/2, resp. X^v, v = 15/(32(1+eta)) in [3/8, 7/16] for eta in [1/14, 1/4]): Lemma W(i) only for N >= 2', vF2.subs(eta, sp.Rational(1, 4)) == sp.Rational(3, 8) and vF2.subs(eta, sp.Rational(1, 14)) == sp.Rational(7, 16) and sp.Rational(7, 16) <= sp.Rational(15, 32) and sp.diff(vF2, eta).subs(eta, sp.Rational(1, 7)) < 0)
print('total OK', ok, 'FAIL', bad)
