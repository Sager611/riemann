#!/usr/bin/env python3
"""Certificate: unit-generated real quadratic orders O = O_Delta, Delta = n^2 - 4N (N = +-1), with Cl(O) = Cl(O)[2]
(Kopp-Lagarias, arXiv 2512.11311, Problem 1.7).  Run:
  uv run --no-project --with numpy --with python-flint python3 cert.py
Part 1  D0 (python-flint arb, rigorous): every 2-torsion order has Delta < D0   [Littlewood78 + class number formula + genus theory]
Part 2  sieve (exact integers, numpy int64): every n with Delta_n < D0, n > NEX, must have no split prime p with 4 p^4 < Delta [Lemma A / Cor. B]
Part 3  exact class groups (reduced indefinite forms, cycles, composition) for all n <= NEX and all sieve survivors
Part 4  final lists; comparison with Kopp-Lagarias Tables 2A/2B (transcribed in KL_TABLES below)
See DERIVATION.md for the mathematics.  All decisions use exact integer arithmetic or arb balls."""
import sys, os, time, math
from math import isqrt, gcd
import numpy as np
import flint
from flint import arb, fmpz
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from forms import class_data, group_structure_2

T0 = time.time()
LOG = []
def out(*a):
    s = ' '.join(str(x) for x in a)
    print(s, flush=True)
    LOG.append(s)

def primes_upto(m):
    s = np.ones(m + 1, dtype=bool); s[:2] = False
    for i in range(2, isqrt(m) + 1):
        if s[i]: s[i * i::i] = False
    return [int(p) for p in np.nonzero(s)[0]]

# ---------------------------------------------------------------- Part 1: D0
flint.ctx.prec = 200
THREE = arb(3)
def g_lower(D, c):
    """rigorous ball for g_c(D) = c sqrt(D) / (log(D+4) * max(3, log log D)); callers use its lower endpoint."""
    D = arb(D)
    ll = D.log().log()
    if ll > THREE: m = ll
    elif ll < THREE: m = THREE
    else: m = arb(ll.upper()) if ll.upper() > 3 else THREE      # ball straddles 3: larger denominator
    return c * D.sqrt() / ((D + 4).log() * m)

PR = primes_upto(200000)
def P_k(k):                     # prod_{i<=k} p_i/(p_i-1), exact rational -> arb
    num = 1; den = 1
    for p in PR[:k]:
        num *= p; den *= (p - 1)
    return arb(num) / arb(den)
def primorial(k):
    r = 1
    for p in PR[:k]: r *= p
    return r

def Dplus(T, c):
    """integer D >= 16 with lower(g_c(D)) > upper(T), smallest up to bisection; g_c is increasing on [16, oo) (DERIVATION 3.3)."""
    lo, hi = 16, 16
    while not (g_lower(hi, c).lower() > T.upper()):
        hi *= 2
    if g_lower(lo, c).lower() > T.upper():
        return lo
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if g_lower(mid, c).lower() > T.upper(): hi = mid
        else: lo = mid
    return hi

KMAX = 60
STEP_K = 15
def D0_of(c, verbose=False):
    """D0(c): if |L(1,chi)| >= c/max(3, log log q) for every primitive chi, every 2-torsion order has Delta < D0(c).
    k = omega(Delta); a 2-torsion order needs g_c(Delta) < 2^k P_k and Delta >= primorial(k) (DERIVATION 3.2-3.4)."""
    D0 = 16
    rows = []
    for k in range(0, KMAX + 1):
        T = arb(2) ** k * P_k(k)
        Dk = Dplus(T, c)
        Pi = primorial(k)
        excluded = Pi >= Dk
        rows.append((k, Dk, Pi, excluded))
        if verbose and (k <= 16 or k % 10 == 0):
            out(f'k={k:3d}  2^k P_k={float(T.mid()):.6g}  D_k+={Dk:.6e}  primorial={float(Pi):.6e}  ' + ('excluded (primorial >= D_k+)' if excluded else 'possible: Delta < D_k+'))
    # base of the induction: the smallest K >= STEP_K with every k in [K, KMAX] excluded
    K = None
    for k in range(STEP_K, KMAX + 1):
        if all(r[3] for r in rows[k:]):
            K = k; break
    assert K is not None and K <= KMAX - 5, 'no induction base'
    D0 = max([16] + [r[1] for r in rows[:K] if not r[3]])
    return D0, K

out('== Part 1: D0 ==')
LCONST = arb('0.11')            # Littlewood78: |L(1,chi)| >= 0.11 / max(3, log log q), every primitive chi mod q >= 3
D0, KB = D0_of(LCONST, verbose=True)
# induction step (DERIVATION 3.4): for k >= 15, g(Pi_{k+1})/g(Pi_k) >= sqrt(p) * X/(2X+1) * log X/log(2X+1) with p = p_{k+1} >= 53, X = log Pi_k >= log Pi_15
X = arb(primorial(STEP_K)).log()
step = arb(PR[STEP_K]).sqrt() * (X / (2 * X + 1)) * (X.log() / (2 * X + 1).log())
need = arb(2) * PR[STEP_K] / (PR[STEP_K] - 1)
out(f'induction step constant at k={STEP_K} (p={PR[STEP_K]}): {step.str(5)}  needs > 2p/(p-1) = {need.str(5)}')
assert step.lower() > need.upper()
out(f'induction base k={KB}: every k in [{KB}, {KMAX}] excluded; k > {KMAX} by induction')
out('D0 =', D0, f'= {D0:.6e}', '(every 2-torsion unit-generated order has Delta < D0)')
D_SIEVE = max(D0, 10 ** 18)
# robustness: smallest constant c (3 decimals, rounded up) with D0(c) <= D_SIEVE
lo, hi = 1, 110                  # in units of 1/1000
assert D0_of(arb(hi) / 1000)[0] <= D_SIEVE
while hi - lo > 1:
    mid = (lo + hi) // 2
    try:
        ok = D0_of(arb(mid) / 1000)[0] <= D_SIEVE
    except AssertionError:
        ok = False
    if ok: hi = mid
    else: lo = mid
out(f'robustness: any bound |L(1,chi)| >= c/max(3, log log q) with c >= {hi/1000:.3f} gives D0(c) <= {D_SIEVE:.0e} = sieve range (Littlewood78 has c = 0.11)')

# ---------------------------------------------------------------- Part 2: sieve
NEX = 1000                      # exact class groups for all n <= NEX; sieve for n > NEX (Delta > 10^6 >= 100, Lemma A applies)
out('== Part 2: split-prime sieve, n > %d, Delta_n < D_SIEVE = %.3e ==' % (NEX, D_SIEVE))
SP = [p for p in PR if 4 * p ** 4 < D_SIEVE + 10] + [next(p for p in PR if 4 * p ** 4 >= D_SIEVE + 10)]
TAB = {}
def allowed(p, N):
    key = (p, N)
    if key not in TAB:
        r = np.arange(p, dtype=np.int64)
        v = (r * r - 4 * N) % p
        leg = np.array([0 if int(x) == 0 else (1 if pow(int(x), (p - 1) // 2, p) == 1 else -1) for x in v])
        TAB[key] = (leg != 1)
    return TAB[key]

def kron_split(D, p):
    """True iff the Kronecker symbol (D/p) = 1, i.e. p splits in O_D."""
    if p == 2:
        return D % 8 == 1
    r = D % p
    return r != 0 and pow(r, (p - 1) // 2, p) == 1

def full_ok(n, N):
    D = n * n - 4 * N
    for p in PR:
        if 4 * p ** 4 >= D: return True
        if kron_split(D, p): return False
    raise RuntimeError('prime list too short')

survivors = {}
CH = 4_000_000
for N in (1, -1):
    nmax = isqrt(D_SIEVE + 4 * N)      # n <= nmax covers all Delta_n < D_SIEVE (superset)
    while nmax * nmax - 4 * N < D_SIEVE: nmax += 1
    n0 = NEX + 1
    tot = 0; chs = 0; fulls = []
    while n0 <= nmax:
        n1 = min(nmax + 1, n0 + CH)
        Dmin = n0 * n0 - 4 * N
        ps = [p for p in SP if p > 2 and 4 * p ** 4 < Dmin]     # primes valid for every n in the chunk
        cand = np.arange(n0, n1, dtype=np.int64)
        tot += cand.size
        for p in ps:
            cand = cand[allowed(p, N)[cand % p]]
            if cand.size == 0: break
        chs += cand.size
        for x in cand:
            if full_ok(int(x), N): fulls.append(int(x))
        n0 = n1
    survivors[N] = fulls
    out(f'N={N:+d}: n in [{NEX+1}, {nmax}] ({tot} values): chunk survivors {chs}, full survivors {len(fulls)}, largest {max(fulls) if fulls else None}')
out('sieve time %.1f s' % (time.time() - T0))

# ---------------------------------------------------------------- Part 3: exact class groups
out('== Part 3: exact class groups ==')
maxD = max(max(n * n + 4 for n in survivors[N]) if survivors[N] else 0 for N in (1, -1))
maxD = max(maxD, NEX * NEX + 4)
M = maxD // 4 + 1
spf = np.zeros(M + 1, dtype=np.int32)
for p in primes_upto(isqrt(M) + 1):
    sub = spf[p * p::p]
    sub[sub == 0] = p
    spf[p * p::p] = sub
z = np.nonzero(spf == 0)[0]
spf[z] = z.astype(np.int32)
found = []
survlog = []
nexact = 0
for N in (1, -1):
    ns = [n for n in range(1, NEX + 1) if (n >= 3 if N == 1 else n >= 1)] + survivors[N]
    for n in ns:
        D = n * n - 4 * N
        d = class_data(D, spf)
        nexact += 1
        if n > NEX:
            survlog.append(f'N={N:+d} n={n} Delta={D} hplus={d["hplus"]} h={d["h"]} wide2={d["wide2"]} narrow2={d["narrow2"]}')
        if d['wide2']:
            st = group_structure_2(D, d)
            found.append((D, N, n, d['h'], d['hplus'], d['narrow2'], st))
out(f'exact class groups computed: {nexact} orders; time %.1f s' % (time.time() - T0))
with open(os.path.join(HERE, 'logs', 'survivors.txt'), 'w') as fh:
    fh.write('# sieve survivors n > %d with Delta_n < D_SIEVE = %d (exact class groups); D0 = %d\n' % (NEX, D_SIEVE, D0))
    fh.write('\n'.join(survlog) + '\n')
out('survivors with n > NEX that are 2-torsion:', [x for x in survlog if 'wide2=True' in x])

# ---------------------------------------------------------------- Part 4: lists and comparison
def fcond(D):
    f = 1
    for p, e in __import__('forms').factor(D).items():
        f *= p ** (e // 2)
    D0f = D // (f * f)
    if D0f % 4 not in (0, 1):      # need Delta_0 = 0,1 mod 4
        f //= 2; D0f *= 4
    return f, D0f
def sname(t):
    if t == (): return '1'
    from collections import Counter
    c = Counter(t)
    return ' x '.join((f'(Z/{o})^{m}' if m > 1 else f'Z/{o}') for o, m in sorted(c.items()))
out('== Part 4: wide 2-torsion list (Cl = Cl[2]) ==')
for N in (1, -1):
    L = sorted(x for x in found if x[1] == N)
    out(f'Delta = n^2 {"-" if N == 1 else "+"} 4: {len(L)} orders')
    for D, _, n, h, hp, n2, st in L:
        f, d0 = fcond(D)
        out(f'   n={n:4d} Delta={D:6d} f={f:3d} Delta0={d0:6d} h={h:3d} h+={hp:3d} Cl={sname(st[1])} Cl+={sname(st[0])} one-class-per-genus={n2}')
distinct = sorted(set(x[0] for x in found))
ocpg = sorted(set(x[0] for x in found if x[5]))
out(f'distinct orders with Cl = Cl[2]: {len(distinct)}; with Cl+ = Cl+[2] (one class per genus): {len(ocpg)}; max Delta {max(distinct)}')

KL_A = {  # Table 2A (Delta = n^2 - 4), real orders: Delta -> (Cl, Cl+) as invariant tuples
 5: ((), ()), 12: ((), (2,)), 21: ((), (2,)), 77: ((), (2,)), 437: ((), (2,)), 32: ((), (2,)), 45: ((), (2,)), 117: ((), (2,)),
 **{D: ((2,), (2, 2)) for D in (60, 140, 165, 285, 357, 572, 957, 1085, 2397, 96, 252, 192, 525, 320)},
 **{D: ((2,), (4,)) for D in (221, 1517, 725)},
 **{D: ((2, 2), (2, 2, 2)) for D in (780, 1020, 1365, 1932, 2805, 4485, 5180, 7917, 8645, 480, 672, 1760, 2912, 1440, 2112)},
 **{D: ((2, 2), (2, 4)) for D in (3965, 7565)},
 **{D: ((2, 2, 2), (2, 2, 2, 2)) for D in (4620, 12540, 26565, 3360, 7392, 14880, 19040, 23712, 27552, 6720)},
 68640: ((2, 2, 2, 2), (2, 2, 2, 2, 2)),
}
KL_B = {  # Table 2B (Delta = n^2 + 4): Cl = Cl+
 **{D: () for D in (5, 8, 13, 29, 53, 173, 293, 20, 68, 125)},
 **{D: (2,) for D in (40, 85, 104, 365, 488, 533, 629, 965, 1448, 1685, 1853, 2813, 260, 200, 845)},
 **{D: (2, 2) for D in (680, 1160, 2120, 2405, 3485, 3848, 5480, 10205, 16133)},
 **{D: (2, 2, 2) for D in (8840, 21320, 32045)},
}
ours_A = {x[0]: (x[6][1], x[6][0]) for x in found if x[1] == 1}
ours_B = {x[0]: x[6][0] for x in found if x[1] == -1}
okA = ours_A == KL_A
okB = ours_B == KL_B and all(x[6][0] == x[6][1] for x in found if x[1] == -1)
out(f'comparison with Kopp-Lagarias Table 2A ({len(KL_A)} real entries): {"IDENTICAL" if okA else "DIFFERENT"}; Table 2B ({len(KL_B)} entries): {"IDENTICAL" if okB else "DIFFERENT"}')
if not okA:
    out('   A only ours:', sorted(set(ours_A) - set(KL_A)), 'only KL:', sorted(set(KL_A) - set(ours_A)), 'struct diff:', [D for D in ours_A if D in KL_A and ours_A[D] != KL_A[D]])
if not okB:
    out('   B only ours:', sorted(set(ours_B) - set(KL_B)), 'only KL:', sorted(set(KL_B) - set(ours_B)))
out('total time %.1f s' % (time.time() - T0))
with open(os.path.join(HERE, 'logs', 'cert_output.txt'), 'w') as fh:
    fh.write('\n'.join(LOG) + '\n')
