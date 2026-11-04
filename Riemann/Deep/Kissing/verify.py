#!/usr/bin/env python3
"""verify.py -- exact check, Python standard library only, of the kissing bound

        tau_38 >= 624020        (previous record 611522)

The configuration is A. Kravatskiy's 611522-point kissing configuration in R^38 = R^24 + R^14
(github.com/alexlegeartis/KissingNumbers, verifications/improved/dim38-joint-latitudes, commit b814435)
with its equator -- Q. Lin's 6855 lines of the norm-48 shell of the Leech lattice, i.e. 13710 points
(q,0)/sqrt(48) -- replaced by the 13104 lines in equator_lines_13104.txt (26208 points).
Count: 611522 - 2*6855 + 2*13104 = 624020.

Coordinates: the integer Leech lattice L of minimum norm 32 built from the Golay code in golay_basis.txt,
the same lattice and the same coordinates as Kravatskiy's caps and Lin's equator.

Checked here (integers only):
  1. the Golay code (weight enumerator 1 + 759 + 2576 + 759 + 1);
  2. L = {x : common parity e, ((x-e)/2 mod 2) in Golay, sum(x) = 4e mod 8}: the minimal vectors generate a
     sublattice of determinant 2^36 = [Z^24 : L], hence all of L; its Gram matrix is 0 mod 8 with diagonal
     0 mod 16, so every norm in L is a multiple of 16; with Golay minimum weight 8 this gives min(L) = 32
     (lemma in leech_min below); the 196560 vectors of norm 32 are generated;
  3. M = {x in L : F x = 0 mod 3} (F = the 6 x 24 matrix in F3.txt) contains no vector of norm 32;
  4. the 13104 equator vectors are distinct lines, lie in L, have norm 48 and lie in M.
     Hence for distinct lines v, w the vectors v - w and v + w are nonzero in M, so |v -+ w|^2 >= 48 and
     |v.w| <= 24 = |v|^2/2.  For every u of norm 32, v -+ u is nonzero in L, so |v.u| <= 24 as well.
  5. independent brute force: every one of the 85850856 pairs has |v.w| <= 24 (exact big-integer convolution).
Compatibility with the rest of Kravatskiy's configuration (see NOTE.md): an equator point (q,0)/sqrt(48)
against a cap (u, sqrt(2) w)/sqrt(48), u of norm 32: cosine q.u/48 <= 1/2 by step 4; against a latitude point
(body of squared length 1/4): cosine <= 1/2 by Cauchy-Schwarz; against an axis point (0, y): cosine 0.
Kravatskiy's own verifiers (refill/verify.py, refill/verify_coordinates.py), with only the equator file name
and line count changed, also pass on this equator: see replay_kravatskiy.py.
Run: python3 verify.py   (about 1-3 minutes)
"""
import sys, os, time, hashlib
from array import array
HERE = os.path.dirname(os.path.abspath(__file__))
def fail(msg):
    print('FAILED:', msg); sys.exit(1)
def need(c, msg):
    if not c: fail(msg)

t0 = time.time()
# ---- 1. Golay code -------------------------------------------------------------------------------------
rows = [l.strip() for l in open(os.path.join(HERE, 'golay_basis.txt')) if l.strip()]
need(len(rows) == 12 and all(len(r) == 24 and set(r) <= {'0', '1'} for r in rows), 'Golay basis format')
words = [0]
for r in rows:
    b = sum(1 << i for i, c in enumerate(r) if c == '1')
    words += [w ^ b for w in words]
code = set(words)
need(len(code) == 4096, 'Golay rank')
wt = {}
for w in words:
    k = bin(w).count('1'); wt[k] = wt.get(k, 0) + 1
need(wt == {0: 1, 8: 759, 12: 2576, 16: 759, 24: 1}, 'Golay weight enumerator')

def in_leech(x):
    e = x[0] & 1
    if any((c & 1) != e for c in x): return False
    r = 0
    for i, c in enumerate(x):
        if ((c - e) // 2) & 1: r |= 1 << i
    return r in code and sum(x) % 8 == 4 * e

# ---- 2. Leech lattice: minimal vectors, Z-basis, norms mod 16, minimum 32 ------------------------------
mins = []
for a in range(24):
    for b in range(a + 1, 24):
        for s in (1, -1):
            v = [0] * 24; v[a] = 4; v[b] = 4 * s; mins.append(tuple(v))
for w in words:
    if bin(w).count('1') != 8: continue
    p = [i for i in range(24) if w >> i & 1]
    for m in range(128):
        if bin(m).count('1') % 2: continue
        v = [0] * 24; v[p[0]] = 2
        for j in range(7): v[p[j + 1]] = 2 - 4 * ((m >> j) & 1)
        mins.append(tuple(v))
for w in words:
    base = [1 - 2 * ((w >> i) & 1) for i in range(24)]
    for a in range(24):
        v = base[:]; v[a] *= -3
        if v[0] > 0: mins.append(tuple(v))
need(len(mins) == len(set(mins)) == 98280, 'Leech minimal lines')
need(all(sum(c * c for c in v) == 32 and in_leech(v) for v in mins), 'minimal vectors: norm 32 and membership')

def hnf(vecs):
    basis = [None] * 24
    for r in vecs:
        r = list(r)
        for i in range(24):
            if r[i] == 0: continue
            if basis[i] is None:
                if r[i] < 0: r = [-x for x in r]
                basis[i] = r; break
            b = basis[i]
            while r[i] != 0:
                q = b[i] // r[i]; b = [x - q * y for x, y in zip(b, r)]; b, r = r, b
            if b[i] < 0: b = [-x for x in b]
            basis[i] = b
            if not any(r): break
    return basis
B = hnf(mins[::37] + mins[1::101])
need(all(b is not None for b in B), 'minimal vectors span R^24')
det = 1
for i in range(24): det *= B[i][i]
need(det == 2 ** 36, 'minimal vectors generate all of L (det 2^36 = [Z^24 : L])')
need(all(in_leech(b) for b in B), 'basis in L')
for i in range(24):
    for j in range(24):
        g = sum(x * y for x, y in zip(B[i], B[j]))
        need(g % 8 == 0 and (i != j or g % 16 == 0), 'Gram matrix of L: inner products 0 mod 8, norms 0 mod 16')
# Lemma (minimum 32): let x in L, x != 0. Odd x: |x|^2 >= 24 and |x|^2 = 0 mod 16, so >= 32. Even x with
# nonzero Golay residue: at least 8 coordinates = 2 mod 4, so >= 8*4 = 32. Even x with zero residue: all
# coordinates = 0 mod 4; norm 16 would mean a single entry +-4, but then sum(x) = 4 mod 8, excluded.
print('Golay code and Leech lattice: ok (%.1fs)' % (time.time() - t0))

# ---- 3. the sublattice M --------------------------------------------------------------------------------
F = [list(map(int, l.split())) for l in open(os.path.join(HERE, 'F3.txt')) if l.strip()]
need(len(F) == 6 and all(len(r) == 24 for r in F), 'F3.txt format')
def in_M(x):
    return all(sum(a * b for a, b in zip(f, x)) % 3 == 0 for f in F)
need(not any(in_M(u) for u in mins), 'M contains a vector of norm 32')
print('M = {x in L : Fx = 0 mod 3} has no vector of norm 32 (all 196560 checked) (%.1fs)' % (time.time() - t0))

# ---- 4. the equator ------------------------------------------------------------------------------------
path = os.path.join(HERE, 'equator_lines_13104.txt')
E = [tuple(map(int, l.split())) for l in open(path) if l.strip()]
print('equator file sha256', hashlib.sha256(open(path, 'rb').read()).hexdigest())
need(len(E) == 13104 and all(len(v) == 24 for v in E), 'equator count/format')
need(all(sum(c * c for c in v) == 48 for v in E), 'equator norms 48')
need(all(in_leech(v) for v in E), 'equator vectors in L')
need(all(in_M(v) for v in E), 'equator vectors in M')
need(all(next(c for c in v if c) > 0 for v in E), 'canonical sign (first nonzero entry positive)')
need(len(set(E)) == len(E), 'equator lines distinct')
print('13104 equator lines: norm 48, in L, in M, distinct  =>  |v.w| <= 24 for all pairs (%.1fs)' % (time.time() - t0))

# ---- 5. brute force over all pairs (exact) ---------------------------------------------------------------
# lane(24j+23) of A_q * P = sum_i (q_i+8)(v_ji+8) < 2^16, no carries: 16-bit lanes.
n = len(E)
buf = array('H', [0]) * (24 * n)
for j, v in enumerate(E):
    for i, c in enumerate(v): buf[24 * j + i] = c + 8
if sys.byteorder != 'little': buf.byteswap()
P = int.from_bytes(buf.tobytes(), 'little')
S = [sum(v) for v in E]
T = [8 * s + 1536 for s in S]
mx = -10**9; mn = 10**9; nbytes = 2 * (24 * n + 24)
for a, q in enumerate(E):
    A = 0
    for i, c in enumerate(q): A |= (c + 8) << (16 * (23 - i))
    lanes = array('H'); lanes.frombytes((A * P).to_bytes(nbytes, 'little'))
    if sys.byteorder != 'little': lanes.byteswap()
    L = lanes[23:24 * n:24]
    sq = 8 * S[a]
    D = [x - t - sq for x, t in zip(L, T)]
    need(D[a] == 48, 'self product')
    D[a] = 0
    hi = max(D); lo = min(D)
    if hi > mx: mx = hi
    if lo < mn: mn = lo
need(mx <= 24 and mn >= -24, 'pairwise |v.w| <= 24')
print('brute force: %d pairs, max v.w = %d, min v.w = %d (%.1fs)' % (n * (n - 1) // 2, mx, mn, time.time() - t0))

# ---- 6. count -----------------------------------------------------------------------------------------
base, old_eq = 611522, 6855
total = base - 2 * old_eq + 2 * n
print('caps 589560 + equator %d + axes 1596 + latitude 6656 = %d' % (2 * n, 589560 + 2 * n + 1596 + 6656))
need(total == 589560 + 2 * n + 1596 + 6656 == 624020, 'count')
print('VERIFIED: tau_38 >= %d  (Kravatskiy 611522 with the 6855-line equator replaced; +%d)' % (total, total - base))

