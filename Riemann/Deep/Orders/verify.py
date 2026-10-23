#!/usr/bin/env python3
"""verify.py - checks the classification of the 24 operation orders
(open question of arXiv 2610.10095, l. 921-925).

Rules (TeX l. 110-129, 151-153): start from the set {2}; move m applies the operation
sigma[(m-1) % 4] to two numbers already made (possibly equal); a difference must be
positive, a quotient a positive integer; the result joins the set.  A sequence is
perfect for N if after N-1 moves every number 1..N is made and every move made a new
number in [1, N].

Checks (all exact):
  A  exhaustive search, all 24 orders, N <= 60, cross-checked by brute force for N <= 14;
  B  the classification table CLAIM agrees with A;
  C  symbolic check of every patch in patches.json (generic in q, P = 4q + r):
     schedule, identities, coverage, operand witnesses, threshold q*;
  D  numeric check of every patch for q* <= q <= q* + 8 (full sequences, rule checker);
  E  the recursive construction (Theorem 2 of NOTE.md) yields a checked perfect
     sequence for every N0 <= N <= 60, for every N <= 400, and for some large N;
  F  the early dead ends used in the impossibility proofs (Proposition 1 of NOTE.md).
Usage: python3 verify.py   (about 6 s, peak memory about 100 MB)
"""
import itertools, json, os, sys, time

sys.setrecursionlimit(10000)
HERE = os.path.dirname(os.path.abspath(__file__))
ORDERS = [''.join(p) for p in itertools.permutations('+x-/')]  # 'x' = times, '/' = divide

# Claimed classification: order -> (finite set of sporadic N, N0 or None)
CLAIM = {
    '+x-/': ([], 9), '+/x-': ([], 9), 'x+-/': ([], 9), '+x/-': ([], 10), 'x+/-': ([], 8),
    '/+x-': ([2, 3, 4], 9), '+/-x': ([4, 6, 10], 12), 'x/+-': ([4, 5, 6, 8, 9, 10], 12),
    'x/-+': ([4, 5, 6], 12), '/+-x': ([2, 3, 4], 12), '/x+-': ([2, 4, 5, 6, 8, 9, 10], 12),
    '/x-+': ([2, 4, 5], None), '/-+x': ([2], None), '/-x+': ([2], None),
}
for o in ORDERS:
    CLAIM.setdefault(o, ([], None))      # the remaining ten orders: never perfect

def claimed(o, N):
    F, N0 = CLAIM[o]
    return N in F or (N0 is not None and N >= N0)

# ---------------------------------------------------------------- rule checker
def check(order, moves, N):
    S = {2}
    for m, (x, op, y, res) in enumerate(moves, 1):
        if op != order[(m - 1) % 4]:
            return 'move %d: operation %s, schedule says %s' % (m, op, order[(m - 1) % 4])
        if x not in S or y not in S:
            return 'move %d: operand not made' % m
        if op == '+': v = x + y
        elif op == 'x': v = x * y
        elif op == '-': v = x - y
        else:
            if y <= 0 or x % y: return 'move %d: inexact division' % m
            v = x // y
        if v != res or not 1 <= v <= N or v in S:
            return 'move %d: result %s wrong, out of range or repeated' % (m, v)
        S.add(v)
    if len(moves) != N - 1 or S != set(range(1, N + 1)):
        return 'final set is not [1,N]'
    return None

# ---------------------------------------------------------------- exact search
def search(order, N):
    """Exhaustive DFS over results; returns a perfect sequence (with operands) or None.
    Pruning (valid): every remaining division must make a new number <= N//2."""
    if N < 2:
        return None
    full = (1 << (N + 1)) - 2
    low = (1 << (N // 2 + 1)) - 2
    divs = [0] * (N + 2)
    for m in range(N - 1, 0, -1):
        divs[m] = divs[m + 1] + (order[(m - 1) % 4] == '/')
    dead = set(); path = []
    def cands(op, S):
        el = [k for k in range(1, N + 1) if S >> k & 1]
        res = 0
        if op == '+':
            for a in el: res |= S << a
        elif op == '-':
            for a in el: res |= S >> a
        elif op == 'x':
            for i, a in enumerate(el):
                for b in el[i:]:
                    if a * b > N: break
                    res |= 1 << (a * b)
        else:
            for b in el:
                if b > 1:
                    for a in el:
                        if a % b == 0: res |= 1 << (a // b)
        res &= full & ~S
        return [k for k in range(N, 0, -1) if res >> k & 1]
    def rec(S, m):
        if S == full: return True
        if S in dead: return False
        if divs[m] > bin(full & ~S & low).count('1'):
            dead.add(S); return False
        for t in cands(order[(m - 1) % 4], S):
            path.append(t)
            if rec(S | 1 << t, m + 1): return True
            path.pop()
        dead.add(S); return False
    if not rec(1 << 2, 1):
        return None
    S = [2]; moves = []
    for m, t in enumerate(path, 1):
        op = order[(m - 1) % 4]
        pair = None
        for x in S:
            y = {'+': t - x, '-': x - t, 'x': t // x if t % x == 0 else 0,
                 '/': x // t if x % t == 0 else 0}[op]
            if y in S:
                pair = (x, y)
                break
        if pair is None:
            raise RuntimeError('no operands for move %d' % m)
        moves.append((pair[0], op, pair[1], t)); S.append(t)
    return moves

def brute(order, N):
    """Plain DFS, no memo, no pruning (cross-check for small N)."""
    def rec(S, m):
        if len(S) == N: return True
        op = order[(m - 1) % 4]; L = sorted(S); new = set()
        for a in L:
            for b in L:
                v = {'+': a + b, 'x': a * b, '-': a - b, '/': a // b if a % b == 0 else 0}[op]
                if 1 <= v <= N and v not in S: new.add(v)
        for v in sorted(new):
            S.add(v)
            if rec(S, m + 1): return True
            S.remove(v)
        return False
    return N >= 2 and rec({2}, 1)

# ---------------------------------------------------------------- symbolic patches
# a value (a, b) means a*q + b; P = 4q + r.  le() is the order for all large q.
def le(u, v): return u[0] < v[0] or (u[0] == v[0] and u[1] <= v[1])
def thr_le(u, v):
    """least q0 >= 0 with u(q) <= v(q) for every q >= q0, or None"""
    (a, b), (a2, b2) = u, v
    if a == a2: return 0 if b <= b2 else None
    if a > a2: return None
    return max(0, -(-(b - b2) // (a2 - a)))
def thr_in(v, lo, hi):
    s, t = thr_le(lo, v), thr_le(v, hi)
    return None if s is None or t is None else max(s, t)

class Shape:
    """Regions of the extension [1,P] -> [1,3P+c] with bulk rounds j = 3..J."""
    def __init__(self, order, r, c, d):
        self.order, self.r, self.c, self.d = order, r, c, d
        self.even = r % 2 == 0
        self.P = (4, r); T = (12, 3 * r); self.T = T
        H = (6, 3 * r // 2) if self.even else (6, (3 * r - 1) // 2)   # H = floor(3P/2)
        i0 = (r - 1) % 4                         # move P has operation order[i0]
        p = 0
        while order[(i0 + p) % 4] != 'x': p += 1
        self.head_ops = [order[(i0 + k) % 4] for k in range(p + 12)]
        ntail = c + 4 * d - p + (0 if self.even else 2)
        ix = order.index('x')
        self.tail_ops = [order[(ix + k) % 4] for k in range(ntail)] if ntail >= 0 else None
        e = 1 if self.even else 3                # middle block = [H-2, H+e+3d]
        self.XB = ((4, r + 1 + d), (6, H[1] - 3))
        self.YB = ((6, H[1] + e + 1 + 3 * d), (12, 3 * r - 8))
        self.Rh = [(12, 3 * r - 6)] if self.even else [(12, 3 * r - 7)]
        self.pool = ([(4, r + i) for i in range(1, d + 1)]
                     + [(6, H[1] + k) for k in range(-2, e + 3 * d + 1)]
                     + [(12, 3 * r + k) for k in range(-7, c + 1)])
        # J >= 2  <=>  P >= 6 + 2d (P even), P >= 7 + 2d (P odd)
        self.qreg = max(0, -(-((6 if self.even else 7) + 2 * d - r) // 4))
    def witness(self, v, done, bulk):
        if v in done: return 0
        best = thr_in(v, (0, 1), self.P)
        if bulk:
            for lo, hi in (self.XB, self.YB):
                t = thr_in(v, lo, hi)
                if t is not None and (best is None or t < best): best = t
        return best

def check_patch(order, r, c, entry):
    """Symbolic proof obligations of one patch; returns q* (raises on failure)."""
    d = entry['d']; sh = Shape(order, r, c, d)
    mv = [(op, (ta, tb), (ua, ub), (wa, wb)) for op, ta, tb, ua, ub, wa, wb in entry['moves']]
    assert sh.tail_ops is not None
    ops = sh.head_ops + sh.tail_ops
    assert [m[0] for m in mv] == ops, 'schedule'
    assert sorted(m[1] for m in mv) == sorted(sh.pool), 'pool coverage'
    nh = len(sh.head_ops); done = set(); q0 = sh.qreg
    for i, (op, t, u, w) in enumerate(mv):
        if op == '+': assert (u[0] + w[0], u[1] + w[1]) == t
        elif op == '-': assert (u[0] - w[0], u[1] - w[1]) == t
        elif op == 'x': assert u[0] == 0 and (u[1] * w[0], u[1] * w[1]) == t
        else: assert w[0] == 0 and w[1] >= 2 and (t[0] * w[1], t[1] * w[1]) == u
        for x in (u, w):
            th = sh.witness(x, done, i >= nh)
            assert th is not None, 'operand %s of move %d has no witness' % (x, i)
            q0 = max(q0, th)
        done.add(t)
        if i == nh - 1: assert set(sh.Rh) <= done, 'head misses R_h'
    if nh == len(mv): assert set(sh.Rh) <= done
    return q0, sh

def bulk(order, P, d):
    """Bulk rounds j = 3..J: x: 3(P-j) = C_j; +: C_j+1; -: C_j-1; /: floor(3P/2)-j."""
    ix = order.index('x'); rot = [order[(ix + k) % 4] for k in range(4)]
    J = P // 2 - 1 - d if P % 2 == 0 else (P - 3) // 2 - d
    out = []
    for j in range(3, J + 1):
        C = 3 * (P - j)
        num = 3 * P - 2 * j if P % 2 == 0 else 3 * P - 1 - 2 * j
        for op in rot:
            out.append({'x': (3, 'x', P - j, C), '+': (C, '+', 1, C + 1),
                        '-': (C, '-', 1, C - 1), '/': (num, '/', 2, num // 2)}[op])
    return out

def patch_moves(order, r, c, entry, q):
    sh = Shape(order, r, c, entry['d']); nh = len(sh.head_ops)
    ev = lambda a, b: a * q + b
    mv = [(ev(ua, ub), op, ev(wa, wb), ev(ta, tb)) for op, ta, tb, ua, ub, wa, wb in entry['moves']]
    return mv[:nh], mv[nh:]

# ---------------------------------------------------------------- construction
DATA = json.load(open(os.path.join(HERE, 'patches.json')))
_cache = {}
def construct(order, N):
    """Perfect sequence for N >= N0: exact search below N1, else N = 3P + c (Theorem 2)."""
    key = (order, N)
    if key in _cache: return _cache[key]
    D = DATA[order]; r = D['r']
    if N < D['N1']:
        res = search(order, N)
        if res is None:
            raise RuntimeError('no perfect sequence for N = %d' % N)
    else:
        c = next(int(k) for k in D['patches'] if (N - int(k)) % 12 == (3 * r) % 12)
        e = D['patches'][str(c)]; P = (N - c) // 3; q = (P - r) // 4
        assert P % 4 == r and q >= e['qstar'] and P >= D['N0']
        hm, tm = patch_moves(order, r, c, e, q)
        res = construct(order, P) + hm + bulk(order, P, e['d']) + tm
    if N <= 5000:
        _cache[key] = res
    return res

def new_after_prefixes(order, nmoves):
    """All sets reachable from {2} by nmoves legal moves, each making a new number
    (no upper bound on the numbers)."""
    sets = [frozenset([2])]
    for m in range(1, nmoves + 1):
        op = order[(m - 1) % 4]; nxt = set()
        for S in sets:
            for a in S:
                for b in S:
                    v = {'+': a + b, 'x': a * b, '-': a - b, '/': a // b if a % b == 0 else 0}[op]
                    if v >= 1 and v not in S: nxt.add(S | {v})
        sets = list(nxt)
    return sets

def main():
    t0 = time.time(); NS = 60; fails = []
    # A + B
    table = {}
    for o in ORDERS:
        perf = []
        for N in range(1, NS + 1):
            s = search(o, N)
            if s is not None:
                if check(o, s, N): fails.append(('A-cert', o, N))
                perf.append(N)
            if N <= 14 and brute(o, N) != (s is not None): fails.append(('A-brute', o, N))
            if claimed(o, N) != (s is not None): fails.append(('B', o, N))
        table[o] = perf
    tA = time.time() - t0
    # C + D
    npatch = 0; qmax = 0
    for o, D in DATA.items():
        r = D['r']
        assert sorted((3 * r + int(c)) % 12 for c in D['patches']) == list(range(12)), 'residues'
        assert (D['N0'] == CLAIM[o][1]) and D['sporadic'] == CLAIM[o][0]
        for c, e in D['patches'].items():
            c = int(c); npatch += 1
            q0, sh = check_patch(o, r, c, e)
            if q0 != e['qstar']: fails.append(('C', o, c))
            qmax = max(qmax, q0)
            for q in range(q0, q0 + 9):
                P = 4 * q + r
                if P < D['N0']: continue
                hm, tm = patch_moves(o, r, c, e, q)
                err = check(o, construct(o, P) + hm + bulk(o, P, e['d']) + tm, 3 * P + c)
                if err: fails.append(('D', o, c, q, err))
    tC = time.time() - t0 - tA
    # E
    big = [1000, 4321, 10000, 12345, 30001]
    for o, D in DATA.items():
        for N in list(range(D['N0'], 401)) + big:
            err = check(o, construct(o, N), N)
            if err: fails.append(('E', o, N, err))
        _cache.clear()
    tE = time.time() - t0 - tA - tC
    # F  dead ends: first operation '-' (no legal move 1); '-' second (no new number at move 2);
    #    order '/x-+': every legal 4-move prefix gives {1,2,3,4,s}, s in 5..8, and move 5 is stuck
    for o in ORDERS:
        if o[0] == '-' and new_after_prefixes(o, 1): fails.append(('F1', o))
        if o[1] == '-' and new_after_prefixes(o, 2): fails.append(('F2', o))
    pre = new_after_prefixes('/x-+', 4)
    if sorted(sorted(S) for S in pre) != [[1, 2, 3, 4, s] for s in range(5, 9)]: fails.append(('F3a',))
    if new_after_prefixes('/x-+', 5): fails.append(('F3b',))
    print('A/B  perfect N <= %d (exhaustive; brute-force cross-check N <= 14):' % NS)
    for o in ORDERS:
        p = table[o]
        if not p: s = 'none'
        elif p[-1] == NS:
            k = len(p) - 1
            while k > 0 and p[k - 1] == p[k] - 1: k -= 1
            s = ', '.join(map(str, p[:k])) + (', ' if k else '') + '%d..%d' % (p[k], NS)
        else: s = ', '.join(map(str, p))
        print('   %s  %s' % (o.replace('x', '*'), s))
    print('C/D  %d patches checked symbolically (max q* = %d) and numerically for q*..q*+8' % (npatch, qmax))
    print('E    constructions checked for N0..400 and N in %s, all 11 orders' % big)
    print('F    dead ends of Proposition 1 confirmed (move 1 / move 2 / move 5 of /x-+)')
    print('N1 (start of the induction): ' + ', '.join('%s %d' % (o.replace('x', '*'), D['N1']) for o, D in DATA.items()))
    print('time A/B %.1fs  C/D %.1fs  E %.1fs' % (tA, tC, tE))
    if fails:
        print('FAILURES:', fails[:20]); sys.exit(1)
    print('ALL CHECKS PASSED')

if __name__ == '__main__':
    main()

