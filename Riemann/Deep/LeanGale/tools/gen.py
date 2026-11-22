#!/usr/bin/env python3
"""Lean data for an explicit finite point set (Riemann/Deep/LeanGale/Data/<NAME>.lean): integer coordinates
(x 10^8), D^2 and the critical squared distance dc^2 (x 10^16), a colouring whose same-colour pairs are at
squared distance <= dc^2, and a lower-bound certificate for "every K-colouring has a same-colour pair at
squared distance >= dc^2":
  --cert:    spread masks of G(>= dc^2) (bit 8u of sp[v]: u adjacent to v) and a refutation token list for
             Riemann/Deep/LeanGale/Cert.lean (forward checking, colour-symmetry cuts):
               ('b', v, k) split c v = k | c v != k;  ('u', v, k) dom v = {k}: assign;
               ('c', v, a, b) domains symmetric in a, b: w.l.o.g. c v != b;  ('w', u) dom u empty: close
  --labels:  the vertex labels of a family member (P_i, Q_j, A_ij, B_ij, C_ij), checked against H(p, q).
Standard library only.  usage: gen.py POINTS K DC2 NAME OUT (--cert | --labels M1)"""
import os, sys, math
from fractions import Fraction
from collections import Counter
sys.setrecursionlimit(1000000)

def load(path):
    pts = []
    for line in open(path):
        line = line.split('#')[0].strip()
        if not line: continue
        fr = [Fraction(s) for s in line.replace(',', ' ').split()]
        assert all((x * 10**8).denominator == 1 for x in fr)
        pts.append([int(x * 10**8) for x in fr])
    return pts

def sqmat(P):
    N = len(P)
    return [[sum((a - b) ** 2 for a, b in zip(P[i], P[j])) for j in range(N)] for i in range(N)]

class Refuter:
    def __init__(self, N, K, nbrs):
        self.N, self.K = N, K
        self.nb = [sorted(x) for x in nbrs]
        self.deg = [len(x) for x in self.nb]
    def assign(self, dom, v, k):
        d = list(dom); d[v] = 1 << k
        for u in self.nb[v]: d[u] &= ~(1 << k)
        return d
    def symm(self, dom, a, b):
        return all(((x >> a) & 1) == ((x >> b) & 1) for x in dom)
    def refute(self, dom, assigned, toks):
        for u in range(self.N):
            if dom[u] == 0:
                toks.append(('w', u)); return
        best = None
        for u in range(self.N):
            if assigned[u]: continue
            key = (bin(dom[u]).count('1'), -self.deg[u], u)
            if best is None or key < best: best = key
        if best is None: raise ValueError('a proper colouring exists')
        v = best[2]
        a2 = list(assigned); a2[v] = True
        if dom[v] & (dom[v] - 1) == 0:
            k = dom[v].bit_length() - 1
            toks.append(('u', v, k))
            self.refute(self.assign(dom, v, k), a2, toks); return
        dom = list(dom)
        cols = [k for k in range(self.K) if dom[v] >> k & 1]
        for b in cols:
            for a in cols:
                if a >= b: break
                if (dom[v] >> a) & 1 and (dom[v] >> b) & 1 and self.symm(dom, a, b):
                    toks.append(('c', v, a, b)); dom[v] &= ~(1 << b); break
        self.branch(dom, a2, v, toks)
    def branch(self, dom, a2, v, toks):
        if dom[v] == 0:
            toks.append(('w', v)); return
        if dom[v] & (dom[v] - 1) == 0:
            k = dom[v].bit_length() - 1
            toks.append(('u', v, k))
            self.refute(self.assign(dom, v, k), a2, toks); return
        k = (dom[v] & -dom[v]).bit_length() - 1
        toks.append(('b', v, k))
        self.refute(self.assign(dom, v, k), a2, toks)
        d2 = list(dom); d2[v] &= ~(1 << k)
        self.branch(d2, a2, v, toks)

def refutation(N, K, nbrs):
    R = Refuter(N, K, nbrs)
    toks = []
    R.refute([(1 << K) - 1] * N, [False] * N, toks)
    return toks

def check_swar(N, K, sp, toks):
    """re-implementation of LeanGale.Cert.run (one integer per state, 8 bits per vertex)"""
    F = (1 << (8 * N)) - 1
    s0 = 0
    for u in range(N): s0 |= ((1 << K) - 1) << (8 * u)
    st = [s0]
    for t in toks:
        if not st: return False
        s = st.pop()
        if t[0] == 'w':
            if not (t[1] < N and (s >> (8 * t[1])) % 256 == 0): return False
        elif t[0] == 'b':
            _, v, k = t
            if not (v < N and k < 8): return False
            a = ((s & (F ^ (sp[v] << k))) & (F ^ (255 << (8 * v)))) | (1 << (8 * v + k))
            st.append(s ^ (s & (1 << (8 * v + k)))); st.append(a)
        elif t[0] == 'u':
            _, v, k = t
            if not (v < N and k < 8 and (s >> (8 * v)) % 256 == 1 << k): return False
            st.append(((s & (F ^ (sp[v] << k))) & (F ^ (255 << (8 * v)))) | (1 << (8 * v + k)))
        else:
            _, v, a, b = t
            if not (v < N and a != b and a < K and b < K): return False
            if not all(((s >> (8 * u + a)) & 1) == ((s >> (8 * u + b)) & 1) for u in range(N)): return False
            st.append(s ^ (s & (1 << (8 * v + b))))
    return not st

def colouring(N, K, nbrs, budget=10**7):
    nb = [sorted(x) for x in nbrs]
    col = [-1] * N
    cnt = [0]
    def rec():
        cnt[0] += 1
        if cnt[0] > budget: return None
        best = None
        for u in range(N):
            if col[u] >= 0: continue
            used = {col[w] for w in nb[u] if col[w] >= 0}
            key = (-len(used), -len(nb[u]))
            if best is None or key < best[0]: best = (key, u, used)
        if best is None: return True
        _, u, used = best
        for c in range(K):
            if c in used: continue
            col[u] = c
            r = rec()
            if r: return True
            if r is None: return None
            col[u] = -1
        return False
    return col if rec() else None

def labels(P, m1):
    """labels (orbit, i, j) of a family member: P_i (second block 0), Q_j (first block 0), A/B/C by the
    signs of the projections on the simplex directions"""
    def unit(v):
        nr = math.sqrt(sum(c * c for c in v)); return [c / nr for c in v]
    Pi = [k for k, s in enumerate(P) if all(c == 0 for c in s[m1:])]
    Qj = [k for k, s in enumerate(P) if all(c == 0 for c in s[:m1])]
    U = [unit(P[k][:m1]) for k in Pi]; W = [unit(P[k][m1:]) for k in Qj]
    lab = {}
    for i, k in enumerate(Pi): lab[k] = ('P', i, 0)
    for j, k in enumerate(Qj): lab[k] = ('Q', 0, j)
    for k, s in enumerate(P):
        if k in lab: continue
        x = unit(s[:m1]); y = unit(s[m1:])
        di = [sum(a * b for a, b in zip(x, u)) for u in U]; dj = [sum(a * b for a, b in zip(y, w)) for w in W]
        i = max(range(len(U)), key=lambda t: abs(di[t])); j = max(range(len(W)), key=lambda t: abs(dj[t]))
        si = di[i] > 0; sj = dj[j] > 0
        lab[k] = ('A' if si and not sj else 'B' if sj and not si else 'C', i, j)
    assert len(set(lab.values())) == len(P)
    return lab, len(Pi), len(Qj)

EDGES = [('P','P','n','-'), ('Q','Q','-','n'), ('P','A','n','-'), ('Q','B','-','n'), ('A','A','n','e'), ('A','A','n','n'),
         ('B','B','e','n'), ('B','B','n','n'), ('P','C','n','-'), ('A','C','n','e'), ('A','C','n','n'), ('Q','C','-','n'),
         ('B','C','e','n'), ('B','C','n','n'), ('C','C','n','n'), ('P','B','e','-'), ('Q','A','-','e')]
def isH(x, y):
    ORB = 'PQABC'
    if ORB.index(x[0]) > ORB.index(y[0]): x, y = y, x
    r1 = '-' if 'Q' in (x[0], y[0]) else ('e' if x[1] == y[1] else 'n')
    r2 = '-' if 'P' in (x[0], y[0]) else ('e' if x[2] == y[2] else 'n')
    return (x[0], y[0], r1, r2) in EDGES

def lean_list(items, per_line=10, indent='  '):
    lines = []
    for i in range(0, len(items), per_line):
        lines.append(indent + ', '.join(items[i:i + per_line]))
    return '[\n' + ',\n'.join(lines) + ']'

def lean_tok(t):
    if t[0] == 'b': return '.br %d %d' % (t[1], t[2])
    if t[0] == 'u': return '.unit %d %d' % (t[1], t[2])
    if t[0] == 'c': return '.cut %d %d %d' % (t[1], t[2], t[3])
    return '.wipe %d' % t[1]

def main():
    path, K, dc2, name, out = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], sys.argv[5]
    P = load(path); N = len(P); n = len(P[0]); sq = sqmat(P)
    D2 = max(max(r) for r in sq)
    iD = next((i, j) for i in range(N) for j in range(N) if sq[i][j] == D2)
    assert all(sq[i][j] > 0 for i in range(N) for j in range(N) if i != j)
    nb_ge = [set(j for j in range(N) if j != i and sq[i][j] >= dc2) for i in range(N)]
    nb_gt = [set(j for j in range(N) if j != i and sq[i][j] > dc2) for i in range(N)]
    col = colouring(N, K, nb_gt)
    assert col and all(sq[i][j] <= dc2 for i in range(N) for j in range(N) if col[i] == col[j])
    L = ['import Riemann.Deep.LeanGale.Cert', '',
         '/-!', '# LeanGale: data for %s (generated by tools/gen.py)' % os.path.basename(path), '',
         '%d points in ℝ^%d (coordinates × 10⁸), D² and the critical squared distance dc² (both × 10¹⁶), a %d-colouring' % (N, n, K),
         'whose same-colour pairs have squared distance ≤ dc², and the data of the lower bound (see the end).', '-/', '',
         'namespace LeanGale.%s' % name, '', 'open LeanGale.Cert', '',
         '/-- the points, coordinates × 10⁸ -/',
         'def pts : List (List ℤ) := ' + lean_list(['[' + ', '.join(str(x) for x in p) + ']' for p in P], per_line=1), '',
         '/-- D² × 10¹⁶: the largest squared distance, attained by the points %d and %d -/' % iD,
         'def D2 : ℤ := %d' % D2, '', 'def iD : ℕ := %d' % iD[0], '', 'def jD : ℕ := %d' % iD[1], '',
         '/-- dc² × 10¹⁶: the critical squared distance -/', 'def dc2 : ℤ := %d' % dc2, '',
         '/-- a colouring whose same-colour pairs are at squared distance ≤ dc² -/',
         'def col : List ℕ := ' + lean_list([str(c) for c in col], per_line=30)]
    if '--cert' in sys.argv:
        toks = refutation(N, K, nb_ge)
        sp = [sum(1 << (8 * u) for u in nb_ge[v]) for v in range(N)]
        assert check_swar(N, K, sp, toks)
        print('tokens', len(toks), Counter(t[0] for t in toks))
        L += ['', '/-- spread masks of G(≥ dc²): bit 8u of entry v set iff the points u, v are at squared distance ≥ dc² -/',
              'def sp : List ℕ := ' + lean_list([str(m) for m in sp], per_line=1)]
        chunk = 1000; nch = 0
        for i in range(0, len(toks), chunk):
            L += ['', '/-- refutation certificate, part %d -/' % nch,
                  'def toks%d : List Tok := ' % nch + lean_list([lean_tok(t) for t in toks[i:i + chunk]], per_line=12)]
            nch += 1
        L += ['', '/-- the refutation certificate (%d tokens) -/' % len(toks), 'def toks : List Tok :=',
              '  ' + ' ++ '.join('toks%d' % i for i in range(nch))]
    if '--labels' in sys.argv:
        m1 = int(sys.argv[sys.argv.index('--labels') + 1])
        lab, p, q = labels(P, m1)
        inv = {v: k for k, v in lab.items()}
        for x in lab.values():
            for y in lab.values():
                if isH(x, y): assert sq[inv[x]][inv[y]] >= dc2, (x, y)
        print('labels: H(%d,%d) maps into G(>= dc2); min over H edges = dc2:' % (p, q),
              min(sq[inv[x]][inv[y]] for x in lab.values() for y in lab.values() if isH(x, y)) == dc2)
        def row(o):
            return [str(inv[(o, i, j)]) for i in range(p) for j in range(q)]
        L += ['', '/-- the indices of P_i, Q_j (H(%d, %d) labels) -/' % (p, q),
              'def idxP : List ℕ := [%s]' % ', '.join(str(inv[('P', i, 0)]) for i in range(p)),
              'def idxQ : List ℕ := [%s]' % ', '.join(str(inv[('Q', 0, j)]) for j in range(q)),
              '/-- the indices of A_ij, B_ij, C_ij (row-major, entry %d·i + j) -/' % q,
              'def idxA : List ℕ := [%s]' % ', '.join(row('A')),
              'def idxB : List ℕ := [%s]' % ', '.join(row('B')),
              'def idxC : List ℕ := [%s]' % ', '.join(row('C'))]
    L += ['', 'end LeanGale.%s' % name, '']
    open(out, 'w').write('\n'.join(L))
    print('wrote', out, 'N', N, 'n', n, 'K', K, 'D2', D2, 'dc2', dc2)

if __name__ == '__main__':
    main()
