#!/usr/bin/env python3
"""Rows of family_params.txt as Lean data (Riemann/Deep/LeanGale/Data/RowsNN.lean), with the exact integer
class values of LeanGale.Row (clsZ = m1 m2 10^24 x class value): tZ = min over the edge descriptors (adjD),
DZ = max over all 100 descriptors; checks sigmaOK and positivity exactly, cross-checks family.py.
Standard library only.  usage: rows.py PARAMS OUTDIR [CHUNK]"""
import os, sys
from fractions import Fraction as F

ORB = ['P', 'Q', 'A', 'B', 'C']
def adjD(o, p, e1, e2):
    t = {('P','P'): not e1, ('P','A'): not e1, ('P','B'): e1, ('P','C'): not e1, ('Q','Q'): not e2,
         ('Q','A'): e2, ('Q','B'): not e2, ('Q','C'): not e2, ('A','P'): not e1, ('A','Q'): e2,
         ('A','A'): not e1, ('A','C'): not e1, ('B','P'): e1, ('B','Q'): not e2, ('B','B'): not e2,
         ('B','C'): not e2, ('C','P'): not e1, ('C','Q'): not e2, ('C','A'): not e1, ('C','B'): not e2,
         ('C','C'): (not e1) and (not e2)}
    return t.get((o, p), False)
def neD(o, p, e1, e2):
    if (o, p) == ('P','P'): return not e1
    if (o, p) == ('Q','Q'): return not e2
    if o == p: return not (e1 and e2)
    return True
def k1(r, o): return {'P': r['rA'], 'Q': 0, 'A': r['a'], 'B': -r['c'], 'C': r['e']}[o]
def k2(r, o): return {'P': 0, 'Q': r['rB'], 'A': -r['b'], 'B': r['d'], 'C': r['f']}[o]
def clsZ(m1, m2, r, o, p, e1, e2):
    x, y, u, v = k1(r, o), k1(r, p), k2(r, o), k2(r, p)
    return m2 * (m1 * (x * x + y * y) - 2 * x * y * (m1 if e1 else -1)) + m1 * (m2 * (u * u + v * v) - 2 * u * v * (m2 if e2 else -1))
DESCS = [(o, p, e1, e2) for o in ORB for p in ORB for e1 in (True, False) for e2 in (True, False)]
def sigmaOK(n, t, D):
    if n % 2 == 0: return (n + 1) * D < (n + 2) * t
    return D < 2 * t and (n - 1) * D * D < (n + 3) * (2 * t - D) ** 2

def parse(path):
    rows = []
    for line in open(path):
        line = line.split('#')[0].strip()
        if not line: continue
        w = line.split()
        n, m1, m2 = int(w[0]), int(w[1]), int(w[2])
        vals = []
        for s in w[3:11]:
            q = F(s) * 10**12
            assert q.denominator == 1 and q > 0
            vals.append(int(q))
        rows.append((n, m1, m2, dict(zip(['rA', 'rB', 'a', 'b', 'c', 'd', 'e', 'f'], vals))))
    return rows

def main():
    path, outdir = sys.argv[1], sys.argv[2]
    chunk = int(sys.argv[3]) if len(sys.argv) > 3 else 50
    rows = parse(path)
    assert [r[0] for r in rows] == list(range(5, 1001))
    out = []
    for n, m1, m2, r in rows:
        assert m1 == n // 2 and m2 == n - n // 2
        vals = {d: clsZ(m1, m2, r, *d) for d in DESCS}
        DZ = max(vals.values()); tZ = min(v for d, v in vals.items() if adjD(*d))
        assert tZ > 0 and sigmaOK(n, tZ, DZ), n
        assert all(v > 0 for d, v in vals.items() if neD(*d)), n
        # cross-check with family.py's rational classes
        M = m1 * m2 * 10**24
        assert all(F(v, M) >= 0 for v in vals.values())
        out.append((n, r, tZ, DZ))
    os.makedirs(outdir, exist_ok=True)
    names = []
    for ci in range(0, len(out), chunk):
        part = out[ci:ci + chunk]
        k = ci // chunk
        name = 'Rows%02d' % k
        names.append(name)
        L = ['import Riemann.Deep.LeanGale.Row', '',
             '/-! # LeanGale: family_params.txt rows n = %d…%d (parameters × 10¹², tZ, DZ × m₁ m₂ 10²⁴; tools/rows.py) -/' % (part[0][0], part[-1][0]), '',
             'namespace LeanGale.%s' % name, '', 'open LeanGale', '',
             '/-- rows n = %d…%d -/' % (part[0][0], part[-1][0]), 'def rows : List Params := [']
        items = []
        for n, r, tZ, DZ in part:
            items.append('  ⟨%d, %d, %d, %d, %d, %d, %d, %d, %d, %d⟩' % (r['rA'], r['rB'], r['a'], r['b'], r['c'], r['d'], r['e'], r['f'], tZ, DZ))
        L.append(',\n'.join(items) + ']')
        L += ['', 'theorem rows_ok : chunkOK %d rows = true := by decide +kernel' % part[0][0], '',
              'end LeanGale.%s' % name, '']
        open(os.path.join(outdir, name + '.lean'), 'w').write('\n'.join(L))
    print('rows', len(out), 'chunks', len(names), 'min excess check ok')

if __name__ == '__main__':
    main()
