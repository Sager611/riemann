#!/usr/bin/env python3
"""Independent check that the Lean data files equal the source files:
Data/P3, P4, P7, P8 .lean (pts) = the point files x 10^8; Data/RowsNN.lean (rows) = family_params.txt x 10^12,
in order n = 5...1000; D2 = the largest squared distance and dc2 = the stated critical value.
Standard library only.  usage: python3 checkdata.py   (from anywhere)"""
import os, re, glob
from fractions import Fraction
here = os.path.dirname(os.path.abspath(__file__))
lg = os.path.join(here, '..')
deep = os.path.join(here, '..', '..')

def load_pts(path):
    out = []
    for line in open(path):
        line = line.split('#')[0].strip()
        if not line: continue
        out.append([int(Fraction(s) * 10**8) for s in line.replace(',', ' ').split()])
    return out

def lean_def(src, name):
    m = re.search(r'def %s : [^:=]+ :=\s*(.*?)(?=\n\n|\n/--|\nend )' % re.escape(name), src, re.S)
    return m.group(1)

def lean_pts(path):
    src = open(path).read()
    body = lean_def(src, 'pts')
    return [[int(x) for x in row.split(',')] for row in re.findall(r'\[([-0-9, ]+)\]', body)]

ok = True
for name, f, dc2 in [('P3', 'Gale/points.txt', 7932058815501876), ('P4', 'Gale45/points4.txt', 8399124363373578),
                     ('P7', 'Gale7/points7.txt', 8900856812318928), ('P8', 'Gale7/points8.txt', 9023994457279939)]:
    a = load_pts(os.path.join(deep, f)); b = lean_pts(os.path.join(lg, 'Data', name + '.lean'))
    src = open(os.path.join(lg, 'Data', name + '.lean')).read()
    D2 = int(re.search(r'def D2 : ℤ := (\d+)', src).group(1)); dc = int(re.search(r'def dc2 : ℤ := (\d+)', src).group(1))
    trueD2 = max(sum((x - y) ** 2 for x, y in zip(p, q)) for p in a for q in a)
    same = a == b and D2 == trueD2 and dc == dc2
    ok &= same
    print('%s: %d points from %s, coordinates equal: %s, D2 = max squared distance: %s, dc2 = %d: %s' % (name, len(a), f, a == b, D2 == trueD2, dc2, dc == dc2))

rows = []
for line in open(os.path.join(deep, 'Gale7', 'family_params.txt')):
    line = line.split('#')[0].strip()
    if not line: continue
    w = line.split()
    rows.append((int(w[0]), [int(Fraction(s) * 10**12) for s in w[3:11]]))
lean_rows = []
for f in sorted(glob.glob(os.path.join(lg, 'Data', 'Rows*.lean'))):
    src = open(f).read()
    start = int(re.search(r'chunkOK (\d+) rows', src).group(1))
    for i, m in enumerate(re.findall(r'⟨([0-9, ]+)⟩', src)):
        v = [int(x) for x in m.split(',')]
        lean_rows.append((start + i, v[:8]))
same = rows == lean_rows
ok &= same
print('family: %d rows of family_params.txt (n = %d...%d), Lean rows equal (n and the 8 parameters x 10^12): %s' % (len(rows), rows[0][0], rows[-1][0], same))
print('ALL OK' if ok else 'MISMATCH')
