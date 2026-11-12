#!/usr/bin/env python3
"""gl2ceil: independent coverage check of a certificate log (cells tile the core box exactly, pairwise disjoint).
usage: python3 covcheck.py logs/ceil_<tag>.json"""
import json, sys
from fractions import Fraction as F
d = json.load(open(sys.argv[1])); HL = F(d['HL']).limit_denominator(2**30); LAMi = round(F(d['LAM'])/HL)
cells = []
for Rc, dd, ia, ib, how in d['cells']:
    Rc, dd = F(Rc), F(dd); cells.append((Rc - dd, Rc + dd, ia, ib))
R0 = min(c[0] for c in cells); R1 = max(c[1] for c in cells)
area = sum((b - a)*(jb - ja) for a, b, ja, jb in cells)
assert all(0 <= ja < jb <= LAMi for _, _, ja, jb in cells)
ok_area = area == (R1 - R0)*LAMi
cells.sort(); over = 0
for i, (a, b, ja, jb) in enumerate(cells):
    for a2, b2, ja2, jb2 in cells[i+1:]:
        if a2 >= b: break
        if ja2 < jb and ja < jb2: over += 1
print(f"{sys.argv[1]}: {len(cells)} cells, R in [{float(R0)}, {float(R1)}], lam-index in [0, {LAMi}], area match {ok_area}, overlapping pairs {over}, certified flag {d['certified']}")
print("COVER OK" if ok_area and over == 0 and d['certified'] else "COVER FAIL")
