#!/usr/bin/env python3
"""odlyzko31 check: fields of minimal |d_K| per (degree, signature) and minimal abelian fields (LMFDB nf_fields API,
data/lmfdb_min.json, fetched 6 Nov 2026) must satisfy the certified bounds: rd_K > bound(n, r1) for 31/32 (every K),
and for abelian K also the B bound. Prints the margin rd_K / bound."""
import json, math
d = json.load(open('certified.json'))
R = {(r['w'], r['n'], r['r1']): r['cert_rd'] for r in d['results'] if r['kind'] == 'table'}
fields = json.load(open('data/lmfdb_min.json'))
bad = 0; seen = set()
for deg, r2, ab, label, disc in sorted(fields, key=lambda t: (t[0], t[1], t[2])):
    if label is None or (label, ab) in seen: continue
    seen.add((label, ab))
    r1 = deg - 2 * r2
    rd = disc ** (1.0 / deg)
    ct = R.get(('ct', deg, r1)); cl = R.get(('classical', deg, r1)); abb = R.get(('ab', deg, r1))
    ok = rd > ct and (not ab or rd > abb)
    bad += 0 if ok else 1
    print(f"{label:22s} n={deg:2d} r1={r1:2d} {'abelian ' if ab else 'min     '} rd={rd:9.4f}  31/32 bound {ct:7.3f} (margin {rd/ct:.4f})"
          + (f"  B bound {abb:7.3f} (margin {rd/abb:.4f})" if ab else '') + f"  classical {cl:7.3f}  {'OK' if ok else 'VIOLATION'}")
print('known-field check:', 'ALL OK' if bad == 0 else f'{bad} VIOLATIONS')
