#!/usr/bin/env python3
"""exp48 cmp.py: compare the scan list with ARSX's data (2605.26020; github.com/erick-ross/research-data,
class-groups-small-exp, Tables 1, 2, 4, 8) restricted to |D| <= X, and report EKN's counts for fundamental D.
Run: python3 cmp.py SCANFILE X     (downloads the four tables to /tmp/exp48_ross if absent; prints their sha256)"""
import sys, re, os, hashlib, urllib.request
scan, X = sys.argv[1], int(sys.argv[2])
cache = '/tmp/exp48_ross'; os.makedirs(cache, exist_ok=True)
URL = 'https://raw.githubusercontent.com/erick-ross/research-data/master/class-groups-small-exp/Table-%d.txt'
ross = {}
for E in (1, 2, 4, 8):
    fn = os.path.join(cache, 'Table-%d.txt' % E)
    if not os.path.exists(fn): urllib.request.urlretrieve(URL % E, fn)
    t = open(fn).read()
    print('Table-%d.txt sha256 %s' % (E, hashlib.sha256(t.encode()).hexdigest()[:16]))
    i = t.index('Discriminant list:'); lst = t[t.index('[', i) + 1: t.index(']', i)]
    ross[E] = sorted(int(x) for x in lst.split(','))
    m = re.search(r'Length = (\d+)', t[i:]); assert int(m.group(1)) == len(ross[E])
mine = {1: [], 2: [], 4: [], 8: []}; fund = {1: [], 2: [], 4: [], 8: []}
for line in open(scan):
    D, f, h, E, n2, n4 = map(int, line.split())
    mine[E].append(D)
    if f: fund[E].append(D)
allok = True
for E in (1, 2, 4, 8):
    r = set(D for D in ross[E] if -D <= X); m = set(mine[E])
    ok = r == m; allok &= ok
    print(f'E={E}: scan {len(m)}, ARSX(<=X) {len(r)}, equal={ok}; only scan: {sorted(m - r)[:5]}; only ARSX: {sorted(r - m)[:5]}; fundamental in scan: {len(fund[E])}, largest |D| {(-min(fund[E])) if fund[E] else None}; largest |D| overall {-min(m) if m else None}')
print('ALL EQUAL' if allok else 'MISMATCH')
