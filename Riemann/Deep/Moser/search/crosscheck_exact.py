#!/usr/bin/env python3
"""Cross-check leaf conditions with Keller's own exact preflight (scripts/upper/exact.py,
Geometry.valid), patched to the certificate's offsets H (same side directions T).
Usage: crosscheck_exact.py REPO_SCRIPTS_DIR CERT_DIR"""
import gzip, json, lzma, sys
from fractions import Fraction as Q
sys.path.insert(0, sys.argv[1] + "/upper")
sys.path.insert(0, sys.argv[1])
import exact  # pyright: ignore[reportMissingImports]
cert = sys.argv[2]
h = json.load(open(cert + "/header.json"))
assert [Q(t) for t in h["t"]] == exact.T, "cross-check needs Keller's side directions"
exact.H[:] = [Q(x) for x in h["h"]]
g = exact.Geometry(h["specs"])
assert (g.kR, g.kL) == (h["kR"], h["kL"])
tree = json.load(gzip.open(cert + "/tree.json.gz", "rt"))
info = {}
sys.setrecursionlimit(1000000)


def walk(n, assigned, rect):
    if n[0] == "l":
        info[n[1]] = (assigned, rect)
    elif n[0] == "p":
        assert len(n) == 4
        for r in (0, 1):
            walk(n[2 + r], assigned | {(n[1], r)}, rect)
    else:
        for cell, c in n[1]:
            walk(c, assigned, cell)


walk(tree, frozenset(), (0, g.kR, 0, g.kL))
man = json.load(open(cert + "/manifest.json"))
floor = g.n ** 2
count = 0
for name in man["files"]:
    with lzma.open(cert + "/" + name, "rt") as f:
        for line in f:
            idx, den, dk, vals = json.loads(line)
            keys, acc = [], 0
            for d in dk:
                acc += d
                keys.append(acc)
            assigned, rect = info.pop(idx)
            for k in keys:
                if k > floor:
                    m, r = divmod(k - floor - 1, 2)
                    assert (m, r) in assigned, "escape row not on branch"
            rows = g.integer_rows(keys)
            assert g.valid(vals, rows, rect), "leaf %d fails Keller's exact preflight" % idx
            count += 1
assert not info, "leaves without multipliers"
print("CROSSCHECK OK: %d leaves pass Keller's exact.Geometry.valid with H = %s" % (count, h["h"]))

