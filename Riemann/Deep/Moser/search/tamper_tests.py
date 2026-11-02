#!/usr/bin/env python3
"""Negative tests: tampered certificates must be rejected by verify_upper.py."""
import gzip, hashlib, json, lzma, shutil, subprocess, sys
from fractions import Fraction
from pathlib import Path

SRC = Path(sys.argv[1])
VER = Path(sys.argv[2])
WORK = Path(sys.argv[3])


def clone(name):
    d = WORK / name
    if d.exists():
        shutil.rmtree(d)
    shutil.copytree(SRC, d)
    return d


def rehash(d):
    man = json.loads((d / "manifest.json").read_text())
    man["sha256"] = {n: hashlib.sha256((d / n).read_bytes()).hexdigest() for n in man["sha256"]}
    (d / "manifest.json").write_text(json.dumps(man, sort_keys=True, indent=1) + "\n")


def edit_header(d, fn):
    h = json.loads((d / "header.json").read_text())
    fn(h)
    (d / "header.json").write_text(json.dumps(h, sort_keys=True, indent=1) + "\n")


def edit_tree(d, fn):
    t = json.load(gzip.open(d / "tree.json.gz", "rt"))
    fn(t)
    with gzip.open(d / "tree.json.gz", "wt") as f:
        json.dump(t, f)


def edit_first_leaf(d, fn):
    man = json.loads((d / "manifest.json").read_text())
    name = man["files"][0]
    lines = lzma.decompress((d / name).read_bytes()).decode().splitlines()
    rec = json.loads(lines[0])
    fn(rec)
    lines[0] = json.dumps(rec, separators=(",", ":"))
    (d / name).write_bytes(lzma.compress(("\n".join(lines) + "\n").encode()))


def first(t, kind):
    stack = [t]
    while stack:
        n = stack.pop()
        if n[0] == kind:
            return n
        if n[0] == "p":
            stack.extend(n[2:])
        elif n[0] == "s":
            stack.extend(c for _, c in n[1])
    raise ValueError("no %r node in tree" % kind)


def run(d, timeout=900):
    p = subprocess.run([sys.executable, str(VER), str(d)], capture_output=True, text=True, timeout=timeout)
    return p.returncode, (p.stdout + p.stderr).strip().splitlines()[-1]


tests = {}
def shrink_h2(h):
    h["h"][2] = str(Fraction(h["h"][2]) * 99 / 100)


d = clone("h_smaller"); edit_header(d, shrink_h2); rehash(d); tests["h2 * 0.99"] = d
d = clone("t_changed"); edit_header(d, lambda h: h["t"].__setitem__(2, "7/17")); rehash(d); tests["t2 = 7/17"] = d
d = clone("p_child"); edit_tree(d, lambda t: first(t, "p").pop()); rehash(d); tests["drop escape branch"] = d
d = clone("s_cover"); edit_tree(d, lambda t: first(t, "s")[1].pop()); rehash(d); tests["drop split cell"] = d


def drop_escapes(rec):
    keys: list[int] = []
    acc = 0
    for x in rec[2]:
        acc += x
        keys.append(acc)
    floor = 63 * 63
    keep = [(k, v) for k, v in zip(keys, rec[3]) if k <= floor]
    ks = [k for k, _ in keep]
    rec[2] = [ks[0]] + [b - a for a, b in zip(ks, ks[1:])]
    rec[3] = [v for _, v in keep]


d = clone("no_escape"); edit_first_leaf(d, drop_escapes); rehash(d); tests["leaf without escape rows"] = d


def foreign_escape(rec):
    rec[2].append(1)  # append key (last key + 1): an escape row of an unassumed placement
    rec[3].append(1)


d = clone("foreign"); edit_first_leaf(d, foreign_escape); rehash(d); tests["unassumed escape row"] = d
d = clone("nohash"); edit_header(d, lambda h: h.__setitem__("kR", h["kR"])); (d / "header.json").write_text((d / "header.json").read_text() + " "); tests["sha mismatch"] = d
ok = True
for name, d in tests.items():
    code, last = run(d)
    good = code != 0 and "REJECT" in last
    ok &= good
    print(("PASS" if good else "FAIL"), name, "->", last)
    shutil.rmtree(d)
print("ALL PASS" if ok else "SOME FAIL")

