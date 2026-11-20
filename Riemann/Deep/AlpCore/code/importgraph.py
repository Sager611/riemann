#!/usr/bin/env python3
# Build the import graph of the Alpoge archive (streamed; header lines only).
import os, sys, array, re
SRC = "/Users/adriansager/playground/riemann/.scratch/alp/lean/src"
OUT = "/Users/adriansager/playground/riemann/.scratch/alpcore"
files = [l.strip()[2:] for l in open(OUT + "/files.txt") if l.strip()]
mods = [f[:-5].replace("/", ".") for f in files]
idx = {m: i for i, m in enumerate(mods)}
src_e, dst_e = array.array("I"), array.array("I")
ext = {}
nlines = array.array("I")
imp_re = re.compile(r"^import\s+(.*)$")
for i, f in enumerate(files):
    n = 0
    with open(os.path.join(SRC, f), "r", encoding="utf-8", errors="replace") as fh:
        in_hdr = True; in_comment = False
        for line in fh:
            n += 1
            if not in_hdr:
                continue
            s = line.strip()
            if in_comment:
                if "-/" in s: in_comment = False
                continue
            if s.startswith("/-"):
                if "-/" not in s[2:]: in_comment = True
                continue
            if not s or s.startswith("--") or s.startswith("module") or s.startswith("prelude"):
                continue
            m = imp_re.match(s)
            if m:
                for t in m.group(1).split():
                    if t in ("public", "meta", "all"): continue
                    j = idx.get(t)
                    if j is None:
                        r = t.split(".")[0]; ext[r] = ext.get(r, 0) + 1
                    else:
                        src_e.append(i); dst_e.append(j)
            else:
                in_hdr = False
    nlines.append(n)
with open(OUT + "/mods.txt", "w") as fo:
    for m in mods: fo.write(m + "\n")
with open(OUT + "/edges.bin", "wb") as fo:
    src_e.tofile(fo); dst_e.tofile(fo)
with open(OUT + "/nlines.bin", "wb") as fo:
    nlines.tofile(fo)
print("files", len(files), "edges", len(src_e), "total lines", sum(nlines))
print("external roots", sorted(ext.items(), key=lambda x: -x[1])[:10])
