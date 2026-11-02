#!/usr/bin/env python3
"""Rewrite a certificate in compact form: leaf multipliers in homogeneous form
(den = 1, numerators divided by their gcd; only ratios matter for the checks),
keys sorted and delta-encoded, files xz-compressed.
Line format of multipliers-NNN.jsonl.xz: [leaf_index, den, key_deltas, numerators].
Usage: compact_cert.py SRC DST"""
import gzip, hashlib, json, lzma, shutil, sys
from math import gcd
from pathlib import Path

src, dst = Path(sys.argv[1]), Path(sys.argv[2])
dst.mkdir(parents=True, exist_ok=False)
man = json.loads((src / "manifest.json").read_text())
shutil.copy(src / "header.json", dst / "header.json")
shutil.copy(src / "tree.json.gz", dst / "tree.json.gz")
files = []
batch, fileno = [], 0


def flush():
    global batch, fileno
    if not batch:
        return
    name = "multipliers-%03d.jsonl.xz" % fileno
    data = ("\n".join(json.dumps(x, separators=(",", ":")) for x in batch) + "\n").encode()
    (dst / name).write_bytes(lzma.compress(data, preset=9 | lzma.PRESET_EXTREME))
    files.append(name)
    batch, fileno = [], fileno + 1


for name in man["files"]:
    opener = lzma.open if name.endswith(".xz") else gzip.open
    with opener(src / name, "rt") as f:
        for line in f:
            rec = json.loads(line)
            if len(rec) == 4:
                idx, den, dk, vals = rec
                keys, acc = [], 0
                for d in dk:
                    acc += d
                    keys.append(acc)
                pairs = list(zip(keys, vals))
            else:
                idx, den, pairs = rec
            g = 0
            for _, n in pairs:
                g = gcd(g, n)
            pairs = sorted((k, n // g) for k, n in pairs)
            keys = [k for k, _ in pairs]
            batch.append([idx, 1, [keys[0]] + [b - a for a, b in zip(keys, keys[1:])], [n for _, n in pairs]])
            if len(batch) >= 2048:
                flush()
flush()
manifest = {"schema": 2, "format": "compact", "files": files,
            "sha256": {n: hashlib.sha256((dst / n).read_bytes()).hexdigest() for n in ["header.json", "tree.json.gz", *files]}}
(dst / "manifest.json").write_text(json.dumps(manifest, sort_keys=True, indent=1) + "\n")
print("compacted", src, "->", dst, len(files), "files")

