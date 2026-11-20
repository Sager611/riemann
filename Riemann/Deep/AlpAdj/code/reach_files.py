#!/usr/bin/env python3
# BFS (alpdep edges) from given names; write the source files of reached THEOREMS (kind T) to stdout (one per line).
import sys, numpy as np
OUT = '/Users/adriansager/playground/riemann/.scratch/alpdep'
eoff = np.load(OUT + '/eoff.npy', mmap_mode='r'); elen = np.load(OUT + '/elen.npy', mmap_mode='r')
store = np.memmap(OUT + '/estore.bin', dtype=np.int32, mode='r')
targets = set(sys.argv[1:]); ids = []
with open(OUT + '/names.tsv', encoding='utf-8') as f:
    for line in f:
        p = line.split('\t', 2)
        if p[1] in targets: ids.append(int(p[0]))
N = len(eoff); seen = np.zeros(N, dtype=np.uint8); stack = list(ids)
for i in stack: seen[i] = 1
while stack:
    u = stack.pop(); o = int(eoff[u]); l_ = int(elen[u])
    if o < 0: continue
    for v in store[o:o + l_].tolist():
        if 0 <= v < N and not seen[v]: seen[v] = 1; stack.append(v)
fidx = set()
with open(OUT + '/names.tsv', encoding='utf-8') as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if seen[int(p[0])] and p[2] == 'T': fidx.add(int(p[3]))
with open('/Users/adriansager/playground/riemann/.scratch/alpcore/files.txt', encoding='utf-8') as f:
    for k, line in enumerate(f):
        if k in fidx: print(line.strip())
