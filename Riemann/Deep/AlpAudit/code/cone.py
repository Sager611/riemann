#!/usr/bin/env python3
# Textual import cone of src/Main.lean (BFS over 'import' headers; streams one file at a time).
import os, sys, re, collections
SRC = sys.argv[1]
OUT = sys.argv[2]
imp_re = re.compile(r'^import\s+(.*)$')
def header_imports(path):
    mods = []
    in_block = 0
    with open(path, encoding='utf-8') as f:
        for line in f:
            s = line.strip()
            if in_block:
                if '-/' in s: in_block = 0
                continue
            if not s or s.startswith('--'): continue
            if s.startswith('/-'):
                if '-/' not in s[2:]: in_block = 1
                continue
            m = imp_re.match(s)
            if m:
                mods += m.group(1).split('--')[0].split()
                continue
            break
    return mods
start = ['Main']
seen = set(); ext = collections.Counter(); order = []
edges = 0
q = collections.deque(start)
while q:
    m = q.popleft()
    if m in seen: continue
    seen.add(m)
    p = os.path.join(SRC, *m.split('.')) + '.lean'
    if not os.path.exists(p):
        ext[m] += 1
        continue
    order.append(p[len(SRC)+1:])
    for d in header_imports(p):
        edges += 1
        if d not in seen: q.append(d)
with open(os.path.join(OUT, 'cone_files.txt'), 'w') as f:
    f.write('\n'.join(sorted(order)) + '\n')
extmods = sorted(m for m in seen if m not in set(x[:-5].replace('/', '.') for x in order))
with open(os.path.join(OUT, 'cone_external_imports.txt'), 'w') as f:
    f.write('\n'.join(extmods) + '\n')
tops = collections.Counter(x.split('/')[0] for x in order)
tot_lines = 0; tot_bytes = 0
for x in order:
    p = os.path.join(SRC, x); tot_bytes += os.path.getsize(p)
    with open(p, 'rb') as f: tot_lines += sum(1 for _ in f)
print('cone project files', len(order), 'lines', tot_lines, 'bytes', tot_bytes, 'import edges', edges)
print('by top dir', dict(tops))
print('external modules imported from cone', len(extmods), 'top roots', collections.Counter(m.split('.')[0] for m in extmods))
