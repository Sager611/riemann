import os, re, sys, json, collections
# usage: closure.py ROOTMOD srcdir1 [srcdir2 ...]
roots = sys.argv[1].split(',')
srcdirs = sys.argv[2:]
def find(mod):
    rel = mod.replace('.', '/') + '.lean'
    for d in srcdirs:
        p = os.path.join(d, rel)
        if os.path.exists(p): return p
    return None
imp_re = re.compile(r'^\s*(?:public\s+|private\s+|meta\s+)*import\s+(?:all\s+)?([A-Za-z0-9_.«»]+)', re.M)
seen = {}; ext = collections.Counter(); missing = []
stack = list(roots)
while stack:
    m = stack.pop()
    if m in seen: continue
    p = find(m)
    seen[m] = p
    if p is None:
        missing.append(m); continue
    src = open(p, encoding='utf-8', errors='replace').read()
    # only the header (imports before first non-import command) matters, but scanning all is fine
    for im in imp_re.findall(src):
        top = im.split('.')[0]
        if find(im): stack.append(im)
        else: ext[top] += 1; seen.setdefault(im, None)
local = [m for m, p in seen.items() if p]
print('local modules in closure:', len(local))
print('external top-level packages:', dict(ext))
unres = [m for m in missing if m.split('.')[0] not in ('Mathlib','Batteries','Aesop','Qq','Lean','Std','Init','ProofWidgets','Plausible','ImportGraph','LeanSearchClient')]
print('unresolved non-std imports:', len(unres), unres[:20])
json.dump(sorted(local), open(os.environ.get('OUT','closure.json'),'w'))
