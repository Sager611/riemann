
import os, re, sys
root = '.lake/packages'
src_roots = {'Mathlib': '.lake/packages/mathlib', 'Batteries': '.lake/packages/batteries', 'Aesop': '.lake/packages/aesop', 'Qq': '.lake/packages/Qq', 'ProofWidgets': '.lake/packages/proofwidgets', 'Plausible': '.lake/packages/plausible', 'LeanSearchClient': '.lake/packages/LeanSearchClient', 'ImportGraph': '.lake/packages/importGraph'}
lean_prefix = os.popen('lean --print-prefix').read().strip()
def find_src(mod):
    top = mod.split('.')[0]
    if top in src_roots:
        p = os.path.join(src_roots[top], *mod.split('.')) + '.lean'
        if os.path.exists(p): return p, os.path.join(src_roots[top], '.lake/build/lib/lean', *mod.split('.'))
    # core
    p = os.path.join(lean_prefix, 'src/lean', *mod.split('.')) + '.lean'
    return (p if os.path.exists(p) else None), os.path.join(lean_prefix, 'lib/lean', *mod.split('.'))
imp_re = re.compile(r'^\s*(?:public\s+)?(?:meta\s+)?(?:private\s+)?import\s+(?:all\s+)?([A-Za-z0-9_.]+)')
seen = {}
stack = [sys.argv[1]]
while stack:
    m = stack.pop()
    if m in seen: continue
    src, base = find_src(m)
    seen[m] = base
    if src is None: continue
    with open(src) as f:
        for line in f:
            s = line.strip()
            if s.startswith('/-!') or s.startswith('/--') : break
            mm = imp_re.match(line)
            if mm: stack.append(mm.group(1))
tot = {'olean':0,'olean.server':0,'olean.private':0}
cnt = 0
for m, base in seen.items():
    for ext in tot:
        p = base + '.' + ext
        if os.path.exists(p):
            tot[ext] += os.path.getsize(p)
    cnt += 1
print(sys.argv[1], 'modules', cnt, {k: round(v/2**20,1) for k,v in tot.items()})

