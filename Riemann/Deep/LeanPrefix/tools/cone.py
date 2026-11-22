#!/usr/bin/env python3
"""Import cone of one or more modules: module count and .olean bytes (public part).
usage (repo root): python3 Riemann/Deep/LeanPrefix/tools/cone.py Mod1 [Mod2 ...]"""
import os, re, sys
roots = {'Mathlib': '.lake/packages/mathlib', 'Batteries': '.lake/packages/batteries', 'Aesop': '.lake/packages/aesop',
         'Qq': '.lake/packages/Qq', 'ProofWidgets': '.lake/packages/proofwidgets', 'Plausible': '.lake/packages/plausible',
         'LeanSearchClient': '.lake/packages/LeanSearchClient', 'ImportGraph': '.lake/packages/importGraph'}
prefix = os.path.expanduser('~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2')
def locate(mod):
    parts = mod.split('.')
    if parts[0] in roots:
        r = roots[parts[0]]
        return os.path.join(r, *parts) + '.lean', os.path.join(r, '.lake/build/lib/lean', *parts)
    return os.path.join(prefix, 'src/lean', *parts) + '.lean', os.path.join(prefix, 'lib/lean', *parts)
imp = re.compile(r'^\s*(?:public\s+)?(?:meta\s+)?(?:private\s+)?import\s+(?:all\s+)?([A-Za-z0-9_.]+)')
seen, stack = {}, list(sys.argv[1:])
while stack:
    m = stack.pop()
    if m in seen:
        continue
    src, base = locate(m)
    seen[m] = base
    if not os.path.exists(src):
        continue
    with open(src) as f:
        for line in f:
            if line.lstrip().startswith(('/-!', '/--')):
                break
            mm = imp.match(line)
            if mm:
                stack.append(mm.group(1))
tot = {e: 0 for e in ('olean', 'olean.server', 'olean.private')}
for base in seen.values():
    for e in tot:
        if os.path.exists(base + '.' + e):
            tot[e] += os.path.getsize(base + '.' + e)
print(' + '.join(sys.argv[1:]), '| modules', len(seen), '|', ', '.join('%s %.0f MB' % (e, v / 2**20) for e, v in tot.items()))

