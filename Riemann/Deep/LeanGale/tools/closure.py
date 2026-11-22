# Import closure of Mathlib modules (from source headers) and the summed sizes of their .olean parts.
# usage: python3 closure.py Mathlib.X ...
import os, re, sys, json
root = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../../../../.lake/packages')
pk = {}
for p in os.listdir(root):
    pk[p] = os.path.join(root, p)
src_roots = []  # (pkg, srcdir)
for p, d in pk.items():
    src_roots.append(d)
def find_src(mod):
    rel = mod.replace('.', '/') + '.lean'
    for d in src_roots:
        f = os.path.join(d, rel)
        if os.path.exists(f):
            return f, d
    return None, None
def olean_size(mod, d):
    rel = mod.replace('.', '/')
    base = os.path.join(d, '.lake/build/lib/lean', rel)
    s = {}
    for ext in ['.olean', '.olean.server', '.olean.private']:
        f = base + ext
        s[ext] = os.path.getsize(f) if os.path.exists(f) else 0
    return s
imp_re = re.compile(r'^(?:public\s+|private\s+)?(?:meta\s+)?import\s+(?:all\s+)?([A-Za-z0-9_.\u00ab\u00bb]+)', re.M)
cache = {}
def imports(mod):
    if mod in cache: return cache[mod]
    f, d = find_src(mod)
    if f is None:
        cache[mod] = (None, [])
        return cache[mod]
    txt = open(f, encoding='utf-8').read()
    # header only: stop at first non-import, non-comment line roughly
    head = txt.split('\n/-!')[0] if '/-!' in txt else txt[:20000]
    ims = imp_re.findall(head)
    cache[mod] = (d, ims)
    return cache[mod]
def closure(mod):
    seen = set(); stack = [mod]; tot = {'.olean':0, '.olean.server':0, '.olean.private':0}; n = 0; missing = 0
    while stack:
        m = stack.pop()
        if m in seen: continue
        seen.add(m)
        d, ims = imports(m)
        if d is None:
            missing += 1
            continue
        n += 1
        s = olean_size(m, d)
        for k in tot: tot[k] += s[k]
        stack.extend(ims)
    return n, missing, tot
for mod in sys.argv[1:]:
    n, miss, tot = closure(mod)
    print(mod, 'mathlib-ish modules', n, 'missing(Init/Lean/Std)', miss, {k: round(v/2**20,1) for k,v in tot.items()})
