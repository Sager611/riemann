#!/usr/bin/env python3
# alpdep method 2 (independent of build_graph.py): backward closure of a seed name set by rg whole-word matching
# on SHORT names (over-approximates Lean name resolution), enclosing declaration taken from alpcore's decls.txt.
import sys, subprocess, re, os, time
SRC = '/Users/adriansager/playground/riemann/.scratch/alp/lean/src'
DECLS = '/Users/adriansager/playground/riemann/.scratch/alpcore/decls.txt'
seeds = sys.argv[1].split(','); roots = set(sys.argv[2].split(',')); maxlev = int(sys.argv[3]); cap = int(sys.argv[4])
KWD = re.compile(r"^(?:(?:@\[[^\]]*\]\s*)|(?:private|protected|noncomputable|partial|unsafe|nonrec|scoped|local)\s+)*(theorem|lemma|def|abbrev|structure|class|inductive|instance|opaque|axiom)\s+([^\s:(\[{]+)")
seen = set(seeds); frontier = list(seeds); lev = 0; t0 = time.time()
users_of = {}
while frontier and lev < maxlev:
    lev += 1
    args = ['rg', '-j1', '-n', '-o', '-w', '--no-heading']
    for n in frontier: args += ['-e', n]
    out = subprocess.run(args + ['.'], cwd=SRC, capture_output=True, text=True).stdout
    hits = {}
    for ln in out.splitlines():
        p = ln.split(':')
        if len(p) < 3: continue
        hits.setdefault(p[0], set()).add((int(p[1]), p[-1]))
    dl = {}
    with open(DECLS, encoding='utf-8', errors='replace') as fh:
        for ln in fh:
            f, _, rest = ln.partition(':')
            if f in hits:
                l_, _, txt = rest.partition(':')
                m = KWD.match(txt)
                if m: dl.setdefault(f, []).append((int(l_), m.group(2).split('.')[-1], m.group(2)))
    new = []
    for f, hs in hits.items():
        ds = sorted(dl.get(f, []))
        for (l_, w) in hs:
            enc = None
            for d in ds:
                if d[0] <= l_: enc = d
                else: break
            if enc is None: continue
            short = enc[1]
            if short == w and enc[0] == l_: continue
            users_of.setdefault(w, set()).add('%s:%d %s' % (f[2:], enc[0], enc[2]))
            if short not in seen:
                seen.add(short); new.append(short)
    hitroots = sorted(roots & seen)
    print('level', lev, 'frontier', len(frontier), 'hitfiles', len(hits), 'new', len(new), 'total', len(seen), 'roots', hitroots, 'secs', round(time.time() - t0), flush=True)
    print('  new:', ' '.join(sorted(new)[:200]), flush=True)
    if hitroots or len(seen) > cap: break
    frontier = new
print('USERS')
for w in sorted(users_of):
    for u in sorted(users_of[w]): print(' ', w, '<-', u)
