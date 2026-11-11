import os, re, difflib
exec(open(__import__('os').path.join(__import__('os').path.dirname(__file__), 'cmpoai.py')).read().split("A = files(ours_root)")[0])
A = files(ours_root); B = files(trk_root)
import collections
kinds = collections.Counter(); out = []
for k in sorted(set(A) & set(B)):
    a = norm(open(A[k], encoding='utf-8').read(), 'ours').split('\n'); b = norm(open(B[k], encoding='utf-8').read(), 'trk').split('\n')
    if a == b: continue
    d = [l for l in difflib.unified_diff(b, a, lineterm='', n=0) if (l.startswith('+') or l.startswith('-')) and not l.startswith('+++') and not l.startswith('---')]
    plus = [l[1:] for l in d if l.startswith('+')]; minus = [l[1:] for l in d if l.startswith('-')]
    # classify: brace fix if removing whitespace and parens makes them equal
    def squash(x): return re.sub(r'[\s()]', '', ''.join(x))
    if squash(plus) == squash(minus): kind = 'parens/whitespace only (set-literal {X}->{(X)})'
    elif any('all_goals' in l for l in plus): kind = 'congr -> all_goals'
    elif any('ULift.ext' in l for l in plus): kind = 'ULift.ext'
    else: kind = 'OTHER'
    kinds[kind] += 1
    out.append((k, kind, len(plus), len(minus), plus[:2], minus[:2]))
for k, v in kinds.items(): print(v, k)
for k, kind, np, nm, p, m in out:
    if kind == 'OTHER' or kind != 'parens/whitespace only (set-literal {X}->{(X)})':
        print(k, kind, '+%d -%d' % (np, nm)); print('    +', [x[:120] for x in p]); print('    -', [x[:120] for x in m])
