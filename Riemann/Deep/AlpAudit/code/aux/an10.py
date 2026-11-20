import re
S = '/Users/adriansager/playground/riemann/.scratch/alp/lean/src/'
dn = re.compile(r'\b(def|theorem|lemma|abbrev|instance|structure|class|inductive|opaque|axiom|alias)\s+(_root_\.)?([A-Za-z0-9_.]*\.)?(RiemannHypothesis|riemannZeta|LFunction|IsPrimitive|NeZero|re|DirichletCharacter|Complex)\b')
out = []
for l in open('/Users/adriansager/playground/riemann/.scratch/alpaudit/scan3.txt', encoding='utf-8', errors='replace'):
    p, ln, t = l.rstrip('\n').split(':', 2)
    m = dn.search(t.split('--')[0])
    if not m: continue
    lines = open(S + p[4:], encoding='utf-8').read().split('\n')[:int(ln)]
    st = []
    for x in lines:
        y = x.strip().split()
        if not y: continue
        if y[0] == 'namespace' and len(y) > 1: st.append(y[1])
        elif y[0] == 'end' and len(y) > 1 and st and st[-1] == y[1]: st.pop()
        elif y[0] == 'end' and len(y) > 1 and st and '.'.join(st).endswith(y[1]):
            k = len(y[1].split('.')); st = st[:-1] if k == 1 else st
    full = '.'.join(st + [(m.group(3) or '') + m.group(4)]) if not m.group(2) else (m.group(3) or '') + m.group(4)
    out.append((full, p[4:] + ':' + ln))
bad = [o for o in out if o[0].split('.')[0] in ('Complex', 'DirichletCharacter') or '.' not in o[0]]
print(len(out), 'decls; in root/Complex/DirichletCharacter namespace:', bad)
import collections
print(collections.Counter(o[0].rsplit('.', 1)[0] if '.' in o[0] else '<root>' for o in out).most_common(30))
