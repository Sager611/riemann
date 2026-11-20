import re, collections
D = '/Users/adriansager/playground/riemann/.scratch/alpaudit/'
S = '/Users/adriansager/playground/riemann/.scratch/alp/lean/src/'
print('== exports')
ex = collections.Counter()
for l in open(D + 'cls_cmd.tsv'):
    loc, txt = l.rstrip('\n').split('\t', 1)
    if re.match(r'^\s*export\b', txt):
        fn, ln = loc.split(':')[:2]
        lines = open(S + fn, encoding='utf-8').read().split('\n')
        ns = [x for x in lines[:int(ln)] if re.match(r'^\s*(namespace|end)\b', x)]
        depth = []
        for x in ns:
            x = x.strip().split()
            if x[0] == 'namespace': depth.append(x[1])
            elif depth and len(x) > 1 and depth[-1].endswith(x[1]): depth.pop()
            elif depth and len(x) == 1: pass
        ex[('.'.join(depth) or '<root>', txt.strip()[:90])] += 1
for (ns, t), v in ex.most_common(): print(v, ns, '|', t)
print('== inductives: nested candidates')
ind = [l.rstrip('\n').split('\t', 1) for l in open(D + 'cls_ind.tsv')]
kinds = collections.Counter(); nested = []; mutual = []
for loc, txt in ind:
    fn, ln = loc.split(':')[:2]; ln = int(ln)
    k = txt.split()[0]; kinds[k] += 1
    lines = open(S + fn, encoding='utf-8').read().split('\n')
    if k == 'mutual': mutual.append(loc); continue
    m = re.match(r'^\s*(?:class\s+)?inductive\s+([^\s({\[:]+)', txt)
    if not m: continue
    name = m.group(1).split('.')[-1]
    block = [lines[ln - 1]]
    for x in lines[ln:ln + 40]:
        if x.strip() == '' or re.match(r'^(theorem|lemma|def|instance|end|namespace|section|open|variable|noncomputable|abbrev|@\[|/--|deriving|structure|inductive|class)\b', x): break
        block.append(x)
    b = ' '.join(block)
    if re.search(r'\b(List|Array|Option|Prod|Sum|Vector|Multiset|Finset|Fin\s*\d*\s*→|Subtype|Sigma|PSigma)\b[^|]*\b' + re.escape(name) + r'\b', b):
        nested.append((loc, b[:200]))
print('kinds', kinds, 'mutual blocks', mutual[:10])
print('nested candidates', len(nested))
for n in nested[:12]: print(n)
