import json, unicodedata, collections
A = '/Users/adriansager/playground/riemann/.scratch/alp/lean/'
W = '/Users/adriansager/playground/riemann/.scratch/alpaudit/web/'
p = {x['name']: (x['url'], x['rev']) for x in json.load(open(A + 'lake-manifest.json'))['packages']}
m = {x['name']: (x['url'], x['rev']) for x in json.load(open(W + 'M_lake-manifest.json'))['packages']}
for k in sorted(set(p) | set(m)):
    print('%-17s project %s | mathlib-pin %s | %s' % (k, p.get(k, ('-', '-'))[1][:12], m.get(k, ('-', '-'))[1][:12], 'SAME' if p.get(k) == m.get(k) else ('n/a (root dep)' if k == 'mathlib' else 'DIFF')))
for f in ['src/Main.lean', 'src/Statements.lean', 'src/Supplements/Kernels/MainTheorems.lean', 'src/Supplements/Kernels/MainTheoremsAll.lean']:
    s = open(A + f, encoding='utf-8').read()
    c = collections.Counter(ch for ch in s if ord(ch) > 126 or (ord(ch) < 32 and ch not in '\n'))
    print(f, {('U+%04X %s' % (ord(ch), unicodedata.name(ch, '?'))): n for ch, n in c.items()})
