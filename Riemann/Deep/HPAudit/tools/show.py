import re, sys, collections
f = sys.argv[1]; cats = sys.argv[2].split(',')
txt = open(f).read().split('\n')
cur = None; out = {}
for l in txt:
    m = re.match(r'== (\S+): (\d+)', l)
    if m: cur = m.group(1); out[cur] = []; continue
    if cur: out[cur].append(l)
lim = int(sys.argv[3]) if len(sys.argv) > 3 else 25
for k in cats:
    if k.startswith('~'):
        k = k[1:]
        c = collections.Counter()
        for l in out.get(k, []):
            t = l.split(': ', 1)[1] if ': ' in l else l
            c[re.sub(r'\s+', ' ', t)[:100]] += 1
        print('##', k, 'distinct', len(c), 'total', len(out.get(k, [])))
        for kk, v in c.most_common(lim): print('  ', v, kk)
        continue
    print('##', k, len(out.get(k, [])))
    for l in out.get(k, [])[:lim]: print('  ', l[:260])
