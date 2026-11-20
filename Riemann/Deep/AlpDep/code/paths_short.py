import sys
L = open('Riemann/Deep/AlpDep/logs/analyze2.log').read().split('\n')
cur = None; paths = {}
for ln in L:
    if ln.startswith('TARGET'):
        p = ln.split(); cur = p[1]
        if 'NOT' in ln: cur = None
        else: paths[cur] = []
    elif cur and ln.startswith('    ') and '[' in ln:
        nm = ln.strip().split(' [')[0]; paths[cur].append(nm.split('.')[-1])
for t in sys.argv[1].split(','):
    print(t.split('.')[-1], len(paths[t]) - 1, ' > '.join(paths[t][1:]))
