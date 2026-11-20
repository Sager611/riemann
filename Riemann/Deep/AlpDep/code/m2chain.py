import re, collections
L = open('/Users/adriansager/playground/riemann/Riemann/Deep/AlpDep/logs/rg_closure.log').read().split('\n')
i = L.index('USERS') if 'USERS' in L else None
par = {}
lev = {}
cur = 0
for ln in L[:i]:
    m = re.match(r'level (\d+) ', ln)
    if m: cur = int(m.group(1))
    if ln.startswith('  new:'):
        for n in ln[7:].split(): lev.setdefault(n, cur)
lev['DivIn'] = 0
for ln in L[i+1:]:
    m = re.match(r'\s+(\S+) <- (\S+) (\S+)', ln)
    if not m: continue
    w, loc, full = m.groups(); u = full.split('.')[-1]
    if u in lev and w in lev and lev[w] == lev[u] - 1:
        par.setdefault(u, (w, loc, full))
n = 'everyCharacter_of_row'; out = []
while n in par:
    w, loc, full = par[n]; out.append('%s (%s, L%d) <- mentions %s' % (full, loc, lev[n], w)); n = w
print('\n'.join(out))
