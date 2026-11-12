"""Summarise certTW_*.txt / certmid / certc logs: kappa, r, N, D, operators, cells, ALL_OK, peak MB, seconds, max C/deltaD."""
import re, glob, os
import os
HERE = os.path.dirname(os.path.abspath(__file__))
os.chdir(HERE)
rows = []
for f in sorted(glob.glob('logs/certTW_r*_*.txt')):
    if 'fail' in f.lower(): continue
    t = open(f).read()
    m = re.search(r'kappa=(\S+) N=(\d+) D=(\d+) gam=(\d+) Nw=(\d+) prec=(\d+)', t)
    res = re.search(r'RESULT kappa=\S+: ALL_OK=(\w+) covered=(\w+) ; peak (\d+) MB, (\d+)s', t)
    fl = re.search(r'lam/k\^2 = ([0-9.]+)', t)
    cells = re.search(r'(\d+) cells for', t)
    ops = len(re.findall(r'chol=', t))
    rat = []
    for dd, cc in re.findall(r'deltaD=\[([0-9.e+-]+) [^\]]*\] C=\[([0-9.e+-]+)', t):
        try: rat.append(float(cc) / float(dd))
        except: pass
    rr = re.search(r'certTW_r(\d+)_', f).group(1)
    rows.append((int(rr), float(eval(m.group(1).replace('/', '/'))) if m else 0, f, m, res, fl, cells, ops, max(rat) if rat else None))
rows.sort()
print('| r | kappa | float sup_v lam/k^2 | N | D | prec | operators (v-cells) | max C/delta_D | ALL_OK | peak MB | s |')
print('|---|---|---|---|---|---|---|---|---|---|---|')
for r, k, f, m, res, fl, cells, ops, mr in rows:
    if m is None or res is None:
        print('| %d | (incomplete: %s) |' % (r, f)); continue
    print('| %d | %s | %s | %s | %s | %s | %d (%s) | %.2f | %s | %s | %s |' % (r, m.group(1), fl.group(1) if fl else '?', m.group(2), m.group(3), m.group(6), ops, cells.group(1) if cells else '?', mr or 0, res.group(1) if res else 'NO RESULT', res.group(3) if res else '?', res.group(4) if res else '?'))
