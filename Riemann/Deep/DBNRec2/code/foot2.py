#!/usr/bin/env python3
# dbnrec2 foot2.py (log summary only): per row bb CPU seconds (sum of the DONE lines), memcap peaks and max RSS of bb runs, certify footprints.
import glob, re, collections
L = 'Riemann/Deep/DBNRec2/logs/'
rows = collections.defaultdict(lambda: dict(n=0, cpu=0, peak=0, rss=0))
for f in sorted(glob.glob(L + 'r*_*.bb.out')):
    r = int(re.search(r'/r(\d+)_', f).group(1)); t = open(f).read()
    m = re.search(r'DONE K .*\| (\d+) s', t); d = rows[r]; d['n'] += 1; d['cpu'] += int(m.group(1)) if m else 0
    for p in re.findall(r'memcap: peak (\d+) MB', t): d['peak'] = max(d['peak'], int(p))
    for p in re.findall(r'(\d+)  maximum resident set size', t): d['rss'] = max(d['rss'], int(p)//2**20)
for r in sorted(rows): d = rows[r]; print('row %d: %d shells (in DBNRec2/logs), bb CPU %d s, memcap peak %d MB, max RSS %d MB' % (r, d['n'], d['cpu'], d['peak'], d['rss']))
for f in sorted(glob.glob(L + 'cert_r*.err')):
    t = open(f).read(); m = re.search(r'memcap: peak (\d+) MB', t); rs = re.search(r'(\d+)  maximum resident set size', t); rl = re.search(r'([\d.]+) real', t)
    print(f.split('/')[-1], 'memcap', m.group(1) if m else '?', 'MB, max RSS', int(rs.group(1))//2**20 if rs else '?', 'MB,', rl.group(1) if rl else '?', 's')
