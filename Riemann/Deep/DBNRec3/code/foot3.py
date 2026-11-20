#!/usr/bin/env python3
# dbnrec3 foot3.py (log summary only): per row: shells run by tranche 3, bb CPU seconds (DONE lines), memcap peaks and max RSS of bb runs; certify footprints.
import re, glob, os
L = 'Riemann/Deep/DBNRec3/logs/'
rows = {}
for f in sorted(glob.glob(L + 'r*_*.bb.out')):
    tag = os.path.basename(f)[:-7]; r = int(tag[1:].split('_')[0]); j = int(tag.split('_')[1])
    if r == 8 and j <= 5: continue   # tranche 2's shells (copied)
    t = open(f).read(); d = rows.setdefault(r, dict(n=0, cpu=0, peak=0, rss=0, ppeak=0))
    m = re.search(r'^DONE .*\| (\d+) s$', t, re.M); d['n'] += 1; d['cpu'] += int(m.group(1)) if m else 0
    d['peak'] = max([d['peak']] + [int(x) for x in re.findall(r'memcap: peak (\d+) MB', t)])
    d['rss'] = max([d['rss']] + [int(x)//1048576 for x in re.findall(r'(\d+) +maximum resident', t)])
    p = open(L + tag + '.prep.json').read(); d['ppeak'] = max([d['ppeak']] + [int(x) for x in re.findall(r'memcap: peak (\d+) MB', p)])
for r in sorted(rows): d = rows[r]; print('row %d: %d shells (tranche 3), bb CPU %d s, bb memcap peak %d MB, bb max RSS %d MB, prep_end memcap peak %d MB' % (r, d['n'], d['cpu'], d['peak'], d['rss'], d['ppeak']))
for f in sorted(glob.glob(L + 'cert_r*.err')):
    t = open(f).read(); pk = re.findall(r'memcap: peak (\d+) MB', t); rs = re.findall(r'(\d+) +maximum resident', t); rl = re.findall(r'([0-9.]+) real', t)
    print('%s: memcap peak %s MB, max RSS %s MB, %s s' % (os.path.basename(f)[:-4], pk[-1] if pk else '-', int(rs[-1])//1048576 if rs else '-', rl[-1] if rl else '-'))
