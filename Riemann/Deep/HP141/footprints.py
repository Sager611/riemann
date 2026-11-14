# footprints.py -- per-run footprint table from every HP141 log written by .scratch/hp141/run.sh
# (memcap top-MEM peak, /usr/bin/time -l max RSS, wall seconds, exit code); output: footprints.md
import os, re, glob
HERE = os.path.dirname(os.path.abspath(__file__))
rows = []
for p in sorted(glob.glob(os.path.join(HERE, '*', 'logs', '*'))):
    if not os.path.isfile(p) or p.endswith(('.json', '.gz', '.md', '.py')): continue
    t = open(p, errors='ignore').read()
    if 'memcap: peak' not in t: continue
    mc = re.search(r'memcap: peak (\d+) MB([^\n]*)', t)
    rss = re.search(r'(\d+)\s+maximum resident set size', t)
    real = re.search(r'([\d.]+) real', t)
    rc = re.search(r'== end \S+ rc=(\d+)', t)
    rows.append((os.path.relpath(p, HERE), mc.group(1) + (mc.group(2).strip() and ' ' + mc.group(2).strip()), '%.0f' % (int(rss.group(1)) / 2**20) if rss else '-',
                 real.group(1) if real else '-', rc.group(1) if rc else '-'))
out = ['| log | memcap peak MB (top MEM) | max RSS MB (time -l) | wall s | rc |', '|---|---|---|---|---|']
out += ['| %s | %s | %s | %s | %s |' % r for r in rows]
open(os.path.join(HERE, 'footprints.md'), 'w').write('\n'.join(out) + '\n')
print('\n'.join(out))
