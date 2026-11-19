# cost.py (bothachain): wall and CPU seconds and peak RSS per chain point from the /usr/bin/time -l records in .scratch/bothachain/ccm/runs
# (scan = rscan2 premise, g = gstate, o/d = prep + LDL of the odd / deflated certificate; dA = shared defl of a pair).
import os, re, json
H = os.path.dirname(os.path.abspath(__file__)); RUNS = os.path.join(H, '..', '..', '..', '.scratch', 'bothachain', 'ccm', 'runs')
def tm(f):
    if not os.path.exists(f): return None
    s = open(f).read(); r = re.search(r'([0-9.]+) real\s+([0-9.]+) user', s); m = re.search(r'(\d+)\s+maximum resident', s)
    return (float(r.group(1)), float(r.group(2)), int(m.group(1))//2**20) if r and m else None
R = [d for d in map(json.loads, open(os.path.join(H, 'logs', 'chain.jsonl'))) if d.get('ok') and not d.get('splice')][3:]
print('| c | n (odd) | G | bits | NG | wall s (scan, g, odd prep+LDL, defl prep+LDL) | CPU s total | wall s total | peak RSS MB |')
print('|---|---|---|---|---|---|---|---|---|')
for d in R:
    t = d['c'].replace('/', 'o'); parts = []
    for k in ('scan%s' % t, 'g%s' % t, 'ch%so.prep' % t, 'ch%so.slt' % t, 'ch%sd.prep' % t, 'ch%sd.slt' % t, 'ch%sdA.prep' % t, 'ch%sdA.slt' % t):
        parts.append(tm(os.path.join(RUNS, k + '.time')))
    w = [p[0] if p else 0 for p in parts]; u = [p[1] if p else 0 for p in parts]; mx = max([p[2] for p in parts if p] or [0])
    print('| %s | %s | %s | %s | %s | %.0f, %.0f, %.0f, %.0f | %.0f | %.0f | %d |' % (d['c'], d['odd'].get('n'), d['odd'].get('G'), d['bits'], d['NG'],
          w[0], w[1], w[2] + w[3], w[4] + w[5] + w[6] + w[7], sum(u), sum(w), mx))
