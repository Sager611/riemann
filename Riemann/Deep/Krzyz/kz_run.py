#!/usr/bin/env python3
"""Run kz_check on sampled production leaves.
usage: kz_run.py <samples-dir-or-file> <n> <out.jsonl> <selection-spec> [limit]
Leaves are taken from the author's exported samples (n=7: returned/batches/*.samples.jsonl.gz;
n=5,6: *_spread_leaves.jsonl).  Mode 'indep' ignores the logged (contracted) mass box and
uses [0, T_hi]; only if that fails within the node limit is 'logged' mode tried."""
import sys, os, json, gzip, glob, random, time
import kz_check as K

def load(src, n):
    rows = []
    if src.endswith('.jsonl'):
        for line in open(src):
            r = json.loads(line); r['kind'] = 'leaf'; rows.append(r)
    else:
        for p in sorted(glob.glob(src + '/*.samples.jsonl.gz')):
            for line in gzip.open(p, 'rt'):
                r = json.loads(line)
                if r['kind'] == 'leaf':
                    rows.append(r)
    return rows

def select(rows, spec, seed=20261104):
    rng = random.Random(seed)
    chosen = {}
    def add(r, tag):
        key = (r.get('attempt', ''), r['path'])
        if key not in chosen:
            chosen[key] = dict(r, tags=[tag] + ([r['tag']] if 'tag' in r else []))
        elif tag not in chosen[key]['tags']:
            chosen[key]['tags'].append(tag)
    for part in spec.split(','):
        kind, _, num = part.partition(':'); num = int(num)
        if kind == 'random':
            for r in rng.sample(rows, min(num, len(rows))): add(r, 'random')
        elif kind == 'deepest':
            for r in sorted(rows, key=lambda r: -r['depth'])[:num]: add(r, 'deepest')
        elif kind.startswith('reason='):
            reason = kind.split('=', 1)[1]
            pool = [r for r in rows if r['reason'] == reason]
            for r in rng.sample(pool, min(num, len(pool))): add(r, 'reason')
        elif kind == 'paths':
            pass
    return list(chosen.values())

if __name__ == '__main__':
    src, n, out, spec = sys.argv[1], int(sys.argv[2]), sys.argv[3], sys.argv[4]
    limit = int(sys.argv[5]) if len(sys.argv) > 5 else 2000
    extra = set()
    if len(sys.argv) > 6:
        extra = set(json.load(open(sys.argv[6])))   # list of paths to force-include
    rows = load(src, n)
    sel = select(rows, spec)
    if extra:
        have = {r['path'] for r in sel}
        for r in rows:
            if r['path'] in extra and r['path'] not in have:
                sel.append(dict(r, tags=['forced'])); have.add(r['path'])
    sl = os.environ.get('SLICE')
    if sl:
        a, b = (int(x) for x in sl.split(':')); sel = sel[a:b]
    t0 = time.time(); done = 0; fails = 0
    with open(out, 'w') as f:
        for r in sel:
            c = r['cell']; s = time.time()
            ok, nodes, reasons, md = K.check_cell(n, c['lo'], c['hi'], None, None, limit=limit)
            mode = 'indep'
            if not ok:
                ok, nodes2, reasons2, md2 = K.check_cell(n, c['lo'], c['hi'], c['ml'], c['mh'], limit=limit)
                mode = 'logged'; nodes, reasons, md = nodes2, reasons2, md2
            fails += (not ok); done += 1
            rec = dict(attempt=r.get('attempt'), path=r['path'], node=r.get('node'), depth=r['depth'], prod_reason=r['reason'],
                       tags=r['tags'], passed=ok, mode=mode, nodes=nodes, sub_depth=md, reasons=reasons, sec=round(time.time() - s, 2))
            f.write(json.dumps(rec) + chr(10)); f.flush()
            if done % 20 == 0:
                print(json.dumps(dict(done=done, of=len(sel), fails=fails, sec=round(time.time() - t0, 1))), flush=True)
    print(json.dumps(dict(done=done, fails=fails, sec=round(time.time() - t0, 1))))
