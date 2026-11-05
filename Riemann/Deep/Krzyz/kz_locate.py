#!/usr/bin/env python3
"""Locate the production leaves (complete n=6 reduced log) that contain genuine solutions of the
reduced moment system (KKT points of 5-atom configurations, found numerically).  Such a leaf can
only be closed soundly by a value test (|a_6| = e^{-t}C < 2/e); a 'dual', 'moment_contraction',
'mass_sum' or 'chart_moment_residual' closure there would expose an unsound contraction/test.
usage: unzip -p logs.zip .../m4_n6_full.jsonl.gz | gunzip -c | kz_locate.py 6 kkt.json out.json"""
import sys, json, math
try:
    import orjson; _loads = orjson.loads
except Exception:
    _loads = json.loads

def inside(pt, lo, hi, eps=1e-12):
    return all(lo[j] - eps <= pt[j] <= hi[j] + eps for j in range(len(pt)))

def main():
    n = int(sys.argv[1]); pts = json.load(open(sys.argv[2])); out = sys.argv[3]; N = n - 1
    P = [s['theta'] + [s['q']] for s in pts]
    tracked = {}; found = []
    for line in sys.stdin:
        if line.startswith('{"kind":"header"'):
            h = json.loads(line)
            for rt in h['roots']:
                c = rt['cell']; ids = [i for i, pt in enumerate(P) if inside(pt, c['lo'], c['hi'])]
                if ids: tracked[rt['path']] = (c, ids)
            continue
        if not line.startswith('{"p":"'):
            continue
        p = line[6:line.index('"', 6)]
        if p not in tracked:
            continue
        cell, ids = tracked.pop(p); rec = _loads(line)
        if 'close' in rec:
            for i in ids:
                lam = pts[i]['lam']
                if not all(cell['ml'][j] - 1e-9 <= lam[j] <= cell['mh'][j] + 1e-9 for j in range(N)):
                    found.append(dict(point=i, an=pts[i]['an_direct'], path=p, depth=len(p) - 1, reason='MASS_OUTSIDE_LEAF_INPUT', cell=cell))
                found.append(dict(point=i, an=pts[i]['an_direct'], path=p, depth=len(p) - 1, reason=rec['close'], cell=cell))
            continue
        par = rec['parent']; kind, k, mid = rec['split']
        ids2 = [i for i in ids if inside(P[i], par['lo'], par['hi'])]
        for i in ids2:
            lam = pts[i]['lam']
            if not all(par['ml'][j] - 1e-9 <= lam[j] <= par['mh'][j] + 1e-9 for j in range(N)):
                found.append(dict(point=i, an=pts[i]['an_direct'], path=p, depth=len(p) - 1, reason='MASS_DROPPED_BY_CONTRACTION', cell=cell, parent=par))
        if len(ids2) != len(ids):
            for i in ids:
                if i not in ids2:
                    found.append(dict(point=i, an=pts[i]['an_direct'], path=p, depth=len(p) - 1, reason='DROPPED_BY_CONTRACTION', cell=cell, parent=par))
        a = dict(lo=list(par['lo']), hi=list(par['hi']), ml=list(par['ml']), mh=list(par['mh']))
        b = dict(lo=list(par['lo']), hi=list(par['hi']), ml=list(par['ml']), mh=list(par['mh']))
        if kind == 'parameter':
            a['hi'][k] = mid; b['lo'][k] = mid
        else:
            a['mh'][k] = mid; b['ml'][k] = mid
        for child, c in ((p + '0', a), (p + '1', b)):
            idc = [i for i in ids2 if inside(P[i], c['lo'], c['hi'])]
            if idc: tracked[child] = (c, idc)
    json.dump(dict(found=found, untracked_left=len(tracked)), open(out, 'w'), indent=1)
    for f in found:
        print(json.dumps(dict(point=f['point'], an=round(f['an'], 9), depth=f['depth'], reason=f['reason'], path=f['path'])))

if __name__ == '__main__':
    main()
