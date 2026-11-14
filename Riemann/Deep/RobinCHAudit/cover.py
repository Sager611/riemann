#!/usr/bin/env python3
# cover.py (robinchaudit) -- box bookkeeping of RobinCH's own runs (boxes.json from dump.py): consecutive passing boxes
# share endpoints, a failed box is retried from the same left end, first box starts at L0, last ends at log X (exact arb).
#   uv run --no-project --with python-flint python3 cover.py boxes.json
import sys, json
from flint import arb, ctx
ctx.prec = 256
def A(me): return arb(me[0])*arb(2)**me[1]
D = json.load(open(sys.argv[1])); B = D['boxes']; LN10 = arb(10).log()
for s in D['summary']:
    if s['X'] == 'probe': continue
    rows = B[s['first']:s['last']]; ok = [r for r in rows if r['ok']]
    L0 = (7 if s['kind'] == 'N' else 10)*LN10; L1 = arb(s['X']).log()
    bad = 0; cur = A(rows[0]['La'])
    for r in rows:
        La, Lb = A(r['La']), A(r['Lb'])
        if not (La <= cur): bad += 1                     # a box must start at or before the covered frontier
        if r['ok']: cur = Lb
    wmin = min(float((A(r['Lb']) - A(r['La'])).mid()) for r in ok)
    print(f"{s['kind']} b0={s['b0']} to {s['X']}: {len(ok)} passing / {len(rows)} boxes; first La <= log 10^{7 if s['kind']=='N' else 10}: "
          f"{bool(A(ok[0]['La']) <= L0)}; frontier gaps: {bad}; last Lb >= log X: {bool(cur >= L1)} "
          f"(log10 end = {float((cur/LN10).mid()):.9f}); min width {wmin:.3g}; max k {max(r['k'] for r in ok)}; min m {min(r['m'] for r in ok):.4g}")
