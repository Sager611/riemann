# FodgeAudit summarize.py: aggregates logs/grid*.jsonl (falsify.py rows) into max ratios per check.
import json, glob
rows = {}; errs = []
for f in sorted(glob.glob('logs/grid*.jsonl')):
    for l in open(f):
        r = json.loads(l)
        if 'error' in r: errs.append(r); continue
        if 'case' in r: rows[json.dumps(r['case'])] = r
rows = list(rows.values())
keys = ['E1', 'N1', 'N2', 'S1', 'S2', 'S3', 'B1', 'M1', 'D1', 'C3']
nst = sum(min(r['nstat'], 40) if r['nstat'] > 40 else r['nstat'] for r in rows)
print(f'{len(rows)} distinct cases, {len(errs)} errors; stationary j tested ~{nst}, nonstationary ~{41*len(rows)}, exceptional ~{2*len(rows)}')
print(f'lambda range {min(r["lam"] for r in rows):.2e} .. {max(r["lam"] for r in rows):.3f}; max nstat {max(r["nstat"] for r in rows)}')
for k in keys:
    have = [r for r in rows if r[k][1] is not None or r[k][0] > 0]
    b = max(rows, key=lambda r: r[k][0])
    print(f'{k}: max ratio {b[k][0]:.6f} over {len(have)} cases; at case {b["case"]} lam={b["lam"]:.4g} info={str(b[k][1])[:70]}')
over = [(k, r['case'], r[k][0]) for r in rows for k in keys if r[k][0] > 1]
print('ratios > 1:', over[:10] if over else 'none')
