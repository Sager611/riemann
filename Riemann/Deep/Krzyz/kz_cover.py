#!/usr/bin/env python3
"""Independent streaming coverage audit of a complete reduced-search node log (native-v1 format),
e.g.  unzip -p logs.zip .../m4_n6_full.jsonl.gz | gunzip -c | kz_cover.py 6 leaves_out.jsonl
Checks (own code, not the author's auditor):
 * header root contains the normalized domain: theta_1 in [0,pi], theta_j in [0,2pi], q in [1,q_cap]
   with T_n(q_cap) > n*T_hi+n-1, masses [0,T_hi];
 * every node is a pending child (each path consumed exactly once);
 * logged contracted parent is a subset of its input box;
 * the parameter part of every contraction is implied by ordering/conjugation/q>=1 (exact rationals
   with a 200-bit bracket of pi), so no admissible (theta,q) is dropped;
 * splits: strict interior binary64 midpoint, children = the two halves of the contracted parent;
 * nothing pending at the end; nodes = 2*leaves-1; per-reason counts.
Leaf input cells are sampled (reservoir + all of rare reasons + deepest) for kz_check."""
import sys, json, random, heapq
try:
    import orjson
    _loads = orjson.loads
except Exception:
    _loads = json.loads
from fractions import Fraction as Fr
from flint import arb, ctx
ctx.prec = 200
_pi = arb.pi()
def _fr(x):
    m, e = x.mid().man_exp(); m, e = int(m), int(e)
    return Fr(m * 2 ** e) if e >= 0 else Fr(m, 2 ** (-e))
PI_LO, PI_HI = _fr(_pi.lower()), _fr(_pi.upper())
TWO_PI_LO, TWO_PI_HI = 2 * PI_LO, 2 * PI_HI
THI = {5: 5.84, 6: 6.99, 7: 8.15}

def main():
    n = int(sys.argv[1]); out = sys.argv[2]; N = n - 1
    rng = random.Random(7)
    R = 3000; reservoir = []; seen_leaves = 0
    rare = {}; deepest = []; value_res = []; value_seen = 0
    pending = {}; orphans = {}
    nodes = leaves = 0; counts = {}; split_kinds = {}; problems = []; maxdepth = 0; maxpend = 0
    def problem(msg):
        if len(problems) < 50: problems.append(msg)
    def cell_t(c):
        return (tuple(c['lo']), tuple(c['hi']), tuple(c['ml']), tuple(c['mh']))
    def process(p, rec, inp):
        nonlocal nodes, leaves, maxdepth, seen_leaves, value_seen
        nodes += 1
        d = len(p) - 1; maxdepth = max(maxdepth, d)
        lo, hi, ml, mh = inp
        if 'close' in rec:
            leaves += 1; why = rec['close']; counts[why] = counts.get(why, 0) + 1
            item = dict(path=p, depth=d, reason=why, cell=dict(lo=list(lo), hi=list(hi), ml=list(ml), mh=list(mh)))
            seen_leaves += 1
            if len(reservoir) < R: reservoir.append(item)
            else:
                k = rng.randrange(seen_leaves)
                if k < R: reservoir[k] = item
            if why == 'chart_value':
                value_seen += 1
                if len(value_res) < 400: value_res.append(item)
                else:
                    k = rng.randrange(value_seen)
                    if k < 400: value_res[k] = item
            elif why not in ('dual', 'chart_moment_residual'):
                lst = rare.setdefault(why, [])
                if len(lst) < 60: lst.append(item)
            if len(deepest) < 40: heapq.heappush(deepest, (d, p, item))
            elif d > deepest[0][0]: heapq.heapreplace(deepest, (d, p, item))
            return
        par = rec['parent']; plo, phi, pml, pmh = par['lo'], par['hi'], par['ml'], par['mh']
        # subset
        for j in range(n):
            if not (plo[j] >= lo[j] and phi[j] <= hi[j] and plo[j] <= phi[j]): problem(('not subset param', p, j))
        for j in range(N):
            if not (pml[j] >= ml[j] and pmh[j] <= mh[j] and pml[j] <= pmh[j]): problem(('not subset mass', p, j))
        # parameter contraction justified by ordering / conjugation / q >= 1
        dlo = list(lo[:N]); dhi = list(hi[:N])
        for j in range(1, N): dlo[j] = max(dlo[j], dlo[j - 1])
        for j in range(N - 2, -1, -1): dhi[j] = min(dhi[j], dhi[j + 1])
        empty = any(dlo[j] > dhi[j] for j in range(N)) or (dlo[0] + dlo[N - 1] > 6.28 and Fr(dlo[0]) + Fr(dlo[N - 1]) > TWO_PI_HI)
        if not empty:
            for j in range(N):
                if plo[j] > dlo[j]: problem(('lo contracted beyond ordering', p, j))
                if phi[j] < dhi[j]:
                    ok = False
                    if j in (0, N - 1):
                        other = dlo[N - 1] if j == 0 else dlo[0]
                        ok = Fr(phi[j]) >= TWO_PI_HI - Fr(other)
                    if not ok: problem(('hi contracted beyond ordering/conjugation', p, j, phi[j], dhi[j]))
            if plo[N] > max(1.0, lo[N]): problem(('q lo contracted', p))
            if phi[N] < hi[N]: problem(('q hi contracted', p))
        kind, k, mid = rec['split']; split_kinds[kind] = split_kinds.get(kind, 0) + 1
        a = [list(plo), list(phi), list(pml), list(pmh)]; b = [list(plo), list(phi), list(pml), list(pmh)]
        if kind == 'parameter':
            if not (plo[k] < mid < phi[k]): problem(('bad param midpoint', p))
            a[1][k] = mid; b[0][k] = mid
        elif kind == 'mass':
            if not (pml[k] < mid < pmh[k]): problem(('bad mass midpoint', p))
            a[3][k] = mid; b[2][k] = mid
        else:
            problem(('unknown split', p))
        for child, c in ((p + '0', a), (p + '1', b)):
            if child in pending: problem(('duplicate child', child))
            pending[child] = tuple(tuple(x) for x in c)
    header = None
    footer = None
    for line in sys.stdin:
        if line.startswith('{"kind":"header"'):
            header = json.loads(line)
            assert header['n'] == n
            for rt in header['roots']:
                c = rt['cell']; lo, hi, ml, mh = c['lo'], c['hi'], c['ml'], c['mh']
                ok = all(x == 0 for x in lo[:N]) and lo[N] == 1 and Fr(hi[0]) >= PI_HI and all(Fr(h) >= TWO_PI_HI for h in hi[1:N])
                qc = Fr(hi[N]); a_, b_ = Fr(1), qc
                for _ in range(n - 1): a_, b_ = b_, 2 * qc * b_ - a_
                ok = ok and b_ > n * Fr(THI[n]) + n - 1
                ok = ok and all(x == 0 for x in ml) and all(Fr(x) >= Fr(THI[n]) - Fr(1, 1000) for x in mh)
                print('root', rt['path'], 'contains normalized domain:', ok, flush=True)
                if not ok: problem(('root', rt))
                pending[rt['path']] = cell_t(c)
            continue
        rec = _loads(line)
        if rec.get('kind') == 'footer':
            footer = rec
            continue
        p = rec['p']
        if p in pending:
            process(p, rec, pending.pop(p))
            # release orphans that became processable
            while orphans:
                done = [q for q in orphans if q in pending]
                if not done: break
                for q in done: process(q, orphans.pop(q), pending.pop(q))
        else:
            if p in orphans: problem(('duplicate node', p))
            orphans[p] = rec
        if len(pending) > maxpend: maxpend = len(pending)
    res = dict(n=n, nodes=nodes, leaves=leaves, nodes_eq_2leaves_minus_1=(nodes == 2 * leaves - 1), pending_left=len(pending),
               orphans_left=len(orphans), counts=counts, split_kinds=split_kinds, maxdepth=maxdepth, max_pending=maxpend, footer=footer,
               problems=problems)
    print(json.dumps(res), flush=True)
    with open(out, 'w') as f:
        for tag, items in (('random', reservoir), ('chart_value', value_res), ('deepest', [x[2] for x in deepest])):
            for it in items: f.write(json.dumps(dict(it, tag=tag)) + chr(10))
        for why, lst in rare.items():
            for it in lst: f.write(json.dumps(dict(it, tag='rare')) + chr(10))

if __name__ == '__main__':
    main()
