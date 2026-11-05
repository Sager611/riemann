#!/usr/bin/env python3
"""Check the n=7, N=6 batch chain from the compact METADATA.json (no bulk logs needed):
 * rebuild each batch's output frontier file (KRZYZ_NATIVE_FRONTIER_V1, C %a hex floats) from the
   JSON summary and compare its SHA-256 with the committed '.frontier' hash and with the next
   batch's input_hash; the first input must be the canonical root frontier;
 * per batch: nodes = 2*leaves + |F_out| - |F_in|; totals; final frontier empty;
 * every frontier cell lies inside the canonical root.
usage: kz_chain.py METADATA.json"""
import sys, json, hashlib

def chex(x):
    if x == 0.0:
        return '-0x0p+0' if str(x).startswith('-') else '0x0p+0'
    s = float(x).hex()                      # e.g. 0x1.4cccccccccccdp+0
    sign = ''
    if s[0] == '-':
        sign, s = '-', s[1:]
    mant, exp = s[2:].split('p')
    if '.' in mant:
        a, b = mant.split('.'); b = b.rstrip('0'); mant = a + ('.' + b if b else '')
    e = int(exp)
    return f"{sign}0x{mant}p{'+' if e >= 0 else '-'}{abs(e)}"

def frontier_bytes(n, cells):
    out = [f'KRZYZ_NATIVE_FRONTIER_V1 {n} {len(cells)}']
    for c in cells:
        out.append(c['path'])
        for key in ('lo', 'hi', 'ml', 'mh'):
            out.append(''.join(chex(v) + ' ' for v in c['cell'][key]))
    return (chr(10).join(out) + chr(10)).encode()

if __name__ == '__main__':
    m = json.load(open(sys.argv[1])); n = m['n']; N = n - 1
    PI_UP, TWO_PI_UP = float.fromhex('0x1.921fb54442d19p+1'), float.fromhex('0x1.921fb54442d19p+2')
    root = dict(path='r', cell=dict(lo=[0.0] * N + [1.0], hi=[PI_UP] + [TWO_PI_UP] * (N - 1) + [1.3], ml=[0.0] * N, mh=[8.15] * N))
    prev_hash = hashlib.sha256(frontier_bytes(n, [root])).hexdigest(); prev_count = 1
    tot_nodes = tot_leaves = 0; problems = []
    for i, b in enumerate(m['batches']):
        cm, sm = b['commit'], b['summary']
        if cm['input_hash'] != prev_hash: problems.append(('input hash mismatch', i))
        F = sm['frontier']; fb = frontier_bytes(n, F); h = hashlib.sha256(fb).hexdigest()
        if h != cm['hashes']['.frontier']: problems.append(('frontier hash mismatch', i, h[:12], cm['hashes']['.frontier'][:12]))
        if len(F) != cm['frontier']: problems.append(('frontier count', i))
        if sm['nodes'] != 2 * sm['leaves'] + len(F) - prev_count: problems.append(('node identity', i, sm['nodes'], sm['leaves'], len(F), prev_count))
        if sm['unresolved'] or cm['unresolved'] or sm['error']: problems.append(('unresolved/error', i))
        if sum(sm['counts'].values()) != sm['leaves']: problems.append(('reason counts', i))
        for c in F:
            cl = c['cell']
            ok = all(0 <= cl['lo'][j] <= cl['hi'][j] for j in range(N)) and cl['hi'][0] <= PI_UP and all(cl['hi'][j] <= TWO_PI_UP for j in range(N))
            ok = ok and 1 <= cl['lo'][N] <= cl['hi'][N] <= 1.3 and all(0 <= cl['ml'][j] <= cl['mh'][j] <= 8.15 for j in range(N))
            if not ok: problems.append(('frontier cell outside root', i, c['path']))
        tot_nodes += sm['nodes']; tot_leaves += sm['leaves']; prev_hash = h; prev_count = len(F)
    res = dict(batches=len(m['batches']), total_nodes=tot_nodes, total_leaves=tot_leaves, final_frontier=prev_count,
               totals_match=(tot_nodes == m['nodes'] and tot_leaves == m['leaves']), problems=problems[:20])
    print(json.dumps(res))
