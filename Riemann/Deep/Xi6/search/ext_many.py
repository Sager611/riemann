#!/usr/bin/env python3
"""ext_many.py - generate many solutions of a symmetric instance by local search (ls2) and try to extend each by t vertices.
usage: ext_many.py <inst path (without .hyp/.pkl)> <shiftmode 0|1> <nseeds> <t> [seed0]"""
import sys, os, subprocess, pickle, time, json
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
from ext import extend
from sym import write_colouring, check_fast
ip, sm, ns, t = sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]); s0 = int(sys.argv[5]) if len(sys.argv) > 5 else 1
n, C, m = pickle.load(open(ip + '.pkl', 'rb'))
seen = set(); stats = []
for seed in range(s0, s0 + ns):
    o = subprocess.run([os.path.join(HERE, 'ls2'), ip + '.hyp', '30', str(seed), str(sm), str(m)], capture_output=True, text=True)
    if not o.stdout.startswith('SOLVED'): print('seed', seed, 'LS fail', o.stdout.split('\n')[0]); continue
    col = list(map(int, o.stdout.split('\n')[1].split()[1:]))
    if sm: M = [[-1 if a == b else (col[C[a][b] // 3] + C[a][b] % 3) % 3 for b in range(n)] for a in range(n)]
    else: M = [[-1 if a == b else col[C[a][b]] for b in range(n)] for a in range(n)]
    key = tuple(col)
    if key in seen: continue
    seen.add(key)
    M2, nb = extend(M, t)
    stats.append(nb)
    print('seed', seed, 'clauses', nb, 'ext', 'SAT' if M2 else 'UNSAT', flush=True)
    if M2:
        ok = check_fast(M2, 6); print('check', ok)
        fn = os.path.join(HERE, 'sol', f'ext_{os.path.basename(ip)}_s{seed}_t{t}.txt')
        write_colouring(fn, M2, note=f'extension of {os.path.basename(ip)} LS seed {seed} by {t}')
        print('WROTE', fn); break
print('distinct bases', len(seen), 'min/max clauses', min(stats) if stats else None, max(stats) if stats else None)

