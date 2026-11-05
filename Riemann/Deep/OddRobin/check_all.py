#!/usr/bin/env python3
# check_all.py -- reruns every certificate behind OddRobin/REPORT.md (one core, ~5 min, < 0.5 GB):
#   cd Riemann/Deep/OddRobin && uv run --no-project --with python-flint --with mpmath --with numpy python3 check_all.py
# (1) gap.py: no odd n with N < n < c0 violates (eps-cost enumeration, arb); dfs_gap.py: the same by exhaustive branch-and-bound;
# (2) oddca.py 5e8: every odd CA number from c0 = 67N to x_1 = 5e8 satisfies the odd inequality (arb, exact products);
# (3) transfer.py: odd CA numbers with x_1 in [5e8, X_R] from the RobinCH Robin certificate plus the odd correction,
#     for the three zero inputs of RobinCH (b0 = 0.87338196661320752, 7/8, none).
import subprocess, sys, time, os
here = os.path.dirname(os.path.abspath(__file__))
t0 = time.time(); ok = True
runs = [['gap.py'], ['dfs_gap.py'], ['oddca.py', '5e8'],
        ['transfer.py', '5e8', '1.43e47', '0.87338196661320752'],
        ['transfer.py', '5e8', '1.42e47', '7/8'],
        ['transfer.py', '5e8', '6.46e46', '1']]
for r in runs:
    p = subprocess.run([sys.executable] + r, cwd=here, capture_output=True, text=True)
    lines = [l for l in p.stdout.splitlines() if l.startswith(('PASS', 'FAIL', 'ALL', 'SOME', 'last', 'largest'))]
    print(f'== {" ".join(r)}: exit {p.returncode}, {time.time()-t0:.0f}s')
    for l in lines[-4:]: print('   ' + l)
    if p.returncode != 0:
        ok = False; print(p.stdout[-2000:], p.stderr[-2000:])
print('ALL PASS' if ok else 'SOME CHECK FAILED', f'total {time.time()-t0:.0f}s')
sys.exit(0 if ok else 1)
