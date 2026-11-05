#!/usr/bin/env python3
"""verify.py -- re-runs every numeric step of PROOF.md, one single-threaded process at a time.

  python3 verify.py                all proof steps A-G (about 1.5 h on one core)
  python3 verify.py --sanity       also the numerical cross-checks S1-S4 against exact values (m = 203)
  python3 verify.py --no-small-m   skip step G (odd m <= 201: ../FasiHigham/arb_cert.py diag 201, ~10 min)

Needs Python 3 with python-flint >= 0.9, mpmath, sympy:  pip install python-flint mpmath sympy
Writes tab203.json (Good-lemma tables) and lg203.pkl (LG error-bound grid) next to this file.
Each step runs as a child process with a wall-clock limit; the step passes only if its success
marker appears in the output.  Exit code 0 iff every step passed.
"""
import os, re, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__))
PY = sys.executable
ENV = dict(os.environ)
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS', 'NUMEXPR_NUM_THREADS'):
    ENV[k] = '1'
NU0, TAB, LG = '203.5', 'tab203.json', 'lg203.pkl'
def fails0(o): return re.search(r'boxes \d+ fails 0 ', o) is not None
def arb_ok(o):
    pos = re.findall(r"^(\d+) \('POS', \d+\)", o, re.M)
    return sorted(int(x) for x in pos) == list(range(3, 202, 2))
STEPS = [
    ('A  identities (sympy / exact / mpmath)', ['identities.py'], 900, lambda o: 'IDENTITIES ALL OK' in o),
    ('B  c_n > 0 for n <= 12 and all real m >= 203', ['smalln.py', '203', '12'], 900, lambda o: 'ALL small-n positive: True' in o),
    ('C  Good-lemma tables at nu0 = 203.5', ['mktab.py', NU0, TAB], 1800, lambda o: o.count(' boxes ') >= 6),
    ('D  LG error-bound grid', ['assembly.py', NU0, TAB, LG, '0', '0', '1'], 1800, lambda o: 'lg grid built' in o),
]
for (a, b) in [(0, 12), (12, 24), (24, 36), (36, 48), (48, 60)]:
    STEPS.append(('E  tau-sweep [%d, %d], n >= 8' % (a, b), ['assembly.py', NU0, TAB, LG, str(a), str(b), '0.02'], 1800, fails0))
STEPS.append(('F  tau >= 60', ['large_tau.py', NU0, TAB, LG, '60'], 1800, lambda o: 'LARGE-TAU ALL OK' in o))
STEPS.append(('G  odd 3 <= m <= 201, all n (Arb)', [os.path.join('..', 'FasiHigham', 'arb_cert.py'), 'diag', '201'], 1800, arb_ok))
SANITY = [
    ('S1 lower bound <= exact c_n, m = 203, n <= 4000', ['sanity.py', '203', '4000', '3'], 1800, lambda o: 'violations 0 ' in o),
    ('S2 Good bound <= exact J_l (quadrature), m = 203', ['check_good2.py', '203', '40:0,150:0,300:0,450:0,420:1,500:1,620:1,900:1,900:2,1300:3,2400:4,2400:5', '1600'], 1800,
     lambda o: o.count(' OK') == 12 and 'VIOLATION' not in o),
    ('S3 I_L and Neg bounds >= numeric, m = 203', ['check_neg_IL.py', '203', '20,407,2000', '400'], 1800, lambda o: o.count('OK') == 6 and 'VIOL' not in o),
    ('S4 crude bounds >= numeric, m = 203', ['check_bounds.py', '203', '420:1,500:2,900:3,1300:5,2000:2,2400:9,3000:1'], 1800, lambda o: 'VIOLATION' not in o and ' OK' in o),
]
def main():
    steps = list(STEPS)
    if '--no-small-m' in sys.argv: steps = [s for s in steps if not s[0].startswith('G')]
    if '--sanity' in sys.argv: steps += SANITY
    for f in (TAB, LG):
        p = os.path.join(HERE, f)
        if os.path.exists(p): os.remove(p)
    allok = True; t00 = time.time()
    for (name, args, tmax, ok) in steps:
        t0 = time.time()
        try:
            r = subprocess.run([PY] + args, cwd=HERE, env=ENV, capture_output=True, text=True, timeout=tmax)
            out = r.stdout + r.stderr; good = (r.returncode == 0 or name.startswith('A')) and ok(out)
        except subprocess.TimeoutExpired as e:
            out = (e.stdout or b'').decode() if isinstance(e.stdout, bytes) else (e.stdout or ''); good = False
            out += '\n[timeout after %d s]' % tmax
        allok = allok and good
        tail = [l for l in out.strip().splitlines() if l.strip()][-3:]
        print('%-4s %-52s %6.0fs' % ('PASS' if good else 'FAIL', name, time.time() - t0), flush=True)
        for l in tail: print('       | ' + l[:150], flush=True)
        if not good and not name.startswith('S'):
            print('       stopping: a proof step failed', flush=True); break
    print('ALL STEPS PASSED (%.0f s)' % (time.time() - t00) if allok else 'VERIFICATION FAILED', flush=True)
    sys.exit(0 if allok else 1)
if __name__ == '__main__':
    main()
