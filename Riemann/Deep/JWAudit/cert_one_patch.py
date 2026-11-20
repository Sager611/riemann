# JWAudit: make a single-cell, verbose copy of SavvaK/certify6.py (never edits the original).
# Usage: python3 cert_one_patch.py OUTDIR ; then in OUTDIR:
#   JW_CELL=917 uv run --no-project --with python-flint --with sympy python3 cert_one.py 80 2.3 1000 1.33
# prints the certified (a, t, w, delta) box endpoints, the X_j enclosures and C6/Cfin/budget/Cloc of that cell (logs/cert_one_cell917.log).
import sys, os
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'SavvaK', 'certify6.py')).read()
src = src.replace("for i in range(NEPS):\n    stack", "import os\nfor i in [int(os.environ['JW_CELL'])]:\n    stack", 1)
src = src.replace("    a, t, w, de = br\n", "    a, t, w, de = br\n    print('JWBOX cell', i, 'eps', el, eh, 'a,t,w,delta =', [(x.lower().str(12, radius=False), x.upper().str(12, radius=False)) for x in br])\n    print('JWX', [x.str(15) for x in X])\n", 1)
src = src.replace("    Cfin = up(C6 * (LAM * KR * sb + omega) ** 6 / 720)\n", "    Cfin = up(C6 * (LAM * KR * sb + omega) ** 6 / 720)\n    print('JWC6 C6', C6.str(8), 'Cfin', Cfin.str(8), 'B/n', b.str(10), 'tot6', tot6.str(8), 'budget', totr.str(6), 'Cloc', Cloc.str(8))\n", 1)
assert src.count('JWBOX') == 1 and src.count('JWC6') == 1 and src.count('JW_CELL') == 1
out = sys.argv[1] if len(sys.argv) > 1 else '/tmp'
os.makedirs(out, exist_ok=True)
open(os.path.join(out, 'cert_one.py'), 'w').write(src)
print('wrote', os.path.join(out, 'cert_one.py'))
