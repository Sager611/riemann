#!/bin/bash
# Rerun the four HSW-dependent certificates with 9.4926 in place of HSW's misprinted 9.3675 (copies in /tmp/stc2).
set -u
R=/Users/adriansager/playground/riemann
MC="bash $R/.scratch/memcap.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
rm -rf /tmp/stc2; mkdir -p /tmp/stc2/logs
for d in Ingham3 Li200k RobinHeight RobinCH; do cp -R $R/Riemann/Deep/$d /tmp/stc2/; done
python3 - <<'EOF'
import re, glob
for p in glob.glob('/tmp/stc2/*/*.py'):
    s = open(p).read(); s0 = s
    s = s.replace("arb('9.3675')", "arb('9.4926')").replace('9.3675', '9.4926').replace('r(93675, 10000)', 'r(94926, 10000)')
    if s != s0: open(p, 'w').write(s); print('patched', p)
EOF
rg -n '9\.3675|93675' /tmp/stc2 --glob '*.py' && echo "LEFTOVER 9.3675" || echo "no 9.3675 left in .py"
run() { local name=$1 lim=$2; shift 2; echo "== $name start $(date -u +%H:%M:%S)"; /usr/bin/time -l $MC 290 $lim "$@" > /tmp/stc2/logs/$name.log 2>&1; echo "rc=$? end $(date -u +%H:%M:%S)"; rg 'memcap:|maximum resident' /tmp/stc2/logs/$name.log; }
cd /tmp/stc2/Li200k && run Li200k 300 uv run --no-project --with python-flint python3 verify.py
cd /tmp/stc2/RobinCH && run RobinCH_check_all2 900 uv run --no-project --with python-flint --with mpmath python3 check_all2.py
cd /tmp/stc2/Ingham3 && run Ingham3 600 uv run --no-project --with python-flint --with mpmath --with sympy python3 cert.py
cd /tmp/stc2/RobinHeight && run RobinHeight_check_all 900 uv run --no-project --with python-flint --with mpmath python3 check_all.py
cd /tmp/stc2/RobinHeight && run RobinHeight_certify 600 uv run --no-project --with python-flint --with mpmath python3 certify.py 7/8 0.87338196661320752 1
echo "== all done $(date -u +%H:%M:%S)"
