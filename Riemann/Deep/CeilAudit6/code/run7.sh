#!/bin/bash
cd /Users/adriansager/playground/riemann/.scratch/ceilaudit6
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
UV="uv run --no-project --with numpy --with scipy python3"
for tc in 1000000 1000013.7 1500000 4500000 1033777.07 1046436.56 1046516.56; do
  $UV wc_vec.py 2,3,5,7,11 $tc 100 16 2>&1 | python3 -c "
import json,sys
d=json.loads(sys.stdin.read().strip().splitlines()[-1])
print('tc', d['tau_c'], ' '.join(f\"[{m['vec']}{m['j']} k={m['kappa']} f10={m['f10']} pk={m['peak_tau']} out={m['out_tau1']}]\" for m in d['modes'] if m['vec']=='v'))"
done
