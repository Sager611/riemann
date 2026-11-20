#!/bin/bash
# ceilaudit6 run6: quasimode test for window-model angles (frequency leakage of top singular vectors beyond tau1).
cd /Users/adriansager/playground/riemann/.scratch/ceilaudit6
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
UV="uv run --no-project --with numpy --with scipy python3"
for spec in "2,3,5,7,11 1046476.56" "2,3,5,7,11 1026625.04" "2,3,5,7,11 1044028.94" "2,3,5,7,11 3077708.3" "2,3,5,7,11 1040929.1" "2,3,5,7,11 1042479.02" "2,3,5,7,11 2447.58" "2,3,5,7 35053.2"; do
  set -- $spec
  $UV wc_vec.py $1 $2 100 16 2>&1 | python3 -c "
import json,sys
d=json.loads(sys.stdin.read().strip().splitlines()[-1])
print('tc', d['tau_c'], ' '.join(f\"[{m['vec']}{m['j']} k={m['kappa']} f10={m['f10']} f25={m['f25']} out={m['out_tau1']} pk={m['peak_tau']} xm={m['x_mean']} xdeep={m['x_gt_half_depth']}]\" for m in d['modes'] if m['j']<2))"
done
/usr/bin/time -l $UV wc_vec.py 2,3,5,7,11 1046476.56 100 16 2>&1 >/dev/null | rg "maximum resident|real"
echo RUN6 DONE $(date -u +%H:%M:%S)
