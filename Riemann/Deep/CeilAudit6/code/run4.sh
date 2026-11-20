#!/bin/bash
# ceilaudit6 run4: N13 artefact check of kappa_top ~ 0.16 at thin pockets: vary window width/depth and centre offset.
cd /Users/adriansager/playground/riemann/.scratch/ceilaudit6
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
UV="uv run --no-project --with numpy --with scipy python3"
P=2,3,5,7,11
for spec in "1046476.56 100 16" "1046476.56 200 16" "1046476.56 100 32" "1046516.56 100 16" "1046436.56 100 16" "1046476.56 50 16" "1033777.07 200 16" "1040929.1 200 16" "1000000 100 16" "1500000 200 16" "4500000 100 16"; do
  set -- $spec
  out=$($UV wc_copy.py $P $1 $2 $3 2>/dev/null | python3 -c "import json,sys; d=json.load(sys.stdin); print('N',d['N'],'Nsc',d['N_sc_window'],'sumk2',d['sum_k2'],'top',d['top'][:5])")
  echo "tc=$1 tmax=$2 depth=$3 $out"
done
/usr/bin/time -l $UV wc_copy.py $P 1046476.56 200 32 2>&1 >/dev/null | rg "maximum resident|real"
echo RUN4 DONE $(date -u +%H:%M:%S)
