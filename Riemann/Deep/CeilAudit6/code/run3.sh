#!/bin/bash
# ceilaudit6 run3: N12 lattice test (lat.py) and N13 thin-pocket test (ACLGH2 window model). One process at a time.
cd /Users/adriansager/playground/riemann/.scratch/ceilaudit6
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
UV="uv run --no-project --with numpy --with scipy --with mpmath python3"
echo "== N12 lat"; /usr/bin/time -l $UV lat.py 3000 24 2>&1 | rg -v "^ +[0-9]+  (page|block|vol|invol|signals|messages|instructions|cycles|peak|average|swaps)"
for spec in "2,3,5,7,11 1000000 1100000" "2,3,5,7,11 3000000 3300000" "2,3,5,7 20000 45000"; do
  set -- $spec; P=$1; tag=$(echo $P | tr -d ,)_$2
  echo "== N13 pockets $spec"; /usr/bin/time -l $UV pk_dump.py $P $2 $3 > logs/pk_$tag.json 2> logs/pk_$tag.time; rg "real|maximum" logs/pk_$tag.time
  sel=$($UV pk_select.py logs/pk_$tag.json 2>> logs/run3_sel.txt); tail -1 logs/run3_sel.txt
  for x in $sel; do tc=${x%%:*}; A=${x##*:}
    out=$($UV wc_copy.py $P $tc 100 16 2>/dev/null)
    echo "tc=$tc A=$A $out" | cut -c1-420
  done
done
echo RUN3 DONE $(date -u +%H:%M:%S)
