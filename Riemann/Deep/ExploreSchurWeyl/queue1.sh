#!/bin/bash
# queue (explore50): waits for colgen N=10, then runs one process at a time
cd /Users/adriansager/playground/riemann/Riemann/Deep/ExploreSchurWeyl
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
while pgrep -f 'colgen.py 10 P3' > /dev/null; do sleep 5; done
MC=../../../.scratch/memcap.sh
bash $MC 290 600 uv run --no-project --with numpy --with scipy python3 relax_lp.py 10 320 >> logs/relax_lp.log 2>&1
for a in "3 24 tri" "4 48 tri" "4 48 pair" "4 64 tri" "4 64 pair" "6 72 tri" "6 72 pair" "8 96 tri" "8 96 pair" "10 120 tri" "10 120 pair" "12 144 tri" "12 144 pair" "16 192 tri" "16 192 pair"; do
  set -- $a
  bash $MC 290 900 uv run --no-project --with numpy --with scipy python3 tri_sym.py $1 $2 $3 >> logs/tri_sym.log 2>&1
done
echo QUEUE_DONE >> logs/tri_sym.log
