#!/bin/bash
# needs the bundle: curl -s -L -o .scratch/fodgeaudit/code.zip https://zenodo.org/api/records/24325329/files/code_exact_stechkin_v31.zip/content; unzip into .scratch/fodgeaudit/code
# sequential: one compute process at a time, each under memcap
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
MC=../../../.scratch/memcap.sh; C2=/Users/adriansager/playground/riemann/.scratch/fodgeaudit/code/appendix_c2
/usr/bin/time -l bash $MC 290 300 uv run --no-project --with python-flint python3 consts2.py $C2 > logs/consts2.log 2> logs/consts2.err
/usr/bin/time -l bash $MC 290 900 uv run --no-project --with mpmath python3 bookkeeping.py blocks > logs/bookkeeping_blocks.log 2> logs/bookkeeping_blocks.err
/usr/bin/time -l bash $MC 290 900 uv run --no-project --with mpmath --with numpy python3 bookkeeping.py c4 100 > logs/bookkeeping_c4.log 2> logs/bookkeeping_c4.err
/usr/bin/time -l bash $MC 290 900 uv run --no-project --with mpmath python3 c5test.py 5 200 > logs/c5test.log 2> logs/c5test.err
(cd $C2 && /usr/bin/time -l bash /Users/adriansager/playground/riemann/.scratch/memcap.sh 290 900 uv run --no-project --with python-flint --with numpy python3 propC1p_v6.py 50.0 50.05) > logs/c7_rerun_50.log 2> logs/c7_rerun_50.err
echo ALLDONE >> logs/c7_rerun_50.log
