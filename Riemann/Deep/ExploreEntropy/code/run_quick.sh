#!/bin/bash
# run_quick.sh (explore32): the fast checks, logs to ../logs (each < 50 s, < 160 MB peak).
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
U="uv run --no-project --with numpy --with scipy --with mpmath --with python-flint python3"
cd "$(dirname "$0")"
L=../logs
for s in t1_identities t2_burnol t3c_scale; do /usr/bin/time -l perl -e 'alarm 900; exec @ARGV' $U $s.py > $L/$s.log 2>&1; done
/usr/bin/time -l perl -e 'alarm 900; exec @ARGV' $U t3_windows.py 40 > $L/t3_windows_N40.log 2>&1
/usr/bin/time -l perl -e 'alarm 1200; exec @ARGV' $U t3b_thresholds.py > $L/t3b_thresholds.log 2>&1
for p in x0 A; do /usr/bin/time -l perl -e 'alarm 900; exec @ARGV' $U t4_pf.py $p > $L/t4_$p.log 2>&1; done
/usr/bin/time -l perl -e 'alarm 900; exec @ARGV' $U t4_pf.py P 3 5 10 20 > $L/t4_P_profiles.log 2>&1
