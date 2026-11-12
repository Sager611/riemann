#!/bin/bash
# run_eig.sh X N Q0 : export sparse class data for admissible n in [N,2N), then top eigenvalues (one process at a time)
X=$1; N=$2; Q0=$3; P=/tmp/rc5e_${X}_${N}
perl -e 'alarm 900; exec @ARGV' ./fam $X 1 $N $((2*N)) $P > /dev/null || { echo "fam failed X=$X N=$N"; exit 1; }
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 /usr/bin/time -l perl -e 'alarm 1200; exec @ARGV' uv run --no-project --with numpy --with scipy python3 eigs.py $P $X $N $Q0 2> /tmp/rc5t.txt
echo "#   mem(eigs) $(rg 'maximum resident' /tmp/rc5t.txt | awk '{printf "%.0f MB", $1/1048576}') $(rg ' real' /tmp/rc5t.txt | awk '{print $1" s"}')"
rm -f $P.trip $P.cls $P.gen /tmp/rc5t.txt
