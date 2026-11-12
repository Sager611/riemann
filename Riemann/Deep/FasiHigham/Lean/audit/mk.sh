#!/bin/bash
# usage: bash .scratch/fasilean/mk.sh <file.lean> [olean-out]   (run from repo root)
f=$1; o=$2; name=$(basename $f .lean)
LP=.scratch/fasilean/olean:$(cat .scratch/fasilean/leanpath.txt)
args="-j1"; [ -n "$o" ] && { mkdir -p $(dirname $o); args="$args -o $o"; }
t0=$(date +%s)
LEAN_PATH=$LP OMP_NUM_THREADS=1 bash .scratch/memcap.sh 290 1800 /usr/bin/time -l ~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean $args $f > .scratch/fasilean/log_$name.txt 2> .scratch/fasilean/err_$name.txt
rc=$?
fp=$(rg -o '([0-9]+)  peak memory footprint' -r '$1' .scratch/fasilean/err_$name.txt | head -1)
pk=$(rg -o 'memcap: .*' .scratch/fasilean/err_$name.txt)
echo "$(date -u +%H:%M:%S) $name rc=$rc $(( $(date +%s) - t0 ))s footprint=$(( ${fp:-0} / 1048576 ))MiB $pk" | tee -a .scratch/fasilean/peaks.log
head -c 3000 .scratch/fasilean/log_$name.txt

