#!/bin/bash
# mon.sh OUT CMD...: run CMD under /usr/bin/time with a 0.5 s RSS poller (kill above 270 MB); one process at a time
OUT=$1; shift
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
/usr/bin/time -l perl -e 'alarm 900; exec @ARGV' "$@" > $OUT.json 2> $OUT.time &
TP=$!; PEAK=0
while kill -0 $TP 2>/dev/null; do
  for P in $(pgrep -f 'python3 vlow.py'); do
    R=$(ps -o rss= -p $P 2>/dev/null | tr -d ' '); [ -z "$R" ] && continue
    [ $R -gt $PEAK ] && PEAK=$R
    if [ $R -gt 276480 ]; then kill $P; echo "KILLED rss=$R KB" >> $OUT.time; fi
  done
  sleep 0.5
done
echo "polled_peak_KB=$PEAK" >> $OUT.time
