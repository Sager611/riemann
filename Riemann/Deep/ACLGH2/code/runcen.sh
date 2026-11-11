#!/bin/bash
# run tag args...: one single-threaded run of centroid.py under alarm 1800 s and /usr/bin/time -l; RSS polled every 0.2 s,
# the job is killed above 270 MB (cap 300 MB); kills are logged in ../logs/<tag>.kill.
tag=$1; shift
cd "$(dirname "$0")"
L=../logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
rm -f $L/$tag.kill
/usr/bin/time -l perl -e 'alarm 1800; exec @ARGV' uv run --no-project --with numpy --with scipy python3 centroid.py "$@" > $L/$tag.json 2> $L/$tag.time &
tp=$!
peak=0
while kill -0 $tp 2>/dev/null; do
  for pid in $(pgrep -f "python3 centroid.py"); do
    m=$(ps -o rss= -p $pid 2>/dev/null | tr -d ' ')
    if [ -n "$m" ]; then
      [ "$m" -gt "$peak" ] && peak=$m
      if [ "$m" -gt 270000 ]; then echo "KILL rss=$m KB" >> $L/$tag.kill; kill $pid; fi
    fi
  done
  sleep 0.2
done
echo "polled_peak_rss_KB=$peak" >> $L/$tag.kill.tmp; cat $L/$tag.kill.tmp >> $L/$tag.time; rm -f $L/$tag.kill.tmp
[ -s $L/$tag.kill ] && cat $L/$tag.kill
awk '/maximum resident set size/{printf "maxrss_MB=%.0f ", $1/1048576} /real/{printf "real_s=%s ", $1} /polled/{print $0}' $L/$tag.time
