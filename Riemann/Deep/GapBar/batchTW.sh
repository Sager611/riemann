#!/bin/bash
# sequential twisted-kernel barrier certificates; watchdog kills any job above 290 MB RSS
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
run() {  # kappa N r vlo vhi shrink [D]
  tag="r$3_$1"
  perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l uv run --no-project --with python-flint --with numpy --with scipy python3 runcerttw.py "$@" > logs/certTW_$tag.txt 2>&1 &
  W=$!
  while kill -0 $W 2>/dev/null; do
    P=$(pgrep -f "python3 runcerttw.py $1 $2 $3" | head -1)
    if [ -n "$P" ]; then R=$(ps -o rss= -p $P 2>/dev/null | tr -d ' '); if [ -n "$R" ] && [ "$R" -gt 296960 ]; then kill $P; echo "KILLED $tag rss=$R" >> logs/batchTW.log; fi; fi
    sleep 2
  done
  echo "$(date -u +%H:%M:%S) $tag: $(rg -o 'RESULT.*' logs/certTW_$tag.txt)" >> logs/batchTW.log
}
run 3.30 72 1 0.8 1.7 0.3
run 3.35 72 2 0.8 1.7 0.5
run 3.5 72 2 0.8 1.7 0.7
run 3.5 72 1 0.8 1.7 0.7
run 4 48 2 0.8 1.7 0.7 76
run 4 72 1 0.8 1.7 0.7
run 5 40 2 0.8 1.7 0.7 76
run 3.5 48 3 0.8 1.8 0.7 66
echo "$(date -u +%H:%M:%S) DONE" >> logs/batchTW.log
