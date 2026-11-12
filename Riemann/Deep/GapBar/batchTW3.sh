#!/bin/bash
# common grid {3.335, 3.35, 3.5, 4, 5} for r <= 3; watchdog kills any job above 280 MB RSS (1 s sampling)
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
run() {  # kappa N r vlo vhi shrink [D]
  tag="r$3_$1"
  perl -e 'alarm 900; exec @ARGV' /usr/bin/time -l uv run --no-project --with python-flint --with numpy --with scipy python3 runcerttw.py "$@" > logs/certTW_$tag.txt 2>&1 &
  W=$!
  while kill -0 $W 2>/dev/null; do
    P=$(pgrep -f "python3 runcerttw.py $1 $2 $3" | head -1)
    if [ -n "$P" ]; then R=$(ps -o rss= -p $P 2>/dev/null | tr -d ' '); if [ -n "$R" ] && [ "$R" -gt 286720 ]; then kill $P; echo "KILLED $tag rss=$R" >> logs/batchTW.log; fi; fi
    sleep 1
  done
  echo "$(date -u +%H:%M:%S) $tag: $(rg -o 'RESULT.*' logs/certTW_$tag.txt)" >> logs/batchTW.log
}
run 3.35 72 1 0.8 1.7 0.7
run 5 72 1 0.8 1.7 0.7
run 3.35 48 3 0.8 1.8 0.7 66
run 4 48 3 0.8 1.8 0.7 76
run 5 40 3 0.8 1.8 0.7 84
echo "$(date -u +%H:%M:%S) DONE3" >> logs/batchTW.log
