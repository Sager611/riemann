#!/bin/bash
# usage: bash .scratch/memcap.sh LIMIT_MB TIMEOUT_S cmd args...
# Runs cmd; about every 3 s reads top's MEM (physical footprint, the cap's measure) for cmd and all its descendants;
# kills the whole tree if any one process exceeds LIMIT_MB or after TIMEOUT_S. Prints the peak to stderr; returns cmd's exit code.
lim=$1; to=$2; shift 2
"$@" &
root=$!
t0=$(date +%s); peak=0; killed=""
tree() { local p; for p in "$@"; do echo $p; tree $(pgrep -P $p); done; }
while kill -0 $root 2>/dev/null; do
  pids=$(tree $root | sort -u)
  args=""; for p in $pids; do args="$args -pid $p"; done
  for m in $(top -l 1 -stats pid,mem $args 2>/dev/null | awk '$1 ~ /^[0-9]+$/ {print $2}'); do
    v=$(echo "$m" | tr -d '+-' | awk '/K$/{printf "%d", ($0+0)/1024; next} /M$/{printf "%d", $0+0; next} /G$/{printf "%d", ($0+0)*1024; next} {printf "%d", ($0+0)/1048576}')
    [ "$v" -gt "$peak" ] && peak=$v
  done
  [ "$peak" -gt "$lim" ] && killed="MEM>${lim}MB"
  [ $(( $(date +%s) - t0 )) -gt "$to" ] && killed="TIME>${to}s"
  if [ -n "$killed" ]; then for p in $pids; do kill -9 $p 2>/dev/null; done; break; fi
  sleep 2
done
wait $root 2>/dev/null; rc=$?
echo "memcap: peak ${peak} MB${killed:+ KILLED $killed} rc=$rc" >&2
exit $rc
