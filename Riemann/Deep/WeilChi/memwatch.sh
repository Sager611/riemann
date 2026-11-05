#!/bin/sh
# memwatch.sh: kill any cert_chi/galerkin python process started by OUR queue scripts (ancestor queueNN.sh) above 1.45 GB RSS
# Fixed 2026-11-05: anc() used the caller's variable p (sh has no local variables), so it killed the ancestor queue script (fix from ../WeilZeta).
anc() { q=$1; i=0; while [ "$q" -gt 1 ] && [ $i -lt 8 ]; do a=$(ps -o args= -p $q 2>/dev/null); case "$a" in *queue[0-9]*.sh*) return 0;; esac; q=$(ps -o ppid= -p $q 2>/dev/null | tr -d ' '); [ -z "$q" ] && return 1; i=$((i+1)); done; return 1; }
while true; do
  ps -Ao pid,rss,args | awk '/python3 (cert_chi|galerkin_chi|ef_check)/ && !/uv run/ && !/time -l/ && $2 > 1450000 {print $1, $2}' | while read pid rss; do
    if anc $pid; then echo "$(date -u) killed $pid rss $rss" >> /Users/adriansager/playground/riemann/Riemann/Deep/WeilChi/logs/memwatch.log; kill $pid; fi; done
  pgrep -f 'queue[0-9]+.sh' > /dev/null || exit 0
  sleep 2
done
