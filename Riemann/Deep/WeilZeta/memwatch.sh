#!/bin/sh
# memwatch.sh: kill any cert_zeta/galerkin python process of OUR queue scripts (ancestor queueNN.sh in WeilZeta) above 1.45 GB RSS.
# Fixed 2026-11-05 21:44: anc() used the caller's variable p (sh has no local variables), so the first version killed the ancestor
# queue script instead of the python process (logs/memwatch.log). The copy in ../WeilChi has the same bug.
anc() { q=$1; i=0; while [ "$q" -gt 1 ] && [ $i -lt 8 ]; do a=$(ps -o args= -p $q 2>/dev/null); case "$a" in *queue[0-9]*.sh*) return 0;; esac; q=$(ps -o ppid= -p $q 2>/dev/null | tr -d ' '); [ -z "$q" ] && return 1; i=$((i+1)); done; return 1; }
while true; do
  ps -Ao pid,rss,args | awk '/python3 (cert_zeta|galerkin_zeta)/ && !/uv run/ && !/time -l/ && $2 > 1450000 {print $1, $2}' | while read pid rss; do
    if anc $pid; then echo "$(date -u) killed $pid rss $rss" >> /Users/adriansager/playground/riemann/Riemann/Deep/WeilZeta/logs/memwatch.log; kill $pid; fi; done
  pgrep -f 'queue[0-9]+.sh' > /dev/null || exit 0
  sleep 2
done
