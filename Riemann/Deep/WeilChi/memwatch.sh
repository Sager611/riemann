#!/bin/sh
# memwatch.sh: kill any cert_chi/galerkin python process started by OUR queue scripts (ancestor queueNN.sh) above 1.45 GB RSS
anc() { p=$1; i=0; while [ "$p" -gt 1 ] && [ $i -lt 8 ]; do a=$(ps -o args= -p $p 2>/dev/null); case "$a" in *queue[0-9]*.sh*) return 0;; esac; p=$(ps -o ppid= -p $p 2>/dev/null | tr -d ' '); [ -z "$p" ] && return 1; i=$((i+1)); done; return 1; }
while true; do
  ps -Ao pid,rss,args | awk '/python3 (cert_chi|galerkin_chi|ef_check)/ && $2 > 1450000 {print $1, $2}' | while read p r; do
    if anc $p; then echo "$(date -u) killed $p rss $r" >> /Users/adriansager/playground/riemann/Riemann/Deep/WeilChi/logs/memwatch.log; kill $p; fi; done
  pgrep -f 'queue[0-9]+.sh' > /dev/null || exit 0
  sleep 2
done
