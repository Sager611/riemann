#!/bin/sh
# rsslog.sh: log RSS (MB) of the running cert_zeta.py process every 5 s (to see which phase peaks); exits after 60 s without one
idle=0
while [ $idle -lt 12 ]; do
  l=$(ps -Ao pid,etime,rss,args | awk '/python3 cert_zeta/ && !/awk/ && !/uv run/ && !/time -l/ {print $1, $2, int($3/1024) " MB", $6, $7}')
  if [ -n "$l" ]; then echo "$(date -u +%H:%M:%S) $l" >> logs/rss.log; idle=0; else idle=$((idle+1)); fi
  sleep 5
done
