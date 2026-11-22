#!/bin/bash
# PROPOSAL, not in use (LeanLi agent, 22 Nov 2026): .scratch/memcap.sh with the cap on the physical
# memory footprint (Adrian's check `top -l 1 -pid <pid> -stats mem`) instead of RSS (`ps -o rss`, which also
# counts the clean, file-backed, mmapped .olean pages). Same usage: bash memcap_fp.sh CAP_MB LIMIT_S cmd [args...]
cap=$1; lim=$2; shift 2
d=$(cd "$(dirname "$0")" && pwd)
fpm=/tmp/fpmem_$(id -u)
[ -x "$fpm" ] || cc -O2 -o "$fpm" "$d/fpmem.c" || exit 2
"$@" &
root=$!
peak=0; peakrss=0; start=$SECONDS; status=""
tree() { echo "$1"; for c in $(pgrep -P "$1" 2>/dev/null); do tree "$c"; done; }
while kill -0 "$root" 2>/dev/null; do
  read -r mb rss < <("$fpm" $(tree "$root"))
  [ "${mb:-0}" -gt "$peak" ] && peak=$mb
  [ "${rss:-0}" -gt "$peakrss" ] && peakrss=$rss
  if [ "${mb:-0}" -gt "$cap" ]; then status="memcap_fp: killed at footprint $mb MB > $cap MB"; break; fi
  if [ $((SECONDS - start)) -gt "$lim" ]; then status="memcap_fp: killed after $lim s"; break; fi
  sleep 0.2
done
if [ -n "$status" ]; then
  for p in $(tree "$root"); do kill -9 "$p" 2>/dev/null; done
  wait "$root" 2>/dev/null
  echo "$status (peak footprint $peak MB, peak RSS $peakrss MB)" >&2
  exit 137
fi
wait "$root"; rc=$?
echo "memcap_fp: peak footprint $peak MB, peak RSS $peakrss MB, rc $rc, $((SECONDS - start)) s" >&2
exit $rc
