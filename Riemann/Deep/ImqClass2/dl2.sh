#!/bin/bash
# dl2.sh K0 K1: as dl.sh (ImqClass), but the four files of a slab are fetched over two connections at a time
# (cl3mod8 with cl7mod8, then cl4mod16 with cl8mod16), each resumed with HTTP ranges until its size equals the
# server's Content-Length; then mjx runs on the slab (one mjx at a time, overlapping with the next slab's
# download) and the slab's files are deleted. At most two slabs on disk. Writes dl.log (per file: size, seconds,
# tries, sha256) and mjx.log. Stops at the first failure, or before a slab when the file STOP exists.
set -u
cd "$(dirname "$0")"
URL='https://beta.lmfdb.org/NumberField/QuadraticImaginaryClassGroups'
mkdir -p dl
pid=""; pk=""
finish() { # wait for the running mjx (slab pk), delete its files
  if [ -n "$pid" ]; then
    if ! wait $pid; then echo "FAIL mjx slab $pk" >> dl.log; exit 1; fi
    rm -f dl/cl3mod8.$pk.gz dl/cl7mod8.$pk.gz dl/cl4mod16.$pk.gz dl/cl8mod16.$pk.gz
    echo "slab $pk done $(date +%H:%M:%S)" >> dl.log
    pid=""
  fi
}
getf() { # getf C K: download dl/C.K.gz completely; result line in dl/res.C.K; status 0 on success
  local c=$1 k=$2; local f=dl/$c.$k.gz q="$URL?filenamebase=$c&k=$k&Fetch=fetch" hdr=dl/hdr.$c.$k
  curl -s -I -L -b human=1 -o $hdr "$q"
  local exp; exp=$(rg -i '^content-length' $hdr | tail -1 | tr -d '\r' | awk '{print $2}')
  rm -f $hdr
  if [ -z "$exp" ]; then echo "FAIL no content-length $c.$k" > dl/res.$c.$k; return 1; fi
  local t0; t0=$(date +%s); local tries=0 sz
  while :; do
    sz=$(stat -f %z "$f" 2>/dev/null || echo 0)
    [ "$sz" -eq "$exp" ] && break
    [ "$sz" -gt "$exp" ] && rm -f "$f"
    tries=$((tries+1)); if [ $tries -gt 40 ]; then echo "FAIL download $c.$k" > dl/res.$c.$k; return 1; fi
    curl -s -L -b human=1 -C - --max-time 900 -o "$f" "$q"
    [ $tries -gt 1 ] && sleep 3
  done
  echo "$c.$k size $exp secs $(( $(date +%s) - t0 )) tries $tries sha256 $(shasum -a 256 "$f" | awk '{print $1}')" > dl/res.$c.$k
  return 0
}
for ((k=$1; k<$2; k++)); do
  [ -f STOP ] && break
  getf cl3mod8 $k & a=$!; getf cl7mod8 $k & b=$!
  wait $a; ra=$?; wait $b; rb=$?
  getf cl4mod16 $k & a=$!; getf cl8mod16 $k & b=$!
  wait $a; rc=$?; wait $b; rd=$?
  for c in cl3mod8 cl7mod8 cl4mod16 cl8mod16; do cat dl/res.$c.$k >> dl.log 2>/dev/null; rm -f dl/res.$c.$k; done
  if [ $ra -ne 0 ] || [ $rb -ne 0 ] || [ $rc -ne 0 ] || [ $rd -ne 0 ]; then finish; exit 1; fi
  finish
  nice -n 10 perl -e 'alarm 900; exec @ARGV' ./mjx state.bin $k dl >> mjx.log 2>&1 &
  pid=$!; pk=$k
done
finish
echo "ALL DONE $1..$2 $(date +%H:%M:%S)" >> dl.log
