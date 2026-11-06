#!/bin/bash
# dl.sh K0 K1 : for slabs k = K0..K1-1, download the 4 Mosunov-Jacobson files of slab k from the LMFDB
# (sequentially, one file at a time, resuming with HTTP ranges) and run mjx on them (one mjx at a time,
# overlapping with the download of the next slab), then delete them.
# Writes dl.log (per file: size, seconds, sha256) and mjx.log. Stops at the first failure.
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
for ((k=$1; k<$2; k++)); do
  [ -f STOP ] && break   # ImqClass2: graceful stop
  for c in cl3mod8 cl7mod8 cl4mod16 cl8mod16; do
    f=dl/$c.$k.gz
    q="$URL?filenamebase=$c&k=$k&Fetch=fetch"
    curl -s -I -L -b human=1 -o dl/hdr.txt "$q"
    exp=$(rg -i '^content-length' dl/hdr.txt | tail -1 | tr -d '\r' | awk '{print $2}')
    rm -f dl/hdr.txt
    if [ -z "$exp" ]; then echo "FAIL no content-length $c.$k" >> dl.log; finish; exit 1; fi
    t0=$(date +%s); tries=0
    while :; do
      sz=$(stat -f %z "$f" 2>/dev/null || echo 0)
      [ "$sz" -eq "$exp" ] && break
      [ "$sz" -gt "$exp" ] && rm -f "$f"
      tries=$((tries+1)); if [ $tries -gt 40 ]; then echo "FAIL download $c.$k" >> dl.log; finish; exit 1; fi
      curl -s -L -b human=1 -C - --max-time 900 -o "$f" "$q"
      [ $tries -gt 1 ] && sleep 3
    done
    t1=$(date +%s)
    echo "$c.$k size $exp secs $((t1-t0)) tries $tries sha256 $(shasum -a 256 "$f" | awk '{print $1}')" >> dl.log
  done
  finish
  nice -n 10 perl -e 'alarm 900; exec @ARGV' ./mjx state.bin $k dl >> mjx.log 2>&1 &
  pid=$!; pk=$k
done
finish
echo "ALL DONE $1..$2 $(date +%H:%M:%S)" >> dl.log
