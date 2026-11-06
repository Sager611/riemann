#!/bin/bash
# runprod.sh N CHUNKFILE: the production sieve, one chunk at a time (lines "name lo hi", bottom-up), each under a
# hard alarm of 1800 s, single-threaded. A chunk whose log already says DONE is skipped; a failed chunk is noted in
# runs/prod.status and the next one is started (failed chunks are split and rerun by hand). Stops before the next
# chunk when the file STOPSIEVE exists.
cd "$(dirname "$0")"
N=$1; F=$2
echo "start N=$N $(date +%H:%M:%S)" >> runs/prod.status
while read -r name lo hi; do
  [ -z "$name" ] && continue
  [ -f STOPSIEVE ] && { echo "stopped before $name $(date +%H:%M:%S)" >> runs/prod.status; break; }
  if [ -f runs/p_$name.log ] && rg -q '^DONE' runs/p_$name.log; then continue; fi
  OMP_NUM_THREADS=1 perl -e 'alarm 1800; exec @ARGV' ./isieve2 $N $lo $hi runs/p_$name
  rc=$?
  if [ $rc -ne 0 ]; then echo "FAIL $name [$lo,$hi] exit $rc $(date +%H:%M:%S)" >> runs/prod.status; else echo "chunk $name [$lo,$hi] end $(date +%H:%M:%S)" >> runs/prod.status; fi
done < "$F"
echo "end $(date +%H:%M:%S)" >> runs/prod.status
