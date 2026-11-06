#!/bin/bash
# dlrun.sh K1: stream the Mosunov-Jacobson slabs from the state's next slab up to K1 with dl2.sh, restarting after
# failures (at most 40 restarts; the failing slab's files are deleted and downloaded again). Stops when the file
# STOP exists (dl.sh checks it before each slab). Run under perl alarm.
cd "$(dirname "$0")"
K1=$1
nk() { if [ -f state.bin ]; then od -A n -t u8 -j 16 -N 8 state.bin | tr -d ' '; else echo 0; fi; }
for ((r=0; r<40; r++)); do
  k=$(nk)
  [ "$k" -ge "$K1" ] && break
  [ -f STOP ] && break
  echo "dlrun start k=$k try $r $(date +%H:%M:%S)" >> dl.log
  ./dl2.sh $k $K1 && break
  k=$(nk)
  rm -f dl/cl3mod8.$k.gz dl/cl7mod8.$k.gz dl/cl4mod16.$k.gz dl/cl8mod16.$k.gz
  sleep 30
done
echo "dlrun end next_k=$(nk) $(date +%H:%M:%S)" >> dl.log
