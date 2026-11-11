#!/bin/bash
# batch_cut.sh: the cut (Nystrom) model for S = {inf,2,..}, t = 1, R = 5, 20, npw = 5, float32; one run at a time.
cd "$(dirname "$0")"
export SP_F32=1
for nm in 1024; do ./run.sh p2_n$nm 2 $nm 5,20 5; ./run.sh p23_n$nm 2,3 $nm 5,20 5; done
for P in 2,3,5 2,3,5,7 2,3,5,7,11; do
  tag=p$(echo $P | tr -d ,)
  for nm in 128 256 512 1024; do
    [ -s ../logs/${tag}_n$nm.json ] && continue
    echo "== $tag $nm $(date -u +%H:%M:%S)"; SP_EIGH=500 ./run.sh ${tag}_n$nm $P $nm 5,20 5
  done
done
echo BATCH_DONE $(date -u +%H:%M:%S)
