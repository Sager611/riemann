#!/bin/bash
# batch_C.sh: counting function of the Laurent model vs tmax (eigenvalues only, float32, N <= 6010).
cd "$(dirname "$0")"
for P in 2 2,3 2,3,5 2,3,5,7 2,3,5,7,11; do
  tag=C_p$(echo $P | tr -d ,)
  for tm in 400 800 1180; do ./runc.sh ${tag}_t${tm}_d16 $P $tm 16; done
done
for P in 2,3,5,7,11,13 2,3,5,7,11,13,17,19,23 2,3,5,7,11,13,17,19,23,29,31,37,41,43,47; do
  top=${P##*,}
  for tm in 400 785; do ./runc.sh C_P${top}_t${tm}_d24 $P $tm 24; done
done
echo BATCH_DONE $(date -u +%H:%M:%S)
