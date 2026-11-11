#!/bin/bash
# batch_L2.sh: Laurent model checks: depth 20/45 (tmax 100), tmax 400 (f32), and t = 1.5 (tmax 100, 200).
cd "$(dirname "$0")"
for P in 2,3,5,7 2,3,5,7,11; do
  tag=L_p$(echo $P | tr -d ,)
  ./runl.sh ${tag}_t100_d20 $P 100 20 5,20
  ./runl.sh ${tag}_t100_d45 $P 100 45 5,20
done
LA_F32=1 ./runl.sh L_p2357_t400_d30 2,3,5,7 400 30 5,20
for P in - 2 2,3 2,3,5 2,3,5,7 2,3,5,7,11; do
  tag=L_p$(echo $P | tr -d ,-); [ "$P" = "-" ] && tag=L_arch
  for tm in 100 200; do ./runl.sh ${tag}_t${tm}_d30_s1.5 $P $tm 30 5,20 1.5; done
done
echo BATCH_DONE $(date -u +%H:%M:%S)
