#!/bin/bash
# run_det.sh X Nlist... : detector-sequence family sums (one process, hard alarm)
X=$1; shift
V=$(python3 -c "import math;print(max(2,int(math.sqrt($X)/2)))")
SPEC=/tmp/rc5spec_$X.txt; : > $SPEC
for N in "$@"; do
  echo "$N 0.8333 $V 0 1 0" >> $SPEC   # main: sigma=5/6, V=sqrt(X)/2, pole killed
  echo "$N 0.8333 $V 0 0 0" >> $SPEC   # same, pole present (no ell-term)
  echo "$N 0.55 $V 0 1 0" >> $SPEC     # sigma near 1/2
  echo "$N 0.55 2 0 1 0" >> $SPEC      # no mollifier (a_V = 1 on admissible n)
  echo "$N 0.8333 $V 5 1 0" >> $SPEC   # t' = 5 (complex)
  echo "$N 0.8333 $V 0 1 1" >> $SPEC   # d/dt' b at t'=0
  echo "$N 0.9 2 0 0 0" >> $SPEC       # pole present, no mollifier
done
/usr/bin/time -l perl -e 'alarm 1700; exec @ARGV' ./fam $X 0 < $SPEC 2> /tmp/rc5time_$X.txt
rg 'maximum resident|real' /tmp/rc5time_$X.txt | sed 's/^ */# /'
rm -f $SPEC /tmp/rc5time_$X.txt
