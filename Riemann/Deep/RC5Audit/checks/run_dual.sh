#!/bin/bash
# run_dual.sh X Nlist... : the dual form Q'(b) of Corollary P over all orders (Lemma P blocks), detector sequences
X=$1; shift
V=$(python3 -c "import math;print(max(2,int(math.sqrt($X)/2)))")
SPEC=/tmp/rc5dspec_$X.txt; : > $SPEC
for N in "$@"; do
  echo "$N 0.8333 $V 0 1 0" >> $SPEC
  echo "$N 0.8333 $V 0 0 0" >> $SPEC
  echo "$N 0.55 $V 0 1 0" >> $SPEC
  echo "$N 0.8333 $V 0 1 1" >> $SPEC
done
/usr/bin/time -l perl -e 'alarm 1700; exec @ARGV' ./fam $X 2 < $SPEC 2> /tmp/rc5dtime_$X.txt
rg 'maximum resident|real' /tmp/rc5dtime_$X.txt | sed 's/^ */# /'
rm -f $SPEC /tmp/rc5dtime_$X.txt
