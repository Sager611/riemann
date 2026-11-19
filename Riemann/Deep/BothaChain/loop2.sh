#!/bin/bash
# loop2.sh (bothachain): as loop.sh for pairs (PAIR=1): the time estimate of the next call is the sum of the last two records' seconds.
# usage: loop2.sh <stop_epoch>
cd "$(dirname "$0")"; S=../../../.scratch/bothachain; export PAIR=1
while true; do
  [ -e $S/STOP ] && { echo "loop2: STOP file $(date -u +%T)" >> logs/chain.err; break; }
  c1=$(cat $S/c1)
  st=$(python3 -c "
import json; from fractions import Fraction as F
r=[json.loads(l) for l in open('logs/chain.jsonl')]; ok=[d for d in r if d.get('ok')]
print('done' if F(ok[-1]['c']) >= F('$c1') else int(sum(d.get('sec') or 900 for d in r[-2:])))")
  [ "$st" = done ] && { echo "loop2: reached $c1 $(date -u +%T)" >> logs/chain.err; break; }
  [ $(( $(date +%s) + st )) -gt $1 ] && { echo "loop2: stop time $(date -u +%T)" >> logs/chain.err; break; }
  bash chain.sh $c1 1800 || { echo "loop2: chain.sh failed $(date -u +%T)" >> logs/chain.err; break; }
done
