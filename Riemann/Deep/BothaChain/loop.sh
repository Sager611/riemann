#!/bin/bash
# loop.sh (bothachain): call chain.sh (one memcap'd step, <= 1800 s) until the file .scratch/bothachain/STOP exists, the next step
# would end after the stop time (arg 1, epoch s; estimate = last record's seconds), the target c1 (file .scratch/bothachain/c1) is reached,
# or chain.sh fails.   usage: loop.sh <stop_epoch>
cd "$(dirname "$0")"; S=../../../.scratch/bothachain
while true; do
  [ -e $S/STOP ] && { echo "loop: STOP file $(date -u +%T)" >> logs/chain.err; break; }
  c1=$(cat $S/c1)
  st=$(python3 -c "
import json; from fractions import Fraction as F
r=[json.loads(l) for l in open('logs/chain.jsonl')]; ok=[d for d in r if d.get('ok')]
print('done' if F(ok[-1]['c']) >= F('$c1') else int(r[-1].get('sec') or 900))")
  [ "$st" = done ] && { echo "loop: reached $c1 $(date -u +%T)" >> logs/chain.err; break; }
  [ $(( $(date +%s) + st )) -gt $1 ] && { echo "loop: stop time $(date -u +%T)" >> logs/chain.err; break; }
  bash chain.sh $c1 1800 || { echo "loop: chain.sh failed $(date -u +%T)" >> logs/chain.err; break; }
done
