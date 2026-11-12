#!/bin/bash
# multistart runs, sequential (one compute process at a time)
for st in GRC uncond; do for k in krzQ rand1 rand2 Rq3; do
  /usr/bin/time -l perl -e 'alarm 900; exec @ARGV' uv run --no-project --with numpy python3 joint.py ms $st $k 4 >> logs/multistart.log 2>> logs/multistart.time
done; done
echo DONE >> logs/multistart.log
