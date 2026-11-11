#!/bin/bash
# usage: mk.sh Name [Name...]  -> builds .scratch/lilean/olean/.../Name.olean sequentially, logs to .scratch/lilean/Name.log
cd /Users/adriansager/playground/riemann
O=/Users/adriansager/playground/riemann/.scratch/lilean/olean
for n in "$@"; do
  lake env sh -c "LEAN_PATH=\"$O:\$LEAN_PATH\" perl -e 'alarm 900; exec @ARGV' lean -o $O/Riemann/Deep/LiRate/Lean/$n.olean Riemann/Deep/LiRate/Lean/$n.lean" > .scratch/lilean/$n.log 2>&1
  echo "EXIT $?" >> .scratch/lilean/$n.log
done
