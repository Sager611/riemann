#!/bin/bash
# usage: mk.sh <Module>  (compiles Riemann/Deep/NBLp/Lean/<Module>.lean; oleans to .scratch/nblplean/olean)
cd /Users/adriansager/playground/riemann
O=$PWD/.scratch/nblplean/olean
mkdir -p $O/Riemann/Deep/NBLp/Lean
n=$1
L=.scratch/nblplean/log_$n.txt
M=.scratch/nblplean/mem_$n.txt
: > $M
perl -e 'alarm 1500; exec @ARGV' lake env sh -c "LEAN_PATH=\"$O:\$LEAN_PATH\" exec lean -o $O/Riemann/Deep/NBLp/Lean/$n.olean Riemann/Deep/NBLp/Lean/$n.lean" > $L 2>&1 &
P=$!
while kill -0 $P 2>/dev/null; do
  for q in $(pgrep -f "lean -o $O/Riemann/Deep/NBLp/Lean/$n.olean"); do
    echo "$(date +%T) $q $(top -l 1 -pid $q -stats mem | tail -1)" >> $M
  done
  sleep 8
done
wait $P
echo "EXIT $?" >> $L
echo "peak: $(awk '{print $3}' $M | sort -h | tail -1)" >> $L

