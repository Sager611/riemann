#!/bin/bash
# repro.sh -- re-run the paper's 8-point certificate search (second window) in checkpointed chunks
cd "$(dirname "$0")"
T=../src/anc/certificate_line/kernel_lower.bin
P=second_window.txt
prev=""
for k in 1 2 3 4 5 6; do
  S=$(date +%s)
  if [ -z "$prev" ]; then
    perl -e 'alarm 1790; exec @ARGV' nice -n 19 ./verify8cp $P $T 150000000 cp_$k.txt > out_$k.json 2> err_$k.txt
  else
    perl -e 'alarm 1790; exec @ARGV' nice -n 19 ./verify8cp $P $T 150000000 cp_$k.txt $prev > out_$k.json 2> err_$k.txt
  fi
  rc=$?
  E=$(date +%s)
  echo "chunk $k exit $rc time $((E-S)) s $(cat out_$k.json)" >> repro.log
  if [ $rc -eq 7 ]; then prev=cp_$k.txt; continue; fi
  break
done
echo "done $(date -u +%H:%M:%S)" >> repro.log
