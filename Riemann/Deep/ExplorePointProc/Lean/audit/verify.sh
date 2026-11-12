#!/bin/bash
# pplean final verification: compile every module in dependency order (one at a time), lean -j1 directly,
# under memcap 290 MB / 1800 s, with /usr/bin/time -l (peak memory footprint = the cap's measure); writes oleans.
cd /Users/adriansager/playground/riemann
O=/Users/adriansager/playground/riemann/.scratch/pplean/olean
LP="$O:$(cat .scratch/pplean/leanpath.txt)"
LEANBIN=/Users/adriansager/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean
out=.scratch/pplean/verify.log; : > $out
for m in FinSum Model NoGood Theorem1 PL Herm Cert Liminf Final Check; do
  t0=$(date +%s)
  LEAN_PATH="$LP" bash .scratch/memcap.sh 290 1800 /usr/bin/time -l $LEANBIN -j1 -o $O/Riemann/Deep/ExplorePointProc/Lean/$m.olean Riemann/Deep/ExplorePointProc/Lean/$m.lean > .scratch/pplean/v_$m.log 2>&1
  rc=$?
  pk=$(rg -o '[0-9]+  peak memory footprint' .scratch/pplean/v_$m.log | awk '{printf "%d", $1/1048576}')
  nw=$(rg -c '^Riemann.*warning' .scratch/pplean/v_$m.log); [ -z "$nw" ] && nw=0
  ne=$(rg -c '^Riemann.*error' .scratch/pplean/v_$m.log); [ -z "$ne" ] && ne=0
  mc=$(rg -o 'memcap: peak [0-9]+ MB' .scratch/pplean/v_$m.log)
  ln=$(wc -l < Riemann/Deep/ExplorePointProc/Lean/$m.lean)
  echo "$m rc=$rc $(( $(date +%s) - t0 ))s footprint=$pk MiB, $mc, warnings=$nw errors=$ne lines=$ln" >> $out
done
cat $out
