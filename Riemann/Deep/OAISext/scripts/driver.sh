#!/bin/bash
# oaisext build driver: chunked, resumable lake build; one lean process; single-threaded.
# Policy: LEAN_NUM_THREADS=1 (lake and lean single-threaded; measured); inherited background QoS (taskpolicy -b) + nice 19
# (orchestrator order 22:58 UTC: no re-clamp).
ROOT=/Users/adriansager/playground/riemann
cd $ROOT || exit 1
S=.scratch/oaisext
L=Riemann/Deep/OAISext/logs
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 LEAN_NUM_THREADS=1
DEADLINE=$(date -j -u -f '%Y-%m-%d %H:%M' '2026-11-12 11:33' +%s)
ts() { date -u +%FT%TZ; }
count() { find .lake/build/lib/lean/Riemann/Deep/QRH .lake/build/lib/lean/Riemann/Deep/OAISext -name '*.olean' 2>/dev/null | wc -l | tr -d ' '; }
memmb() { top -l 1 -pid $1 -stats mem 2>/dev/null | tail -1 | tr -d '+-' | awk '{v=$1; u=substr(v,length(v)); n=substr(v,1,length(v)-1); if(u=="K") print int(n/1024); else if(u=="M") print int(n); else if(u=="G") print int(n*1024); else print 0}'; }
dffree() { df -g / | awk 'NR==2{print $4}'; }
trunc_ir() { # truncate C/setup files of finished modules (never delete)
  find .lake/build/ir/Riemann -type f \( -name '*.c' -o -name '*.setup.json' \) -size +0 -mmin +5 -print0 2>/dev/null | xargs -0 -n 200 perl -e 'for (@ARGV) { truncate($_, 0) }'; }
waitchange() { # wait for targets.txt change, GO file, STOP file or deadline
  local m0=$(stat -f %m $S/targets.txt)
  while :; do
    [ -f $S/STOP ] && return 1
    [ -f $S/GO ] && { rm -f $S/GO; return 0; }
    [ $(stat -f %m $S/targets.txt) != $m0 ] && return 0
    [ $(date +%s) -ge $DEADLINE ] && return 1
    sleep 30
  done; }
n=$(ls $L/chunk_*.log 2>/dev/null | wc -l | tr -d ' ')
echo "$(ts) driver start pid $$ built=$(count)" >> $L/build.log
while :; do
  [ -f $S/STOP ] && { echo "$(ts) STOP file: driver exits" >> $L/build.log; exit 0; }
  [ $(date +%s) -ge $DEADLINE ] && { echo "$(ts) deadline: driver exits" >> $L/build.log; exit 0; }
  [ "$(dffree)" -lt 10 ] && { echo "$(ts) disk $(dffree)G < 10G: truncating ir .c/.setup.json" >> $L/build.log; trunc_ir; }
  n=$((n+1)); C=$(printf '%03d' $n)
  T=$(tr '\n' ' ' < $S/targets.txt)
  b0=$(count); t0=$(date +%s)
  echo "$(ts) chunk $C start built=$b0 targets=$(wc -l < $S/targets.txt | tr -d ' ') df=$(dffree)G" >> $L/build.log
  rm -f $S/pgid
  PGF=$ROOT/$S/pgid perl -e 'setpgrp(0,0); open(F,">",$ENV{PGF}); print F $$; close F; alarm 1800; exec @ARGV' lake build -v $T > $L/chunk_$C.log 2>&1 &
  P=$!
  sleep 2; PG=$(cat $S/pgid 2>/dev/null)
  maxl=0; maxmem=0; maxmod=""; maxlake=0; warned=""
  while kill -0 $P 2>/dev/null; do
    sleep 30
    [ -z "$PG" ] && PG=$(cat $S/pgid 2>/dev/null)
    lp=$(ps -Ao pid=,pgid=,comm= | awk -v g=$PG '$2==g && $3 ~ /\/lean$/ {print $1}')
    nl=$(echo $lp | wc -w | tr -d ' '); [ $nl -gt $maxl ] && maxl=$nl
    lk=$(ps -Ao pid=,pgid=,comm= | awk -v g=$PG '$2==g && $3 ~ /\/lake$/ {print $1}' | head -1)
    [ -n "$lk" ] && { km=$(memmb $lk); [ "$km" -gt $maxlake ] && maxlake=$km; }
    for p in $lp; do
      m=$(memmb $p); mod=$(ps -o args= -p $p | tr ' ' '\n' | rg -m1 '\.lean$' | sed "s#$ROOT/##")
      fl=$(top -l 1 -pid $p -stats faults 2>/dev/null | tail -1 | tr -d ' +-'); pr=$(ps -o pri=,ni= -p $p | tr -s ' ')
      echo "$(ts) $C lean $p mem ${m}M faults $fl pri/ni $pr $mod" >> $L/samples.log
      [ "$m" -gt $maxmem ] && { maxmem=$m; maxmod=$mod; }
      if [ "$m" -gt 600 ]; then
        echo "$(ts) MEMORY STOP: lean $p ${m}M > 600M on $mod; killing group $PG" >> $L/build.log
        kill -TERM -- -$PG 2>/dev/null; sleep 3; kill -KILL -- -$PG 2>/dev/null
        echo "$mod ${m}M" > $S/MEMSTOP; touch $S/STOP; exit 2
      fi
      if [ "$m" -gt 450 ] && [[ "$warned" != *" $p "* ]]; then echo "$(ts) WARN lean $p ${m}M > 450M on $mod" >> $L/build.log; warned="$warned $p "; fi
    done
    [ -f $S/STOP ] && { echo "$(ts) STOP file during chunk: killing group $PG" >> $L/build.log; kill -TERM -- -$PG 2>/dev/null; sleep 3; kill -KILL -- -$PG 2>/dev/null; exit 0; }
  done
  wait $P; rc=$?
  [ -n "$PG" ] && { kill -TERM -- -$PG 2>/dev/null; sleep 2; kill -KILL -- -$PG 2>/dev/null; }
  b1=$(count); dt=$(( $(date +%s) - t0 ))
  echo "$(ts) chunk $C end rc=$rc dt=${dt}s built=$b1 (+$((b1-b0))) maxlean=$maxl maxleanmem=${maxmem}M ($maxmod) maxlakemem=${maxlake}M df=$(dffree)G" >> $L/build.log
  if [ $rc -eq 0 ]; then
    echo "$(ts) BUILD OK for current targets; waiting for new targets/GO" >> $L/build.log
    waitchange || { echo "$(ts) driver exits (STOP/deadline)" >> $L/build.log; exit 0; }
  elif [ $rc -ne 142 ]; then
    rg -n 'error' $L/chunk_$C.log | head -20 >> $L/build.log
    echo "$(ts) BUILD ERROR rc=$rc; waiting for fix (targets change or GO)" >> $L/build.log
    waitchange || { echo "$(ts) driver exits (STOP/deadline)" >> $L/build.log; exit 0; }
  fi
done
