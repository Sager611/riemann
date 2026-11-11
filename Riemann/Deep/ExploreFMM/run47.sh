#!/bin/bash
# c = 47, one sector (run47.sh even|odd eps bits nu...; needs ./slt, ./rscan and scan/scan_zeta_47_3683_rig.json from ../WeilZeta/logs): lean prep (Python), Cauchy sums in 3 chunks, structured LDL; one process at a time, each < 1790 s
cd "$(dirname "$0")"
export WZ_RSCAN=$PWD/rscan WZ_SCAN=$PWD/scan OMP_NUM_THREADS=1 LEAN=1
SEC=$1; EPS=$2; P=$3; shift 3
b=runs/c47$SEC
log=runs/c47$SEC.log
echo "start $(date -u)" > $log
[ -s $b.idx ] || /usr/bin/time -l perl -e 'alarm 1790; exec @ARGV' uv run --no-project --with mpmath --with numpy --with python-flint python3 prep_slt.py 47 $SEC $EPS $P $b > $b.prep.json 2> $b.prep.time
echo "prep done $(date -u)" >> $log; cat $b.prep.json >> $log
n=$(head -1 $b.idx | cut -d' ' -f1)
c1=$((n/3)); c2=$((2*n/3))
for rg in "0 $c1" "$c1 $c2" "$c2 $n"; do
  part=$b.F_$(echo $rg | tr ' ' '_')
  [ -s $part.done ] || { SLT_RANGE="$rg" SLT_FOUT=$part /usr/bin/time -l perl -e 'alarm 1790; exec @ARGV' ./slt $b $P > /dev/null 2> $part.time && echo ok > $part.done; }
  echo "sums $rg done $(date -u)" >> $log; tail -2 $part.time | head -1 >> $log
done
parts=$(ls $b.F_* | rg -v 'time|done' | tr '
' ',' | sed 's/,$//')
SLT_POINT=1 SLT_FIN=$parts /usr/bin/time -l perl -e 'alarm 1790; exec @ARGV' ./slt $b $P -g $b.gen "$@" > $b.slt.json 2> $b.slt.time
echo "ldl done $(date -u)" >> $log; cat $b.slt.json >> $log; rg 'maximum res|real' $b.slt.time >> $log
