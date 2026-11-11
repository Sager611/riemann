#!/bin/bash
# validation runs: prep (Python, Chekmachev's node sets as in WeilZeta) + slt (structured LDL in Arb); one process at a time
cd "$(dirname "$0")"
export WZ_RSCAN=$PWD/rscan WZ_SCAN=$PWD/scan OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
run() { # c sector eps bits margin
  b=runs/c$1$2
  [ -s $b.idx ] || /usr/bin/time -l perl -e 'alarm 1700; exec @ARGV' uv run --no-project --with mpmath --with numpy --with python-flint python3 prep_slt.py $1 $2 $3 $4 $b > $b.prep.json 2> $b.prep.time
  lo=$(python3 -c "print('%.15e' % ($5*(1-1e-9)))"); hi=$(python3 -c "print('%.15e' % ($5*(1+1e-9)))")
  /usr/bin/time -l perl -e 'alarm 1700; exec @ARGV' ./slt $b $4 -g $b.gen 0 > $b.slt.json 2> $b.slt.time
  /usr/bin/time -l perl -e 'alarm 1700; exec @ARGV' ./slt $b $4 -l $b.gen -b $lo $hi ${BIS:-12} > $b.bis.json 2> $b.bis.time
  echo "== c=$1 $2 eps=$3 bits=$4 WeilZeta margin $5"; cat $b.prep.json $b.slt.json $b.bis.json
  rg -o '[0-9]+  maximum resident set size' $b.prep.time $b.slt.time $b.bis.time; rg -o 'Cauchy sums.*|total.*' $b.slt.time
}
for spec in "$@"; do run $spec; done
