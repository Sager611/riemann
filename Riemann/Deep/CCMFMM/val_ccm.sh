#!/bin/bash
# validation (ccmfmm): reproduce ../WeilSimple's L_2 (defl) and L_odd certificates at c = 10, 16, 20 (logs/int.jsonl: same NG, bits, K2, eps)
# in the particle form: point-mode LDL at nu = m(1 - 1e-9) must pass and at m(1 + 1e-9) must fail, m = WeilSimple's margin.
cd "$(dirname "$0")"
R=../../../.scratch/ccmfmm/runs; L=logs/val
mkdir -p $L
nu() { python3 -c "m=$1; print('%.12e %.12e' % (m*(1-1e-9), m*(1+1e-9)))"; }
for spec in "10 80 320 201 7.348097e-40 9.55135345279e-37 1.349542e-41 7.19952806309e-40" \
            "16 120 512 427 4.085821e-71 1.12391899385e-67 4.580865e-73 2.93457950161e-71" \
            "20 160 576 667 2.841006e-92 1.62050936675e-88 2.497138e-94 6.60901632398e-92"; do
  set -- $spec; c=$1; NG=$2; bits=$3; k2=$4; e2=$5; m2=$6; eo=$7; mo=$8
  [ -s $R/g$c.json ] || /usr/bin/time -l perl -e 'alarm 1790; exec @ARGV' uv run --no-project --with mpmath --with numpy --with python-flint python3 gstate.py $c $NG $bits $R/g$c.json 0 > $L/g$c.out 2> $L/g$c.time
  CCM_G=$R/g$c.json ./one.sh v${c}d $c defl $e2 $bits $k2 $(nu $m2) > $L/v${c}d.out 2>&1
  ./one.sh v${c}o $c odd $eo $bits $k2 $(nu $mo) > $L/v${c}o.out 2>&1
  cp $R/v${c}d.prep.json $R/v${c}d.slt.json $R/v${c}d.prep.time $R/v${c}d.slt.time $R/v${c}o.prep.json $R/v${c}o.slt.json $R/v${c}o.prep.time $R/v${c}o.slt.time $L/ 2>/dev/null
  echo "done $c $(date -u +%H:%M:%S)"
done
