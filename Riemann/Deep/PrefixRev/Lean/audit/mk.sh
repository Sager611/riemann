#!/bin/bash
# usage: bash .scratch/prefixlean/mk.sh Module [timeout] [lines]  (compiles Riemann/Deep/PrefixRev/Lean/Module.lean to an olean, under memcap)
m=$1; d=Riemann/Deep/PrefixRev/Lean; o=.scratch/prefixlean/olean/$d
LP=.scratch/prefixlean/olean:$(cat .scratch/prefixlean/leanpath.txt)
L=~/.elan/toolchains/leanprover--lean4---v4.35.0-rc2/bin/lean
mkdir -p $o .scratch/prefixlean/logs
LEAN_PATH=$LP bash .scratch/memcap.sh 290 ${2:-1800} /usr/bin/time -l $L -j1 -o $o/$m.olean $d/$m.lean > .scratch/prefixlean/logs/$m.log 2>&1
rc=$?
pk=$(rg "peak memory footprint|memcap|real" .scratch/prefixlean/logs/$m.log | tr -s ' ' | tr '\n' ' ')
echo "$(date -u +%H:%M) $m rc=$rc $pk" | tee -a .scratch/prefixlean/peaks.log
rg -v "peak memory|memcap|real|resident|instructions|cycles|page|swaps|block|messages|signals|context|footprint|^\s*[0-9]+\s+[a-z ]+$" .scratch/prefixlean/logs/$m.log | head -${3:-40}
