#!/bin/bash
# chain.sh (bothaaudit): continue CCMFMM's certified chain from c = 30 (its U(30) = 3.27353749418e-150) with its scripts, copied unchanged
# into ccm/ (only prep_ccm.py's search path for ../WeilZeta patched); step o1 heuristics from logs/chain_int.jsonl (our scan, scaled 0.6).
# Whole tree under memcap 290 MB; T = total seconds.   usage: chain.sh <c1> <T>
cd "$(dirname "$0")"; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export CCM_RUNS=$(cd ../../..; pwd)/.scratch/bothaaudit/ccm/runs CCM_BIN=$(cd ../../..; pwd)/.scratch/bothaaudit/ccm
bash ../../../.scratch/memcap.sh 290 $2 uv run --no-project --with mpmath --with numpy --with python-flint python3 ccm/chain_ccm.py logs/chain_int.jsonl logs/chain.jsonl 30 $1 >> logs/chain.out 2>> logs/chain.err
