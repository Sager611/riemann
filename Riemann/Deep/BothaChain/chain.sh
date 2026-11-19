#!/bin/bash
# chain.sh (bothachain): continue bothaaudit's chain (logs/chain.jsonl starts with its 3 records, last ok c = 385/12, U = 1.72477599362e-161)
# with ccm/ (CCMFMM's scripts via bothaaudit; driver changes listed in ccm/chain_ccm.py's header). One ok step per call (CCM_MAXSTEPS=1);
# whole process tree under memcap 290 MB and T <= 1800 s.   usage: chain.sh <c1> <T> [int.jsonl] [chain.jsonl] [c0]
cd "$(dirname "$0")"; export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
R=$(cd ../../..; pwd)/.scratch/bothachain/ccm
export CCM_RUNS=${CCM_RUNS:-$R/runs} CCM_BIN=$R WZ_RSCAN2=$R/rscan2 CCM_MAXSTEPS=${CCM_MAXSTEPS:-1}
echo "start $(date -u +%H:%M:%S) c1=$1" >> logs/chain.err
if [ "$PAIR" = 1 ]; then   # one pair of steps per call, shared defl certificate (ccm/chain_pair.py)
  bash ../../../.scratch/memcap.sh 290 $2 uv run --no-project --with mpmath --with numpy --with python-flint python3 ccm/chain_pair.py ${3:-logs/chain_int.jsonl} ${4:-logs/chain.jsonl} $1 >> logs/chain.out 2>> logs/chain.err
else
  bash ../../../.scratch/memcap.sh 290 $2 uv run --no-project --with mpmath --with numpy --with python-flint python3 ccm/chain_ccm.py ${3:-logs/chain_int.jsonl} ${4:-logs/chain.jsonl} ${5:-30} $1 >> logs/chain.out 2>> logs/chain.err
fi
rc=$?; echo "end $(date -u +%H:%M:%S) rc=$rc" >> logs/chain.err; exit $rc
