#!/bin/bash
# dbnrec run_row.sh ROW T u0 u1 u2 ... : for each shell [u_j, u_{j+1}]: prep_end.py (Arb enclosures of Q^+) then bb_rec (thr = T + s u_j), sequentially,
# each process under memcap 290 MB / 1800 s, single thread. Logs: Riemann/Deep/DBNRec3/logs/r<ROW>_<j>.{prep.json,bb.out}, unres dump .unres.txt.
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
ARC=.scratch/dbnrec3/dl/Proof_and_Replay__Lambda_Peer_Review__verification__arcs__ACTUAL_67_ROW_ARC_INPUTS.json
AFF=.scratch/dbnrec3/dl/Proof_and_Replay__Lambda_Peer_Review__verification__referee__ALL_67_ORIGINAL_AFFINE_NATIVE.json
L=Riemann/Deep/DBNRec3/logs; BB=.scratch/dbnrec3/bb_rec
row=$1; T=$2; shift 2
us=("$@"); n=${#us[@]}
for ((j=0; j<n-1; j++)); do
  a=${us[$j]}; b=${us[$((j+1))]}; tag=r${row}_$(printf %02d $j)
  [ -f $L/$tag.bb.out ] && rg -q '^CERT' $L/$tag.bb.out && continue
  [ -n "$SHELL_DEADLINE" ] && [ $(date +%s) -gt $SHELL_DEADLINE ] && { echo "$(date -u +%H:%M:%S) $tag skipped: shell deadline" >> $L/progress_r${row}.txt; break; }
  bash .scratch/memcap.sh 290 600 uv run --no-project --with mpmath --with python-flint python3 Riemann/Deep/DBNRec3/code/prep_end.py $ARC $AFF $row $a $b $T .scratch/dbnrec3/$tag.dat > $L/$tag.prep.json 2> .scratch/dbnrec3/$tag.prep.err
  tail -1 .scratch/dbnrec3/$tag.prep.err >> $L/$tag.prep.json
  /usr/bin/time -l bash .scratch/memcap.sh 290 1800 $BB .scratch/dbnrec3/$tag.dat 48 400 0 0 0.001 $L/$tag.unres.txt > $L/$tag.bb.out 2> .scratch/dbnrec3/$tag.bb.err
  rg 'memcap|maximum resident|peak memory' .scratch/dbnrec3/$tag.bb.err >> $L/$tag.bb.out
  echo "$(date -u +%H:%M:%S) $tag [$a,$b] $(rg -o 'CERT.*' $L/$tag.bb.out) $(rg -o 'unres [0-9]+ \| maxUB3' $L/$tag.bb.out) $(rg -o '[0-9]+ s$' $L/$tag.bb.out | tail -1)" >> $L/progress_r${row}.txt
done
echo "$(date -u +%H:%M:%S) row $row done" >> $L/progress_r${row}.txt
