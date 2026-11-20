#!/bin/bash
# run_certs.sh (explore67): lower-bound certificates (cert_simple_b.py, Chekmachev lattice minorant + Woodbury + Cholesky, rigorous scan),
# one process at a time, each under .scratch/memcap.sh 290 MB; appends to logs/certs.jsonl. Run from the repo root: bash <this> <listfile>
# List lines: label | env | c sector eps bits
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
export WZ_RSCAN=$PWD/.scratch/explore67/rscan WZ_SCAN=$PWD/.scratch/explore67/scan
mkdir -p .scratch/explore67/scan
D=Riemann/Deep/ExploreWigner
while IFS='|' read -r label envx args; do
  [ -z "$label" ] && continue
  ( env $envx /usr/bin/time -l bash .scratch/memcap.sh 290 1500 uv run --no-project --with python-flint --with numpy --with mpmath python3 $D/code/cert_simple_b.py $args > .scratch/explore67/out.json ) 2> .scratch/explore67/err.txt
  L=$(tail -1 .scratch/explore67/out.json); M=$(rg -o 'memcap: peak [0-9]+ MB' .scratch/explore67/err.txt | head -1)
  R=$(rg 'maximum resident' .scratch/explore67/err.txt | awk '{print int($1/1048576)}')
  [ -z "$L" ] && L='null'
  [ -z "$R" ] && R=0
  echo "{\"label\": \"$label\", \"memcap\": \"$M\", \"rss_mb\": $R, \"res\": $L}" >> $D/logs/certs.jsonl
  echo "$label: $(echo "$L" | rg -o '"rigorous": [a-z]+|"margin": "[^"]*"|"times": \[[^]]*\]' | tr '\n' ' ') $M rss $R"
done < "$1"
