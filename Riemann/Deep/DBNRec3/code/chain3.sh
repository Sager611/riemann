#!/bin/bash
# dbnrec3 chain3.sh (= dbnrec2 chain2.sh, paths swapped) TAG CERT.json ... : exact chain with the given certified rows (tranche 1's chain(): mkrec source -> make_flow_rec.py -> mkrec verifier -> verify_arcs_rec.py),
# every compute step under memcap 290 MB; logs in Riemann/Deep/DBNRec3/logs/*_TAG.*
cd /Users/adriansager/playground/riemann
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1
L=Riemann/Deep/DBNRec3/logs; C=Riemann/Deep/DBNRec3/code; W=.scratch/dbnrec3/ws
UVF="uv run --no-project --with mpmath --with python-flint"
UVP="uv run --no-project --python 3.12 --with python-flint==0.8.0"
log() { echo "$(date -u +%H:%M:%S) $*" >> $L/master.log; }
tag=$1; shift
$UVF python3 $C/mkrec.py source "$@" > $L/mkrec_$tag.json 2>&1; log "mkrec source $tag rc=$?: $(cut -c1-300 $L/mkrec_$tag.json)"
rm -f $W/dbnrec/FLOW_REC.json.gz
/usr/bin/time -l bash .scratch/memcap.sh 290 1800 $UVP python3 $C/make_flow_rec.py --workspace $W --out $W/dbnrec/FLOW_REC.json.gz > $L/make_flow_$tag.log 2>&1
log "make_flow $tag rc=$? $(rg -o '\"endpoint_decimal\": [0-9.e-]+' $L/make_flow_$tag.log)"
cp $W/dbnrec/FLOW_REC_SUMMARY.json $L/FLOW_REC_SUMMARY_$tag.json
$UVF python3 $C/mkrec.py verifier > $L/mkrec_ver_$tag.json 2>&1; log "mkrec verifier $tag rc=$?: $(cut -c1-200 $L/mkrec_ver_$tag.json)"
(cd .scratch/dbnrec3/verify && /usr/bin/time -l bash ../../memcap.sh 290 1800 $UVP python3 verify_arcs_rec.py --workspace ../ws > ../../../$L/verify_$tag.log 2>&1)
log "verify $tag rc=$? $(rg -o '\"PASS\": [a-z]+' $L/verify_$tag.log)"
cp .scratch/dbnrec3/verify/EXACT_ARC_VERIFICATION_FRESH_V1.json $L/EXACT_ARC_VERIFICATION_$tag.json 2>/dev/null
ls -la .scratch/dbnrec3/verify/ >> $L/verify_$tag.log
