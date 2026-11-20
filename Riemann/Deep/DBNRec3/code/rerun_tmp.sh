#!/bin/bash
# dbnrec3 rerun_tmp.sh (= dbnrec2 rerun_tmp.sh adapted): reproduce chain G (Λ ≤ 0.107410439495; rows 1-3 tranche 1, 4-7 tranche 2, 8-12 tranche 3) from scratch in /tmp/dbnv3; writes nothing in the repo.
# 1 fetch the 4 inputs with Riemann/Deep/DBNHigh/rzip.py (Zenodo 23752103, HTTP range) and check the sha256 pins of REPORT $Setup -> /tmp/dbnv3/{dl,ws}
# 2 code copy -> /tmp/dbnv3/code, every repo-relative or absolute repo path sed-ed to /tmp/dbnv3 (hard check: none left)
# 3 certs -> /tmp/dbnv3/certs   4 mkrec.py source (expect source sha256 0045c6f68e1f1fbe30c1c7b0b5ef646f9c167cfca1c6c4119c9dab6b3977fc97; regenerated make_flow_rec.py = the repo's)
# 5 make_flow_rec.py   6 mkrec.py verifier, diff against the repo's verify_arcs_rec.py   7 verify_arcs_rec.py   8 compare with chain G's logs
# Every python step: bash .scratch/memcap.sh 290 <s> (verify 1500 s) under /usr/bin/time -l, single thread, one at a time. Runtime ~8 min, disk ~0.45 GB.
# Not reproducible by design: sha256 of the gzip files (flow, witness trace) and hence of the verifier, which pins the flow sha256: Python's gzip header
# carries the write time. Checked instead: the flow summary (all fields but seconds), the verifier text apart from its flow pin, the uncompressed transcript.
# usage: bash Riemann/Deep/DBNRec3/code/rerun_tmp.sh   (any cwd; replaces /tmp/dbnv3 only if it carries this script's .dbnv3 marker)
set -uo pipefail
REPO=$(cd "$(dirname "$0")/../../../.." && pwd); V=/tmp/dbnv3; MC="$REPO/.scratch/memcap.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1 PYTHONDONTWRITEBYTECODE=1 RZIP_CACHE=$V/cd.bin
UVF="uv run --no-project --with mpmath --with python-flint"
UVP="uv run --no-project --python 3.12 --with python-flint==0.8.0"
say() { echo "$(date -u +%H:%M:%S) $*"; }
die() { say "FAIL: $*"; exit 1; }
mc() {   # mc SECONDS TAG STDOUT_FILE cmd ... : one compute step under memcap 290 MB; stderr (incl. memcap and time -l) -> logs/TAG.err
  local s=$1 tag=$2 out=$3; shift 3
  /usr/bin/time -l bash "$MC" 290 "$s" "$@" > "$out" 2> "$V/logs/$tag.err"; local rc=$?
  local pk rss rt
  pk=$(rg -o 'memcap: peak [0-9]+ MB( KILLED [^ ]+)?' "$V/logs/$tag.err" | tail -1)
  rss=$(rg -o '[0-9]+ +maximum resident set size' "$V/logs/$tag.err" | tail -1 | awk '{printf "%d MB", $1/1048576}')
  rt=$(rg -o '[0-9.]+ real' "$V/logs/$tag.err" | tail -1)
  say "  $tag: rc=$rc | $pk | max RSS $rss | $rt"; return $rc; }
chk() { local h; h=$(shasum -a 256 "$1" | cut -c1-64); [ "$h" = "$2" ] && say "  sha256 OK $h $(basename "$1")" || die "sha256 of $(basename "$1"): $h, pin $2"; }
[ -f "$MC" ] && [ -f "$REPO/Riemann/Deep/DBNHigh/rzip.py" ] && [ -f "$REPO/Riemann/Deep/DBNRec/logs/cert_r1.json" ] && [ -f "$REPO/Riemann/Deep/DBNRec2/logs/cert_r4.json" ] || die "repo not found at $REPO"
[ -e /tmp/lambda_publication_flint_runtime ] && die "/tmp/lambda_publication_flint_runtime exists: the verifier puts it first on sys.path (it would shadow uv's python-flint)"
if [ -e $V ]; then [ -f $V/.dbnv3 ] || die "$V exists without the .dbnv3 marker: move it away"; rm -rf $V; fi
mkdir -p $V/dl $V/ws/lambda_round14/material/m8_late $V/ws/lambda_round15/flow_main_v2 $V/code $V/certs $V/verify $V/logs && touch $V/.dbnv3 $V/.start || die "mkdir"
cd $V; T0=$(date +%s); say "repo $REPO -> $V"

say "1. fetch inputs (rzip.py, Zenodo 23752103 Proof_and_Replay.zip) and check the sha256 pins"
P=Proof_and_Replay/Lambda_Peer_Review; D=$V/dl/Proof_and_Replay__Lambda_Peer_Review
mc 600 rzip $V/logs/rzip.out python3 "$REPO/Riemann/Deep/DBNHigh/rzip.py" get $V/dl $P/inputs/workspace/lambda_round14/material/m8_late/M8_LATE_PLAN.json \
  $P/inputs/workspace/lambda_round15/flow_main_v2/MAIN_CLOCK_SOURCE.json $P/verification/arcs/ACTUAL_67_ROW_ARC_INPUTS.json \
  $P/verification/referee/ALL_67_ORIGINAL_AFFINE_NATIVE.json || die "rzip"
chk ${D}__inputs__workspace__lambda_round14__material__m8_late__M8_LATE_PLAN.json ec346b7c1136eda38021eb6323033b3b48ed21e6f7458eeeac7cc9579db2497a
chk ${D}__inputs__workspace__lambda_round15__flow_main_v2__MAIN_CLOCK_SOURCE.json b59905d7b994d3424ba14f9347be75a2429ae4d4b0de2c2404971ffda1d4d536
chk ${D}__verification__arcs__ACTUAL_67_ROW_ARC_INPUTS.json 29902c66cabd632cf4ddadef7d9e23ed094f5b6b0f118a2ad76f8f53eed1a5cb
chk ${D}__verification__referee__ALL_67_ORIGINAL_AFFINE_NATIVE.json 61bde51a2789609b691fcb5878ee1c66b215cbb2a221c26312c7a036c780f8bf
cp ${D}__inputs__workspace__lambda_round14__material__m8_late__M8_LATE_PLAN.json $V/ws/lambda_round14/material/m8_late/M8_LATE_PLAN.json
cp ${D}__inputs__workspace__lambda_round15__flow_main_v2__MAIN_CLOCK_SOURCE.json $V/ws/lambda_round15/flow_main_v2/MAIN_CLOCK_SOURCE.json

say "2. code copy -> $V/code, repo paths sed-ed to $V"
cp "$REPO"/Riemann/Deep/DBNRec3/code/* $V/code/ && rm -f $V/code/rerun_tmp.sh && cp "$MC" $V/memcap.sh || die "copy"
for f in $V/code/*; do sed -i '' -e "s#$REPO#$V#g" -e 's#/Users/adriansager/playground/riemann#/tmp/dbnv3#g' -e 's#\.scratch/dbnrec3#/tmp/dbnv3#g' \
  -e 's#Riemann/Deep/DBNRec3#/tmp/dbnv3#g' -e 's#Riemann/Deep/DBNRec2/logs#/tmp/dbnv3/certs#g' -e 's#Riemann/Deep/DBNRec2#/tmp/dbnv3/dbnrec2#g' \
  -e 's#Riemann/Deep/DBNRec/logs#/tmp/dbnv3/certs#g' -e 's#Riemann/Deep/DBNHigh#/tmp/dbnv3/dbnhigh#g' \
  -e 's#Riemann/Deep/CeilAudit5/logs#/tmp/dbnv3/ceilaudit5#g' -e 's#\.\./\.\./\.\./\$L/#$L/#g' -e 's#\.\./\.\./memcap\.sh#/tmp/dbnv3/memcap.sh#g' \
  -e 's#\.scratch/memcap\.sh#/tmp/dbnv3/memcap.sh#g' "$f"; done
left=$(rg -n 'Riemann/|\.scratch|/Users/|\.\./\.\./' $V/code); [ -z "$left" ] || die "repo paths left in the code copy: $left"
odd=$(rg -n -o '/tmp/[A-Za-z0-9_.-]+' $V/code | rg -v '/tmp/dbnv3|/tmp/lambda_publication_flint_runtime'); [ -z "$odd" ] || die "unexpected /tmp paths: $odd"
for f in make_flow_b2.py verify_arcs_141.py make_flow_rec.py verify_arcs_rec.py; do cmp -s $V/code/$f "$REPO/Riemann/Deep/DBNRec3/code/$f" || die "$f changed by sed"; done
say "  no repo path left; producer/verifier templates byte-identical to the repo's (their only absolute path, /tmp/lambda_publication_flint_runtime, is Raymond's and absent)"

say "3. certs -> $V/certs (tranche 1: rows 1-3, tranche 2: rows 4-7, tranche 3: rows 8-12)"
for r in 1 2 3; do cp "$REPO/Riemann/Deep/DBNRec/logs/cert_r$r.json" $V/certs/; done
for r in 4 5 6 7; do cp "$REPO/Riemann/Deep/DBNRec2/logs/cert_r$r.json" $V/certs/; done
for r in 8 9 10 11 12; do cp "$REPO/Riemann/Deep/DBNRec3/logs/cert_r$r.json" $V/certs/; done
(cd $V/certs && shasum -a 256 cert_r*.json | sed 's/^/  /')

say "4. mkrec.py source (Raymond's MAIN_CLOCK_SOURCE.json with first_hi/last_hi of rows 1-12 replaced)"
mc 300 mkrec_source $V/logs/mkrec_source.json $UVF python3 $V/code/mkrec.py source $(for r in 1 2 3 4 5 6 7 8 9 10 11 12; do echo -n "$V/certs/cert_r$r.json "; done) || die "mkrec source"
HS=$(shasum -a 256 $V/ws/dbnrec/MAIN_CLOCK_SOURCE_REC.json | cut -c1-64)
[ "$HS" = 0045c6f68e1f1fbe30c1c7b0b5ef646f9c167cfca1c6c4119c9dab6b3977fc97 ] && say "  source sha256 $HS (expected 0045c6f68e1f1fbe30c1c7b0b5ef646f9c167cfca1c6c4119c9dab6b3977fc97): OK" || die "source sha256 $HS (expected 0045c6f68e1f1fbe30c1c7b0b5ef646f9c167cfca1c6c4119c9dab6b3977fc97)"
cmp -s $V/code/make_flow_rec.py "$REPO/Riemann/Deep/DBNRec3/code/make_flow_rec.py" && say "  regenerated make_flow_rec.py identical to the repo's" || die "make_flow_rec.py differs from the repo's"

say "5. make_flow_rec.py (theta141 start, step 1/2048000, exact FLINT rationals)"
mc 900 make_flow $V/logs/make_flow.out $UVP python3 $V/code/make_flow_rec.py --workspace $V/ws --out $V/ws/dbnrec/FLOW_REC.json.gz || die "make_flow"
say "  $(rg -o '"endpoint_decimal": [0-9.e-]+' $V/ws/dbnrec/FLOW_REC_SUMMARY.json) | flow sha256 $(shasum -a 256 $V/ws/dbnrec/FLOW_REC.json.gz | cut -c1-64) (chain G's run: b484ea2a...d7d09; gzip mtime differs)"

say "6. mkrec.py verifier; diff against the repo's chain-G verifier"
mc 300 mkrec_verifier $V/logs/mkrec_verifier.json $UVF python3 $V/code/mkrec.py verifier || die "mkrec verifier"
NEWV=$V/verify/verify_arcs_rec.py; OLDV="$REPO/Riemann/Deep/DBNRec3/code/verify_arcs_rec.py"
cmp -s $NEWV $V/code/verify_arcs_rec.py || die "mkrec wrote two different verifiers"
say "  regenerated: sha256 $(shasum -a 256 $NEWV | cut -c1-64) | $(rg -o "'flow':\('[^']+','[0-9a-f]{64}'\)" $NEWV)"
say "  repo (G):    sha256 $(shasum -a 256 "$OLDV" | cut -c1-64) | $(rg -o "'flow':\('[^']+','[0-9a-f]{64}'\)" "$OLDV")"
diff "$OLDV" $NEWV > $V/logs/verifier.diff; say "  diff repo vs regenerated: $(rg -c '^[<>]' $V/logs/verifier.diff || echo 0) changed lines"; cut -c1-200 $V/logs/verifier.diff | sed 's/^/    /'
norm() { sed -E "s#('flow':\('dbnrec/FLOW_REC\.json\.gz',')[0-9a-f]{64}#\1<FLOW>#" "$1"; }
diff <(norm "$OLDV") <(norm $NEWV) > /dev/null && say "  verifiers identical apart from the flow pin" || die "verifier differs beyond the flow pin"

say "7. verify_arcs_rec.py (exact arc verifier; memcap 290 MB / 1500 s)"
(cd $V/verify && mc 1500 verify $V/logs/verify.out $UVP python3 verify_arcs_rec.py --workspace ../ws) || die "verify (see $V/logs/verify.out, verify.err)"

say "8. compare with chain G's logs"
cat > $V/compare.py <<'EOF'
import json, sys
new_v, old_v, new_s, old_s, mk = sys.argv[1:6]; ok = True
def bad(m):
    global ok; ok = False; print('  FAIL:', m)
N, O = json.load(open(new_v)), json.load(open(old_v)); w = N['all_arc_witnesses']; e = N['exact_full_cleanup_endpoint']
print('  PASS:', N['PASS'], '| arcs', N['arcs_checked'], '| exact checks', N['exact_checks'], '| rows covered', N['all_rows_covered'])
print('  endpoint', e, '= %.17f' % (int(e.split('/')[0])/int(e.split('/')[1])), '| headline', N['headline'])
print('  uncompressed canonical transcript sha256', w['uncompressed_canonical_transcript_sha256'], '(expected 22facbf6...d667c)')
if N['PASS'] is not True: bad('PASS is not true')
if w['uncompressed_canonical_transcript_sha256'] != '22facbf693085de8875c6b6e9ed2eafaea492dacee1c4d259530c50ee98d667c': bad('transcript sha256')
d = [k for k in sorted(set(N) | set(O)) if k not in ('seconds', 'checker_sha256', 'input_sha256', 'all_arc_witnesses') and N.get(k) != O.get(k)]
if {k: v for k, v in N['input_sha256'].items() if k != 'flow'} != {k: v for k, v in O['input_sha256'].items() if k != 'flow'}: d.append('input_sha256 plan/source')
if (w['path'], w['uncompressed_canonical_transcript_sha256']) != (O['all_arc_witnesses']['path'], O['all_arc_witnesses']['uncompressed_canonical_transcript_sha256']): d.append('witness path/transcript')
if d: bad('verification report differs from chain G in: ' + ', '.join(d))
else: print('  verification report = chain G (logs/EXACT_ARC_VERIFICATION_G.json) in every field but seconds, checker sha256, flow sha256, witness gzip bytes/sha256')
NS, OS = json.load(open(new_s)), json.load(open(old_s))
d2 = [k for k in sorted(set(NS) | set(OS)) if k != 'seconds' and NS.get(k) != OS.get(k)]
if d2: bad('flow summary differs in: ' + ', '.join(d2))
else: print('  flow summary = chain G (logs/FLOW_REC_SUMMARY_G.json) in every field but seconds: steps', NS['steps'], ', headline', NS['headline'])
M = json.load(open(mk)); rows = [r['row'] for r in M['rows'] if not r.get('skipped')]
if rows != list(range(1, 13)): bad('mkrec replaced rows %s' % rows)
else: print('  mkrec: rows 1-12 replaced (dA, dB):', ', '.join('%d (%.6f, %.6f)' % (r['row'], r['dA_float'], r['dB_float']) for r in M['rows']))
print('  RESULT:', 'PASS' if ok else 'FAIL'); sys.exit(0 if ok else 1)
EOF
mc 300 compare $V/logs/compare.out python3 $V/compare.py $V/verify/EXACT_ARC_VERIFICATION_FRESH_V1.json "$REPO/Riemann/Deep/DBNRec3/logs/EXACT_ARC_VERIFICATION_G.json" \
  $V/ws/dbnrec/FLOW_REC_SUMMARY.json "$REPO/Riemann/Deep/DBNRec3/logs/FLOW_REC_SUMMARY_G.json" $V/logs/mkrec_source.json; rc=$?; cat $V/logs/compare.out
new=$(find "$REPO/Riemann/Deep/DBNRec3" "$REPO/Riemann/Deep/DBNRec2" "$REPO/Riemann/Deep/DBNRec" "$REPO/Riemann/Deep/DBNHigh" "$MC" -newer $V/.start -print)
[ -z "$new" ] && say "repo: nothing written in DBNRec3/, DBNRec2/, DBNRec/, DBNHigh/, memcap.sh" || say "WARN: repo files newer than the start: $new"
say "done in $(( $(date +%s) - T0 )) s; $V holds $(du -sh $V | cut -f1) (witness trace $(du -h $V/verify/ALL_221010_ARC_INEQUALITY_WITNESSES.jsonl.gz 2>/dev/null | cut -f1)); delete with: rm -rf $V"
exit $rc
