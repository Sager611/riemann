#!/bin/bash
# post.sh -- certification runs after the scan, strictly one process at a time, each with a hard wall-clock limit
cd "$(dirname "$0")"
H=$(pwd); A=$H/repo/ramanujans-inequality-main/certify
P=${PYTHON:-python3}
( cd $A && $H/to.sh 600 $P coverage_check.py $H/out/scan_mine_0.txt ) > out/coverage.log 2>&1
( cd $A/arb && $H/to.sh 1800 $P transfer_arb.py $H/out/scan_mine_0.txt ) > out/transfer.log 2>&1
./to.sh 600 $P cell_detail.py > out/cell_detail.log 2>&1
./to.sh 900 $P zfr.py > out/zfr_selftest.log 2>&1
( cd $A && $H/to.sh 900 $P smoothing_const.py ) > out/smoothing_const.log 2>&1
( cd $A/arb && $H/to.sh 1800 $P gammaint_cert.py ) > out/gammaint_cert.log 2>&1
./to.sh 600 $P zfr_lemma5.py > out/zfr_lemma5.log 2>&1
./to.sh 900 $P check_w3.py > out/check_w3.log 2>&1
./to.sh 3500 $P zfr_run1.py 6 1e-4 out/zfr_table_check.json > out/zfr_run_check.log 2>&1
echo ALLDONE >> out/post.done

