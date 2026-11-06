# Rerun recipe (imqclass, N = 5000)

Work directory: any scratch directory W (we used .scratch/imqclass). Single core; peak disk about 1 GB; peak memory
about 700 MB (mjx on slab 0) and 300 MB (isieve). No pip, no venv, no Docker.

    cd W
    R=<repo>/Riemann/Deep/ImqClass
    cc -O2 -o mjx $R/mjx.c -lz
    cc -O2 -o isieve $R/isieve.c -lm
    cp $R/dl.sh $R/run5000.sh . && chmod +x dl.sh run5000.sh && mkdir -p runs

1. Tables (about 60 min at 10–17 MB/s; LMFDB serves the files from beta.lmfdb.org behind a cookie gate, human=1; the
   server cuts some transfers after about 40 s, dl.sh resumes them with HTTP ranges):

        perl -e 'alarm 3600; exec @ARGV' ./dl.sh 0 56      # rerun with the next K0 if the alarm stops it
        ./mjx --report state.bin final                      # final_perh.txt, final_sha.txt

   mjx stops with FAIL on any inconsistency (unexpected d, missing fundamental d, invariant product or chain, gzip error,
   wrong slab order). The rows N <= 5000 of final_sha.txt must equal the table in REPORT.md; the rows N > 5000 are partial
   (d < 56·2^28 only). To regenerate a list: add a filter to mjx (h <= N) or reuse merge_update's text lines.

2. Thresholds (python-flint; a few seconds):

        uv run --no-project --with python-flint python3 $R/threshold.py 2000 5000 10000

3. Sieve (about 50 min on one Apple M4 core; four chunks; fallback statements N = 4000, 4500 need chunks a–c only):

        ./run5000.sh            # runs/n5000_{a,b,c,d}.log, .t4 (empty), runs/n5000.status

   Each log ends with "DONE N=5000 [LO,HI]: lanes L = certF + cert1 + t3 (...)"; t4 must be 0 (otherwise decide the
   listed d with t4.gp).

4. Checks: hk_compare.py (with class_group_orders.txt from people.kth.se/~kurlberg/class_group_data/), regcheck.py and
   regx.gp (needs lists/imq_h2000.txt.xz of cln-imq-2k), wincheck.py (needs a window dump: echo "LO HI path" > dump.cfg
   before mjx processes that slab). Outputs as in tests.txt.
