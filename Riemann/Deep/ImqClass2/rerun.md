# Rerun recipe (imqclass2, N = 10^4)

Work directory: .scratch/imqclass2 of the repo (tablechecks.sh uses paths relative to it). One core for the sieve plus the downloader
(curl and the mjx extraction bursts); peak disk about 1 GB (two slabs of LMFDB files); peak memory about 700 MB
(mjx on slab 0) and about 450 MB (isieve2 near the top of the range). No pip, no venv, no Docker.

    cd W
    R=<repo>/Riemann/Deep/ImqClass2
    cc -O2 -o mjx $R/mjx.c -lz
    cc -O3 -mcpu=apple-m4 -o isieve2 $R/isieve2.c -lm      # the production build (sha256 of source and binary: sieve_logs/build.sha)
    cc -O2 -o isieve_old $R/isieve.c -lm                     # ImqClass's sieve, only for the cross-check on W1
    cp $R/*.sh $R/*.py $R/chunks_order*.txt $R/dump.cfg $R/t4.gp . && chmod +x *.sh && mkdir -p dl dumps runs cln
    curl -s -L -o cln/postproc.gp https://raw.githubusercontent.com/iwaokimura/cln-imq-2k/v1.4/qsieve/postproc.gp

1. Tables k = 0..110 (d < 111*2^28 = 29,796,335,616; 444 files, about 53 GB streamed, at most two slabs on disk; 2 h 51 min
   on 6 Nov 2026 at 1.9-7 MB/s; dlrun restarts after server failures):

        perl -e 'alarm 19000; exec @ARGV' ./dlrun2.sh 111
        ./tablechecks.sh          # mjx --report (final_sha.txt, final_perh.txt) and the table checks -> checks.txt

   mjx stops with FAIL on any inconsistency (unexpected d, missing fundamental d, invariant product or chain, gzip
   error, wrong slab order). tablechecks.sh needs class_group_orders.txt (people.kth.se/~kurlberg/class_group_data/),
   PARI/GP, runs/p_x5.t4 and runs/x5.found (step 3) and the dump dumps/tail.txt written by mjx (dump.cfg).

2. Thresholds (python-flint; seconds):

        uv run --no-project --with python-flint python3 $R/threshold.py 6000 7000 7500 8000 8500 9000 9500 10000

3. Sieve, N = 10^4 on [2.8*10^10, 7.627*10^11] in 38 chunks, one at a time, each under alarm 1800 (about 3.7 h on one
   shared Apple M4 efficiency core):

        ./runprod.sh 10000 chunks_order.txt     # c03..c09 (then stopped by hand; the rest is in the next files)
        ./runprod.sh 10000 chunks_order2.txt    # c02, c01, x1, x2, x3, c10
        ./runprod.sh 10000 chunks_order3.txt    # x4, x5, c11..c33
        awk '{print $1}' runs/p_x5.t4 > runs/x5.t4d && D4=runs/x5.t4d N=10000 gp -q -f t4.gp < /dev/null > runs/x5.found
        cat chunks_order*.txt | sort -u -k2,2n > chunks_all.txt
        python3 assemble.py chunks_all.txt 10000 28000000000 thresholds.txt

   Each log ends with "DONE N=10000 [LO,HI]: lanes L = certF + cert1 + t3 (...)"; assemble.py checks every chunk
   (counts, no failures, lane count = number of d = 3, 7, 11, 15, 4, 8 mod 16 in the chunk) and the coverage.
   Only chunk x5 has tier-4 entries (33); t4.gp decides them (one field: d = 28,154,696,707, h = 9819).

4. Sieve tests: tests.txt (commands and outputs; wincheck2.py, wincheck.py, t4.gp on W1).
