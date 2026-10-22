# EVIDENCE2: certified paired-probe floors (dbn5) and Planat's run
Every number in REPORT6.md is produced by the scripts here. Python with numpy and python-flint (the .scratch/dbn4 venv). Scratch and caches go to $DBN5_WORK (default .scratch/dbn5).

**Files**
- floors_dbn5.json: one record per distinct (cell, height), with a floor for each of the 10 shifts.
- certs_dbn5.json.gz: per (cell, height, shift) the 7 frequency subranges, each with [s1, s2, C_top, c, K, floor].
- classical_dbn5.json.gz and classical_dbn5_AUDIT.json: history and audit (Planat's unchanged code, U = 0.147373586543538).
- 07_planat_run_dbn5.txt: the producer and auditor run log.
- 08_tests.txt: interval arithmetic, toy brute force and the Euler bound.
- 09_mass_check.txt: bin masses against direct sums.
- 10_recheck.txt: an independent re-verification of a sample of the consumed floors.

**Steps** (repo root; $PY = .scratch/dbn4/venv/bin/python)
1. `cd Riemann/Deep/DBN/EVIDENCE2/scripts`, then run `$PY run5.py I 3 OUT_I.jsonl` for I = 0, 1, 2. About 30 CPU-min each. They write shard files under out/ in the working directory.
2. Merge the shards with merge5.py, run in the folder that holds out/.
3. Run `$PY s7_new_floor_run.py`. It rebuilds the sources from .scratch/dbn4v/evidence (hash-checked), then runs Planat's producer and auditor.
4. Run `DBN5_E2=.. $PY check5.py 48` for the sample recheck; tests5.py and tests5_mass.py for the tests.

**Premises still assumed:** A2, A3, and the analytic part of A1 (Raymond Foundations §4–5, §8.1 and (9.1); the positivity and monotonicity lemma for head coefficients; the envelope (8.1)–(8.3) used for s ≥ 8/5).
