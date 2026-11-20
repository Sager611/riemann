RECORD: Λ ≤ 0.107410439495 (12 of 67 rows certified: rows 1–12 in full; exact arc verifier PASS), 6.05·10⁻⁵ below tranche 2's 0.107470962867, by tranche 1's reviewed lemma and pipeline (plus a reviewed, output-identical two-group prefilter, 1.4–1.5× faster) on row 8 (completed) and rows 9–12.
# dbnrec3: record certificate tranche 3 (N2): near-band head credits on row 8 (rest) and rows 9–15, 17–19 of Raymond's chain
Orca agent dbnrec3, 2026-11-20, from 14:52 UTC. No Lean; edits only in Riemann/Deep/DBNRec3/ and .scratch/dbnrec3/. Tranches 1–2 (Riemann/Deep/DBNRec/, DBNRec2/) read-only.

## Setup 14:52–14:58 UTC
- code/ = DBNRec2/code with path strings swapped (.scratch/dbnrec2 → .scratch/dbnrec3, Riemann/Deep/DBNRec2 → Riemann/Deep/DBNRec3); changed lines: certify.py 1, check.py 3, mkrec.py 3, run_row.sh 8, summarize.py 1; others byte-identical. Tranche-2 orchestration scripts not copied (master*.sh, regen_r4.sh, equiv_hook.sh, final.sh, chain2.sh, credits2.py, foot2.py, rerun_tmp.sh); new: master3.sh, regen_r8.sh, equiv_hook3.sh, bb_rec3.c, final3.sh, chain3.sh, credits3.py, rerun_tmp.sh.
- Inputs refetched (rzip.py, Zenodo 23752103, RZIP_CACHE in .scratch/dbnrec3): sha256 M8_LATE_PLAN ec346b7c… (pin), MAIN_CLOCK_SOURCE b59905d7… (pin), ACTUAL_67_ROW_ARC_INPUTS 29902c66… (= tranche 2), ALL_67_ORIGINAL_AFFINE_NATIVE 61bde51a… (= tranche 2); rzip memcap 13 MB, max RSS 24 MB.
- bb build: cc -O2 bb_rec2.c (source byte-identical to tranche 2's) → .scratch/dbnrec3/bb_rec sha256 b6864793… (tranche 2's binary ae084653… not kept; binary hash differs, source and flags the same).
- regen_r8.sh: 6 of 6 prep.json lines of tranche 2's r8_00–05 IDENTICAL (also pins the two unhashed inputs); their bb.out/prep.json/unres.txt copied into logs/.
## Speed-up 14:59–15:03 UTC (within 20 min)
- bb_rec3.c = bb_rec2.c + exact two-group prefilter ANY2 (31-line diff): bitset NP[q][b1] of group-2 bins where candidate q survives gi = 1, 2 (same float ops); 11-bit tops of the group-1 and group-2 anchor phases (span 2^53 − 1 < 2^54: two bins each); ANY2 bit 0 → skip anchor with n1 += nl, as the original.
- reviewer: pass 1 GAP (only: printed lines carry elapsed CPU seconds, so whole-output byte identity is literal-false; speed not shown); pass 2 with the claim qualified (identical apart from elapsed-time fields, stripped as in tranche 2's test): **CORRECT**.
- Hook (code/equiv_hook3.sh, armed as .scratch/dbnrec3/hook.sh, run by master3.sh before row 9): old vs new on r8_05 with K = 6e8 (identical outputs apart from seconds, speed from bb's own CPU seconds), new on the full shell r8_00 vs tranche 2's bb.out + dump, mkrec.py path check (rows 1–8 of tranches 1–2 → chain F's source 6321ae69…), adoption by atomic rename only if identical and faster.
## Queue (code/master3.sh, from 14:56:58 UTC; logs/master.log)
- check.py rows 8–15, 17–19 (logs/check.out): I6 5377388.879 > 0; every row: mode retained, no jet, M = Σw = AFFINE mass, (1 − β)M = retained_main_mass, M c₀ − main_lo > 0; β = 1/512, rows 11–14 β = 1/128.
- Row 8 shells 06–10 at T = 0.913695 (tranche 2's shell list …, 0.09, 0.11, 0.13, 0.15, 0.17, 0.1868), bb_rec2 build: r8_06 149 s, r8_07 137 s CPU, no cell above threshold.
- Equivalence (hook 15:11:21–15:13:46 UTC, logs/equiv.txt): r8_05 with K = 6e8: bb_rec2 build 15 s, bb_rec3 build 10 s (bb CPU), outputs IDENTICAL apart from seconds (7 lines incl. 4 progress lines); new build on the full shell r8_00 (thr 0.913118): 90 s against tranche 2's 125 s, bb.out (seconds stripped) and the 350-cell unresolved dump IDENTICAL to tranche 2's; memcap peak 40 MB, max RSS 42 MB. mkrec.py with tranches 1–2's certs r1–r8 rebuilds chain F's source sha256 6321ae69…9aab (path check). **Adopted** 15:13:46 (atomic rename; .scratch/dbnrec3/bb_rec = bb_rec3 build, sha256 a36207e0…): every bb run from row 9 on uses it; row 8's shells 06–10 ran on the bb_rec2 build (b6864793…).
### Row 8 completed (t ∈ [1/50, 9/400], s = 0.36199570; T = 0.913695 as tranche 2)
- Shells 06–10 (u 0.09 → 0.1868), bb_rec2 build: no cell above threshold; bb CPU 149, 137, 127, 105, 57 s; memcap peak 39 MB. Shells 00–05: tranche 2's (shell 0: 350 cells, S₀ = 0.913692678).
- certify.py (logs/cert_r8.json; memcap 15 MB, max RSS 36 MB): A′ = 0.943522063, B′ = 0.913695000; far u_F = 0.1868: a_F = 0.941003, b_F = 0.911176. Raymond's A = 0.999627216, B = 0.967229000: **dA = 0.056105, dB = 0.053534** (tranche 2's partial row 8: 0.022604, 0.019954).
- Timeline (master.log): row 8 shells 06–10 14:58–15:11; hook 15:11–15:14; row 9 15:13:53–15:39:39, row 10 –16:03:38, row 11 –16:28:09, row 12 –16:59:27 (≈ 24–31 min per row at load 5–16); 16:59:28 new-row deadline → master3 done → chain G.
### Rows 9–12 (dorow: shell 0 at T₀ from plan_row.py, T from its dump, outer shells; bb_rec3 build; memcap peak ≤ 40 MB, max RSS ≤ 42 MB per bb run; certify memcap ≤ 15 MB, max RSS ≤ 36 MB)
- Row 9 (t ∈ [9/400, 1/40], s = 0.36905060): shell 0 at T₀ = 0.937000: 327 unresolved cells; binding cell anchor 3244584625 (u = 0.009756), Arb UB 0.94071100 (float 0.940711011), S₀ = 0.937110544; T = 0.937113; shells 1–11 (to u_F = 0.1946): no cell above threshold. bb CPU 1214 s (96+82+127+112+147+132+113+105+102+96+93+9); certified (logs/cert_r9.json): A′ = 0.968061928, B′ = 0.937113000; far a_F = 0.965001, b_F = 0.934052. Raymond's A = 1.028942841, B = 0.995313672: **dA = 0.060881, dB = 0.058201**.
- Row 10 (t ∈ [1/40, 11/400], s = 0.37654481): shell 0 at T₀ = 0.963992: 322 unresolved cells; binding cell anchor 3244584625 (u = 0.009756), Arb UB 0.96775515 (float 0.967755153), S₀ = 0.964081573; T = 0.964084; shells 1–11 (to u_F = 0.1993): no cell above threshold. bb CPU 1115 s (95+81+132+109+142+118+91+90+83+79+75+20); certified (logs/cert_r10.json): A′ = 0.996157403, B′ = 0.964084000; far a_F = 0.992781, b_F = 0.960707. Raymond's A = 1.059829839, B = 1.024978855: **dA = 0.063672, dB = 0.060895**.
- Row 11 (t ∈ [11/400, 3/100], s = 0.38266597): shell 0 at T₀ = 0.994391: 274 unresolved cells; every dumped cell has Arb UB − s log((x′ − r − ε − πt₊/4)/X) ≤ T₀, so S₀ = T₀ = 0.994391000 (worst: none); T = 0.994393; shells 1–12 (to u_F = 0.2149): no cell above threshold. bb CPU 1107 s (85+73+104+98+124+114+101+86+83+74+77+80+8); certified (logs/cert_r11.json): A′ = 1.028159892, B′ = 0.994393000; far a_F = 1.024045, b_F = 0.990278. Raymond's A = 1.098120775, B = 1.061361222: **dA = 0.069961, dB = 0.066968**.
- Row 12 (t ∈ [3/100, 13/400], s = 0.39025816): shell 0 at T₀ = 1.019367: 253 unresolved cells; every dumped cell has Arb UB − s log((x′ − r − ε − πt₊/4)/X) ≤ T₀, so S₀ = T₀ = 1.019367000 (worst: none); T = 1.019369; shells 1–12 (to u_F = 0.2191): no cell above threshold. bb CPU 1204 s (90+84+120+104+136+117+94+99+89+82+80+91+18); certified (logs/cert_r12.json): A′ = 1.054227648, B′ = 1.019369000; far a_F = 1.049699, b_F = 1.014841. Raymond's A = 1.126986096, B = 1.089028215: **dA = 0.072758, dB = 0.069659**.
- reviewer on the certify.py diff (tranche 2 → 3, inline): **CORRECT** ("changes only L and S_; no mathematics or certification logic changes"). LEMMA.md, bb_rec2.c and every other reviewed file unchanged; bb_rec3.c is a new file (reviewed above).

## Credit table (logs/credits.json: exact A′, B′; rows 1–3 tranche 1, rows 4–7 tranche 2, read-only certs)
| row | source | A′ | B′ | Raymond A | Raymond B | dA | dB | shells | u_F |
|---|---|---|---|---|---|---|---|---|---|
| 1 | tranche 1 | 0.798070454 | 0.774112724 | 0.841700097 | 0.815596449 | 0.043630 | 0.041484 | 10 | 333/2000 |
| 2 | tranche 1 | 0.815478225 | 0.790748000 | 0.860827810 | 0.833960859 | 0.045350 | 0.043213 | 10 | 1623/10000 |
| 3 | tranche 1 | 0.838799632 | 0.813139000 | 0.886329286 | 0.858447493 | 0.047530 | 0.045308 | 10 | 419/2500 |
| 4 | tranche 2 | 0.860699881 | 0.834173000 | 0.910052483 | 0.881251491 | 0.049353 | 0.047078 | 11 | 857/5000 |
| 5 | tranche 2 | 0.879826297 | 0.852524000 | 0.931199801 | 0.901548972 | 0.051374 | 0.049025 | 11 | 439/2500 |
| 6 | tranche 2 | 0.900470237 | 0.872323000 | 0.954090179 | 0.923512974 | 0.053620 | 0.051190 | 11 | 1803/10000 |
| 7 | tranche 2 | 0.925425090 | 0.896308000 | 0.981247488 | 0.949621090 | 0.055822 | 0.053313 | 11 | 921/5000 |
| 8 | tranche 3 | 0.943522063 | 0.913695000 | 0.999627216 | 0.967229000 | 0.056105 | 0.053534 | 11 | 467/2500 |
| 9 | tranche 3 | 0.968061928 | 0.937113000 | 1.028942841 | 0.995313672 | 0.060881 | 0.058201 | 12 | 973/5000 |
| 10 | tranche 3 | 0.996157403 | 0.964084000 | 1.059829839 | 1.024978855 | 0.063672 | 0.060895 | 12 | 1993/10000 |
| 11 | tranche 3 | 1.028159892 | 0.994393000 | 1.098120775 | 1.061361222 | 0.069961 | 0.066968 | 13 | 2149/10000 |
| 12 | tranche 3 | 1.054227648 | 1.019369000 | 1.126986096 | 1.089028215 | 0.072758 | 0.069659 | 13 | 2191/10000 |
- Rows 13–15, 17–19: not reached (the queue's new-row deadline 16:52 UTC passed during row 12); Raymond's fields kept.

## Exact chain with rows 1–12 (chain G, the record; code/final3.sh → code/chain3.sh, started 16:59:28 UTC when the queue ended at the new-row deadline, before the 17:07 stop time)
- mkrec.py: rows 1–3 from tranche 1's cert_r1–3.json, rows 4–7 from tranche 2's cert_r4–7.json, rows 8–12 from logs/cert_r8–12.json; all 12 rows replaced (dA, dB as in the table), every other field and row asserted identical; source sha256 0045c6f68e1f1fbe30c1c7b0b5ef646f9c167cfca1c6c4119c9dab6b3977fc97.
- make_flow_rec.py: 219,875 arcs, endpoint T + q_T/2 = 8509931755964727070902736459/79228162514264337593543950336 = 0.10741043949406991 (T = 1759/16384 = 219875/2048000), headline ceiling 21482087899/200000000000 = **0.107410439495**; startup gap 9.31·10⁻¹³; 104 s (wall 107 s), memcap 23 MB, max RSS 42 MB; flow sha256 b484ea2acb771e320309fc2e00b73c5f0217c0e2be82b39a262894ffa6fd7d09.
- verify_arcs_rec.py: **"PASS": true**, 219,875 arcs, 5,770,634 exact checks, all 67 rows covered, endpoint and headline identical; inputs plan ec346b7c…, source 0045c6f6…, flow b484ea2a…; 328 s (wall 334 s, load ≈ 9), memcap 44 MB, max RSS 60 MB; witness trace 420,842,212 bytes, sha256 de8a79a0c2e1e2c5a9d10717921e13216a5a1436dbcc02e5087f4646d0023f85, **uncompressed canonical transcript sha256 22facbf693085de8875c6b6e9ed2eafaea492dacee1c4d259530c50ee98d667c** (trace deleted; regenerable).
- Gain: 6.052·10⁻⁵ below tranche 2's 0.107470962867 (row 8 second half, rows 9–12; ≈ 1.2–1.4·10⁻⁵ per full row, rising with the row index); 1.051·10⁻⁴ below tranche 1's 0.107515584737; 1.404·10⁻⁴ below 0.107550869585.
- Rerun: code/rerun_tmp.sh (self-contained, /tmp/dbnv3 only: inputs refetched and hash-checked, code copy with every repo path sed-ed and checked, certs copied, every step under .scratch/memcap.sh 290 MB, verifier 1500 s; expects source 0045c6f6…, PASS, transcript 22facbf6…; prints all three); written, not run.

## Footprints (one compute process at a time, single thread, every run under .scratch/memcap.sh 290 MB)
```
row 8: 5 shells (tranche 3), bb CPU 575 s, bb memcap peak 39 MB, bb max RSS 39 MB, prep_end memcap peak 0 MB
row 9: 12 shells (tranche 3), bb CPU 1214 s, bb memcap peak 40 MB, bb max RSS 40 MB, prep_end memcap peak 0 MB
row 10: 12 shells (tranche 3), bb CPU 1115 s, bb memcap peak 40 MB, bb max RSS 40 MB, prep_end memcap peak 0 MB
row 11: 13 shells (tranche 3), bb CPU 1107 s, bb memcap peak 40 MB, bb max RSS 40 MB, prep_end memcap peak 0 MB
row 12: 13 shells (tranche 3), bb CPU 1204 s, bb memcap peak 40 MB, bb max RSS 40 MB, prep_end memcap peak 0 MB
cert_r10: memcap peak 15 MB, max RSS 33 MB, 6.57 s
cert_r11: memcap peak 15 MB, max RSS 34 MB, 6.64 s
cert_r12: memcap peak 15 MB, max RSS 34 MB, 7.44 s
cert_r8: memcap peak 15 MB, max RSS 34 MB, 7.62 s
cert_r9: memcap peak 14 MB, max RSS 33 MB, 8.26 s
```
- prep_end.py (Arb enclosures, ≈ 7–10 s each; 65 runs: 55 shells, 6 regenerations of row 8's shells 00–05, 4 row probes): memcap samples ≤ 1 MB (runs shorter than its sampling); check.py, mkrec.py ≤ 3 s. Equivalence runs: memcap 40 MB, max RSS 42 MB. Rerun script not run.
- No overlap of compute processes: the equivalence test ran as a queue hook between rows; the reviewer runs are LLM agents, not compute; log summaries (rowsum3.py, credits3.py, foot3.py) read logs only.

## Not done (precisely)
- Rows 13–15 and 17–19 keep Raymond's fields (the queue finished row 12 at 16:59 UTC, after the 16:52 new-row deadline); jet rows 0, 16, 44 remain outside the lemma; rows 20–66 gain about half as much per CPU-minute.
- Further speed-ups: residue-class (in-order Diophantine) enumeration of the anchors the group-1 table skips (≈ 80 % of anchors, now the dominant per-anchor cost); a three-group table; lower T₀ on rows where no dumped cell binds (row 11: S₀ = T₀).
- Same rigor ledger as tranches 1–2 (float tables with a-priori margins, libm cos); imported inputs I1–I10 as in DBNRec/REPORT.md.

## Files (run from the repo root)
- code/: tranche 2's code with path strings swapped (unchanged otherwise), plus bb_rec3.c (two-group prefilter), master3.sh (queue), regen_r8.sh, equiv_hook3.sh (equivalence test + adoption), final3.sh (stop rule + chain), chain3.sh (exact chain, = chain2.sh with paths swapped), credits3.py, rowsum3.py, foot3.py (log summaries), rerun_tmp.sh (from-scratch rerun of chain G under /tmp/dbnv3, nothing written in the repo; written for the orchestrator, not run here).
- logs/: r<row>_<shell>.prep.json/.bb.out/.unres.txt (rows 8–12; row 8 shells 00–05 copied from tranche 2), progress_r*.txt, cert_r8–*.json/.err, credits.json, check.out, regen_r8.log, equiv.txt, master.log, mkrec_G.json, make_flow_G.log, FLOW_REC_SUMMARY_G.json, mkrec_ver_G.json, verify_G.log, EXACT_ARC_VERIFICATION_G.json.
- Deleted at the end: downloads (.scratch/dbnrec3/dl, ws), .dat files, binaries, flow file and witness trace (regenerable: regen via prep_end.py, make_flow_rec.py, verify_arcs_rec.py); .scratch/dbnrec3/ keeps STATE.md.

RH relevance: N2
