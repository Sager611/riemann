RECORD: Λ ≤ 0.107515584737 (4 of 67 rows certified: rows 1, 2, 3 in full, row 4 to u = 0.09; exact arc verifier PASS), i.e. 3.53·10⁻⁵ below the certified 0.107550869585, from a proved lemma (Raymond's source with the signed head and the unspent main slope; reviewer CORRECT) and Arb-checked branch-and-bound head sups.
# dbnrec: record certificate tranche 1 (N2) for de Bruijn–Newman Λ via near-band head credits in Raymond's chain
Orca agent dbnrec, 2026-11-20, from 09:22 UTC. No Lean; edits only in Riemann/Deep/DBNRec/ and .scratch/dbnrec/.

## Imported estimates (lesson 344): statement, range, source line
U.tex = Raymond, Zenodo 23752103, Release_Support/reconstruction/baseline_originals/Unconditional_0107965.tex.

| # | imported statement | range used | source |
|---|---|---|---|
| I1 | source identification: J_i ≥ Mc − Re P_{t,c}(χ) − (fully paid analytic + physical debit) | every row, t in its cell, x ≥ X | U.tex l.3329–3332 (Thm src:source l.1011, its proof l.1039–1620) |
| I2 | affine native: N_actual(L) ≤ N(0) + βM(L − L₀)/4, L ≥ L₀; unspent main slope exactly (1 − β)M | retained rows ("R", l.3619ff), all L ≥ L₀ | l.1495–1509, l.3209–3219 |
| I3 | refund keeps β and its slope; final native debit N_i = D_i^sc (rows ≠ 0,16,44); clock/MAIN gains change only first/last head fields | same rows | l.1523–1525, l.3241, l.3266–3274, l.3605–3607 |
| I4 | d_t(n) ≥ 0, nondecreasing in t; triangular recurrence; d_t(p) = log p·b_t(p) | 2 ≤ n ≤ 4096, t ∈ [0, 0.109] | Lemma src:coeff l.1266–1272, l.3311–3320 |
| I5 | head = Re P_{t,c}(χ) = Σ_{n≤4096} ½ d_t(n) n^{−1/2−tc} F(n) cos(τ log n), τ = x/2 + πt/8, F(n) = Σ_j w_j n^{−h_j/2} cos(a_j log n/2) | same | l.1595–1604, l.3301–3325, clock l.3427, l.3543 |
| I6 | X_*/(4π) > N₀², X_* = X − 4 (so Δ = c − c₀ ≥ ¼ log(x/X)) | x ≥ X | l.1090–1093 (rechecked in Arb: check.py) |
| I7 | M_i = main_lo ≤ Mc₀, M = Σ_j w_j (exact; mass and β from ALL_67_ORIGINAL_AFFINE_NATIVE.json, `retained_main_mass` = (1 − β)M) | certified rows | check.py; prep_end.py asserts |
| I8 | Raymond's 67-row flow producer and exact arc verifier, restarted from θ₁₄₁ (exact FLINT rationals, 2^−96 outward rounding) | whole chain | HP141/REPORT.md §DBNB2 (make_flow_b2.py, verify_arcs_141.py; PASS 220,162 arcs) |
| I9 | θ₁₄₁ = 0.87334585610671… (zero-free half-plane start) | flow start | PR 141, via HP141 |
| I10 | dbnhead bb.c: exact uint64 phases (drift < 2^−34 turn, EPSU 2^−24, PHERR), sampled tables with midpoint + Lipschitz bounds; reviewer CORRECT after the PHERR fix | all bb runs | DBNHead/REPORT.md §Far shells (review) |

## What is proved and how (LEMMA.md has the proof; reviewer: GAP on the x′ domain, already satisfied → CORRECT)
- Lemma (L1): for a retained scalar row, J_i(t,x) ≥ M_i − N_i + (1 − β)MΔ − Q_i(t,x), Δ ≥ ¼ log(x/X), Q_i the signed head with actual phases (I1–I3, I6). Hence any rational A′ ≥ B′ > 0 with Q_i − s_i log(x/X) ≤ (1 − θ)A′ + θB′ (all t in the cell, x ≥ X; s_i = (1 − β)M/4) may replace (first_hi, last_hi); Raymond's credits for that row are discarded.
- Correction to the brief's step 3/4 ("lower first_hi/last_hi by credit_i − Raymond's credit_i"): over one time cell the head charge falls by A_i − B_i = 0.026–0.046 (the chord in t, l.3331–3345), as large as the credits. A box sup over the whole cell (dbnhead) certifies only a constant A′ = B′ = sup, which loses ≈ (A − B)/2 per row against the uniform shift that dbnhead's float endpoint assumed. Fix used here: termwise, Q_i(t,x) ≤ (1 − θ)Q⁻(x′) + θQ⁺(x′) ≤ Q⁺(x′) + (1 − θ)D with Q^± the end-of-cell bound functions (upper chord, common lower bound) and D = Σ_n (H⁻_n − H⁺_n) (LEMMA.md Part B). One branch-and-bound run per shell on Q⁺ gives B′ = S, A′ = S + D: the full time decay is kept.
- Exactness: d_t(n) (triangular recurrence), F_i(n) (exact rational nodes) and every enclosure in Arb (python-flint, 160 bits), rounded outward to doubles (prep_end.py); phases exact uint64 from 50-digit mpmath (dbnhead); bb_rec prunes only when the float bound is ≤ thr − 10⁻⁹ (float summation error < 10⁻¹¹: ≤ 4096 terms of size ≤ 1, relative error 2⁻⁵³ each; table entries carry +10⁻¹⁵ and the Lipschitz factor 1.0000001); every unresolved final cell is re-evaluated in Arb at 128 bits (certify.py: upper endpoints of the balls lo·cos, hi·cos, plus the Lipschitz constant times radius + location error), and its slope credit uses log((x′ − r − ε − πt₊/4)/X) in Arb.

## Row ranking (brief step 3; DBNHigh chain.py sensitivities logs/chain.json, ceilaudit5 selfcons credits, cost ∝ X(e^{u_max} − 1))
- Gain per CPU-minute is flat over rows 0–19 (cells of length 0.0025: |dE/dS| 1.7–2.3·10⁻⁴, model gain 1.1–1.6·10⁻⁵ per row, u_max 0.17–0.25) and half as large for rows 20–66 (cells 0.00125, |dE/dS| ≤ 8.5·10⁻⁵, falling to 0 at row 66).
- Rows 0, 16, 44 carry Raymond's joint-jet overlay (conic functional α·scalar + 2·jet, l.3588–3593) and are outside LEMMA.md's scope; skipped.
- Order used: rows 1, 2, 3, … (non-jet, smallest u_max ≈ 0.16–0.17, so the cheapest full certificates), one sequential queue (code/master.sh).
- Cost measured: 3.5·10⁷ anchors/s, ≈ 4–5 min per shell of width 0.02 in u (≈ 1.2·10¹¹ in x), ≈ 40 min per row to u_F ≈ 0.16; the far shells do not get cheaper (the anchor loop dominates).

## Results
### Row 1 (t ∈ [1/400, 1/200], s = (1 − 1/512)M/4 = 0.31934597, M = 15806563055613/12350000000000 = 1.2798836482…)
- 10 shells u ∈ [0, 0.01, 0.02, 0.035, 0.05, 0.07, …, 0.15, 0.1665], thresholds thr_j = T + s·u_j (shell 0: 0.7715, then T = 0.774037 from shell 0's dump); bb_rec 1650 s CPU in total, memcap peak 39 MB each; prep_end (Arb) 8 s, ≤ 45 MB max RSS each.
- Shell 0: 861 unresolved cells of length D/32768; binding cell at u = 0.00976 (x′ − x′_s = 5.88·10¹⁰): float UB 0.777152, Arb UB 0.777228 (old dump format, location error ≤ 10⁻⁴ charged as Lipschitz slack), slope credit s·u = 0.003116 → S₀ = 0.774112724; witness 0.776740 at the same point (bb's qact: the head with the t₊ coefficients and the amplitudes of x = X e^{u₀} = X, at that x′; with the amplitudes of the actual x there, 0.7767242, orchestrator R513). Shells 1–9: no cell above threshold, so S_j = T = 0.774037. D_j = 0.023897 → 0.024033 (rising with L_lo); far range u ≥ 0.1665: a_F = 0.794484, b_F = 0.770451.
- Certified fields (exact, 2⁻⁹⁶ grid, logs/cert_r1.json): A′ = 0.798070454 (from the last shell, S + D₉), B′ = 0.774112724; Raymond's A = 0.841700097, B = 0.815596449 (his arithmetic payments for row 1: C₁ = 0.008310, g₁ = 0.000725, U.tex table l.3621): dA = 0.043630, dB = 0.041484, about 5 times his payments.
### Exact chain with row 1 only (chain A, superseded by chain B; it validated the pipeline at 10:24)
- mkrec.py: MAIN_CLOCK_SOURCE_REC.json = Raymond's source with row 1's first_hi/last_hi replaced (all other fields and rows asserted identical), sha256 b31a069a…2101a.
- make_flow_rec.py (= HP141 make_flow_b2.py, SOURCE path/hash swapped; θ₁₄₁ start, step 1/2048000, 2⁻⁹⁶ outward): 220,142 arcs, endpoint T + q_T/2 = 2130070039411895523492855447311/19807040628566084398385987584000 = 0.10754105468637595, headline ceiling 107541054687/10¹² = 0.107541054687; startup gap 9.31·10⁻¹³ > 5·10⁻¹³; 73 s, memcap 24 MB, max RSS 40 MB.
- verify_arcs_rec.py (= HP141 verify_arcs_141.py with source/flow hashes and the pinned terminal constants of the new flow; the per-row arc-count check generalized to any stopping row; no mathematical check changed): **"PASS": true**, 220,142 arcs, 5,776,508 exact checks, endpoint identical; 251 s, memcap 44 MB, max RSS 60 MB; witness trace sha256 5727d148… (423 MB, deleted; regenerable).
- Gain over the certified 0.107550869585: 9.815·10⁻⁶ (float model dE/dS ≈ 2.31·10⁻⁴ × 0.0426 ≈ 9.8·10⁻⁶: consistent).

### Row 2 (t ∈ [1/200, 3/400], s = 0.32476842)
- Shell 0 at T₀ = 0.788556 (tri⁺ − 0.965·cr*): 729 unresolved cells; binding cell again at u = 0.009756, x′ − x′_s = 5.88·10¹⁰ (the same point as row 1: the phases of the small primes align there), Arb UB 0.79391479 (float 0.793914797); location-aware S₀ = 0.790746348 ≤ T = 0.790748 (T carries a +2·10⁻⁶ planning margin). Shells 1–9 (to u_F = 0.1623): no cell above threshold. bb_rec 1367 s CPU, memcap peak 39 MB.
- Certified: A′ = 0.815478225 (S + D₉), B′ = 0.790748000; Raymond's A = 0.860827810, B = 0.833960859 (payments C₂ = 0.008846, g₂ = 0.000770): dA = 0.045350, dB = 0.043213. D₀ = 0.024595; far a_F = 0.814388, b_F = 0.789657.
- certify.py first failed on row 2 (assert: the float-chosen maximum of the shell values was not the exact maximum, since S_j = thr_j − s·u_j agree to 10⁻¹⁷); fixed by taking exact Arb maxima (row 1's run had passed the same assert); rerun by hand at 10:56 (13 s, 19 MB; it overlapped row 3's bb run: the only time two compute processes ran).

### Row 3 (t ∈ [3/400, 1/100], s = 0.33084008)
- Shell 0 at T₀ = 0.811148: 661 unresolved cells; binding cell once more at u = 0.009756 (anchor 3244584625, x′ − x′_s = 5.88·10¹⁰), Arb UB 0.81636463; S₀ = 0.813136951 ≤ T = 0.813139. Shells 1–9 (to u_F = 0.1676): no cell above threshold. bb_rec 1488 s CPU, memcap peak 39 MB.
- Certified: A′ = 0.838799632, B′ = 0.813139000; Raymond's A = 0.886329286, B = 0.858447493 (C₃ = 0.009464, g₃ = 0.000825): dA = 0.047530, dB = 0.045308. D₀ = 0.025527; far a_F = 0.837515, b_F = 0.811854.

### Row 4, partial (t ∈ [1/100, 1/80], s = 0.33690724; shells to u = 0.09 only, stopped by the shell deadline)
- Shell 0 at T₀ = 0.832387: 586 unresolved cells; binding cell at u = 0.009756 (same point), Arb UB 0.83745754; S₀ = 0.834170669 ≤ T = 0.834173. Shells 1–5 (to u = 0.09): no cell above threshold; shells 6–10 (to u_F = 0.1714) not run.
- Hence the far chord from u = 0.09 binds: A′ = tri⁻ − s·0.09 = 0.886993449, B′ = tri⁺ − s·0.09 = 0.860554154 (tri^± of the last shell, conservative); Raymond's A = 0.910052483, B = 0.881251491: dA = 0.023059, dB = 0.020697 (about half of a full row; the full row would give ≈ 0.047).

### Exact chain with rows 1–4 (chain B, the record)
- mkrec.py: rows 1, 2, 3, 4 replaced (dA = 0.043630, 0.045350, 0.047530, 0.023059; dB = 0.041484, 0.043213, 0.045308, 0.020697), every other field and row identical; sha256 04f16da4…e000.
- make_flow_rec.py: 220,090 arcs, endpoint T + q_T/2 = 212956555508380966927421960397/1980704062856608439838598758400 = 0.10751558473669763 (T = 22009/204800), headline ceiling 107515584737/10¹² = **0.107515584737**; startup gap 9.31·10⁻¹³; 90 s, memcap 24 MB, max RSS 40 MB.
- verify_arcs_rec.py: **"PASS": true**, 220,090 arcs, 5,775,364 exact checks, same endpoint; inputs plan ec346b7c…, source 04f16da4…, flow f767fdf2…; 274 s, memcap 44 MB, max RSS 60 MB; witness trace deleted (regenerable).
- Gain 3.528·10⁻⁵ over 0.107550869585 (rows 1, 2, 3: ≈ 1.0·10⁻⁵ each; partial row 4: ≈ 0.5·10⁻⁵).

## Footprints (one compute process at a time, single thread, every run under .scratch/memcap.sh 290 MB)
| step | runs | wall/CPU | memcap peak | max RSS |
|---|---|---|---|---|
| prep_end.py (Arb enclosures) | 41 | ≈ 8 s each | ≤ 1 MB sampled | 44 MB (test run) |
| bb_rec (branch-and-bound) | 41 shells | 60–274 s each; rows 1/2/3/4: 1650/1367/1488/1074 s | 39 MB | 39 MB |
| certify.py (Arb re-evaluation) | 5 | ≤ 13 s | 19 MB | — |
| make_flow_rec.py | 2 | 73 s, 90 s | 24 MB | 40 MB |
| verify_arcs_rec.py | 2 | 251 s, 274 s | 44 MB | 60 MB |
- One overlap: the manual rerun of certify.py for row 2 (13 s, 19 MB) ran beside row 3's bb run (logs/master.log 10:56, 11:18).

## Premise checks (code/check.py, Arb 200 bits; logs/check.out)
- I6: X − 4 − 4πN₀² = 5377388.879 > 0.
- I7, rows 1–7: retained, no jet; M = Σ_j w_j equals the AFFINE file's mass; (1 − β)M equals retained_main_mass; β = 1/512; Mc₀ − main_lo ∈ [3·10⁻³¹, 1.2·10⁻²⁹] > 0 (main_lo is Mc₀ rounded down on the 2⁻⁹⁶ grid).

## Rigor ledger of the branch-and-bound (what is interval arithmetic and what is float with an a-priori margin)
- Interval/exact: d_t(n), F_i(n), all enclosures lo/hi (Arb 160 bits → outward doubles); x′_s, phase increments (50-digit mpmath → uint64 turns, drift < 2⁻³⁴ turn ≪ EPSU = 2⁻²⁴, start-phase PHERR, dbnhead review); every unresolved final cell (Arb 128 bits, exact x′ from anchor index and offset); slope credits, D_j, far bounds, A′, B′ (Arb/fmpq); producer and verifier (FLINT rationals).
- Float with margin: the level-1–3 tables T_L[g][b] (samples of g_g with midpoint + Lipschitz bound; per-entry float error ≤ 3·10⁻¹⁵·Σ|c| against an added 10⁻¹⁵ + 10⁻⁷·½·lip·dθ ≥ 4·10⁻¹¹·Σ|c|) and the prune sums (≤ 4055 terms, error ≤ 10⁻¹² against the 10⁻⁹ margin in Gneed); libm cos assumed accurate to 1 ulp. Not re-done in interval arithmetic.

## Not done (precisely)
- Rows other than those listed above keep Raymond's fields; jet rows 0, 16, 44 are outside the lemma (their source is α·scalar + 2·jet, l.3588–3593; the jet part needs its own head filter); rows 20–66 (cells of length 0.00125) gain half as much per CPU-minute.
- The bb tables (levels 1–3) and prune sums are float with an a-priori margin (Rigor ledger), not interval arithmetic; cos from libm (1 ulp).
- Imported, not re-proved: Raymond's source identification and debit ledger (I1–I3), the coefficient lemma (I4), the flow producer/verifier (I8, as in HP141), θ₁₄₁ (I9).
- Not tried: narrower shells around the binding cell (row 1: the box over [0, 0.01] costs ≈ 4·10⁻⁴ against the witness), splitting the time cell, a Diophantine enumeration of anchors (the anchor loop is the whole cost: 3.5·10⁷ anchors/s).

## Suggested lines
- RECORDS.md (N2): the certified Λ bound of the verdict line, with the rows used.
- Lesson: a head credit certified as a sup over a whole time cell lowers first_hi and last_hi only to that constant; Raymond's chord falls by 0.026–0.046 per cell, so lowering both by one box credit (the brief's step 3–4 and dbnhead's float endpoint) is not justified. Certify the end-of-cell head (t₊ coefficients) and add the chord difference D at t₋.
- No NOGO line.

## Files (run from the repo root; inputs fetched by Riemann/Deep/DBNHigh/rzip.py into .scratch/dbnrec/dl and .scratch/dbnrec/ws, hash-checked)
- LEMMA.md (lemma, proof with U.tex lines, Part B construction, both reviews).
- code/: prep_end.py (Arb enclosures of the end-of-cell head Q⁺ on a shell, chord difference D, exact uint64 phase data), bb_rec.c (dbnhead bb.c + prune margin 10⁻⁹ + exact anchor index in the dump; cc -O2), run_row.sh (shells of one row, each process under memcap 290 MB), plan_row.py (float helper: thresholds, shell lists), certify.py (Arb re-evaluation of every unresolved cell, S_j, A′, B′ on the 2⁻⁹⁶ grid), check.py (premises I6, I7), mkrec.py (new source file with only the certified rows' first_hi/last_hi replaced; writes make_flow_rec.py and verify_arcs_rec.py as hash/constant swaps), make_flow_rec.py, verify_arcs_rec.py, master.sh (the sequential queue), summarize.py (tables), endpoint_rec.py (float chain estimate, not run). Unchanged copies: bb.c, prep.py, endpoint.py (dbnhead), make_flow_b2.py, verify_arcs_141.py (HP141).
- logs/: credits.json (the certified first_hi/last_hi of rows 1–4, exact, with Raymond's values and dA, dB), r<row>_<shell>.prep.json / .bb.out / .unres.txt (every shell), progress_r<row>.txt, cert_r<row>.json, check.out, mkrec_*.json, make_flow_{A,B}.log, FLOW_REC_SUMMARY_{A,B}.json, verify_{A,B}.log, EXACT_ARC_VERIFICATION_{A,B}.json, master.log.
- Deleted at the end: Zenodo downloads, the flow files (FLOW_REC.json.gz, regenerable by make_flow_rec.py in ≈ 75 s) and the 423 MB witness traces; .scratch/dbnrec/ keeps STATE.md.

RH relevance: N2
