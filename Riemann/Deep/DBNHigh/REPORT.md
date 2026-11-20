NO RESULT: the θ141 chain (Λ ≤ 0.107550869585) is bound by its finite head (labels n ≤ 4096) and its own pair R(0,√q), tied to the main mass by the tight far-field kernel limit; high-order h-derivative designs have no main or far-field term and act only at log n ≥ 10 (whole native debit worth 1.6·10⁻⁴) and on the own pair, so per-row LPs gain ≤ 0.024 in S − O and 2.6·10⁻⁵ in the float endpoint (0.1075244; 5.3·10⁻⁵ if every per-row gain were doubled), and orders ≳ 20–25 at h0 ≈ 4 cannot pay Raymond's Cauchy physical-error charge through the own pair.

# dbnhigh: high-order height-derivative probe designs vs Raymond's certified chain (N2)
Orca agent dbnhigh, 2026-11-20, from 05:05 UTC. Float tests unless marked. Notes appended as the work goes.

## 0. Inputs read
- DBNAudit/REPORT.md: chain = arcs 1/2048000, rate r = 2(S_i(t0) - O_i(q0)) - 1e-12, stop q = 1e-4, bound t + q/2; 98.66% of gain from signed multi-height sources; kernel = boundary inequality R_i(d, sqrt q) <= 4/(d^2+4q), Poisson-extended to 0 <= b <= y.
- CeilAudit3 B1: the witness rows (lattice / rigidified CUE) of explore51's 0.0957 floor fail designs of >= 25 height-derivatives at h0 ~ 4 (float); floor for multi-height phase-uniform data open (single-pair witness 0.071-0.077).
- Files: Zenodo 23752103 Proof_and_Replay.zip (201 MB) read by HTTP range (.scratch/dbnhigh/rzip.py): M8_LATE_PLAN.json (sha256 ec346b7c… OK), MAIN_CLOCK_SOURCE.json (b59905d7… OK), ALL67_OPERATIONAL_SOURCE_TABLE.json, DIRECT_SOURCE_RECONSTRUCTION.json, ALL67_SOURCE_DECIMAL_GUIDE.csv, ALL67_PHYSICAL_DEBIT_MAP.csv, ACTUAL_67_ROW_ARC_INPUTS.json; Unconditional_0107965.tex (adbd147a… OK).

## 1. What the chain needs of a source (TeX §2.3, §3.1, §3.3, App. B)
- Functional: L_i = finite signed sum of symmetric translations (w/2)[S_t(x+a,h)+S_t(x-a,h)], S_t(x,h) = -Im H_t'/H_t(x+ih), plus finite height-derivative ("joint jet") terms (rows 0, 16, 44 use orders 0-4 at p = 3.05 etc.). So h-derivatives of any finite order are inside the formalism.
- Dirichlet-side profile: v(u) = sum_j w_j e^{-h_j u/2} cos(a_j u/2) (u = log n); TeX's F_i(n) = v(log n). A k-th h-derivative design (1/k!)(-d/dh)^k at h0 has v(u) = (u/2)^k e^{-h0 u/2}/k!.
- Source (Thm src:source): J_i(t,x) >= M_i - N_i - (1-θ)A_i - θB_i for all x >= X = 6000000185827, with
  - main M_i = c0 * v(0) (c0 = L0/4 = 6.72302; row 0: v(0) = 1.26098, M_0 = 8.47760);
  - finite head A_i, B_i = (1/2) sum_{n=2}^{4096} d_t(n) n^{-1/2-tL/4} |v(log n)| (triangle sum over all 4095 labels, eq. direct67:raw-heads) minus arithmetic credits (row 0: raw 0.83887, credits 0.00800, first_hi 0.83091);
  - native N_i: labels n > 4096, Riemann-Siegel and analytic errors (row 0: 0.00735).
- Zero side: R(d,b) = L_i[K_h(d,b)] = int_0^inf v(u) cosh(bu/2) cos(ud/2) du; own pair O_i(q) = R(0, sqrt q).
- Kernel (Lemma fnd:poisson): R(d, sqrt q) <= 4/(d^2+4q) for all d and all q in [0, cap_i], and every probe pole outside the zero strip (h^2 > q_max).
- Rate on an arc: r = 2(J_i(t0) - O_i(q0)) - 1e-12 (TeX §4); comparison Q' = -2 - rQ; endpoint T + q_T/2.

## 2. Float model of the chain and of one source row (written from the TeX; scripts chain.py, model.py)
- chain.py (scheme 1 of HP141/DBNB2/float_check.py, rows from ACTUAL_67_ROW_ARC_INPUTS.json, q0 = (2θ141-1)^2+1e-12): endpoint 0.10755086958424454, 220,162 arcs (= DBNB2's float check; certified 0.1075508695842417). 5 s, 9 MB.
- Per-row sensitivity dE/dS_i (raise row i's source by 1e-3 on its whole time interval): -2.33e-4 (row 0) ... -1.74e-4 (row 19; rows 0-19 last 1/400), -8.5e-5 (row 20) ... -1.7e-6 (row 65; rows 20-66 last 1/800); sum over all rows -6.44e-3. So a 1e-4 endpoint gain needs a uniform source gain of 0.0155 in every row (0.1-0.2% of S_i = 7.64...18.37), or 0.43 in row 0 alone.
- Visited (t, q) per row: row 0 t 0-0.0025, q 0.5576-0.5356, S 7.640, O 1.448, r 12.38; row 30 q 0.143-0.138, S 12.21, O 3.47, r 17.5; row 60 q 0.0164-0.0135, S 17.15, O 6.29, r 21.7. O grows 1.45 -> 7.04: the own pair is the largest debit after the main term in every late row.
- Every row has sum_j w_j h_j = 1.998-1.9998, i.e. the far-field kernel limit d^2 R -> 2 sum w h <= 4 is tight in all 67 rows; v(0) = M grows 1.261 -> 3.069 (main 8.48 -> 20.63), so the effective height 2/M falls 1.59 -> 0.65 while the actual heights stay in [1.22, 6].
- model.py reproduces, for Raymond's own rows 0, 16, 30, 60: main_lo exactly (8.478318, 11.415451, 13.772436, 19.310454); the raw finite head (row 0: 0.838902 vs certified raw 0.838870; credits then give first_hi 0.830909); O(cap) and the kernel (row 0 min relative slack 3e-6 to 7e-6, far field d^2 R = 3.999994). The model's native tail sum over 4096 < n <= 690988 is 0.6-0.5x of the certified native debit (0.0049 vs 0.0074; 0.0287 vs 0.0463; 0.0304 vs 0.0622): the model is optimistic for designs living at large n.

## 3. Decomposition of the 0.10755 endpoint into its binding pieces (decomp.py; float chain, every row changed at once)
| change in every row | endpoint | dE |
|---|---|---|
| none (θ141 chain) | 0.1075508696 | - |
| native debit N_i -> 0 (everything at labels n > 4096, RS errors, activation) | 0.1073902654 | -1.61e-4 |
| native debit doubled | 0.1077125113 | +1.62e-4 |
| all arithmetic credits removed | 0.1076622405 | +1.11e-4 |
| finite head (n <= 4096, raw triangle sum) -10% | 0.1067464810 | -8.04e-4 |
| finite head -> 0 | 0.1003211884 | -7.23e-3 |
| main term +1% | 0.1067988137 | -7.52e-4 |
| own pair O_i(q) x0.99 | 0.1073752602 | -1.76e-4 |
| own pair O_i(q) x0.90 | 0.1058422597 | -1.71e-3 |
Binding pieces: the finite head (labels n <= 4096, i.e. log n <= 8.32) and the own pair O_i(q) = R_i(0, sqrt q), both coupled to the main mass v(0) through the kernel (far field sum w h = 2 tight in all rows). The whole native debit is worth 1.6e-4.
Consequence for the lead: a k-th h-derivative at h0 has v(u) = (u/2)^k e^{-h0 u/2}/k!, peaked at log n = 2k/h0; ceilaudit3's designs (h0 ~ 4, k = 20-206) peak at log n = 10-100, i.e. in the native range or beyond the native cutoff log N0 = 13.45 (only the lowest orders, k ~ 20, reach into the head range log n <= 8.32 with their lower tail). They have v(0) = 0 (no main term) and, for k >= 2, no far-field term; their zero-side kernel is concentrated at |d| < ~2 (R_k(d) ~ (1+d^2/h0^2)^{-(k+1)/2}), so they cannot relax the tight kernel points (d ~ 9 at q = cap, far field). Their channels: (a) cancel the native tail of the combined profile v (ceiling 1.6e-4 for removing ALL native debit, impossible in practice), (b) lower O through R(0,y) at the price of band at large n, (c) for k ~ 20, some cancellation in the head's signed sums. The LP below measures (a)-(c) together.

## 4. Per-row design LP (lp.py) — does the high-order class raise the best source S_i - O_i(q)?
LP per row i at its start time t_i and its visited mid q*: maximize F = c0 v(0) - sum_k W_k |v(u_k)| - R(0, sqrt q*) over signed weights (|w| <= 100 in column-normalised units) of a dictionary, s.t. R(d, sqrt q') <= 4/(d^2+4q') (cutting planes from a fine grid: d in [0, 3000], 9 values q' in [0, cap]; final fine violation <= 2e-5 relative unless noted) and far field sum (d^2 R) <= 4. W_k: exact d_t(n) n^{-1/2-t c0}/2 for n <= 150, u-bins for n <= N0 = 690988, and ("tail") a PNT continuation (1/2) e^{u/2} e^{-t u (L0-u)/4} du for log N0 < u < 40.
Dictionaries: "base" = the row's own heights x shifts 0, 0.195, ..., 1.95 (+ Raymond's row as one column "ray"); "jets" = k-th h-derivatives (1/k!)(-d/dh)^k at h0 = 3 (k = 1, 2, 4, ..., 30) and h0 = 4 (k = 1, 2, 4, ..., 40), shift 0.
Pricing of frequencies beyond the native cutoff is essential: without the tail the LP puts the jets beyond log N0 for free and reports row 0 F = 7.69 (O = -3.68), an artifact; at t = 0 the tail is the exact zeta'/zeta series (PNT), for t > 0 a proxy (the RS truncation of Raymond's framework costs more: divisor-type coefficients, plus his native is 1.5-2x the triangle tail), so the model is optimistic for jets.
| row (t) | Raymond's row, model | base+ray LP | + jets LP | jets gain |
|---|---|---|---|---|
| 0 (0) | 6.18889 | 6.19365 | 6.19392 (jets h0 3, 4, k <= 30/40 step 2) | +0.0003 |
| 10 (0.025) | 7.04550 | 7.05926 | 7.06190 | +0.0026 |
| 20 (0.05) | 8.10510 | 8.13300 | 8.13765 | +0.0047 |
| 30 (0.0625) | 8.74140 | 8.76107 | 8.76404 | +0.0030 |
| 40 (0.075) | 9.41610 | 9.44266 | 9.44741 | +0.0048 |
| 50 (0.0875) | 10.13830 | 10.19446 (fine kernel violation 6.9e-3: base not converged) | 10.19937 | (+0.0049) |
| 60 (0.1) | 10.87332 | 10.95032 | 10.96633 | +0.0160 |
Mechanism of the small jet gain (row 60): O falls 6.268 -> 6.224 for native +0.014: at late t the heat factor makes localized negative v at log n ~ 8-11 lower the own pair faster than it raises the band (ratio 2cosh(yu/2) e^{-u/2} e^{tu(L0-u)/4} > 1 for t > ~0.07, u ~ 8-10; < 1 at every u for t = 0).
Rich jets (h0 2, 3 with k = 1..20/30, h0 5, 6 with k <= 40 step 2) on more rows (sweep4.jsonl): row 5 6.60314 -> 6.60701 (+0.0039), row 15 7.56614 -> 7.57030 (+0.0042), row 25 8.43658 -> 8.44270 (+0.0061), row 65 11.34993 -> 11.37366 (+0.0237).
Richer jet dictionaries saturate (row 60: h0 3,4 with k = 1..30 and shifts 0, 0.585 -> 10.96801; h0 2,3,5,6 -> 10.96891, i.e. +0.018-0.019; row 40: 9.44810 / 9.44846, +0.0055-0.006). Memory peaks 199-232 MB (memcap), max RSS <= 255 MB.
Raymond-style native accounting (ENV=1: labels n > 4096 and the tail charged by the per-node envelope sum_l |w_l| |v_l|, as in his per-node native charges and the spectral envelope |F_i(λ)| <= sum c_l e^{-h_l λ/2}; head keeps signed sums):
| row | Raymond's row, model | base+ray | + jets (h0 2,3,4,5; k <= 20/30/40/40) | jets gain |
|---|---|---|---|---|
| 0 | 6.18876 | 6.19204 | 6.19238 | +0.00034 |
| 20 | 8.10391 | 8.12384 | 8.12427 | +0.00043 |
| 40 | 9.29752 | 9.42789 | 9.42966 | +0.0018 |
| 60 | 10.87082 | 10.91549 | 10.92239 | +0.0069 |

## 5. Endpoint with the new class (gain.py: chain rerun with per-row source gains dS_i, linear interpolation in the row index, flat beyond the last measured row)
| case | dS_i | float endpoint | gain |
|---|---|---|---|
| A: jets gain, signed-|v| native (optimistic) | 0.0004 (row 0), 0.0026, 0.0047, 0.0030, 0.0058, 0.0049, 0.0186 (rows 10-60; richest dictionaries for rows 40, 60) | 0.1075284674 | 2.24e-5 |
| B: jets gain, envelope native | 0.00034 ... 0.0069 | 0.1075443924 | 6.5e-6 |
| C: flat 0.006 rows 0-49, 0.02 rows 50-59, 0.03 rows 60-66 | - | 0.1075068601 | 4.40e-5 |
| D: jets gain, signed-|v| native, rich dictionaries where run (adds rows 5, 15, 25, 65) | 0.0003, 0.0039, 0.0026, 0.0042, 0.0047, 0.0061, 0.0030, 0.0058, 0.0049, 0.0186, 0.0237 (rows 0, 5, ..., 30, 40, 50, 60, 65) | 0.1075243970 | 2.65e-5 |
| E: twice D | - | 0.1074979569 | 5.29e-5 |
All below the brief's 1e-4 gate, so step 3 (certification) was not started.
Side observation (not the lead, model-dependent): re-optimising the weights over each row's own heights with a finer shift grid (base vs Raymond's row in the same model) gains 0.003-0.077 per row in the signed-|v| model, but those LP designs saturate the weight bound (|w| = 100) and in Raymond's envelope native accounting much of it disappears (row 0: +0.0033; row 60: +0.045; row 40's model of his own 24-node row is 0.12 below his certificate, so the envelope model is too pessimistic there). Not evidence of slack in the record.
Final reruns (final_lp.jsonl): row 60 envelope 10.91549 -> 10.92239, signed-|v| + jets 10.96891; row 0 signed-|v| 6.19365 -> 6.19392. Every LP design checked (final_lp.jsonl, sweep4.jsonl) has O(q) nondecreasing on [0, cap] (min step +6.9e-5 ... +8.0e-4 over 200 steps), so the arc lemma holds for them; Raymond's Cauchy jet charge (section 6) is not in the LP, so its jet gains are upper bounds.

## 6. Admissibility check (step 1 of the brief)
Sufficient conditions (the ones the chain uses; TeX Lemma fnd:poisson, Thm src:source, eq. src:physical-guard, src:jet-physical, arc lemma) for a row L + sum_k x_k D_k, D_k = (1/k!)(-d/dh)^k at h0 with symmetric shift a:
- (K) boundary kernel for the COMBINED row: R(d, sqrt q) <= 4/(d^2+4q) for all real d and all q in [0, cap], with D_k[K](d,b) = (1/2) sum_{+-} Re[(h0-b+i(d+-a))^{-k-1} + (h0+b+i(d+-a))^{-k-1}];
- (P) h0^2 > cap (poles outside the zero strip; then the Poisson extension gives K_y(d,b) >= (y/2)R(d,b) for 0 <= b <= y);
- (G) guard on the jet's licensed Cauchy circle |ζ-p| <= r: heights in [1.51 - t ell/4, 6.01], |a| <= 39/20;
- (S) source charges for the combined profile v + sum x_k (u/2)^k e^{-h0 u/2} cos(au/2)/k!: main c0 v(0); head (signed sum before modulus, n <= 4096); native and differentiated analytic remainders (spectral envelope, frequencies >= log 4097); the jet physical charge (η/(pr)) sum_r' (r'+1)|b_r'|(p/r)^r';
- (Ar) on each arc: Jlow(t) >= Jlow(t0) (head chords decrease in t) and O(q) nondecreasing in q.
So the designs ARE admissible in form (Raymond's rows 0, 16, 44 already carry jets of orders 4, 4, 8). Consequences (checked with mpmath, admis_check.py: closed forms to 7e-16 and 2e-41):
- v_k(0) = 0 for k >= 1: no main term (TeX: a jet's main term is b_0 c0/p only).
- Far field: lim d^2 D_k[K](d,0) = -2 (k = 1; positive x relaxes the tight far-field limit, as two close nodes do), 0 (k >= 2).
- Own pair: unshifted O_k(q) = (h0-sqrt q)^{-k-1} + (h0+sqrt q)^{-k-1} > 0, increasing in q; shifted O_{k,a}(q) = Re[(h0-sqrt q+ia)^{-k-1} + (h0+sqrt q+ia)^{-k-1}] can have either sign and slope (e.g. h0 = 4, a = 1.95, k = 3).
- Gain channels: the own pair and cancellation inside the combined head/native profile (signed sums before modulus). The LP of section 4 uses the combined profile, any sign, shifted and unshifted jets and k = 1, so it measures all channels at once.
- Physical charge (not in the LP): with p = h0 = 4, r <= 2, δ = e^-15, the Cauchy charge of one order-k term over its own-pair effect is 0.18, 2.6, 36 for k = 20, 25, 30 at y = 0.74 and 0.81 (k = 20, y = 0.5), 5.5 (k = 20, y = 0.12): ceilaudit3's orders 20-44, 60-68, 120-206 at h0 ~ 4 cannot pay for themselves through the own pair; only sharper error bounds (e.g. the exact zeta'/zeta series at t = 0) would remove this charge.
Review: one reviewer subagent, claims pasted inline. Verdict GAP on both claims as first stated: (i) "iff" too strong (now "sufficient"); (ii) missing arc condition Jlow(t) >= Jlow(t0) and the differentiated native/analytic remainders in (S) (added); (iii) the own-pair sign and monotonicity hold only for a = 0 (corrected above); (iv) positive x D_1 relaxes the far field (corrected); (v) head cancellation of the combined profile means a gain does not require x < 0, and Cauchy ratios > 1 exclude only own-pair use (corrected; the LP covers these channels). The reviewer confirmed the kernel-derivative and Dirichlet-profile formulas and the r = 2 Cauchy ratios. Not re-reviewed after the corrections.

## 7. Certificate
None: the float gain (section 5) is below the brief's 1e-4 gate, so no source was certified and the exact arc verifier was not rerun.

## 8. Verdict and suggested lines
- Step 1 (falsify the lead): the designs are admissible in form, but the chain's limiting factor sits elsewhere: the finite head (labels n <= 4096) and the own pair R(0, sqrt q), tied to the main mass by the far-field limit sum w h <= 2 that is tight in all 67 rows. High-order derivatives have no main term, no far-field term (k >= 2), a kernel concentrated at |d| < ~2, and live at log n = 2k/h0 >= 10, where the whole native debit is worth 1.6e-4 of endpoint; orders >~ 20-25 at h0 ~ 4 also cost more in Raymond's Cauchy physical-error charge than they save in own pair.
- Step 2 (float model): the best admissible jets found by per-row LPs add 0.0003 (row 0) to 0.024 (row 65) to S - O; float endpoint 0.1075243970 (gain 2.6e-5; envelope accounting 6.5e-6; every per-row gain doubled 5.3e-5). Below 1e-4: no step 3.
- Why CeilAudit3's exclusion does not transfer: there the designs test a fixed hypothetical zero row (two-sided bound |design - mean| <= Band at the row's coherent mode log n = L/2); a chain source needs a one-sided lower bound on the repulsion felt by the top zero, uniform in x, and a mean-zero design contributes only -Band to it.
- Suggested TREE N2 line: "Raymond-type chains (θ141 restart, Λ <= 0.107550869585) are bound by the finite head (n <= 4096; removing it gives 0.1003) and the own pair (10% less gives 0.10584); high-order height-derivative designs (h0 ~ 4, k >= 20) add 2.6e-5 (float; 5.3e-5 with every per-row gain doubled) [DBNHigh](Riemann/Deep/DBNHigh/REPORT.md)".
- Suggested NOGO line: "dbnhigh (N2): high-order h-derivative probe designs (h0 ~ 4, orders 20-206; ceilaudit3's row-excluding class) cannot lower Raymond's chain by 1e-4: mean-zero (no main term), no far-field term, they act only at log n >= 10 (native debit worth 1.6e-4 in total) and on the own pair (per-row LP gain <= 0.024, float endpoint gain 2.6e-5, 5.3e-5 if doubled), and orders >~ 20-25 fail the Cauchy physical-error charge; the binding pieces are the finite head (n <= 4096) and the own pair [DBNHigh](Riemann/Deep/DBNHigh/REPORT.md)".

RH relevance: N2

## Files
- Riemann/Deep/DBNHigh/: REPORT.md; chain.py (θ141 float chain + per-row sensitivities; args ACTUAL_67_ROW_ARC_INPUTS.json out.json); decomp.py (endpoint decomposition; args ARC_INPUTS DIRECT_SOURCE_RECONSTRUCTION.json out.json); model.py (one-row source model; validation run on rows 0,16,30,60); lp.py (per-row design LP; args ARC_INPUTS chain.json row dict [tail]; env HSET, XB, NCUT, NIT, ENV); gain.py (endpoint with per-row gains); admis_check.py (mpmath checks, Cauchy ratios); rzip.py (HTTP-range reader for Zenodo 23752103 Proof_and_Replay.zip).
- logs/: chain.out, chain.json (lp.py input), decomp.json, gain.json, sweep1.jsonl (rows 10-50), sweep2.jsonl (rich jets rows 40, 60), sweep3.jsonl (envelope accounting), sweep4.jsonl (rich jets rows 5, 15, 25, 65), final_lp.jsonl (final reruns), model_check.out, admis_check.out.
- Inputs (deleted after the run; re-fetch with rzip.py get <dir> <member>): M8_LATE_PLAN.json, MAIN_CLOCK_SOURCE.json, ACTUAL_67_ROW_ARC_INPUTS.json, DIRECT_SOURCE_RECONSTRUCTION.json, ALL67_OPERATIONAL_SOURCE_TABLE.json, ALL67_SOURCE_DECIMAL_GUIDE.csv, ALL67_PHYSICAL_DEBIT_MAP.csv, ARC_VERIFICATION_EXPLANATION.md, Unconditional_0107965.tex.

## Logs (memory, time)
One compute process at a time, OMP_NUM_THREADS=1, every run under .scratch/memcap.sh 290. Peaks: chain.py 9 MB (memcap) / 18 MB max RSS, 166 s; decomp.py 21 MB, 17 s; model.py 107 MB memcap / 124 MB RSS; gain.py 9 MB, 10 s; admis_check.py 35 MB RSS; lp.py runs 20-232 MB memcap, max RSS <= 255 MB, 5-180 s each (one 703 s run). Cap breaches: three early lp.py runs (254 columns, dense kernel grid) were killed by memcap at 382, 361 and 333 MB before cutting planes and smaller grids were introduced; no run finished above the cap. Total compute about 45 min.
