CONFIRMED: a rigorous branch-and-bound sup of the actual x- and t-dependent head (all 4095 labels) puts rows 0, 30, 66 at 6.94%, 9.60% and ≥ 11.46% below the torus max on the binding windows (model 6.58/8.57/10.33%); to 3–4u* the certified credits are 0.05711 (model 0.05738 at κ 1.32), ≥ 0.13748 and ≥ 0.23798; float endpoint 0.1069862 unchanged; the other 64 rows and the shells beyond 3–4u* are still model.
# dbnhead: is the near-band sup of Raymond's finite head really 6–12% below its torus max? (falsification of ceilaudit5's N2 lead)
Orca agent dbnhead, 2026-11-20, 07:38–09:38 UTC. No Lean, no edits outside this folder and .scratch/dbnhead/.

## Object
- Head (Zenodo 23752103, Unconditional_0107965.tex l.1595–1604, common model l.1041–1077): Q_i(t,x) = Σ_{2≤n≤4096} c_n(t,x) cos(τ log n), c_n = (1/2) d_t(n) n^{−1/2−tL/4} F_i(n), τ = x/2 + πt/8, L = log(x/4π); F_i includes the joint-jet overlay (DBNHigh model.v_row); all 4095 labels, composites included. The source charges −Q_i; the worst x maximizes Q_i (ceilaudit5, reviewed).
- Enclosure over the row's time cell [t0,t1] and x ∈ [X, X e^u]: d_t(n) ≥ 0 nondecreasing in t (Lemma src:coeff, l.1267) and n^{−tL/4} decreasing in t, x give c_n ∈ sign F_i(n)·[½d_{t0}n^{−1/2−t1L_hi/4}|F|, ½d_{t1}n^{−1/2−t0L_lo/4}|F|] (relative margin 1e−11). Phase: τ log n = (x′/2) log n, x′ = x + πt/4 ∈ [X + πt0/4, X e^u + πt1/4].
- Bounded function: U(x′) = Σ_g g_g(θ_g(x′)), groups g = primes p (all powers p^k ≤ 4096 together) and composite labels, g_g(θ) = Σ_k max(lo_k cos kθ, hi_k cos kθ) ≥ every head value at x′ for every admissible (t, x).
- Windows u_i = ceilaudit5's selfcons shell optimum at κ = 1.32 (rows 0, 30, 66: u = 0.00726, 0.0109, 0.01061; x up to 1.0073X, 1.011X, 1.0107X; lengths 4.4e10, 6.6e10, 6.4e10).
- Gap = 1 − sup/U0, U0 = mx_pp + composite triangle (the full head's torus max bound; mx_pp = frozen prime-power torus max as in ceilaudit5, row 0 0.836455); the model's gap_pp(u*) is put in the same convention. The credit the chain uses is tri − sup (model: gap_pp(u*), composites at triangle).

## Method (code/prep.py, code/bb.c; rigorous up to stated float margins)
- Phases exact in uint64 turns (mpmath, 50 digits): anchor x′_k = x′_s + kD, D = 4π/log 2 (double), increments rounded to 1 unit (2^−64 turn); total drift < 2^−34 turn, covered by a 2^−24-turn margin in every lookup.
- Four cell levels: D/32, D/256, D/2048 (0.567, 0.0708, 0.00885), then Lipschitz bisection to 1e−3. Levels 1–3: per-group tables T[b] ≥ max g_g over the bin plus the cell's arc (8 samples per 1/1024-turn bin, midpoint+Lipschitz cell bound, deque window max); cell bound = Σ_{top m} T + Σ_rest max g_g (m = 48, 400, all 4055). Final level: Σ_g g_g(mid) + (Σ_g Lip_g ω_g)·len/2 + 1e−9.
- Output: CERT sup ≤ max(thr, max UB of unresolved 1e−3 cells); qact = head with the actual t = t0, x = X coefficients at the same x′ (witness, float).
- Test: row 0 on [X, X + 1.0e8] at thr 0.765 recovers ceilaudit5's fullscan maximum at X + 87571525.10 (UB 0.76970; qact 0.76542; fullscan frozen 0.76787), 1 s.

## Results
| row | t-cell | window (x − X) | U0 (torus max bound) | model gap_pp(u*)/U0 (κ 1.32 / 7) | certified sup | certified gap 1 − sup/U0 | qact (t0) at the sup | credit tri − sup (model gap_pp) | Raymond credit | time, memcap peak |
|---|---|---|---|---|---|---|---|---|---|---|
| 0 | [0, 0.0025] | 4.372e10 (u 0.00726) | 0.837838 | 6.58% / 6.18% | 0.779720 | 6.94% | 0.775349 at X + 2.927e9 | 0.059182 (0.055100) | 0.00799 | 148 s, 39 MB |
| 30 | [0.0625, 0.06375] | 6.576e10 (u 0.0109) | 1.539288 | 8.57% / 8.05% | 1.391447 | 9.60% | 1.385975 at X + 5.882e10 | 0.148577 (0.131886) | 0.02819 | 272 s, 39 MB |
| 66 | [0.1075, 0.10875] | 6.400e10 (u 0.01061) | 2.225970 | 10.33% / 9.77% | ≤ 1.970981 (= thr; no cell reached level 3) | ≥ 11.46% | — | ≥ 0.255000 (0.229967) | 0.02656 | 301 s, 39 MB |

- Column "certified gap" = 1 − sup/U0, U0 = mx_pp + composite triangle (the full head's torus max bound; row 0 0.837838, row 30 1.539288, row 66 2.225970); the model's gap is gap_pp(u*)/mx_pp = (cr* − s·u*)/mx_pp from ceilaudit5's selfcons.json.
- The credit-relevant test (the model's credit keeps composites at triangle, so it needs sup ≤ tri − gap_pp(u*)): row 0 0.779720 ≤ 0.783802 (margin 0.0041), row 30 1.391447 ≤ 1.408134 (margin 0.0167), row 66 1.970981 ≤ 1.996011 (margin ≥ 0.0250). The certified window credits tri − sup exceed even the model's full cr* = gap_pp + s·u* (0.057383, 0.137480, 0.237978).
- Convention trap: ceilaudit5's fullscan quoted 1 − sup_full/mx_pp (full head against the prime-power max). That is harmless in row 0 (composites 0.16% of the head) but wrong in late rows: row 30's composites are 4.8% of the head and sit in phase with the primes at the torus max (all phases 0), so that ratio gives 5.08% for row 30 against the model's 9.00%, a false alarm. Measured against the full torus max, or as a credit against the triangle, row 30 beats the model.
- Rows 0 and 30 attain their sups at distinct isolated points (row 0 at X + 2.927e9, the first quarter; row 30 at 0.89 of the window); the actual t0-value there is 0.0044 (row 0) and 0.0055 (row 30) below the certified bound (cost of the box over the time cell).

## Float endpoint (code/endpoint.py, DBNHigh chain.py)
- credit_i = min(tri − sup_i, cr*_i) on rows 0, 30, 66 (shells beyond u* stay model), cr*_i elsewhere: endpoint 0.10698609692 (κ 1.32) / 0.10702914357 (κ 7), identical to ceilaudit5's model endpoints, because all three certified window credits exceed cr*. Base 0.10755086958; gain 5.65·10⁻⁴ / 5.22·10⁻⁴.

## Far shells and binding point
- far1 (box enclosure, x ∈ [X e^{0.00726}, X e^{0.02904}], 1.33e11, thr 0.78268): certified sup 0.787696 at x/X = 1.0279 (qact 0.783446); 1269 unresolved cells dumped. Location-aware credit min over cells of tri − UB + s·log(x/X) (s = 0.31402) and tri − thr + s·u0 for the pruned rest: 0.055244, set by a cell at u = 0.00976 (UB 0.786722). So the binding point lies just past u*, as the model says, and the box enclosure certifies 0.0552 there: between the model's κ = 7 (0.05381) and κ = 1.32 (0.05738) credits. 418 s, 39 MB.
- The box loses ≈ 0.0044 to the time cell (qact vs UB at the row-0 maximum), because it pairs d_{t1} with n^{−t0L/4}. For primes d_t(p) = log p·b_t(p) exactly, so c_p = ½ log p·p^{−1/2} e^{t log p (log p − L)/4}|F_i(p)| is decreasing in t and L (log p < L): enclosure [c_p(t1, L_hi), c_p(t0, L_lo)] (prep.py tight = 1; composites and prime powers keep the box). Row 0 U0 drops 0.838301 → 0.836001.
- Review (reviewer, bb.c points A–C): partition, pruning, tables, wrap-around and the cached p = 2 key: correct; M1 = 32 ≤ 64 bounds the j-list. GAP: the final-level bound lacked an allowance for the start-phase error (< 2^−30 turn): 2π·2^−30·Σ lip_g ≤ 4.2·10⁻⁸ for all rows run. Fixed in bb.c (PHERR) after the runs above; every certified sup above holds with +1·10⁻⁷, which changes no digit used below.
- tight04 (tight primes, x ∈ [X, X e^{0.02904}] = 4u*, 1.77e11, thr 0.7815 = tri − cr*(κ 1.32) + 0.00002): certified sup 0.785881 at x/X = 1.0279 (qact 0.783449), 848 unresolved cells dumped, none with u ≤ u*. Location-aware credit 0.057109 (cell at u = 0.00976, UB 0.784857; pruned rest tri − thr = 0.057402): 0.5% below the model's κ = 1.32 credit 0.057383, 6% above its κ = 7 credit 0.053813. 473 s, 39 MB.
- Row 30, tight03 (tight primes, x ∈ [X, X e^{0.0327}] = 3u*, 1.99e11, thr 1.402544 = tri − cr*(κ 1.32)): no cell above thr, so the certified credit is ≥ 0.137480 (the model's κ = 1.32 credit) out to 3u*. 0 unresolved, 674 s, 39 MB.
- Row 66, tight03 (x ∈ [X, X e^{0.0318}] = 3u*, 1.94e11, thr 1.988003 = tri − cr*(κ 1.32)): no cell above thr, so the certified credit is ≥ 0.237978 out to 3u*. 0 unresolved, 852 s, 39 MB.
- Float endpoint with row 0 at 0.057109 (to 4u*; beyond: model), rows 30, 66 at ≥ model (to 3u*; beyond: model), model elsewhere (code/endpoint.py, logs/endpoint2.out): 0.10698616 (κ 1.32), 0.10702914 (κ 7); base 0.10755087.

## Assessment
- Falsification test passed: on every tested binding window the rigorous sup of the actual head (all 4095 labels, the x- and t-dependent coefficients, composites and joint jets at their true phases) sits at or below the model's level, never near the torus max. Certified gaps 6.94% / 9.60% / ≥ 11.46% (rows 0 / 30 / 66) against the model's 6.58% / 8.57% / 10.33% (same convention); all ≥ 5%, none < 2%. The extrapolation from log N ≤ 20 to log N ≈ 25–27 holds here.
- Near-far shells: row 0 to 4u* gives a location-aware certified credit of 0.057109, within 0.5% of the model's κ = 1.32 credit and above its κ = 7 credit; the binding point is u ≈ 0.0098 (model: 0.0065–0.0073). Rows 30 and 66 to 3u* are certified at the model's κ = 1.32 credits 0.137480 and 0.237978.
- Not shown: rows other than 0, 30, 66; shells beyond 4u* (row 0) and 3u* (rows 30, 66), still model; exact arithmetic (coefficients from float dt_coeffs with 1e−11 relative margins, float sums, review-fixed phase allowance); the reading of Raymond's source theorem that allows a signed head sup plus the slope s·log(x/X) in place of the triangle charge (ceilaudit5, reviewer CORRECT; not re-audited here).

## Steps of a record certificate (ceilaudit5's (i)–(iii), sharpened)
1. Per row: bb.c with tight enclosures over [X, X e^{u_max}], u_max = (cr − (tri − U0))/s (row 0: ≈ 0.18), location-aware credit min over cells of tri − UB + s·log(x/X). Measured cost: 150–470 s per row-window of 4–18·10¹⁰; the far shells have rising thresholds and prune faster. Estimate 67 rows × 15–40 min ≈ 17–45 h single-threaded, independent per row and shell. Composites and prime powers: split the time cell (box slack ∝ t1 − t0).
2. Exactness: coefficients d_t(n) and F_i(n) in Arb (python-flint) with outward rounding; re-evaluate the final cells (≈ 10³ per row) in interval arithmetic; phases are already exact (uint64 from 50-digit mpmath).
3. Raymond's exact arc verifier with raised sources: first_hi/last_hi lowered by credit_i − Raymond's credit_i, outward rounding, the 67-row chain rerun from θ₁₄₁ (HP141 restart). Float target 0.106986 (κ 1.32) to 0.107029 (κ 7); gain 5.2–5.6·10⁻⁴ over 0.1075509.
4. A written lemma (Raymond l.1505, 1600–1604, 1626–1740) that the source inequality holds with the signed head and the slope (1 − β)M·log(x/X)/4, certified β.

## Suggested lines
- No NOGO line: the near-band head credit survives falsification on rows 0, 30, 66.
- Lesson: compare a full-head sup with the full torus max (or as a credit against the triangle), not with the prime-power max. Late rows' composites (5–9% of the head) sit in phase with the primes at the maximum. The mixed ratio turned row 30's certified 9.6% into an apparent 5.1%.

## Files
- code/prep.py (coefficient enclosures, exact uint64 phase data; args ARC row u gfrac out [u0] [tight]), code/bb.c (branch-and-bound; args data m1 m2 [Kmax] [thr] [hmin] [dumpfile]; cc -O2 bb.c -lm), code/endpoint.py (float chain with certified credits; imports DBNHigh chain.py). Run from the repo root; ARC.json = ACTUAL_67_ROW_ARC_INPUTS.json fetched by Riemann/Deep/DBNHigh/rzip.py; selfcons.json from Riemann/Deep/CeilAudit5/logs.
- logs/: prep_*.json (row data), bb_r0_g0759.out, bb_r30_t138.out, bb_r66_c255.out (binding windows), bb_r30_g10_killed.out (first row-30 attempt, threshold set against mx_pp, stopped), bb_r0_far1.out + unres_r0_far1.txt, bb_r0_tight04.out + unres_r0_tight04.txt (unresolved cells: x′ − x′_s, UB), bb_r30_tight03.out, bb_r66_tight03.out (empty dumps: no unresolved cells), endpoint.out, endpoint2.out, certs*.json.
- Each bb run: one single-threaded process under memcap 290 MB; peak 39–40 MB (footprint 40.7 MB by /usr/bin/time -l); prep 18 MB.

RH relevance: N2
