OPEN: (N2) explore51's phase-uniform row ceiling 0.0957 (and Raymond-source 0.1018), freedom = multi-height designs (≥ 25 height-derivatives of one probe, same triangle-inequality lemma), float: both witnesses violate a valid ζ inequality at t = 0 (lattice row ratio 1.03 at orders 20–44, h₀ = 3.97, up to 2.64; all nine rigidified CUE rows 1.7–2.3), so the class floor is open (single-pair witness ≈ 0.071–0.077); (N3) R400's 'one worst-case zero' Li limit, freedom = explicit zero density, float: the binding step is ≈ 1.6·10¹¹ zeros near T, one zero allows 175 T², a KLN-shaped density with constant ≤ 10 gives 110–117 T² (constant only); (N9) Explicit611's floor 0.5075, already passed by Fodge's explicit B-process for t > 8.7·10¹² (constant only). ROBUST: R449's Davenport–Heilbronn barrier (Λ ≤ 0.2306 unreachable by arguments whose inputs all hold for DH).

# ceilaudit3: falsifying our own barriers on N2, N3, N5, N9, N13

Started 2026-11-20 04:13 UTC. Float tests only; nothing certified. Notes appended as the work goes.

## Barriers chosen
- B1 (N2) explore51's phase-uniform model ceiling 0.0957 = (16/L²)log cosh(L(2θ₁₄₁−1)/4) (ExploreRelax/REPORT.md l.1, 17, 96, 127; TREE N2 "single-height phase-uniform probes"; NOGO l.683).
- B2 (N3) R400: "the Li-range constant (106.25 T²) is limited by one worst-case zero above T near Re s = B, so explicit zero density cannot raise it" (NOGO l.546).
- B3/B4: chosen after B1, B2 (N5 nb2nd sign of C₂; N13 / N9 floors).

## B1 (N2): explore51's phase-uniform row ceiling 0.0957

Step fixed: the model step "a row of off-line pairs placed in one window above T₀ passes the chain's inputs", i.e. the witness of the barrier: the heat flow of cos(πρz) − cosh(πρy₀) (one pair per period P = 2/ρ, pairs at ±y₀ = 2θ₁₄₁ − 1 = 0.7466917 at t = 0) reaches the line at t_B = (16/L²)log cosh(Ly₀/4) = 0.0957314 (L = log(x/4π) = 26.891756 at x = 6·10¹², πρ = L/4).
Class as stated: ExploreRelax l.17 "(i) the strip and (ii) phase-uniform σ_eff > 1 data that the row satisfies"; tested only for single-height probes J(x;h) (l.63–75, margin 3.15); l.1, 96, 127 conclude "every phase-uniform input", "the floor for any phase-uniform design stays 0.0957", "none can go below 0.0957". TREE N2 and NOGO l.683 scope it to single-height probes but conclude "Raymond-type chains have about 5–11% room".
Freedom excluded by the class: linear combinations of probes at several heights (equivalently h-derivatives of one probe). The same lemma (absolute convergence + triangle inequality / Kronecker on Σ Λ(n)w(n)n^{−ix/2}) gives, for every real design v(u) = Σ_h c_h e^{−hu/2} (h > 1):
  ζ: |Σ_h c_h J_ζ(x;h) − πρΣc_h| ≤ Band(v) := ½ Σ_n Λ(n) n^{−1/2}|v(log n)|  for all x (independent phases per n: an upper bound for the torus sup).
  row: Σ_h c_h J_row(x;h) − πρΣc_h = Σ_{m≥1} 2πρ cosh(mπρy₀) v(m u*) cos(mπρx), u* = 2πρ = L/2 (Poisson sum per period), so its sup-norm ≥ A₁ = 2πρ cosh(πρy₀)|v(u*)| = 1017.94|v(u*)|.
  The row is excluded by this datum iff A₁ > Band(v). Single height v = e^{−h(u−u*)/2}: ratio A₁/Band = deficit/Pen of ExploreRelax (max ≈ 0.064 at t = 0, h ≈ 1.1–1.3).
Test (b1_base.py; exact Λ(n) for n ≤ 10⁷, PNT density e^{u/2}du beyond, to u = 250): v(u) = (u/u*)^K e^{−(K/u*+½)(u−u*)} Q(u), deg Q = m, i.e. the h-derivatives of orders K…K+m of one probe at h₀ = 1 + 2K/u*; Q by LP (min Band s.t. v(u*) = 1, HiGHS dual), then Band re-evaluated on all prime powers ≤ 10⁷.
| K | m | h₀ | Band | ratio A₁/Band |
|---|---|---|---|---|
| 0 (h₀ = 2) | 16 | 2.0 | 2584 | 0.394 |
| 20 | 16 | 3.97 | 1199 | 0.849 |
| 60 | 8 | 9.92 | 932.2 | **1.092** |
| 60 | 16 | 9.92 | 697.2 | 1.460 |
| 120 | 16 | 18.85 | 482.2 | 2.111 |
| 190 | 0 | 29.26 | 1016.7 | 1.001 |
| 190 | 16 | 29.26 | 385.0 | 2.644 |
Reading: the barrier's own witness (the periodic row at the strip edge) violates a valid phase-uniform σ_eff > 1 inequality of ζ once the design uses ≈ 60–70 height-derivatives (ratio 1.09; the bare 190th derivative alone gives 1.001). The row's coherent mode at frequency πρ sits at log n = L/2 (n ≈ √(x/4π) = 6.9·10⁵), amplified by cosh(πρy₀) = 75.7; a design localised in log n to width ≈ 2 sees it above the prime mass there.


B1 second witness. The rigidified CUE rows (ExploreRelax analysis3.py config(seed, l_c), 200 pairs per period at ±y₀, the witnesses for "Raymond's sources ≥ 0.1018", hit 0.0993–0.1028) tested with the same designs (b1_cue.py, exact Fourier modes of the periodic configuration, D(x) − M_c = (2π/P)Σ_m v(u_m)Re[e^{2πimx/P}Σ_k e^{u_m y_k/2}e^{−2πimx_k/P}], u_m = 4πm/P; the lattice row reproduces A₁/Band exactly):
| design (K, m), centre u* | lattice | l_c = 8 (seeds 1,2,3) | l_c = 12 | l_c = 16 |
|---|---|---|---|---|
| (60, 8) | 1.092 | 1.044, 0.983, 0.997 | 1.088, 0.989, 1.046 | 1.066, 0.969, 1.045 |
| (120, 16) | 2.111 | 1.848, 1.829, 1.830 | 1.732, 1.789, 1.855 | 1.775, 1.718, 1.781 |
| (190, 16) | 2.644 | 2.311, 2.186, 2.226 | 2.138, 2.140, 2.211 | 2.140, 2.011, 2.091 |
(centres u* ± 1 give smaller ratios, 0.05–1.3.) So all nine rigidified rows also violate a valid ζ inequality at t = 0: the pair structure keeps 80–90 % of the lattice's coherent mode at log n = L/2.
Kernel locality (b1_kernel.py): the design's zero-side kernel K(Δx) = ½∫e^{y₀u/2}v(u)cos(uΔx/2)du has |K(20P)|/|K(0)| ≤ 8·10⁻⁶ for all three designs, so a row window of ≥ 64 pairs suffices; ζ's zeros outside the window add < 10⁻⁵ relative.
Mean terms: v has a zero of order K ≥ 60 at u = 0, so every polynomial-in-h part of the archimedean term (incl. 1/s′, 1/(s′−1)) is annihilated; the ψ part gives an oscillatory integral at frequency x/4 ≈ 1.5·10¹² of an analytic function (negligible).
Periodic-row floor under these designs (row with initial height y excluded iff 2πρ cosh(πρy) > Band): (120,16): y > 0.6355, t_B = 0.0792; (190,16): y > 0.6020, t_B = 0.0742. As the order grows, Band → 0 and every coherent row is excluded (the lattice of real zeros needs Band < 2πρ = 13.4, i.e. width < 0.032 in log n).
Witness left: a single off-line pair among real zeros (ExploreRelax l.42, local equilibrium, hits 0.071–0.077). A single pair's design response is ≤ 2cosh(πρy₀)e^{−πρ} ≈ 0.18 of the band for every localised design (point value vs. local prime mass), so it plausibly passes every phase-uniform datum; not tested along t.
Verdict B1: OPEN (float). The class "single-height probes" (TREE, NOGO) is correct as scoped, but the REPORT's "every phase-uniform σ_eff > 1 input"/"any phase-uniform design stays 0.0957" and the 11 % half of TREE's "so Raymond-type chains … have about 5–11% room" do not follow (the 5 % half, rows passing Raymond's present 67 sources, stands): both witnesses fail multi-height (height-derivative) phase-uniform data that the same lemma proves. Model floor for multi-height phase-uniform arguments: between the single-pair witness ≈ 0.071–0.077 (if it passes along the whole trajectory) and the record 0.10755. Not a Λ bound.
Certificate steps for the exclusion: (1) Arb evaluation of Band(v) = ½ΣΛ(n)n^{−1/2}|v(log n)| with an explicit ψ(x) − x tail bound beyond 10⁷; (2) Arb value of A₁ (closed form) or of the CUE rows' exact Fourier modes; (3) the window-splice bound (row inside ζ's zeros, kernel tails); (4) for a Λ result: new chain sources built from such designs and the single-pair witness's trajectory (open).


## B2 (N3): R400 "the Li-range constant (106.25 T²) is limited by one worst-case zero above T near Re s = B, so explicit zero density cannot raise it"

Step fixed (Li200k verify.py l.10–14, 126–145): every zero above T is charged q_n(γ) = 4 sinh²(n a/(2γ²)) + 4n a²/γ³ (a = θ₁₄₁ − ½, all zeros allowed at β = B), summed with the zero density (1/2π)log(γ/2π) on [T, 100T]; positivity needs Σ_{γ>T} q_n < λ_n's RH-model value ≈ (n/2)(log n − 1 − log 2π + γ_E).
Class of the barrier: "one worst-case zero". Freedom the class misses: the binding quantity is the number of zeros allowed at β ≈ B just above T, not one zero (b2_li.py, mpmath):
- method model (all zeros above T at β = B): c = n/T² = 106.267 (record 106.2557: reproduced to 0.01 %); the sum is carried by ≈ 1.6·10¹¹ zeros in [T, 1.0126T];
- one worst quartet at height T: c = 175.0 (the true one-zero ceiling);
- a KLN-shaped explicit density N(σ, T′) ≤ A·T′^{8(1−σ)/3}(log T′)^{5−2σ} near σ = B (worst placement: all allowed zeros at height T): A = 1 → c = 117.2, A = 10 → 110.4, A = 100 → no gain.
Verdict B2: OPEN (float, constant only): the stated step is wrong; zero-density input (Kadiri–Lumley–Ng type, already proved in the literature) can raise the constant unless its constant at σ ≈ 0.873, T ≈ 3·10¹² exceeds ≈ 70. Not checked: KLN's actual table constants. A certificate needs an explicit N(σ, T) at σ ∈ [0.8, B], T ∈ [T₀, 1.05T₀] and the modified sum in verify.py. RECORDS-level (range constant), not a node move.


B1 review and repairs. Reviewer (verdict GAP, INCOMPLETE; formulas pasted inline): (a) ζ-side bound correct (absolute convergence, triangle inequality, uniform in x); (b) Poisson expansion correct, but ‖Σ A_m cos mφ‖_∞ ≥ A₁ is false in general (only ≥ (π/4)A₁); (c) the test is sufficient, not "iff"; (d) Band uses a density tail, not a certified bound; (e) window splice not bounded in sum; (f) t = 0 suffices to reject the witness, but rejecting it does not settle the class floor. Repairs here:
- (b) the harmonics are explicit: A_m = 2πρ cosh(mπρy₀)v(mu*); b1_harm.py: |A₂ + A₃ + A₄|/A₁ = 7.2·10⁻⁶, 8.9·10⁻⁸, 3.3·10⁻¹⁵ for the three designs, so sup_x|D_row − M| = A₁(1 ± 10⁻⁵); b1_cue.py computes the sup over x from all modes in the design's support.
- mean term: v(u) = O(u^K) at 0, so D[m] = Σ_{j≥K} p_j(−2)^j ∂_h^j m(x, h₀); |∂_h^j m| ≤ j!(2/x)^{j+1}·O(1) (ψ^{(j)}(s′/2) and 1/s′ terms at |s′| ≈ x/2), i.e. < 10⁻⁶⁰⁰ for j ≥ 60, x = 6·10¹².
- (d), (e): open certificate steps (listed above); the float margins (ratio 2.11, 2.64) are large against both.
Lowest order found (b1_base.py, run s2): orders 20–44 at one height h₀ = 3.97 (σ_eff = 2.49, inside Raymond's height range) already give ratio 1.028 for the lattice row; 10–34 at h₀ = 2.49: 0.69. The K = 20/m ≥ 32 and K = 10/m ≥ 32 LPs lose control beyond their point range (Band re-evaluation explodes) and are not used; K = 20, m = 40 hit the 600 s cap.
Context: Raymond's 23752103 sources are already signed multi-height designs (DBNAudit l.66: row 0 has 9 nodes at heights 2.5–6, weights −14.3…16.9), and explore51's raycheck found the rigidified rows pass all 67. So the robust part of the barrier is: rows pass single-height probes and Raymond-size designs (≈ 9 nodes, heights 2.5–6); they fail designs of ≈ 25+ height-derivatives (ratio 1.03 at orders 20–44, 1.09 at 60–68, 2.1–2.6 at 120–206).

## B3 (N9): Explicit611's "floor 0.5075 belongs to the method" (t^{1/6} log t pipeline)
Step fixed: at the binding height t ≈ 2.3·10⁶ the trivial part 2Σ_{n<4K}n^{−1/2} with K ≥ t^{1/3}, r ≥ 4 (Revers/HPY lemma hypotheses). Class: Riemann–Siegel blocks with van der Corput bounds and an absolute-value head. Freedoms: (i) a different inequality in the same step: Fodge's explicit B-process (Lemma C.3) already passes the floor for the large-t constant: Prop. C.1 = 0.3080532 log t + 2.405897 + 506.6 t^{−1/6} (per t^{1/6}) is below 0.5075 log t for t > 8.66·10¹² (float); (ii) exact (Arb) evaluation instead of the head bound at the binding heights. Verdict: OPEN, already passed by Fodge (17 Nov) for t > 8.7·10¹²; constant-level only, no test beyond the crossing.

## B4 (N2): R449 "de Bruijn–Newman arguments blind to the Euler product cannot prove Λ ≤ 0.2306"
Step fixed: none (a counterexample barrier). Class: arguments every input of which also holds for Davenport–Heilbronn's f (functional equation and gamma factor, order and growth, Riemann–von Mangoldt count); Λ_DH > 0.2306 certified (Missiaen, rerun R449). Freedoms the method's lemmas admit outside the class: RH verified to T₀ (DH has 0.8085 + 85.70i), the θ₁₄₁ strip (DH has zeros up to σ* = 1.12), phase-uniform log-derivative data at σ_eff > 1 (f′/f has no convergent Λ-series below σ*); every record uses them. Verdict: ROBUST for its class (a counterexample covers every freedom inside it by construction).

## Not examined
- N5 nb2nd ("N ≤ 5120 cannot fix the sign of C₂"): NB2nd saved only window means of κ(N) (t8_kappa.log), not per-N values; the obvious freedom (regressing out the first zeros' almost-periodic terms cos(γ_k log N) before averaging) needs a rerun of the Gram computations, beyond this budget.
- N13 (CCM "certifying L > 2 is open", weil-simple limits): engineering limits (memory, precision), no stated method step to falsify in the time left.


## Novelty (OPEN items)
- B1: python3 .scratch/lst.py 'Bruijn|Newman|heat flow|phase.uniform' math.NT math.CV math.CA --zenodo '"de Bruijn-Newman" barrier' (2026-10, 2026-11): titles 2611.06245, 2610.12803, 2610.18214, 2610.27343 (all seen by explore51), 2610.05890 (Brascamp–Lieb, unrelated); Zenodo: nothing on dBN barriers (24377750/53/57 "Computing the Riemann Hypothesis" I–III, 19 Nov, Thorin programme, titles only). Our files: ExploreRelax states the overclaim; DBNAudit l.66 notes Raymond's signed multi-height rows; no file tests high-order designs against the rows. An audit correction, not a discovery (the mechanism is band-limiting in log n).
- B2: lst.py 'Li coefficient|Keiper|Li.s criterion|zero.density' math.NT --zenodo '"Li coefficients"': 2610.15989 "Zero density in short strips and short-interval bounds toward Hypothesis S" (unread; a short-interval density near T is exactly the input B2 needs), 2610.37371 (unrelated); Zenodo: no range-with-density item. Our files: ExploreLatticeGas/REPORT.md l.122 already found "one worst-case zero above T gives a Li reach of 174 T²" (here 175.0), which contradicts R400; nobody drew the consequence.

## Suggested NOGO lines
- ceilaudit3 (N2, explore51 refined): the row ceilings 0.0957 (lattice) and 0.1018 (rigidified CUE) hold for single-height probes and Raymond-size signed designs (≈ 9 nodes, heights 2.5–6) only; ≥ 25 height-derivatives of one probe (h₀ ≈ 4, the same triangle-inequality lemma) exclude both witnesses at t = 0 (float: ratio 1.03–2.64 lattice, 1.7–2.3 CUE); the floor for general phase-uniform σ_eff > 1 data is open (single-pair witness ≈ 0.071–0.077) [CeilAudit3](Riemann/Deep/CeilAudit3/REPORT.md)
- ceilaudit3 (N2, R449 checked): DH barrier Λ > 0.2306 binds exactly the arguments whose every input holds for Davenport–Heilbronn (functional equation, gamma factor, growth, zero count); RH to T₀, the θ₁₄₁ strip and σ_eff > 1 log-derivative data lie outside its class [CeilAudit3](Riemann/Deep/CeilAudit3/REPORT.md)
- Correction for R400 (N3): replace "limited by one worst-case zero …, so explicit zero density cannot raise it" by "limited by ≈ 1.6·10¹¹ zeros in [T, 1.013T] all charged at β = B; one zero allows 175 T²; an explicit density N(σ, T′) ≤ A T′^{8(1−σ)/3}(log T′)^{5−2σ} near σ = B gives 117 T² (A = 1), 110 T² (A = 10) (float, ceilaudit3)".

RH relevance: N2 (the model floor for phase-uniform Λ chains is no longer 0.0957; no Λ bound), N3 and N9 (constants only).

## Logs (memory, time)
Peaks (max RSS by /usr/bin/time -l, memcap peak): b1_sieve 74 MB; b1_lp 184 MB; b1_dual 167 MB; b1_base s1 158 MB; b1_cue 163 MB (73 s); b1_kernel, b1_harm < 160 MB; b2_li 19 MB (233 s); b1_base s2 killed at its 600 s cap (175 MB). Two cap breaches, both mine: the first b1_sieve version (int64 overflow in P**k) reached 672 MB max RSS for 2.8 s before memcap sampled it; one b1_base run (K = 0, 51,175 LP columns) was killed by memcap at 326 MB. Total compute ≈ 25 min, one process at a time, OMP_NUM_THREADS = 1.

## Files
Riemann/Deep/CeilAudit3/: REPORT.md; scripts b1_sieve.py (Λ(n), n ≤ 10⁷ → /tmp/ceil3_lambda.npz, run first), b1_lp.py, b1_dual.py (early LP forms), b1_base.py (jets K, m; args 'K list' 'm list' tag), b1_cue.py (lattice and rigidified CUE rows; arg: centre offsets), b1_kernel.py, b1_harm.py, b2_li.py; logs b1_base_s1.json, b1_cue.json, b1_dual_t1.json, b1_lp_2_0,2,4.json. Scratch notes: .scratch/ceilaudit3/STATE.md.

