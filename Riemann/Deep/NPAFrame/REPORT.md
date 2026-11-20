NO IDEA: the inside level-2 gain survives the continuum falsification but nothing is certified — at fixed ε no zero-cost fake exists (rays act on nothing; convergent cluster families always raise a constrained diagonal quartic), off-line columns stay unused (lattices N = 6, 7), and the fixed-θ gain shows no decay up to N = 12: for words with Σ|ξ| ≤ 1 the open-band discretisation (AF's finite programme: Po(N) ≈ 0.68183 + 0.27/N², while explore52/npafals' closed analogue is inflated by ≈ 0.6/N by the single boundary word E|p_N|² = N) and the closed one agree at N = 10 and 12 (6.0 and 6.7·10⁻⁴; suggestive of a positive limit, not proof); no certificate for stationary laws, so goal (2) was not started (node N8)

# NPAFrame — does the level-2 (NPA) gain survive in the continuum? (node N8, npaframe, Orca)

Statement. Falsification-first test of explore52/npafals' level-2 lead in the continuum (stationary laws, bandwidth-one form factor on
the open band, Rudnick–Sarnak-range triple/quartic cumulants on Σ|ξ| ≤ 2 − ε, off-line pairs with weight e^{2π|ξ|y} = T^{|ξ|δ}).
(i) Under the model hypothesis H the continuum question at fixed ε is the N → ∞ limit of the N-periodic problems at fixed θ = 1 − ε/2 (F1: level-2
cumulant data transfer exactly, pair data on the grid; reduction R). (ii) No zero-cost fake exists at fixed ε:
rays act on nothing, convergent cluster families always raise a constrained diagonal quartic (F2; escaping sequences open); priced
off-line columns stay unused (lattices N = 6, 7). (iii) The closed CUE(N) analogue used so far is biased by ONE boundary word,
E|p_N|² = N at ξ = 1, which lies outside the open band (the continuum problem leaves Γ({1}) free): it inflates the bandwidth-one
value by ≈ 0.6/N (0.071 at N = 8); the open band (AF's finite programme) gives Po(N) ≈ 0.68183 + 0.27/N² (fit on N = 8, 12; predicted N = 14 to 10⁻⁶; an extrapolation: AF's
ceiling looks like the N → ∞ value of this model), and its optimal laws carry a Bragg spike of mass ≈ 0.87 at ξ = 1. (iv) Fixed-θ
gains show no decay: at θ = 1/2 (Σ|ξ| ≤ 1) the closed gains rise, the open ones fall from large small-N values (non-monotone), and
the two agree at N = 10 (0.000604 / 0.000596) and N = 12 (0.000667 / 0.000666) within the heuristic precision — suggestive of a
common positive limit, not distinguishing it from a late decay; at θ = 2/3 they close in (N = 9: 0.002142 / 0.001182); at θ = 3/4 they are still apart at N = 8 (open 0.003852, falling from 0.0335; closed 0.001564, rising from 0.000945). (v) Nothing is certified for
stationary laws, so goal (2) was not started; its blocking step is a configuration-wise multi-bandwidth matrix inequality (§7).

## 1. Continuum model (stationary laws, no N)
- Configuration: stationary marked point process X = Σ_j m_j δ_{x_j} on ℝ (on-line zeros, marks m_j ≥ 1) plus off-line pairs at x ± iy
  (y > 0, mean-spacing units; y = (β − 1/2) log T/2π); density 1 counted with multiplicity; objective n₁ = density of simple on-line zeros.
- Linear statistic X_f = Σ_j m_j f(x_j) + Σ_pairs [f(x + iy) + f(x − iy)], f̂ ∈ C_c^∞; a factor of frequency ξ gets weight ≍ e^{2π|ξ|y} = T^{|ξ|δ}.
- Level 1 (P): Bartlett spectrum Γ = |ξ|dξ on the OPEN band (−1, 1) (Montgomery/BGSTB; r̂ ∈ C_c(−1,1)).
- Level 2 (L2(ε)): third and fourth cumulant spectra B, T vanish on {Σξ_j = 0, Σ|ξ_j| ≤ 2 − ε} (NPAFals lemma (C): mixed words = Gaussian
  values + O(1/log T), no RH, any compact subset of the open simplex). For three frequencies Σ|ξ_j| = 2 max|ξ_j|: triples live on max|ξ_j| ≤ 1 − ε/2.
- Finite-N dictionary (the discretisation used by explore52/npafals): N-periodic laws, frequencies ξ = k/N, Γ({k/N}) = E|p_k|²/N²,
  E p_a p_b p̄_{a+b} = N·B(a/N, b/N) + o(N), E p_λ p̄_μ − (Gaussian part, O(N²)) = N·T + o(N); 'level' n = |λ| ↔ Σ|ξ| = 2n/N.
  Fixed ε ↔ fixed θ = 1 − ε/2: inside words at levels ≤ θN. The continuum gain at ε is G_∞(θ) := lim_N G(N, ⌊θN⌋) (if it exists),
  G(N, n₀) = min[P + words of length 3–4 at levels ≤ n₀] − min[P] (closed: pairs k ≤ N; open G_o: pairs k ≤ N − 1, sets Po, To<n₀>);
  the open-simplex gain is lim_{θ→1} G_∞(θ) = sup_θ G_∞(θ) (monotone).

## 2. Structural facts (F1 exact; F2 proved for rays and convergent clusters; F3 a reading, checked at level 1 in §4b)
- (F1) Level-2 constraints transfer EXACTLY between N-periodic laws and stationary laws. For an N-periodic law, rotation-averaged (WLOG: all
  words and n₁ are rotation-invariant; reviewer: state it, cumulants are of the rotated q_k), with uniform random shift
  and fixed-width statistics X_f = Σ m_j f(x_j − τ), Poisson summation gives X_f = N^{-1}Σ_k f̂(k/N) p_k e^{−2πikτ/N}, so
  E X_{f1}⋯X_{fm} = N^{-m} Σ_{k1+⋯+km=0} Π f̂_j(k_j/N) E p_{k1}⋯p_{km}: the cumulant spectral measures of the periodic law are atomic on
  (1/N)ℤ^m with atoms (cumulant part of the word)/N. Hence the periodic law satisfies "B = T = 0 on a region R" iff every finite word
  sampled in R equals its Diaconis–Shahshahani value (the DS value = the Gaussian/pairing value, all rotation-invariant cumulants 0).
  The PAIR constraint does not transfer exactly (atomic Γ = Σ_k E|p_k|²N^{-2}δ_{k/N} vs |ξ|dξ): AF's ceiling needs the r', r'' error
  terms for exactly this reason. So the finite model = continuum level-2 data + sampled bandwidth-one data.
- (F2) No zero-cost fake at fixed ε from rays or convergent cluster families (first review GAP: clusters; repaired, partly). A pair at depth y weighs ≍ e^{2π|ξ|y} per
  factor (periodic model: p_k gets e^{2πikx/N}·2cosh(2πky/N)), so ≍ e^{4πσy} in a word of level σ = Σ|ξ|/2, ≍ e^{4πξy} in E|p_ξ|².
  (a) Rays (isolated pairs, y → ∞, density ν → 0, νe^{4πy} bounded, inserted independently): every constraint of P_open + L2(ε)
  (pair words on [−r, r], r < 1: O(e^{−4π(1−r)y}); words with σ ≤ θ := 1 − ε/2: O(e^{−2πεy})) tends to 0, cross terms with on-line
  points carry smaller exponentials: rays, explore52's free fakes of top-level words, act on nothing at fixed ε (second review: correct).
  (b) Growing clusters (B pairs per cluster, cluster rate ν/B; contributions νB^{m−1}e^{4πσy} to order-m words, reviewer's example):
  for B → ∞ the order-4 terms lead, and the leading effect on the level-σ quartics is Σ_c w_c p_c(ξ₁)p_c(ξ₂)p̄_c(ξ₃)p̄_c(ξ₄), where
  p_c(ξ) = Σ_i m_i e^{2πiξx_i}·2cosh(2πξy_i) (normalised) is an entire function of ξ. Its diagonal part Σ_c w_c |p_c(ξ)|²|p_c(σ − ξ)|² is a constrained
  inside word for EVERY ξ ∈ (0, σ) (σ ≤ θ); it vanishes only if p_c(ξ)p_c(σ − ξ) ≡ 0 on an interval, i.e. p_c ≡ 0 (entire functions
  form an integral domain), and then every word of the cluster vanishes. So a cluster family whose normalised profiles converge
  locally uniformly always raises a constrained diagonal quartic: a priced column, never free; this removes NPAFals (B)'s corner (3,1)
  candidates for that class (on the integer grid p_c may vanish at the needed integers, on an interval it cannot). NOT covered (second
  review, GAP, routine): escaping sequences, e.g. mixed regimes B = e^{ay} with νB³e^{4πθy} = 1, whose profiles e^{2πy(ξ − θ/2)} have no
  locally uniform limit; the coincident-stack example is still diagonally priced (its limit effect is 1 on every level-θ diagonal
  quartic), but a proof for all escape limits is missing. (2,2) words are closed in every limit (Cauchy–Schwarz, partners at the same
  level). The remaining LP question, whether priced off-line columns pay, is tested on lattices in §4 (immunity).
- (F3) The finite top level is a boundary layer. A spectral feature of width w needs coherence length ≍ 1/w; an N-periodic law has
  coherence ≤ N, so its top level |λ| = N (Σ|ξ| = 2) is the continuum boundary layer of width 2/N, where a stationary law may carry
  structure invisible to every constraint on the open simplex (a ray: Γ_ray(ξ) ≍ νe^{4πξy} → 0 on [0, 1) while νe^{4πy} stays
  finite, so its mass escapes to ξ ≥ 1; in the N-periodic model the same ray sits on the top level k = N). Hence finite-N top-level
  gains (L4 − L4o ≈ 0.02 at N ≤ 9) need not have a continuum counterpart, and L4o − P = G(N, N − 1) is not a fixed-ε quantity
  (its cutoff ε = 2/N → 0 sits inside the layer); level-1 instance in §4b: the single boundary word E|p_N|² = N moves the finite
  optimum by ≈ 0.6/N (0.071 at N = 8). Under H the continuum question is whether G_∞(θ) = lim_N G(N, ⌊θN⌋) > 0 for
  some fixed θ < 1.
- Monotonicity: G_∞ is non-decreasing in θ (more constraints), so "survives as ε → 0" is automatic once G_∞(θ) > 0 for one θ < 1;
  the only way the gain dies is G(N, ⌊θN⌋) → 0 for EVERY θ < 1 (npafals' reduction, here with the boundary layer of F3 made explicit).


## 3. Method
- Continuum question reduced (F1–F3) to the fixed-θ gain G(N, ⌊θN⌋) of the N-periodic model, N → ∞. Measured with cg_theta.py /
  cg_law.py / cg_sets.py (copies of NPAFals/colgen_simple.py: continuous positions, RH world, patterns = doubles + at most one triple,
  CUE(N) initial columns, Adam pricing, tolerance 2·10⁻⁵; new sets T<n₀> = all pair words + all words of length 3–4 at levels ≤ n₀,
  Dg<n₀> = diagonal quartics |p_a p_b|² only, Ph<n₀> = phase words (λ ≠ μ) only; cg_law.py also saves the final law).
  Values are master values of feasible laws (upper bounds for each finite min); lb = master − max reduced cost found (heuristic).
- Off-line: npa_cg_lean.py (NPAFals copy, rowspec T<n₀> added in the local npa_cg.py): lattice Z_M, M = 3N, depths {0.05, 0.2, 0.6,
  1.5, 4}, explore52's rays, exhaustive pricing; on vs off + rays on the same lattice = the immunity test at fixed θ.
- Open band: sets Po (pairs k ≤ N − 1, AF's finite programme) and To<n₀> (Po + words of length 3–4 at levels ≤ n₀); cg_lite.py =
  cg_law.py with 1000 instead of 3000 CUE initial columns (memory). Every run single-threaded under .scratch/memcap.sh 290
  (+ /usr/bin/time -l); queues run_q1–q17.sh, one process at a time; memcap peaks ≤ 278 MB for every finished run.


## 4. Evidence (RH world unless stated; logs/ in this folder; * = npafals' logs, same code family)
Fixed-θ gains G(N, ⌊θN⌋) = T_{⌊θN⌋} − P (master values; lb within 2·10⁻⁵ unless stated; P(N) = 0.847708*, 0.815624*, 0.793924*,
0.778268*, 0.766430*, 0.757058*, 0.749498*, 0.743256, –, 0.733730 for N = 3–12; N10_P_s2 (seed 2, 2× batch, 2× steps) 0.743250):
| θ | N = 3 | 4 | 6 | 8 | 9 | 10 | 12 |
|---|---|---|---|---|---|---|---|
| 1/4 | – | – | – | 0.000020 (≈ noise) | – | – | 0.000108 |
| 1/2 | – | 0.000073 | 0.000197 | 0.000378 | – | 0.000596 | 0.000666 |
| 2/3 | 0.000217* | – | 0.000686* | – | 0.001182 | – | – |
| 3/4 | – | 0.000945* | – | 0.001564* | – | – | (killed, memory) |
- Every fixed-θ column rises with N (θ = 1/2 increments +1.2, +1.8, +2.2, +0.7 ·10⁻⁴: rising, possibly starting to flatten at N = 12);
  a boundary layer of O(1) levels (gain depending on N − n₀ only) would make these columns fall: at N = 12, θ = 1/4 sits 9 levels
  below the top and still gains 1.1·10⁻⁴ (5× its N = 8 value).

- Which words carry the closed-model fixed-θ gain (θ = 1/2): N = 8: Dg4 (diagonal quartics |p_a p_b|² only) 0.000064, Ph4 (phase words
  λ ≠ μ only) 0.000324, T4 (all) 0.000378; N = 10: Dg5 0.000068, Ph5 0.000534, T5 0.000596 (logs/N8_Dg4, N8_Ph4, N10_Dg5, N10_Ph5).
  The moduli statistics (type mixing, Cov(|p_a|², |p_b|²)) give a small flat part; the growth is in the triples / off-diagonal words.
- Pricing robustness: N = 10, seed 2, 2× batch, 2× Adam steps: P 0.743250 (seed 1: 0.743256), T5 0.743846 (0.743852): same gain 0.000596.
- Off-line immunity at fixed θ (lattices Z_{3N}, depths {0.05, 0.2, 0.6, 1.5, 4}, explore52's rays incl. 3-pair clusters, exhaustive
  pricing): N = 6, T4 (θ = 2/3): on 0.838057 = off + rays 0.838057 (P 0.837345*); N = 7, T5 (θ = 5/7): on 0.825180 = off + rays 0.825180
  (P 0.824515, reproduces npafals); off-line mass 0, no ray column priced in: the lattice gains (0.000712, 0.000665) are untouched.

### 4b. Open band — the interior-frequency discretisation (reviewer's k = N point; AF's finite programme drops the pair word k = N)
- The closed model's top pair word E|p_N|² = N sits at ξ = 1, outside the open band (Γ({1}) is free in the continuum problem), yet it is load-bearing:
  Po(N) (pairs k ≤ N − 1) = 0.711349, 0.698234, 0.688938, 0.685999, 0.685210, 0.684572, 0.683681, 0.683188 (N = 3, 4, 6, 8, 9, 10, 12, 14) against
  closed P(N) = 0.847708*, 0.815624*, 0.778268*, 0.757058*, 0.749498*, 0.743256, 0.733730, 0.726777 (N·(P − Po) = 0.41, 0.47, 0.54,
  0.57, 0.58, 0.59, 0.60, 0.61). Po(N) ≈ 0.681827 + 0.267/N² (fit on
  N = 8, 12; N = 4, 6, 10 within 3.1·10⁻⁴; it predicted N = 14 as 0.683189, run: 0.683188; AF's N = 256 optimum 0.6818287 within
  2·10⁻⁶; a 3-parameter fit on N = 8, 10, 12 gave 0.6812, refuted by N = 14: extrapolations, not bounds); the closed excess ≈ 0.6/N is a pure boundary effect (F3 seen at level 1). explore52's and
  npafals' closed-model values carry it, and so does the closed table above.
- What the boundary word does (logs/formfactor.log, E|p_k|²/k of the optimal laws): every bandwidth-one optimum has E|p_k|² = k up to
  the last constrained k, then a Bragg-type spike at the first free frequency, E|p_k|²/N² ≈ 0.86–0.88 (open: k = N, i.e. ξ = 1; N = 8:
  7.04N, N = 12: 10.37N; closed: k = N + 1; N = 12: 9.32N), then a depleted gap (E|p_k|²/k ≈ 0.03–0.12 for the next 3–5 k). The
  continuum adversary is thus lattice-like: a Bragg peak of mass ≈ 0.87 at ξ = 1, invisible to the open band; the closed model pushes
  it to 1 + 1/N at a cost ≈ 0.6/N. The level-2-optimal laws keep the spike (N = 12: open To6 10.32N, closed T6 9.25N at k = 13).
- Open fixed-θ gains G_o(N, ⌊θN⌋) = To_{⌊θN⌋} − Po:
| θ | N = 3 | 4 | 6 | 8 | 9 | 10 | 12 |
|---|---|---|---|---|---|---|---|
| 1/4 | – | – | – | 0.000029 (≈ noise) | – | – | 0.000046 |
| 1/2 | – | 0.000869 | 0.000664 | 0.000512 | – | 0.000604 | 0.000667 |
| 2/3 | 0.039604 | – | 0.002900 | – | 0.002142 | – | – |
| 3/4 | – | 0.033534 | – | 0.003852 | – | – | (killed, memory/time) |
- The two discretisations come together: at small N the open gains are large (the inside words sit next to the free top frequency)
  and fall (non-monotonically), the closed gains are small (the top word already binds) and rise; at θ = 1/2 they meet at N = 10
  (open 0.000604, closed 0.000596) and again at N = 12 (open 0.000667, closed 0.000666): two different discretisations agree to
  10⁻⁵ (inside the 2·10⁻⁵ heuristic precision) on a gain of ≈ 6.7·10⁻⁴, still rising slowly. At θ = 2/3 the gap closes (N = 9: open 0.002142, closed 0.001182). At θ = 1/4 both are tiny and rise (open 0.000029, 0.000046;
  closed 0.000020, 0.000108 at N = 8, 12; logs/N8_To2, N12_To3): no decay, not yet met. At θ = 3/4 they are still apart at N = 8
  (closed 0.001564*, open 0.003852); the N = 12 pair (94–95 features) hit 317/294 MB at iteration 0–1 under memcap and was
  killed, and a 1000-CUE-column retry needed ≈ 20 s per iteration (stopped at 0.81 after 14 iterations); a θ = 0.7, N = 10 pair (To7/T7) was stopped likewise (≈ 15 s per iteration under the
  machine load, 0.79 after 15 iterations; logs/*_killed.log); so was the θ = 1/2, N = 14 pair (To7 at 0.805 after 20 iterations,
  18 s per iteration; Po(14) and P(14) finished).

## 5. Checks and review
- Reproduction (brief: the finite runs must reproduce npafals' L4o − P at N = 3–5): N = 4 T3 (= L4o) 0.816571 (lb 0.816552) vs npafals
  0.816569 (logs/N4_T3.log); N = 3: P 0.847708, T2 (= L4o) 0.847925, gain 0.000217 (npafals: 0.847708, 0.847925, 0.000217); N = 5: P 0.793926
  (npafals 0.793924), T4 (= L4o) 0.795886, gain 0.001960 (npafals 0.795885, 0.001961). P(12) reproduced exactly by two codes (cg_theta 0.733730, cg_law 0.733730);
  P(10) seed 1 / seed 2 with 2× batch and 2× Adam steps: 0.743256 / 0.743250; T5(10) seed 2: 0.743846 (seed 1: 0.743852).
- Review 1 (fresh reviewer, formulas F1, F2, R pasted inline): F1 GAP (state that cumulants are of the rotated variables q_k) →
  repaired in the text, the equivalence itself confirmed ('standard equivalence after correction'); F2 GAP (density control misses
  clusters with growing multiplicity) → F2(b) added; R: monotonicity CORRECT, the periodic-to-stationary identification is a MODEL
  HYPOTHESIS (H) (existence of the limit, recovery, control of n₁), the cutoff must be read on the closed region Σ|ξ| ≤ 2 − ε (n ≤ θN,
  as used here), the extra pair word k = N needs a check (open-band runs Po/To, §4). Overall ROUTINE.
- Review 2 (fresh reviewer, repaired F2 pasted inline): (a) CORRECT with the rates now stated; (b) correct for convergent profiles,
  GAP for escaping sequences (text restricted accordingly). Overall ROUTINE.
- Review 3 (fresh reviewer, §4b numbers pasted inline): GAP on all three readings as first worded ('faithful discretisation' needs
  the convergence of H; the endpoint is outside the band, nullity alone is no reason; the 1/N² fit is an extrapolation; 'from above
  and below' overstated: open gains are non-monotone; agreement at N = 10, 12 is within the heuristic precision and cannot exclude a
  late decay) → wording corrected as now in the Statement and §4b. Overall ROUTINE: 'useful numerical diagnostics and conjectures'.

## 6. Novelty
- No discovery is claimed (nothing certified), so the arXiv/Zenodo scan was not run. Our files: explore52 and npafals use the closed
  analogue only (their P(N) carry the ≈ 0.6/N boundary excess of §4b); AF (2608.13637) set up the open programme at N = 256 but
  consider pair data only; npafals (B) proves only corner (3,1) words fakeable on the integer grid, F2(b) removes them for convergent
  cluster families in the continuum.

## 7. Goal (2) gate and the blocking step
- (2) is gated on a certified gain in (1); none exists (§4 values are floats of the periodic model, F1: pair data only on a grid),
  so no degree-4 certificate was attempted. Blocking step, stated precisely: zeta23's frame turns pair words into a configuration-wise
  statement through the exact Poisson–Gabor identity Σ_k φ̂(z − α_k)φ̂(z′ − α_k) = Lφ²^(z − z′) (AF Lemma 'Poisson–Gabor'), which makes
  ‖G̃‖²_HS a pair sum with band-limited kernel; the level-2 analogue needs traces of products of Gram matrices of DIFFERENT
  bandwidths α_j with Σα_j ≤ 2 − ε (a single bandwidth needs kα < 2: AF §7(e)), whose kernels are products K_{α1}(γ−γ′)K_{α2}(γ′−γ″)⋯
  and are evaluated by lemma (C) only after the same exact periodisation; no such multi-bandwidth Gabor identity is known to us, and
  the finite-model certificates are polynomials in p_1..p_N, not trace functions of G̃. Shape of the missing inequality: a lower bound
  for rank P (P = on-line Gram matrix) by mixed traces tr(P K₁ P K₂ ⋯) with band-limiting K_j; with ONE K the spectral moments of
  K^{1/2}PK^{1/2} give Chebyshev–Markov–Stieltjes bounds, but three kernels need bandwidth < 1/2 and four < 1/4, and the rank of the
  compression is then capped at that bandwidth (AF's dimension cap); mixed K_j need an NPA moment matrix of non-commuting words in
  P and the K_j, which is not available. Lemma (C) supplies every needed expectation
  (all mixed moments of order ≤ 4 on the open simplex); what is missing is the configuration-wise matrix inequality itself.
- Size check against the node: in the open band the periodic values are already near the limit (Po(12) = 0.683681; To6(12) = 0.684348,
  To6(8) = 0.689851 at θ = 3/4). Beating AF's ceiling needs P_∞ + G_∞(θ) > 0.6818287: with P_∞ ≈ 0.68183 (the fit that predicted N = 14 to 10⁻⁶)
  the θ = 1/2 gain alone would give ≈ 0.6825 in the model; and a near-optimal certificate would still be needed —
  the best certified bandwidth-one value (0.6736393) is itself 0.0082 below the ceiling.


## 8. Caveat
- Everything numerical is floats of N-periodic models (master values of feasible laws; lb heuristic): evidence, not certificates.
- (H) Identifying the stationary problem with N → ∞ of periodic problems is a model hypothesis (reviewer 1); the two periodic models
  (closed, open) are only consistent with a common positive limit at θ = 1/2, not a proof of it; N ≤ 12.
- F2(b) covers rays and convergent cluster families, not every escaping sequence (reviewer 2); priced off-line columns were tested on
  lattices N = 6, 7 only. F3 is a heuristic reading (coherence length), F1 is exact.
- The RH-world adversary uses patterns 'doubles + at most one triple' (npafals' choice); higher multiplicities untested here.
- Nothing here is about ζ itself beyond lemma (C) (NPAFals) supplying the inside words; no bound for ζ.

RH relevance: N8 (unconditional proportion of simple zeros on the line: whether Rudnick–Sarnak-range words strictly inside the simplex,
the unconditional part of level 2, carry a gain over the bandwidth-one optimum in the continuum, i.e. past AF's ceiling 0.6818287).

## Suggested NOGO line
- npaframe (Orca, continuum falsification of the level-2 lead → N8) NO IDEA: at fixed ε (Rudnick–Sarnak words with Σ|ξ| ≤ 2 − ε, the
  unconditional part) no zero-cost fake exists (rays act on nothing; convergent cluster families always raise a constrained diagonal
  quartic, an entire-function argument that also removes npafals' corner candidates) and off-line columns stay unused (lattices);
  the closed CUE(N) analogue of explore52/npafals is inflated by ≈ 0.6/N by the one boundary word E|p_N|² = N (open band = AF's
  programme: 0.68183 + 0.27/N², a Bragg spike of mass ≈ 0.87 at ξ = 1); fixed-θ gains do not decay (Σ|ξ| ≤ 1: open = closed =
  6.7·10⁻⁴ at N = 12); nothing certified for stationary laws [NPAFrame](Riemann/Deep/NPAFrame/REPORT.md)

## Files
- REPORT.md (this); cg_theta.py, cg_law.py, cg_sets.py (copies of NPAFals/colgen_simple.py; headers list the changes); npa_cg.py,
  npa_cg_lean.py, offline_lp.py (copies from NPAFals, rowspec T<n₀> added in npa_cg.py); colgen_open_copy.py (ExploreNPA copy, unused);
  cg_lite.py (cg_law.py with fewer CUE initial columns, memory); analyze_law.py, formfactor.py (law structure); run_q1–q17.sh
  (queues, one process at a time; superseded tails noted in their headers; killed runs kept as *_killed.log); logs/ (every run:
  .log, .err with memcap peak and /usr/bin/time -l, .npz with duals and the final law; logs/summary.log one line per run).
