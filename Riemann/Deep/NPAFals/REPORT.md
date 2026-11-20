VERDICT: SURVIVES (numerical evidence to N = 9, not a proof): in continuous positions the immune inside gain L4o − P rises at every N = 3–9 (0.0002 → 0.0064, ≈ 0.001 per unit N), its fixed-θ lower bounds rise too, the best fit tends to 0.031 (≈ the full level-2 gain), a fit tending to 0 must peak beyond N ≈ 3000, and no pair cluster (J ≤ 4) fakes a word for free; but explore52's lattice trend (0.0008 → 0.0108) and inside shares (3 → 24 %) were coarse-lattice artifacts, and at N ≤ 9 the inside gain is only 1–25 % of the level-2 gain.
(C): LEMMA PROVED on the open simplex Σα_j < 2 (mixed Rudnick–Sarnak-range words of ζ's zeros, complex, off-line pairs included, equal their Gaussian = CUE values + O(1/log T) without RH; Hughes–Rudnick's argument with unequal supports; reviewer GAP at the prime-cluster count, repaired); NOT at the boundary Σα_j = 2, explore52's top level, where the off-diagonal step fails.

# NPAFals — falsification of explore52's level-2 (NPA) lead, node N8 (npafals, Orca, 19–20 Nov 2026, 22:28–00:10 UTC)

Summary. (A) explore52's numbers came from lattices Z_M whose resolution fell with N (M/N = 12 → 3); in continuous positions
(column generation, objective min E[#simple]/N) the inside gain is 2–3× smaller than the lattice said but grows with N, and so does
the gain on any fixed fraction θ < 1 of the RS range (θ = 2/3: ×3.2 from N = 3 to 6; θ = 3/4: ×1.65 from N = 4 to 8); inside
words stay immune at N = 7; the off-line retention of full L4 at fixed resolution is 100, 100, 98.9, 98.0 % (N = 4–7): falling
slowly, never below the inside share. (B) No free fake from clusters: impossible for every diagonal-closed set (L4c) at every level
(rank-one PSD); in L4 only 'corner' (3,1) words can be candidates (proved), and none is faked by J ≤ 4 pairs at N ≤ 9 (search with
a positive control). (C) The load-bearing input is a theorem for the inside words (no RH); the top level stays open.
Never a bound for ζ: everything in (A)–(B) is the finite CUE(N) analogue.

## (A) The trend of the inside gain and of the off-line retention
Method. colgen_simple.py (copy of explore52's colgen_open.py; change: objective min E[#simple]/N instead of max E[Σm²]/N, CUE(N)
samples as initial columns, warm-started pricing): continuous positions, RH world, multiplicity patterns doubles + at most one triple;
the master value is a feasible law (upper bound for the min), lb = master − max reduced cost found (heuristic lower bound);
the vectorised gradient reproduces the old code's N = 4 L4o value 0.816569 exactly, 3.3× faster (logs/A_N4_L4o.log, A_summary.log). Inside
words are immune (on = off on every lattice tried: explore52 N = 4–6, A3 N = 7), so the RH-world L4o value is taken as the
off-line value (immunity in continuous positions is assumed, not checked).
npa_cg_lean.py (copy of npa_cg.py; change: lean column generation, rays priced as zero-mass columns): lattice Z_M, off-line pairs
at depths {0.05, 0.2, 0.6, 1.5, 4.0} plus explore52's rays (infinitely deep one-, two- and three-pair clusters at the top level).
Reduction (monotonicity). Removing constraints lowers the min, so for every θ < 1 and N: inside gain(N) ≥ G(N, ⌊θN⌋), where
G(N, n₀) is the gain of the words of length ≥ 3 at levels ≤ n₀ (pair words at all levels). Hence the inside gain tends to 0 only
if G(N, ⌊θN⌋) → 0 for EVERY fixed θ < 1, i.e. only if in the continuum limit the level-2 information sits entirely on the RS
boundary Σ|ξ| = 2, which ζ cannot use (the RS range is open). Test at θ = 3/4: G(4, 3) (= L4o at N = 4) against G(8, 6) (= L4oo,
cutoff N − 2, at N = 8).
A1. The lattice trend was an artifact. Continuous positions against explore52's lattices (min E[#simple on-line]/N, RH world):
| N | P continuous | L4o continuous (lb) | inside gain, continuous | explore52 lattice M: P, L4o, inside gain |
|---|---|---|---|---|
| 3 | 0.847708 | 0.847925 (0.847913) | 0.000217 | 96: 0.848059, –, – |
| 4 | 0.815624 (rerun at tol 2·10⁻⁵: 0.815638) | 0.816569 (0.816569) | 0.000945 (0.000931) | 48: 0.817137, 0.817885, 0.000748 |
| 5 | 0.793924 | 0.795885 (0.795884) | 0.001961 | 40: 0.795205, 0.800269, 0.005064 |
| 6 | 0.778268 | 0.780697 (0.780659) | 0.002429 | 24: 0.806701, 0.814086, 0.007385 |
| 7 | 0.766430 | 0.770179 (0.770174) | 0.003749 | 21: 0.824515, 0.835332, 0.010817 |
| 8 | 0.757058 | 0.761720 (0.761703) | 0.004662 | – (explore52 colgen P 0.757057) |
| 9 | 0.749498 | 0.755867 (0.755739; wall-clock stop, max reduced cost 1.3·10⁻⁴) | 0.006369 (≥ 0.006241) | – |
  The lattice restricts the adversary (values are upper bounds) and M/N fell from 12 to 3, inflating P by up to 0.058 (N = 7) and
  L4o more than P: the lattice inside gain grew 14-fold, the true one 4-fold. The full level-2 gain (continuous, N = 4: L4 0.841747,
  gain 0.026123) is 28× the inside gain at N = 4; with the lattice full gains (slight over-estimates: N = 4 lattice 0.0269 vs
  continuous 0.0261) the inside share is 3.6 % (N = 4), ≈ 6 % (5), ≈ 8 % (6), ≈ 13 % (7: lattice L4 − P = 0.0299, A3), ≈ 18 %
  (8: full gain ≈ 0.026 from explore50's RH-world colgen, P − P3 = 0.0201 in Σm² units plus ≈ 0.006 for the quartics), not
  explore52's 2.8 → 16.4 → 23.5 % (N = 4–6).
A2. Growth, not decay. Inside gain g(N) = 0.000217, 0.000945, 0.001961, 0.002429, 0.003749, 0.004662, 0.006369 (N = 3–9):
  increasing at every step, with an odd/even ripple (steps +0.73, +1.02, +0.47, +1.32, +0.91, +1.71 ·10⁻³; the last is the largest). θ = 3/4 test: G(8, 6) = 0.758622 − 0.757058 = 0.001564 (L4oo,
  lb 0.758615) against G(4, 3) = 0.000945: the fixed-θ gain grew by 65 % from N = 4 to 8, so the monotone lower bound
  g(N) ≥ G(N, ⌊3N/4⌋) is itself rising. θ = 2/3 test: G(6, 4) = 0.778954 − 0.778268 = 0.000686 (L4oo at N = 6, lb 0.778945)
  against G(3, 2) = 0.000217 (L4o at N = 3): tripled. N = 3: P 0.847708, L4o 0.847925 (inside gain 0.000217).
  Fits (least squares on N = 3–9, N = 9 at its master value; M1–M4 also by fit_A.py with scipy curve_fit, logs/fit_A.log,
  M5–M7 in the JS sandbox; residuals are dominated by the parity ripple, 2–5·10⁻⁴):
| model | parameters | max abs residual | g(12), g(24), g(∞) |
|---|---|---|---|
| M5 G₁e^{−c/N} | G₁ = 0.0315, c = 14.8 | 3.2·10⁻⁴ | 0.0092, 0.0170, 0.0315 |
| M7 a + bN (no saturation) | b = 0.00099 per unit N | 5.0·10⁻⁴ | 0.0088, 0.0207, – |
| M4 A(N − N₀)e^{−N/L} (→ 0) | N₀ = 2.99, L ≥ 2964 (grid edge) | 5.3·10⁻⁴ | 0.0087, 0.0203, 0 (peak beyond N ≈ 3000) |
| M2 g∞ − A/N + B/N² | g∞ = 0.0143, A = 0.0946, B = 0.158 | 6.2·10⁻⁴ | 0.0075, 0.0106, 0.0143 |
| M6 g∞(1 − N₀/N)² | g∞ = 0.0124, N₀ = 2.99 | 8.2·10⁻⁴ | 0.0070, 0.0095, 0.0124 |
| M3 A N^α | α = 2.27 (curve_fit; log-log: 2.87) | 3.1·10⁻⁴ | 0.0122, 0.058 (exceeds the full gain: not a limit model) |
| M1 g∞ − A/N | g∞ = 0.0076, A = 0.0248 | 1.5·10⁻³ | 0.0055, 0.0066, 0.0076 (misses N = 9 by 1.5·10⁻³) |
  (N = 3–8 only: M5 G₁ = 0.0240, M4 peak at N ≈ 490; N = 4–8: M5 0.0233, M4 peak ≈ 200: each new N pushed a turnover further out
  and the fitted limit up.) The best fit's limit ≈ 0.03 is the full level-2 gain (≈ 0.026–0.03): the data are consistent with the
  top-level share decaying and the inside gain filling the whole gain; they exclude the quickly saturating M1.
  A model tending to 0 must turn over beyond N ≈ 3000 to fit N ≤ 9. What would decide the size of the limit: N = 12 separates
  M2/M6 (0.0070–0.0075) from M4/M5/M7 (0.0087–0.0092) by ≈ 4 × the parity ripple; no finite N excludes a late decay, only a
  proof that G(N, ⌊θN⌋) → G_∞(θ) > 0 for some θ < 1 would. N = 12 needs 153 L4o features (N = 8: 51 features, 207 MB, 6 min;
  N = 9: 71, 203 MB, 14 min at the wall-clock stop): it needs a leaner master than colgen_simple.py.
A3. Immunity and retention at N = 7 (lattice Z₂₁, depths {0.05, 0.2, 0.6, 1.5, 4.0}, 27,937 rays incl. three-pair clusters;
  logs/L_N7_M21_*.log). L4o off-line + rays = 0.835332 = L4o on (explore52): no off-line column prices negative, the inside words
  stay immune at N = 7. L4 on 0.854409 (first N = 7 value; explore52's run hit the cap), off + rays 0.853805, P 0.824515: off-line
  loss 0.000604, retention 98.0 % (N = 6, M = 24, explore52's 8 depths: 88.6 %). The grids differ, so the fixed-grid series below
  decides whether retention falls with N.
A4. Retention at fixed resolution (M = 3N, the same 5 depths, rays; npa_cg_lean.py; logs/L_N*_M*_*.log):
| N, M | P | L4 on | L4 off + rays | off-line loss | retention (L4off − P)/(L4on − P) |
|---|---|---|---|---|---|
| 4, 12 | 0.871546 | 0.894590 | 0.894590 | 0 | 100 % |
| 5, 15 | 0.851336 | 0.876481 | 0.876481 | 0 | 100 % |
| 6, 18 | 0.837345 | 0.864625 | 0.864332 | 0.000293 | 98.9 % |
| 7, 21 | 0.824515 | 0.854409 | 0.853805 | 0.000604 | 98.0 % |
  At fixed resolution the retention falls slowly (≈ 1 point per unit N from N = 5) and the off-line loss grows with N (more room
  for clusters), but explore52's steeper fall (96.1 → 93.1 → 88.6 %) mixes in its finer grids: at N = 4 refining M = 12 → 48
  (and 5 → 12 depths) alone lowers the retention from 100 % to 96.1 %. Retention keeps falling, slowly; its continuous-angle value
  for N ≥ 5 was not computed (the off-line colgen needs depth variables in the pricing). Floor: every loss sits on top-level words,
  so the retention never drops below the inside share (A1), which rises.

## (B) Clusters: which words can J deep pairs fake for free?
Set-up (explore52 §4). A configuration X contributes v_X(λ, μ) = Re(p_λ p̄_μ)/√(z_λz_μ); a ray is a limit of ε_t v_{X_t}, ε_t → 0
(zero mass, zero objective). A FREE FAKE of a word set S is a ray that vanishes on every diagonal word (λ, λ) ∈ S but not on S.
Equivalently, along the sequence Q_W(X) := |v_X(W)| / max_{(λ,λ)∈S} v_X(λ, λ) → ∞ for some off-diagonal W ∈ S.
- B1 (closure, theorem; any N, any number J of pairs, any level). If S is diagonal-closed ((λ, μ) ∈ S ⇒ (λ, λ), (μ, μ) ∈ S; e.g. P,
  Q4, L4c) it has no free fake. Proof: for every X, G_X = (p_λ p̄_μ/√(z_λz_μ))_{λ,μ} is Hermitian PSD of rank one, so Re G_X is real
  PSD and |v_X(λ, μ)| ≤ √(v_X(λ, λ)v_X(μ, μ)) ≤ ½(v_X(λ,λ) + v_X(μ,μ)); the inequality passes to positive multiples and limits, so
  Q_W ≤ 1 for every W ∈ S and every configuration, deep or not. ∎ (Free fakes for L4c are impossible below the top level as well
  as on it; what remains are priced trade-offs.)
- B2 (L4: only corner words). In L4 = {ℓ(λ) + ℓ(μ) ≤ 4, |λ| = |μ| ≤ N} every word has both diagonal partners except the (3,1) words
  W = (a, b, c)|(n), n = a + b + c. If n + min(a, b, c) ≤ N (say c = min), then |p_a p_b p_c p̄_n| = |p_a p_b|·|p_c p_n|
  ≤ ½(|p_a p_b|² + |p_c p_n|²), and (a, b)|(a, b), (c, n)|(c, n) are diagonal words of L4 (levels a + b < n, c + n ≤ N): Q_W ≤ const.
  So free fakes in L4 can only move corner words n + min(a, b, c) > N, i.e. n > 3N/4: on the top level all (3,1) words, below
  it only n = N − r with min part ≥ r + 1 (first cases: (2,2,2)|(6) at N = 7, (2,2,3)|(7) at N = 8, (3,3,3)|(9) at N = 11).
- B3 (necessary conditions for a corner fake; rigorous part). Let ε_t v_{X_t} → v with v(W) ≠ 0, v = 0 on the L4 diagonal, W at
  level n. Write P_k = p_k(X_t) ≈ e^{ks}c_k + o(e^{ks}), c_k = Σ_j e^{−kΔ_j}w_j^k over the deepest group (pairs with s − s_j bounded;
  shallower groups and on-line points are o(e^{ks})), |c_k| ≤ J. Since |ε P_aP_bP_cP_n| ≤ J⁴εe^{2ns}(1 + o(1)) must stay ≍ 1, εe^{2ns}
  is bounded below; so for every m ∈ [n, N] the pair word ε|P_m|² → 0 forces c_m → 0, and every quartic diagonal (d, m − d), m ∈ [n, N],
  forces c_d c_{m−d} → 0. Hence the deepest group's limit power sum vanishes on the run [n, N]; as it satisfies a linear recurrence of
  order ≤ J with nonzero roots and c_0 = J ≠ 0, a run of J zeros is impossible: J ≥ N − n + 2 (J − 1 vanishing top levels, explore52's
  count). It must also satisfy the covering condition (c_d = 0 or c_{m−d} = 0) on all levels [n, N]. Not proved: that these force more
  pairs than the analogue allows (2J ≤ N); zeros off the run are not bounded by the recurrence (1 + (−1)^k vanishes at every odd k),
  and v(W) can come from subleading groups (c_a may vanish too). Reviewer (second pass): these two gaps are real; B3 is stated
  without them. Example of the subleading mechanism being priced: N odd, deepest group {1, −1} (c_k = 0 for odd k) meets both
  conditions at the top level, but then one of a, b, c is odd and |P_bP_c|²/|P_{N−a}|² stays bounded, so W is priced by the
  diagonal words (a, N − a) and (N): no free fake.
- B4 (numerical search, b_cluster_search.py, logs/B_cluster.log). For a configuration of J off-line pairs only (depths s_j ≤ s0,
  angles θ_j; on-line points are negligible in the deep limit) maximise log Q_W = log|v(W)| − log max_D v(D) over (s_j, θ_j) by
  30 random-restart L-BFGS runs (smooth max, temperature 40), for depth caps s0 = 2, 4, 8, 12. A free fake shows as log Q_W growing
  linearly in s0; a priced word stays bounded. Positive control: in the set P (pair words only) the top triple (2,4)|(6) at N = 6,
  J = 2 gives log Q = 3.0, 5.0, 9.0, 13.0 (slope 1: free, explore52's Claim A); in L4 the same triple gives 0.000 (B1: Q ≤ 1).

  Results, max over s0 of log Q (bounded = priced; no word grows with s0 except the positive control):
| N | J | word (L4) | level | log Q at s0 = 2 / 4 / 8 / 12 |
|---|---|---|---|---|
| 6 | 2 | (2,2,2)\|(6), (1,2,3)\|(6), (1,1,4)\|(6) | top | 0.71/0.65/0.71/0.65; 1.68/1.68/1.52/1.52; 1.13/1.13/1.13/0.81 |
| 6 | 3 | same three | top | 1.75/1.75/1.75/1.75; 2.41/2.41/2.35/2.35; 1.88/1.88/1.77/1.77 |
| 7 | 2 | (2,2,2)\|(6); (2,2,3)\|(7) | N − 1; top | −0.53 (all s0); 1.05/0.88/0.88/0.88 |
| 7 | 3 | (2,2,2)\|(6); (2,2,3)\|(7) | N − 1; top | 0.27 (all s0); 2.06/1.49/2.06/1.36 |
| 8 | 4 | (2,2,3)\|(7); (2,3,3)\|(8); (1,3,4)\|(8) | N − 1; top; top | 1.06/1.06/1.06/0.53; 2.79/2.53/2.53/2.53; 2.84/2.85/2.39/2.55 |
| 9 | 4 | (2,3,3)\|(8); (2,2,4)\|(8) | N − 1; N − 1 | 0.80/0.80/0.80/0.79; 1.93/1.93/1.92/1.93 |
  B result: no free fake at any level for J = 2, 3 (N = 6, 7) and J = 4 (N = 8, 9; at N = 8 the largest J the analogue allows),
  below the top level included ((2,2,2)|(6) at N = 7, (2,2,3)|(7) at N = 8, (2,3,3)|(8) and (2,2,4)|(8) at N = 9); L4c and every
  diagonal-closed set: none possible (B1). What the clusters do is price, as rays do (A3: the off-line L4 loss is 0.0006 at N = 7).

## (C) The input: mixed mock-Gaussian lemma (unconditional, complex zeros)
Notation (Hughes–Rudnick math/0208220, l.136–144, 203–236): L = log T/2π, H = T^a (0 < a ≤ 1), w ≥ 0, ∫w = 1, supp ŵ ⊂ [−c, c],
⟨F⟩ = ∫F(τ)w((τ−T)/H)dτ/H. For f̂ ∈ C_c^∞(ℝ) (f̂(u) = ∫f(x)e^{−2πixu}dx) f(z) = ∫f̂(u)e^{2πizu}du is entire, and
|f(x+iy)| ≪_A (1+|x|)^{−A}e^{2πα|y|} if supp f̂ ⊂ (−α, α). Every nontrivial zero is written ρ = 1/2 + iγ_ρ with γ_ρ ∈ ℂ
(off-line pairs ρ, 1 − ρ̄ have ordinates γ, γ̄); N_f(τ) = Σ_ρ f(L(γ_ρ − τ)) with multiplicity (absolutely convergent);
N̄_f(τ) = (1/2π)∫f(L(r−τ))Ω(r)dr + f(L(i/2−τ)) + f(L(−i/2−τ)), Ω(r) = ½Ψ(¼+½ir) + ½Ψ(¼−½ir) − log π (explicit, no zeros).

LEMMA (mixed words, RS simplex). Let m ≥ 1, f̂_1, …, f̂_m ∈ C_c^∞(ℝ), supp f̂_j ⊂ (−α_j, α_j), Σ_j α_j < 2a, X_j := N_{f_j} − N̄_{f_j}.
Then, with no hypothesis on the zeros, ⟨X_1⋯X_m⟩ = Σ_π Π_{{i,j}∈π} C(f_i, f_j) + O(1/log T), C(f, g) := ∫_ℝ |u| f̂(u) ĝ(−u) du,
the sum over perfect matchings π of {1, …, m} (empty for m odd); the O depends on m, a, c, w and the α_j, ‖f̂_j‖_{C¹}.
Dictionary: f̂_j a bump at frequency ±ξ_j gives the analogue's word p_λ p̄_μ (ξ = k/N), C ↔ E|p_k|² = k, Σ_π ↔ δ_{λμ}z_λ
(Diaconis–Shahshahani), Σα_j < 2 ↔ |λ| = |μ| < N (inside words); the top level Σα_j = 2 is not covered (limit only).

Proof. (1) Weil's explicit formula (HR Lemma 2.1, l.208–213; contour integration of ζ'/ζ, no hypothesis on the zeros, complex
g allowed by linearity) with g_j(u) = e^{−iτu}f̂_j(u/log T)/log T, h_j(r) = ∫g_j e^{iru}du = f_j(L(r−τ)) gives the identity
  X_j(τ) = −(1/log T) Σ_{n≥2} Λ(n)n^{−1/2}[f̂_j(ℓ_n)n^{−iτ} + f̂_j(−ℓ_n)n^{iτ}],  ℓ_n = log n/log T,  only n < T^{α_j} occur.
(2) Multiply out (finite sums): Π_j X_j = (−1/log T)^m Σ_{ε∈{±1}^m} Σ_{n_1..n_m} Π_j Λ(n_j)n_j^{−1/2}f̂_j(−ε_jℓ_{n_j}) e^{iτΣε_j log n_j},
and ⟨e^{iτx}⟩ = e^{iTx}ŵ(−Hx/2π).
(3) Off-diagonal terms vanish identically (HR l.437–443 with unequal α_j). Let M = Π_{ε_j=+1}n_j, M' = Π_{ε_j=−1}n_j (empty
product 1), x = log(M/M'). If M > M' then x ≥ log(1 + 1/M') ≥ 1/(2M') ≥ 1/(2√(MM')) > T^{−Σα_j/2}/2, as MM' = Πn_j < T^{Σα_j};
ŵ(−Hx/2π) ≠ 0 needs |x| ≤ 2πcT^{−a}; since Σα_j/2 < a this is impossible for T ≥ T₀. Same for M < M'. Only M = M' survives.
(4) Diagonal terms (M = M', ŵ(0) = 1). By unique factorisation a tuple splits into prime clusters: indices I with n_i = p^{k_i},
Σ_{i∈I,ε_i=+}k_i = Σ_{i∈I,ε_i=−}k_i =: K ≥ 1, weight (log T)^{−|I|}(log p)^{|I|}p^{−K}·Π|f̂_i| (Λ(p^k) = log p). K = 1 forces |I| = 2,
one index of each sign, k = 1: a matched pair. For K ≥ 2 (r = |I| ≥ 2) the balanced exponent vectors with common sum K number
≤ C(K−1, r₊−1)C(K−1, r₋−1) ≤ K^{r−2}, and Σ_{K≥2}K^{r−2}Σ_p (log p)^r p^{−K} ≪_r Σ_p (log p)^r p^{−2} < ∞ (reviewer's repair: the number of
exponent vectors grows with T, so it must be summed), so O((log T)^{−r}) = O(log^{−2}T) per such cluster; dropping the distinct-prime condition only enlarges absolute sums (bounded: matched pairs are
O(1), step 5). Hence the diagonal part = Σ over bijections {ε=+} → {ε=−} of Π_pairs S_{ij} + O(log^{−2}T) (tuples matching two ways
are clusters with K = 2), S_{ij} = (log T)^{−2}Σ_p (log p)²p^{−1}f̂_i(−ℓ_p)f̂_j(ℓ_p); the sign (−1)^m = ((−1)²)^{m/2} cancels.
Summing over ε: each perfect matching π of {1..m} appears with both orientations of each pair: S_{ij} + S_{ji}.
(5) PNT (θ(x) = x + O(xe^{−c√log x}), unconditional, partial summation): (log T)^{−2}Σ_p (log p)²p^{−1}F(ℓ_p) = ∫_0^∞uF(u)du
+ O_F(1/log T) for F ∈ C_c¹. So S_{ij} + S_{ji} = ∫_0^∞u[f̂_i(−u)f̂_j(u) + f̂_i(u)f̂_j(−u)]du + O(1/log T) = C(f_i, f_j) + O(1/log T).
(HR's min(1,|u|) is |u| here: |u| < min(α_i, α_j) < a ≤ 1.) Combining (3)–(5) proves the lemma. ∎
What is used: only (1) Weil's formula, (3) the arithmetic gap |log(M/M')| ≥ 1/(2 min(M, M')), (4) unique factorisation, (5) PNT.
None of them sees where the zeros are; HR's own theorem (l.151–152, 'will not assume the Riemann Hypothesis' l.129–131) is the
case f_1 = ⋯ = f_m, α_j = α < 2a/m; the polarised cube max α_j < 2a/m is strictly smaller than the simplex Σα_j < 2a used here.

Review (two reviewer subagents, fresh, proof pasted inline). Pass 1: GAP at step (4): the cluster bound was stated per fixed
exponent vector, whose number grows with T; repair (the reviewer's): count K^{r−2} vectors and sum over K, as now written in (4);
it confirmed the explicit formula with complex ordinates, the unequal-support off-diagonal step, the covariance and the signs,
and that no RH is used. Pass 2 (repaired text): 'the repaired summation in C(4) resolves the stated counting issue'; its other
objections concerned (B3), now restated as necessary conditions only.
RS96 itself: no legal copy reachable (IAS publication archive behind a Cloudflare challenge; Rudnick's TAU pages 404; Sarnak's
Princeton preprint page does not list it; Duke is paywalled). Secondary source only: arXiv 2002.00595 (a survey of RS96) restates
RS Theorem 1.1 (smoothed n-level correlation R_n(T, f, h), GL_m cuspidal, m ≤ 3 or Hypothesis H, supp f̂ ⊂ Σ|ξ_j| < 2/m, h = ĝ,
g ∈ C_c^∞) with no RH among the hypotheses (its TeX l.199–206); RH is invoked to order the γ's as real numbers (l.105) and for the
unsmoothed C_n(f, T) version (l.529). AF's l.603 cites Hej94/RS96 'under the Riemann hypothesis' for tr G̃³. The lemma above does
not rely on either: it is HR's argument, unconditional, run on the simplex. Where it stops: at Σα_j = 2a (explore52's top level)
step (3) fails, since |log(M/M')| ≥ 1/(2√(MM')) ≍ T^{−a} no longer beats the window resolution T^{−a}, and the O-term is not
uniform as Σα_j → 2a (T₀ grows); in CUE the DS range ends there too. So the lemma covers exactly the inside words (compact supports
in the open simplex), not the top level.

## Suggested NOGO line
- npafals (Orca, falsification of explore52 → N8; the lead survives, these parts die) explore52's growing inside gain (0.0008 → 0.0108, N = 4–7) and inside share (3 → 24 %) were coarse-lattice artifacts (M/N 12 → 3): continuous positions give L4o − P = 0.0002, 0.0009, 0.0020, 0.0024, 0.0037, 0.0047, 0.0064 (N = 3–9), still rising (fixed-θ gains rising, best fit → 0.03), 1–25 % of the level-2 gain; no free fake from J ≤ 4 pair clusters (none possible for diagonal-closed sets, rank-one PSD; in L4 only corner (3,1) words could be faked); the input (mixed RS-range words, unconditional, complex zeros) is proved on the open simplex Σα < 2 (HR's diagonal argument), not at the boundary Σα = 2 where explore52's small-N gain sits; no bound for ζ [NPAFals](Riemann/Deep/NPAFals/REPORT.md)

RH relevance: N8 (unconditional proportion of simple zeros on the line: which correlation data beyond bandwidth one survive off-line zeros; here: only RS-range words strictly inside the simplex are both unconditional (lemma) and immune, and their gain in the analogue is small at N ≤ 9 but growing).

## Files
- colgen_simple.py — copy of ExploreNPA/colgen_open.py (header names the changes: min-simple objective, L4oo/L4c sets, CUE initial
  columns, warm-started vectorised pricing, tolerance 2·10⁻⁵, wall-clock stop); logs/A_N*_*.log, logs/A_summary.log, logs/summary.log.
- npa_cg_lean.py — lean column generation over explore52's npa_cg.py (copied unmodified as npa_cg.py, with offline_lp.py);
  logs/L_N*_M*_*.log, logs/B_validate_N4_M48.log (reproduces explore52's N = 4, M = 48: L4 on 0.844041, off + rays 0.842988).
- b_cluster_search.py — (B) fake-ratio search; logs/B_cluster.log. fit_A.py — fits (the report's fits were run in the JS sandbox).
- run_A.sh, run_Q.sh, run_Q2.sh — sequential queues (one process at a time, every job under .scratch/memcap.sh 290); peaks ≤ 207 MB
  (colgen) and 287 MB (N = 7 L4 off-line, before the view-leak fix in npa_cg_lean.py); logs/A_N7_L4o_oldcode_killed.log and
  logs/A_N9_L4o_prune4000_killed.log are superseded runs (stopped by npafals, slow settings).
