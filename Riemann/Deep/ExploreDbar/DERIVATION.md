# ExploreDbar — derivation notes (explore59, 20 Nov 2026)

> Referee verdict GAP (20 Nov ~09:40 UTC). REPAIR.md supersedes §1 (∂̄ step → explicit Mittag-Leffler corrections), §3 (window,
> trivial-zero kernels, growth, closure) and replaces Lemma R by the open Lemma R′; §2, §4, §6, §7 stand.

## 0. Problem (TREE N5, continuous Nyman–Beurling distance)
- BBLS's D(λ): distance in L²(0,∞) from 1_(0,1) to span{ {θ/x} : λ ≤ θ ≤ 1 }. Mellin–Plancherel on Re s = 1/2
  ({θ/x} ↦ −θ^s ζ(s)/s, 1_(0,1) ↦ 1/s): D(λ)² = inf_G ∫ |1 − ζG|² w dt, w = 1/(2π(1/4+t²)), G ∈ M_λ = closed span of θ^s,
  i.e. G(s) = ∫_[0,L] e^{−us} dν(u), L = log(1/λ). d_N ≥ D(1/N).
- Known: liminf L·D² ≥ Σ_ρ m_ρ²/|ρ|² (Burnol math/0103058, unconditional). BBLS conjectured lim = Σ m_ρ/|ρ|²;
  Burnol (math/0112254 §5) noted that this forces RH + simple zeros. Upper bounds: BCF 1211.5191 (RH + Σ_{γ≤T}|ζ′(ρ)|⁻² ≪ T^{3/2−δ})
  gives d_N² ~ C/log N, hence lim L·D² = C; under RH alone only Balazard–de Roton d_N² ≪ (log log N)^{5/2+ε}/√log N (quoted BCF l.61–64),
  so D² ≪ (log L)^{5/2+ε}/√L. explore42 (ExploreCarleson §5): "an RH rate for D(λ) this way [zero evaluators] needs BCF's arithmetic control".

## 1. The ∂̄ transfer (Hörmander, one variable)
- Write E = 1 − ζG. G entire ⇔ E = 1 at every zero of ζ (with multiplicity), trivial zeros included.
- Hörmander/Cauchy step: start from G₀ = Ψ·Φ/ζ with Ψ, Φ band-limited, Ψ = 0 at the low zeros; the remaining poles of G₀
  (high zeros, trivial zeros) are removed by ∂̄-corrections v with ∂̄v = ∂̄χ·G₀ (χ a cutoff around each pole cluster),
  solved in L²(e^{−ψ}) with a weight whose curvature lives where the poles are; in one variable the minimal-norm correction
  can be written as Cauchy integrals over circles c around clusters:
  R_c(s) = (1/2πi)∮_c Ψ(w)Φ(w)B(s − w)/(ζ(w)(s − w)) dw, B band-limited with B(0) = 1.
- Division isometry (the point of the transfer): v ∈ L²(e^{−ψ}/|ζ|²) ⇔ v/ζ ∈ L²(e^{−ψ}); the datum on c is ∂̄χ/ζ ~ 1/(r·|ζ|_c),
  and v = ζ·(v/ζ) is multiplied back by |ζ| ≍ |ζ|_c near c, so the size of ζ′(ρ) cancels. Only |ζ| on the circle (Cartan
  minimum modulus, unconditional: log|ζ| ≥ −C log γ · log(e/r) on suitable circles) enters, never 1/|ζ′(ρ)| alone.
- Where BCF lose: their V_N = Σ μ(n)(1 − log n/log N)n^{−s} mimics 1/ζ, so its error carries the residues 1/ζ′(ρ) (J_{−1} sums).

## 2. The low part is an explicit band-edge interpolation problem
- V_L := {E(t) = ∫_[0,L] e^{−iut} dμ(u)} ⊂ L²(w); ∫ e^{−i(u−v)t} w dt = e^{−|u−v|/2} (Ornstein–Uhlenbeck), so with a = 1/2
  ‖E‖² = ∫∫ e^{−a|u−v|} dμ dμ̄ and the reproducing kernel is explicit (OU inverse = (a² − ∂²)/(2a) + endpoint masses):
  K(z, z₀) = (1/2a)[(a² + z̄₀²)(1 − e^{−iL(z−z̄₀)})/(i(z − z̄₀)) + (a − iz̄₀) + (a + iz̄₀)e^{−iL(z−z̄₀)}],  z = t − i(σ − 1/2).
  K(γ,γ) = (1/4 + γ²)L + 1, so one on-line zero costs 1/K ≈ 1/(L|ρ|²): Burnol's per-zero constant. Checked: reproducing
  property by quadrature (dbar_U.py check: ⟨k₁,k₂⟩ 18.8397−47.0088i vs K 18.8609−47.0182i, tail 1/4000).
- U_Z(L) := min{‖E‖² : E ∈ V_L, E = 1 on Z} = 1ᵀ[K(z_i,z_j)]⁻¹1.
- m² from the band edge: the constant 1 sits at the spectral endpoint u = 0; in the centred band the data (1, 0, …, 0) of a
  multiplicity-m zero become (1, iL/2, (iL/2)², …), i.e. evaluation at the band edge, where Σ_{k<m} (2k+1) = m² (Legendre
  P_k(1) = 1). So a cluster of m zeros closer than ~1/L costs m²/(L|ρ|²), whatever ζ′ is there.

## 3. Construction (RH assumed; proof sketch — Lemma R below is the one open technical step)
Parameters: ε = L^{−1/2}, L′ = (1 − ε)L, T₁ = L², T₂ = L⁶ (zero density (1/2π)log(T₂/2π) ≤ 6 log L/(2π) ≪ L′/(2π));
the window costs O(1/T₁) = O(L⁻²) = o(1/L). (First draft used T₁ = e^{L/4}, T₂ = e^{L/3}; then C_L below is e^{L/3}L and the
margin εL·5/4 does not beat it, so the polynomial choice is needed.)
1. Low zeros. Z_low = {ρ : |γ| ≤ T₂} (on the line by RH). E_low ∈ V_{L′} the minimal-norm interpolant, E_low = 1 on Z_low
   (Hermite conditions at multiple zeros / at clusters of diameter < c/L). Ψ := 1 − E_low vanishes on Z_low.
   Lemma R (open here): after blocking Z_low into clusters at scale c/L, the normalised evaluators of V_{L′} at the clusters
   (band-edge data) form a Riesz system with constants depending only on c, the density ratio ≤ 1/3 and the cluster size,
   growing at most like e^{O(m)}; then L′·‖E_low‖² → Σ_ρ m_ρ²/|ρ|² by dominated convergence (cluster sizes m ≪ log γ, weights 1/γ²).
2. Window. Φ(t) = F(t/T₁), F = h²/h(0)², h = ∫φ(ξ)e^{ixξ}dξ with φ ≥ 0 an even Gevrey-2 bump on [−½, ½]: F even, type ≤ 1,
   0 ≤ F ≤ 1, F(0) = 1, 1 − F(x) ≪ x², F(x) ≪ exp(−c|x|^{1/2}) (rescale c to 1). G₀ := ΨΦ/ζ has poles only at |γ| > T₂ and at the trivial zeros −2k.
3. Trivial zeros. R_tr(s) = Σ_k r_k e^{−L(s+2k)} β(s+2k)/(s+2k), r_k = ΨΦ(−2k)/ζ′(−2k), β of type ≤ εL/2, β(0) = 1.
   Spectrum of each term ⊂ [L − εL/2, L] ⊂ [0, L]. Since |Ψ(−2k)| ≤ C_L e^{L′(2k+1/2)} with C_L = 1 + ‖μ‖_TV ≪ #{|γ| ≤ L⁶} ≪ L⁷, on the critical line
   |r_k e^{−L(s+2k)}β| ≤ C_L e^{−εL(2k+1/2)/2}/|ζ′(−2k)|; |ζ′(−2k)| = (2k)! ζ(2k+1)/(2(2π)^{2k}) grows factorially.
4. High zeros. For each cluster c above T₂ (circle of radius ≍ 1/log γ chosen by Cartan's lemma so that
   log|ζ| ≥ −C log γ · log log γ on c), R_c as in §1. On the line |Ψ| ≤ C, |Φ| ≤ exp(−(γ/T₁)^{1/2}) ≤ e^{−L²} for γ ≥ L⁶,
   which beats γ^{C log log γ}, so ‖ζ R_high‖² ≪ e^{−L²}.
5. Membership. G := G₀ − R_tr − R_high is entire and |G(s)| ≤ C_L(1+|s|)^{N_L} e^{L max(0,−σ)} (left: 1/ζ decays
   factorially off the trivial zeros, R_tr has type L; right: RH gives |1/ζ| ≤ t^{o(1)} for σ ≥ 1/2 + 1/log log t, maximum
   principle on circles nearer the line). Paley–Wiener–Schwartz: G is the Fourier–Laplace transform of a distribution on
   [0, L]; convolving with an approximate identity (§7: spectrum kept in [0, L − 1]) puts ζG/s in the closed NB span, error → 0 by
   dominated convergence (|ζG/s| ∈ L²).
6. Error. On the line 1 − ζG = (1 − Φ) + Φ E_low + ζR_tr + ζR_high, so
   D(λ)² ≤ (1 + o(1))‖E_low‖² + O(1/T₁) + O(e^{−εL}) ⇒ limsup L·D(λ)² ≤ Σ_ρ m_ρ²/|ρ|².
   With Burnol: RH (+ Lemma R) ⇒ lim_{λ→0} D(λ)² log(1/λ) = Σ_ρ m_ρ²/|ρ|²; no hypothesis on ζ′(ρ), multiple zeros allowed.

## 4. Consistency checks
- Band margin is forced. With L′ > L step 3 would beat Burnol's bound; it fails because then ν_k := |Ψ(−2k)|e^{−L(2k+1/2)}/|ζ′(−2k)|
  is O(1)·e^{(L′−L)(2k+1/2)}: already at ε = 0 the numerics give log₁₀ν_k ∈ [−0.8, +0.2] (L = 10), [−1.6, −0.5] (L = 40), while
  ε = 0.1 gives ≤ 10^{−2.4} (L = 10), ≤ 10^{−6.5} (L = 40), ε = 0.25 ≤ 10^{−12.1} (L = 40) (dbar_pair (b)). The Γ-factor's
  trivial zeros pin the band; they are why V_{L+Δ} cannot be used.
- Off-line zeros. A zero right of the line has z₀ = γ − ib, K(z₀,z₀) → ((1/4 + γ² − b²)/(2b) + 1/2 + b) (bounded), so U ≥ floor
  ≈ 2b/|ρ|²: the construction fails exactly at off-line zeros (the Hardy-space evaluator bound), as it must.
- Finite-L efficiency. V_L ignores ζ's phase e^{−iθ(t)}, which adds local band log(γ/2π) (explore38's diagonal L + log(γ/2π)):
  L·U ≈ 0.047 at L = 7.155 (1000 zeros) against explore38's L·B = 0.0360 ≤ L·D. Both tend to C; a sharp finite-L
  construction needs Sonine-type spaces e^{−iθ}V_L.

## 5. Lemma R — proof idea (not written in full)
- For real t ≠ s: |K(t,s)| ≤ (¼ + s²)·2/|t − s| + 2|s| + 1 and K(t,t) = (¼ + t²)L + 1, K Hermitian, so the normalised entries obey
  |A_ts| ≲ (2/(L|t − s|))·min(|s|/|t|, |t|/|s|) + 2/(L·max(|t|,|s|)).
- Schur test for zeros below T₂ = L⁶ separated by ≥ c/L: the near sum Σ_j 2/(L|γ_i − γ_j|) ≲ (2 log T₂/(πL))·log(LT₂) ≍ (log L)²/L
  and the far sum Σ_j 2/(L|γ_j|) ≲ (log T₂)²/(2πL) ≍ (log L)²/L, so ‖A − I‖ = O((log L)²/L) → 0: Riesz bounds → 1.
- Clusters (gaps < c/L): replace the evaluators of a cluster by its jet functionals (orders < m) at one point; under RH
  m ≤ N(γ + c/L) − N(γ) ≪ log γ/log log γ ≤ C log L/log log L. Inside a cluster the normalised jet Gram depends only on the
  rescaled configuration (compact set) and m, with condition number e^{O(m)} = L^{o(1)}; between clusters at distance ≥ c/L the
  Schur bound above applies with an O(1/c) nearest-neighbour term (take c large).
- Limit: each cluster's cost is the band-edge Christoffel cost of its jet, → m_ρ²/|ρ|² once L resolves it (§2), and is
  ≤ C^m m²/(L|ρ|²) uniformly; Σ_ρ γ^{o(1)} log γ/γ² < ∞ gives dominated convergence: L′·U → Σ_ρ m_ρ²/|ρ|².
- Numerics (dbar_riesz.py, 800 zeros, min gap 0.221, no clusters): Riesz bounds of the normalised evaluator Gram
  [0.072, 2.006] at L = 7.155 (density/Nyquist 0.73), [0.175, 1.861] at 10 (0.52), [0.611, 1.371] at 20 (0.26),
  [0.776, 1.225] at 40 (0.13), [0.923, 1.081] at 80 (0.065): they close up as the density ratio → 0, as the Schur bound says.
- Clusters (dbar_cluster.py): triple at γ₁, pair at γ₅, spacing 10⁻³: L·U equals Σ m²/|ρ|² (0.129046) to 0.1 % at L = 20–80.

## 6. What remains arithmetic (discrete d_N)
- d_N² = D(1/N)² + P(N) with P the integer-dilation penalty (explore38). With §3, RH ⇒ d_N² ~ C/log N ⇔ P(N) = o(1/log N);
  explore38 measured P ≈ 0.064/log²N for ζ (second order, N ≤ 1280) and 0.020/log N for Davenport–Heilbronn (first order).
  So after the ∂̄ step the open RH-alone question for d_N (22647730 remark (ii), 22650745 §Open problems) is exactly an
  estimate of P, which is where the Euler product must enter.

## 7. Step 5 made precise (closure)
- Run §3 with L replaced by L − 1 everywhere (R_tr kernels e^{−(L−1)(s+2k)}β(s+2k)); the cost changes by a factor L/(L − 1) → 1.
  Then |G(s)| ≤ C_L(1+|s|)^{N_L} e^{(L−1)max(0,−σ)}, so by Paley–Wiener–Schwartz G = FL(T), T a distribution of finite order
  supported in [0, L − 1]. For φ_δ ≥ 0 smooth on [0, δ], δ < 1, ∫φ_δ = 1: h_δ = T * φ_δ ∈ C^∞, supp ⊂ [0, L], and
  G_δ = G·φ̂_δ with |φ̂_δ(½ + it)| ≤ 1, φ̂_δ → 1 pointwise. ζG_δ/s = ∫ h_δ(u) ζ(s)e^{−us}/s du is an L²-Bochner integral of
  span elements (u ↦ ζe^{−us}/s is L²-continuous), hence in the closed span, and ‖ζ(G_δ − G)/s‖ → 0 by dominated convergence
  with majorant 2|ζG/s| = 2|1 − E|/|s| ∈ L². (A contraction u ↦ αu instead would need control of G on the line Re s = α/2.)
- Pole removal at a cluster c (steps 3–4) is Mittag-Leffler with kernel B: R_c(s) = Σ_{p∈c} Res_{w=p}[F(w)B(s−w)/(s−w)],
  F = ΨΦ/ζ, B entire, B(0) = 1; principal parts match F's, R_c is holomorphic off c, and for s on the line
  ζR_c = ζ·I_{c′}(s) − ΨΦ(s)·1_{s inside c′}, I_{c′} = (1/2πi)∮_{c′} F(w)B(s−w)/(s−w) dw, so only sup_{c′}|F| enters.
