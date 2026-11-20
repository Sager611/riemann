IDEA: Hörmander-type pole removal (explicit Mittag-Leffler corrections with a band-limited kernel; ∂̄ as heuristic) reduces the continuous Nyman–Beurling rate to an interpolation cost free of 1/ζ′(ρ): under RH, D(λ)² ≤ (1 + o(1))·U_{Z(L⁶)}(L − √L) + O(L^{−3/2}), so RH + Lemma R′ (open tail domination) + the standard contour Lemma C give lim D(λ)² log(1/λ) = Σ_ρ m_ρ²/|ρ|²; referee GAP, repaired 6 of 7 points (REPAIR.md) (node N5)

# ExploreDbar — explore59 (Orca exploration, 20 Nov 2026)

## 1. Transfer
- Tool: Hörmander's weighted L² estimate for ∂̄ in one variable, used to remove poles: E = 1 − ζG must equal 1 on the zero
  divisor of ζ, so one builds G₀ = ΨΦ/ζ from band-limited pieces and kills its remaining poles by ∂̄-corrections ∂̄v = ∂̄χ·G₀
  solved in L²(e^{−ψ}/|ζ|²). The division isometry v ↦ v/ζ moves the weight 1/|ζ|² onto the datum, and v = ζ·(v/ζ) gives it
  back near the pole, so only |ζ| on a circle around a pole cluster enters (Cartan minimum modulus, unconditional), never 1/|ζ′(ρ)|.
  In one variable the minimal correction is a Cauchy integral over the circle with a band-limited kernel (DERIVATION §1, §3).
  After the referee's verdict this heuristic is realised only by explicit Mittag-Leffler corrections (REPAIR (3iii)); the
  weighted ∂̄ estimate and the division isometry are no longer used in any proof step.
- Node: N5, the continuous distance D(λ) of Báez-Duarte–Balazard–Landreau–Saias (dilations θ ∈ [λ,1]; d_N ≥ D(1/N)).
  TREE N5: Burnol's liminf L·D² ≥ Σ m²/|ρ|² (unconditional); upper bounds only under RH + a moment hypothesis on 1/ζ′(ρ)
  (BCF 1211.5191, weakened in Zenodo 22647730, 22650745), or under RH alone with the weak rate (log L)^{5/2+ε}/√L
  (Balazard–de Roton, quoted BCF l.61–64), L = log(1/λ). explore42 found that its zero-evaluator route to an RH rate for D(λ)
  "needs BCF's arithmetic control"; the ∂̄ route below does not.

## 2. Statement (exact status after the referee's GAP verdict; arguments in REPAIR.md, which supersedes DERIVATION §1/§3)
- Low part, exact. V_L = {∫_[0,L] e^{−iut}dμ(u)} ⊂ L²(w dt), w = 1/(2π(1/4+t²)), has Gram e^{−|u−v|/2} (Ornstein–Uhlenbeck), hence
  the explicit kernel K(z,z₀) = (1/4 + z̄₀²)(1 − e^{−iL(z−z̄₀)})/(i(z−z̄₀)) + (1/2 − iz̄₀) + (1/2 + iz̄₀)e^{−iL(z−z̄₀)}, z = t − i(σ − ½),
  K(γ,γ) = (1/4+γ²)L + 1 (referee: correct). U_Z(L) := min{‖E‖² : E ∈ V_L, E = 1 on Z} = 1ᵀK_Z⁻¹1.
- Claim F (proved, unconditional; REPAIR (2)): for every fixed finite zero set, L·U_Z(L) → Σ_Z m_ρ²/|ρ|² — the constant 1 sits at
  the spectral endpoint u = 0 and the shifted Legendre values q_k(0)² = (2k+1)/L sum to m²/L. At finite L the constant depends
  on the rescaled configuration (referee's point (2), accepted); only true multiplicities survive the limit.
- Reduction (proved under RH modulo the standard contour Lemma C; REPAIR (3i)–(3v)): with L′ = L − √L, T₁ = L², T₂ = L⁶,
  G = ΨΦ/ζ − R_tr − R_high, Ψ = 1 − E_low (E_low ∈ V_{L′} the optimal interpolant of 1 on the zeros below L⁶), Φ a one-sided
  Gevrey window, R_tr = Σ_k r_k e^{−(L−1)(s+2k)}/(s+2k) (prefactor from the reproducing kernel, cost ≪ e^{−5√L}), R_high explicit
  Mittag-Leffler corrections with B(x) = ((1 − e^{−ax})/(ax))², a = (L−1)/2 (cost ≪ e^{−cL²}); G is entire with Paley–Wiener–
  Schwartz growth, hence admissible; no ∂̄ solution, trace estimate or condition on ζ′(ρ) is used, and s = 1 needs nothing for
  Burnol's unconstrained D(λ) (subtract G(1)·1 for Nyman's constrained variant). Result:
  D(λ)² ≤ (1 + L^{−1/2})·U_{Z(L⁶)}(L − √L) + O(L^{−3/2}).
- Open: Lemma R′ (tail domination): under RH, L·U_{Z(L⁶)}(L − √L) ≤ L·U_{Z(H)}(L − √L) + ε(H) with ε(H) → 0 uniformly in L. The
  short-interval Schur argument of DERIVATION §5 does not prove it (referee (4), accepted: RH allows ≍ log γ/log log γ zeros in a
  window of length 1/L). Lemma R″ (sup_L L·U < ∞) would give the order D(λ)² ≪ 1/log(1/λ).
- Conditional theorem: RH + Lemma R′ (+ Lemma C) ⇒ lim_{λ→0} D(λ)² log(1/λ) = Σ_ρ m_ρ²/|ρ|² (= 2 + γ − log 4π if all zeros are
  simple); with Burnol's lower bound this is BBLS's conjecture (lim = Σ m_ρ/|ρ|², which Burnol's m² bound shows forces simple
  zeros, math/0112254 §5) in its m²-form, for the continuous problem. Previously known under RH alone: (log L)^{5/2+ε}/√L.
- Where RH enters: the |1/ζ| bounds off the line used for the growth of G, and the location of the zeros below L⁶ (an off-line
  zero there costs its Hardy-space evaluator floor ≈ 2b/|ρ|², b = β − ½).

## 3. Tests (double precision; 40-digit mpmath cross-check of the worst-conditioned case agrees to 4.4·10⁻¹⁰; logs/)
Smallest case and scale (zeros from WangAudit/zeros1250.txt, both signs ±γ; C = 2 + γ_E − log 4π = 0.0461914; C_K = truncated sum):
| zeros K | L = 7.155 | L = 10 | L = 20 | L = 40 | C_K | valid for L > log(γ_K/2π) |
|---|---|---|---|---|---|---|
| 100 | 0.04164 | 0.03999 | 0.04036 | 0.03995 | 0.03997 | 3.63 |
| 400 | 0.04559 | 0.04373 | 0.04397 | 0.04352 | 0.04353 | 4.68 |
| 1000 | 0.04706 | 0.04505 | 0.04521 | 0.04473 | 0.04475 | 5.42 |
(entries L·U_Z(L)). Below the Nyquist threshold the interpolation is impossible (cond ≥ 10¹⁷); above it cond ≤ 2·10⁴ and
L·U → C_K with cross-terms of ±1 % at L = 7–20 (dbar_U.py, peak 272 MB).
- J_{−1}-freeness (dbar_pair.py (a)): split γ₁ = 14.1347 into γ₁ ± η/2. L·U stays bounded for η = 2 … 10⁻⁴ (J_{−1}-type proxy
  |F′(ρ)|⁻² ~ η⁻² → 10⁸): at L = 20, 0.0539 (η = 2) → 0.07358 (η ≤ 10⁻³); at L = 80, 0.0538 → 0.07348. Predicted: two resolved
  zeros 0.053530, one double zero (m² = 4) 0.073526. The crossover is at η ≈ 1/L, as the band-edge law says.
- Trivial zeros pin the band (dbar_pair.py (b)): ν_k = |Ψ(−2k)|e^{−L(2k+½)}/|ζ′(−2k)| for the optimal E_low (K = 400):
  ε = 0: log₁₀ν_k ∈ [−0.8, +0.2] at L = 10, [−1.6, −0.5] at L = 40 (an O(1) penalty: the full band cannot be used);
  ε = 0.1: ≤ −2.4 (L = 10), ≤ −6.5 (L = 40); ε = 0.25: ≤ −12.1 (L = 40), while L′·U stays 0.0434–0.0443. This is the
  consistency check of the sketch: with L′ > L the construction would beat Burnol's lower bound, and it fails exactly here.
- Clusters (dbar_cluster.py): a triple at γ₁ and a pair at γ₅ with spacing h = 10⁻³ cost m² = 9 and 4: L·U / prediction 0.129046
  = 1.0011, 0.9978, 0.9999 at L = 20, 40, 80 (h = 10⁻²: 1.0005, 0.9955, 0.9909, resolving as Lh → 1).
- Evidence for Lemma R′, not a proof (dbar_riesz.py, 800 zeros, no clusters): Riesz bounds of the normalised evaluators [0.072, 2.006] at
  L = 7.155, [0.611, 1.371] at 20, [0.923, 1.081] at 80, closing as density/Nyquist → 0 (Schur bound O((log L)²/L), DERIVATION §5).
- Finite-L efficiency: at L = 7.155 (N = 1280) explore38's lower bound has L·B = 0.0360 ≤ L·D, against L·U ≈ 0.047: V_L ignores
  ζ's phase e^{−iθ(t)}, which adds local band log(γ/2π). Both sides tend to C; the second-order term is not touched here.

## 4. Adversarial
- Davenport–Heilbronn (dbar_dh.py): 106 on-line zeros below 180 (sign changes of the completed Λ, |Im Λ|/|Λ| ≤ 1.2·10⁻¹⁷) plus the
  off-line 0.808517+85.699348i, 0.650830+114.163343i, 0.574356+166.479306i, 0.724258+176.702461i (refined): 106 + 2·4 = 114 = N(180)
  (smooth count 114.5), so the list is complete. On-line part: L·U_on = 0.1859, 0.1849, 0.1820, 0.1832, 0.1836, 0.1835 at
  L = 7.155 … 160 → Σ_on 1/|ρ|² = 0.18362, exactly as for ζ. All zeros: L·U_all = 0.1888 → 0.2240, U_all → 1.40·10⁻³ at L = 160,
  growing like L·(floor 2.5·10⁻⁴), the Hardy-space evaluators of the four right zeros. Moved-pair models (ζ's γ₁ or a DH-like
  zero moved to β = ½ + b): same floors (b = 0.3085 at γ = 14.13: L·U = 0.525 at L = 80).
- Verdict: the mechanism (band-edge m², J_{−1}-free) is RH-blind; it holds verbatim for DH's on-line zeros and for any zero set on
  the line. ζ is told apart from DH only by the zeros' location, i.e. by RH itself. The result is a conditional rate, not a step
  toward RH; it does not touch the integer-dilation penalty d_N² − D(1/N)² (explore38), the one place where ζ and DH differ.

## 5. Novelty (adversarial)
- arXiv API: Ohsawa–Takegoshi, Hörmander, ∂̄, "L2 estimates", Suita, corona × zeta / RH / Nyman: no hits. Monthly listings
  Oct–Nov 2026 (math.NT, CA, CV, SP, PR, FA, CO): 7 "Beurling" titles, none on Nyman–Beurling.
- Zenodo "Nyman-Beurling" (60 records). Closest: 22647730 (7 Sep 2026; discrete d_N, BCF polynomial, hypothesis weakened to
  Σ|ζ′(ρ)|⁻²(1+|γ|)^{−2+δ₀} < ∞) and 22650745 (Aug 2026; J_{−p} family). Both list "any bound of the optimal order under RH alone"
  as open (22647730 remark "open questions" (ii); 22650745 §Open problems item 2, citing Balazard–de Roton as best under RH alone).
  17957671 (Helson-blur) claims an unconditional rate for free dilations, which would imply RH; not examined.
- GXL evidence (3 claims, incl. "under RH alone D(λ)² log(1/λ) → 2 + γ − log 4π"): no support or refutation; only Burnol's lower
  bound and the BBLS-conjecture passage of math/0112254 §5 surface. arXiv web search (approach arxiv): no relevant result.
- Classics read: BBLS 2000 conjecture (via Burnol), Burnol 2002, BCF 2012 TeX, Burnol's arXiv abstracts. Not read: BBLS 2005
  (Ramanujan J.), Balazard–Saias "Notes IV", Burnol's Sonine papers (math/0105120, 0208121, 0203120, 0407443) beyond abstracts:
  residual risk that the continuous RH-alone asymptotic is there. Repo: no RH upper bound for D(λ) in TREE/NOGO/COROLLARIES;
  explore42 expected the opposite.

## 6. What it gives, and limits
- Gives (conditional on RH, sketch): lim L·D(λ)² = Σ m²/|ρ|², the exact constant, replacing (log L)^{5/2+ε}/√L, with no
  moment hypothesis and no simplicity; for the continuous problem the J_{−1} sums of BCF and the two Zenodo papers are an
  artefact of approximants that copy the poles of 1/ζ (BCF's V_N). Unconditional by-product (sketch): D(λ)² is bounded above, up
  to 1 + o(1) and O(L⁻²), by the interpolation cost of the zeros below L⁶.
- Does not give: anything on RH (RH-blind mechanism, §4); the discrete d_N: the correction is not a Dirichlet polynomial, so
  RH ⇒ d_N² ~ C/log N becomes exactly the bound P(N) = d_N² − D(1/N)² = o(1/log N) on explore38's integer-dilation penalty
  (≈ 0.064/log²N for ζ, 0.020/log N for DH): the arithmetic left after the ∂̄ step; the second-order term (V_L misses the θ-phase).
- Open before COROLLARIES/DISCOVERIES: Lemma R′ (tail domination, REPAIR (4)); a written check of the standard contour Lemma C
  (REPAIR §C); then a second referee pass on REPAIR.md and the unread French sources. Suggested next: one Orca agent on Lemma R′
  (candidate dbarnb: multi-scale blocking of RH-allowed clusters, or a counterexample configuration), one reviewer.

## 7. Files (Riemann/Deep/ExploreDbar/)
REPAIR.md (point-by-point answer to the referee, the repaired construction, Lemmas R′ and C); DERIVATION.md (setup, kernel, first construction, consistency checks); dbar_U.py (kernel, reproducing check, U tables; peak 272 MB);
dbar_pair.py (split pair, ν_k, off-line models; 89 MB); dbar_dh.py (DH zeros and U_DH; 23 MB); dbar_mpcheck.py (40-digit check); dbar_riesz.py (Riesz bounds; 246 MB); dbar_cluster.py (m = 2, 3 clusters);
logs/ (all outputs, dh_zeros_on.npy, dh_zeros_off.npy).

## Suggested NOGO line
- explore59 (Orca, exploration: complex analysis: Hörmander's L² method for ∂̄ — pole removal, realised as explicit Mittag-Leffler corrections with a band-limited kernel → N5, continuous Nyman–Beurling distance D(λ)) IDEA, conditional and incomplete: the low part is an explicit Ornstein–Uhlenbeck interpolation (K(γ,γ) = (¼+γ²)L + 1) whose fixed-zero-set limit is m²/(L|ρ|²) (shifted Legendre at the spectral endpoint; split pair η = 2…10⁻⁴: L·U 0.0539 → 0.07358 = the m² = 4 value, J_{−1} proxy 10⁸); under RH (+ standard contour Lemma C) D(λ)² ≤ (1 + L^{−1/2})U_{Z(L⁶)}(L − √L) + O(L^{−3/2}) with no condition on ζ′(ρ) (trivial zeros force the band margin √L: ν_k = O(1) at margin 0); the limit Σ m_ρ²/|ρ|² (BBLS–Burnol, continuous case; was (log L)^{5/2}/√L) still needs the open tail-domination Lemma R′ (referee GAP; 6 of 7 points repaired, REPAIR.md); RH-blind (DH on-line part identical, 0.1836; off-line floor 2.5·10⁻⁴), no effect on d_N (remaining step: explore38's integer-dilation penalty = o(1/log N)) [ExploreDbar](Riemann/Deep/ExploreDbar/REPORT.md)

RH relevance: N5 (conditional rate for the continuous Nyman–Beurling distance under RH; no move toward RH itself).
