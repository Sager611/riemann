NO IDEA: on N5 the Carleson tools act only through Burnol's Sonine-space decomposition D(λ)² ≥ F(λ) + (zero-evaluator form) (equality is his Theorem A if it transfers; F a ζ-free prolate floor that corrects explore38's B: 92 % of its N = 2 gap; brackets 0.2–0.7 % at N ≤ 10), where Carleson's separation constant is the wrong control (λ_min ≈ 1.14σ_min^1.84 → 0 for ζ, CUE and DH alike while the form obeys a multiplicity box law, Burnol's m² at finite N) and an RH rate needs Pavlov A₂ control of |ζ|^{±2} near the line, BCF's input; on N8 Carleson's constant of the shifted zeros is exactly 2ε|ξ′(ρ)|/|ξ(ρ+2ε)|, a GUE gap statistic

# explore42: Carleson measures and interpolating sequences → N5 (N8 side test) (Orca exploration, 12 Nov 2026, 13:01–14:20 UTC)

- Derivation: [DERIVATION.md](DERIVATION.md).
- Code: [code/](code/). Rerun with `uv run --no-project --with numpy --with scipy --with mpmath python3 <file>`, OMP etc. = 1. `s0_cache.py` first; it writes `.scratch/explore42/zeros_cache.npz`.
- Logs: [logs/](logs/).
- Notation (explore38): λ = 1/N, L = log N, C = 2 + γ − log 4π. D(λ) is the continuous NB distance, d_N ≥ D(1/N). Y_ρ are Burnol's zero vectors.
- Φ_N(t) = t log N + 2θ(t) + 6 arctan 2t, so Φ_N′ = log N + log(t/2π) + O(t⁻²). σ is a gap measured in Clark cells, Φ_N′Δγ/2π.

## 1. Bridges (measured before choosing)
| node | Carleson-field reading | verdict |
|---|---|---|
| N5 | explore38's zero Gram (1 − e^{iΔΦ_N})/(iΔγ) is a model-space kernel: Shapiro–Shields/Carleson-measure bounds for it, Nikolski's distance formula for off-line points | chosen |
| N8 | Carleson's constant δ(T) of the points ½ + ε + iγ | side test (§3e): a gap statistic |
| N1 | Carleson embedding of Σ_γ δ_γ in PW_a = the frame bound of Weil's form | covered: explore19, explore41 (Landau count), explore36 |
| Nikolski 1995 | dist(1, BH²)² = 1 − \|B(1)\|² = Beurling's floor | covered: explore13 (Hankel/AAK). Localization from certified d_N is empty: excluding an off-line zero above 3·10¹² needs (2β−1)/\|ρ\|² > d_N² ≈ C/log N, i.e. log N > 5.6·10²³ |

## 2. Transfer and statements
- **S1 (structure; Burnol's Sonine picture in explore38's normalization).** VC_λ ⊂ L²[λ,∞) and C_λ ⊂ L²(0,1), so C_λ ⊂ W_λ = {f ∈ L²(0,1) : Vf = 0 on (0,λ)}, a Sonine-type space. Hence D(λ)² = F(λ) + ‖P_{W⊖C}χ₁‖².
  - F(λ) = dist²(χ₁, W_λ) is ζ-free: F = ⟨q, (P(I−T)P)⁻¹q⟩. Here T is the even prolate concentration operator on (0,λ) with band (0,1), q = CJχ₁ = −2[(1 + log a)Ci(a) + J(a)] with a = 2πv, and P projects onto H₃ = {h ⊥ 1, log(λ/v), log²(λ/v)}.
  - F ≥ ε := ‖(1−Q_λ)Vχ₁‖², and Burnol's B omits even ε.
  - By duality W⊖C contains the closed span of W's evaluators at the zeros (P_WY_ρ). Equality is synthesis: C_λ = {f ∈ W_λ : f vanishes at every zero}.
  - Burnol proves the analogue in his extended Sonine space L_a (math/0203120, Theorem A: for a < 1, L_a = \bar P_a ⊕ Y_a; P_a = co-Poisson sums of g ∈ L¹(a,1/a), a = N^{−1/2}).
  - W_λ is the 3-dimensional extension {F = 0 on (0,a), CF ∈ span{1, log, log²} on (0,a)} of de Branges' K_a (after J and dilation); L_a is a 2-dimensional extension. The transfer is a finite-dimensional bookkeeping step that is not proved here (reviewer: UNCLEAR). The brackets of §3b support it.
  - If it holds, explore38's "D ≈ B within 2.6–4.5 %" (Burnol's guess) is exact once the floor F is added and Y_ρ is replaced by P_WY_ρ.
- **S2 (Carleson/Shapiro–Shields transfer).** With G̃ the normalized Gram and |c̃|² = Σ|c_ρ|²/Φ_N′(γ) ≈ Σ 1/(|ρ|²(L + log(γ/2π))), the transfer gives |c̃|²/λ_max ≤ cᴴH⁻¹c ≤ |c̃|²/λ_min.
  - Carleson: λ_min > 0 needs uniform separation in Clark cells, and λ_max < ∞ needs a box condition.
  - Uniform separation (inf σ > 0), with the uniform density gap L/2π below the Clark density, would make the evaluators a Riesz sequence (standard sufficiency for de Branges spaces with slowly varying phase; not checked for W_λ's exact phase) and give D(λ)² ≤ F + |c̃|²/λ_min ≍ 1/log(1/λ). Under GUE, inf σ = 0.
- **S3 (box law, measured).** The form is controlled by multiplicity at Clark-cell scale, not by separation.
  - A cluster of m unresolved zeros contributes ≈ m²|c|²/Φ′ (Burnol's m², Hilbert-matrix identity Σ(H_m⁻¹)_{ij} = m²), i.e. ≤ κ|c̃|² with κ the box count.
  - Coalescing pairs double the form while λ_min → 0 like 1.6σ² (§3d).
- **S4 (N8).** For ε > B − ½ (any ε > 0 under RH), ∏_ρ b_{ρ+ε} = ξ(s−ε)/ξ(s+ε). So Carleson's constant is δ_n(ε) = (2β_n + 2ε − 1)|ξ′(ρ_n)|/|ξ(ρ_n + 2ε)| (DERIVATION §5).

## 3. Tests
- **(a) Floor, smallest cases** (s14; ε agrees with explore38's power series to 7 digits, s11).
  - F − ε = 8.07·10⁻³ (N = 2), 5.12·10⁻⁴ (3), 4.5·10⁻⁵ (4), 6.0·10⁻⁶ (5), 7.7·10⁻⁹ (10), 8·10⁻¹² (20).
  - Top prolate eigenvalue 0.038 (N = 2), 5.2·10⁻⁵ (N = 10).
- **(b) Brackets for D(1/N)².**
  - Lower bound (s15): ‖P_{R+span y}χ₁‖² with K = 800 zero pairs, plus explore38's diagonal tail. Rigorous for finite K without the tail, and increasing in K: N = 2 gives 0.228874 / 0.229514 / 0.229725 at K = 100 / 400 / 800.
  - Upper bound: Farey relaxation, all a/b ∈ [1/N,1] with b ≤ 100 (s12, reproduces explore38's Bd ≤ 60 values).

| N | ε | F | B (explore38, K = 800 + tail) | lower (F ⊕ zero evaluators) | Farey Bd = 100 | width |
|---|---|---|---|---|---|---|
| 2 | 0.2078267 | 0.2158957 | 0.0136489 | 0.2297253 | 0.2302010 | 0.21 % |
| 3 | 0.0435314 | 0.0440438 | 0.0142912 | 0.0583486 | 0.0585506 | 0.35 % |
| 5 | 0.0043026 | 0.0043085 | 0.0131340 | 0.0174411 | 0.0175596 | 0.68 % |
| 10 | 0.0001488 | 0.0001488 | 0.0116481 | 0.0117969 | 0.0118802 | 0.70 % |

  - At N = 2, Farey − (B + ε) = 8.73·10⁻³. F − ε accounts for 8.07·10⁻³ (92 %), the R–Y coupling for 1.8·10⁻⁴, and K/Bd truncation for the rest.
  - The Farey sequence still falls by 2·10⁻⁵ per step at Bd = 80→100. At N = 10, extrapolations (s13) place D² 0–2·10⁻⁴/L above B + ε, which is consistent with exactness.
- **(c) Two-sided estimate on the head** (s16; K = 400 pairs; B·L in [bounds]).

| N | λ_min, λ_max | B·L | [\|c̃\|²/λ_max, \|c̃\|²/λ_min]·L | B/\|c̃\|² |
|---|---|---|---|---|
| 10 | 0.174, 1.879 | 0.0253 | [0.0132, 0.1424] | 1.021 |
| 1280 | 0.471, 1.498 | 0.0346 | [0.0231, 0.0735] | 0.9988 |
| 10¹⁰ | 0.776, 1.236 | 0.0402 | [0.0325, 0.0518] | 0.9990 |

  - Smallest case {ρ₁, ρ̄₁}: λ = 0.987–1.013 at N = 10, B/|c̃|² = 1.009.
- **(d) Separation versus form** (s16, s18).
  - Over 112 blocks of 200 ζ zeros (N = 1280, T ≤ 2·10⁴): λ_min ∈ [0.010, 0.692], fitted as λ_min ≈ 1.14σ_min^{1.84} (corr 0.997; a two-point pair gives π²σ²/6). Box count n_max ≤ 2 per Clark cell.
  - Two coalescing zeros at ρ₁: B/|c̃|² = 0.97–1.00 for σ ≥ 1, then 1.63 (σ = 0.5), 1.89 (0.25), 1.98 (0.1) and 1.996 (0.001), while λ_min = 0.35, 9.4·10⁻², 1.6·10⁻² and 1.6·10⁻⁶ (N = 10; N = 1280 and 10¹⁰ alike).
  - Inserted near-duplicates (σ = 0.05) in the first 200 pairs:

| N | none | 30 % pairs | all pairs | 30 % triples |
|---|---|---|---|---|
| 1280 | 0.998 | 1.644 | 2.043 | 1.904 |
| 10¹⁰ | 0.999 | 1.607 | 1.988 | 2.213 |

    The values stay close to Σm²/Σm (1.46 / 2 / 2.1) while λ_min falls to 10⁻³–4·10⁻⁸.
  - At N = 10, doubling every zero pushes the density above the Clark density Φ′/2π, and the form jumps to 5.5. Real zeros never do this: their density is log(γ/2π)/2π < Φ′/2π.
- **(e) N8** (s17).
  - S4 holds to relative error ≤ 5·10⁻⁹: the truncated product over 22,491 zeros plus the density tail equals mpmath's 2ε|ξ′|/|ξ| at n = 1, 2, 5, 100 for ε = 0.4 and 0.1.
  - At fixed ε the shifted zeros are never H^∞-interpolating: the mean of δ·(T/2π)^ε is 5.0 (ε = 0.4) and 1.8 (ε = 0.1), with inf δ = 0.053 and 0.397 on the sample.
  - At ε = a/log(T/2π) (T ≈ 5600, 14200; median δ):

| a | ζ | CUE | Poisson |
|---|---|---|---|
| 0.5 | 0.950 | 0.949 | 0.85 |
| 1 | 0.820 | 0.812 | 0.58 |
| 2 | 0.494 | 0.490 | 0.21 |

    ζ's minimum (0.597 / 0.324 / 0.125) sits at CUE's 0.5 % quantile.

## 4. Adversarial
- S1's floor F depends only on the gamma factor; the zero form is explore38's RH-blind kernel geometry. For Davenport–Heilbronn, explore38 found B/Farey = 0.982 with off-line zeros entering as interior evaluators; the Beurling floor shows only at N ≈ 10³³⁴ (explore13).
- Riesz bounds (s16, blocks of 100 on-line DH zeros in [98, 356], density log(5γ/2π)):
  - N = 10: DH λ_min = 0.24 / 0.17 against ζ's 0.44 / 0.37 in the same windows.
  - N = 10¹⁰: DH 0.82 / 0.77 against ζ's 0.88 / 0.86.
  - The bounds follow density and gaps, not RH.
- Over 200-zero blocks at N = 1280, the median λ_min is 0.08–0.35 for ζ against 0.09–0.22 for CUE; Poisson gives ~10⁻⁷, and n_max is 2 for ζ and CUE alike.
  - The box law (S3) is a statement about multiplicities, which GUE and ζ share at the heights tested.
  - An off-line pair with (β−½)Φ′ ≪ 1 is a near-coincident on-line pair (explore41's Kadec stability), so it changes the form by at most the m² factor.
- None of S1–S4 separates ζ from DH or CUE: no RH information.

## 5. What it gives on N5 (TREE N5 line: Burnol's liminf; explore38's B(N); BCF's rate under RH + Σ1/|ζ′(ρ)|² ≪ T^{3/2−δ}; Balazard–de Roton d_N² ≪ (log log N)^{5/2+ε}/√log N under RH, quoted in BCF)
- **No new bound.** The rigorous lower bound improves explore38's B by F, which is ≥ ε: +1.3 % at N = 10, < 10⁻⁵ relative from N = 20 on. This is irrelevant to the node, since d_N is computed exactly for N ≤ 1280.
- **Correction to explore38.** Rigorously, D² ≥ F + c_Wᴴ G_W⁻¹ c_W ≥ explore38's B + ε. With synthesis in W_λ (Burnol's Theorem A, transfer unproved), equality holds: the continuous distance is then fully determined by the gamma factor and the zeros' evaluator geometry, and its second-order term is the zero form's (explore38's −0.1035).
- **Barrier for the Carleson route to an upper bound or rate.**
  - The Riesz lower bound degenerates with min gaps (λ_min ≈ 1.14σ^1.84), but the form does not: it is a box (multiplicity) quantity.
  - The natural tool for non-separated sets is divided-difference / generalized interpolation (Hruščev–Nikolski–Pavlov 1981, Avdonin–Ivanov). It needs, besides box counts, an A₂-type condition on the generating function.
  - At heights γ ≫ N the zero density tends to the Clark density, so that condition becomes Pavlov's A₂ for |ζ(½ + y/Φ′ + it)|², which explore41 found infinite on lines through zeros and growing near the line. This is negative-moment input of BCF's kind (not tested here beyond explore41).
  - Result: an RH rate for D(λ) this way needs the same arithmetic control as BCF. The integer-dilation penalty d_N² − D² (explore38's RH-sensitive part) is a distance inside C_λ and is not seen by any evaluator at the zeros.

## 6. Novelty
- Prior work:
  - Burnol math/0105120 ("the quotient [of a Sonine space by the ζ-subspace] contains vectors attached to the zeros").
  - Burnol math/0203120 Theorem A (completeness and minimality of the zero evaluators in L_a), read in TeX.
  - Burnol math/0103058 (λ → 0 only).
  - Burnol 1106.4751 (ζ(s)/(s−ρ)^k complete and minimal).
  - BCF 1211.5191 intro (Balazard–de Roton bound).
  - Noor 2203.05030 (N^⊥ via de Branges–Rovnyak).
- arXiv API, 13 queries (Nyman/Beurling × Carleson, interpolating, Toeplitz, model space, de Branges, prolate, distance asymptotic; Sonine × zeta; Carleson measure × Dirichlet series; au:Burnol): no source states F(λ), the finite-λ brackets, the box law, or S4.
- S1 is Burnol's theorem; S4 follows from Hadamard's product. F's prolate formula, the brackets and the box law are routine or numerical.
- Nothing for DISCOVERIES.
- Side note for the orchestrator: 2609.11560 (10 Sep 2026, NEW to the repo) claims the local embedding ∫_θ^{θ+1}|P(½+it)|^p dt ≤ C_p‖P‖^p_{ℋ^p} for every 2 < p < ∞. This would be a Carleson-embedding theorem for Dirichlet series that was conjectured false beyond even p. Candidate for an audit; relevant to NBLp-type criteria.

## 7. Caveats
- Double precision throughout; nothing interval-certified.
- Lower-bound values include explore38's diagonal tail beyond K, which is not rigorous. Without it they are rigorous for each finite K, up to rounding.
- F uses Nyström quadrature (M = 40 and 80 agree to 7 digits).
- The identification W_λ ↔ Burnol's L_a holds up to the finite-dimensional bookkeeping of (1−M)³ (W_λ/K_a has dimension 3, his L_a/K_a has 2). It is checked numerically by the brackets, not line by line.
- Reviewer (verdict only, claims inline):
  - Claim 1 (C_λ ⊂ W_λ, the decomposition, W^⊥) and claim 2 (H₃, F formula): CORRECT.
  - Claim 3 (exactness via Burnol's Theorem A): UNCLEAR. It needs the explicit C_λ ↔ P_a identification, multiplicity bookkeeping and the infinite-Gram limit.
  - Claim 4 (S4): CORRECT under the zero-location hypothesis.
  - Novelty: none demonstrated.
- S3 and §5's A₂ reduction are measured or heuristic, not proved.
- Memory, measured with the top MEM watchdog:
  - Kept runs ≤ 195 MB.
  - getrusage maxrss reached 369 MB on s12 at Bd = 100, where the watchdog showed 195 MB.
  - Two first versions were killed by the watchdog and then fixed: s12 at 475 MB (vectorized cotangent chunks of 2·10⁶ → 2.5·10⁵) and s18 at 315 MB (400 → 200 pairs).
  - One process at a time throughout.

## Suggested NOGO line
- explore42 (Orca, exploration: Carleson measures and interpolating sequences — Carleson interpolation, Carleson embedding/Shapiro–Shields, Nikolski distance formulae, model-space Carleson measures → N5, N8) NO IDEA: C_λ sits in a Sonine space W_λ, so D(λ)² ≥ F(λ) + (zero-evaluator form), with equality by Burnol's math/0203120 Thm A if it transfers (synthesis, unproved here; F a ζ-free prolate floor, F − ε = 8.1·10⁻³ at N = 2, which corrects explore38's B; brackets 0.2–0.7 % at N ≤ 10); Shapiro–Shields reads the form's Riesz bounds as min-gap functions (λ_min ≈ 1.14σ^1.84 in Clark cells, ζ ≈ CUE, DH alike), but the form follows a box law (coalescing pairs: ×2 while λ_min → 1.6σ²), so separation does not control the NB distance, and a divided-difference route needs Pavlov A₂ control of |ζ|^{±2} near the line, BCF's input; Nikolski localization from certified d_N is empty (needs log N > 5.6·10²³ above 3·10¹²); N8: Carleson's constant of ρ + ε is 2ε|ξ′(ρ)|/|ξ(ρ+2ε)| (checked to 5·10⁻⁹), a GUE gap statistic (median 0.820 vs CUE 0.812 at ε = 1/log T) [ExploreCarleson](Riemann/Deep/ExploreCarleson/REPORT.md)

RH relevance: N5 (structure and rate of the continuous Nyman–Beurling distance; no bound, no RH information); N8 (none).

## Files
- code/ (one at a time, OMP etc. = 1; watchdog top-MEM peaks ≤ 195 MB; getrusage maxrss ≤ 260 MB except s12 at Bd = 100, 369 MB):
  - s0_cache.py: zeros cache, from explore38.
  - s3_scan.py: explore38's B(N).
  - s11_eps.py: ε.
  - s12_farey.py: lean Farey upper bound, Bd ≤ 100.
  - s13_fit.py: extrapolation.
  - s14_floor.py: F(λ).
  - s15_joint.py: lower bound F ⊕ zero evaluators.
  - s16_riesz.py: Riesz bounds, separation, CUE, Poisson, DH.
  - s17_n8.py: S4 and δ statistics.
  - s18_cluster.py: box law.
- logs/: outputs of each script, named after it.
