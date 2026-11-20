IDEA: Integer dilations span Bagchi's piecewise space W1, so they are never a Marcinkiewicz–Zygmund set for the continuous Nyman–Beurling family, but the defect is free: the integer-dilation penalty splits exactly as roughness R + in-W1 discretization Q with d_N² ≤ (1 + Q/R)·D(1/N)², and Q obeys Burnol's zero formula Q·log N = Σ_ρ|𝓜[roughness](ρ)|² (1–3 %), so Q/R ≈ 1.5/log N and, under RH, BBLS's discrete limit d_N² log N → C reduces to Q = o(R) (node N5).

# explore61 — sampling discretization → N5 integer-dilation penalty (Orca exploration, 20 Nov 2026, 11:44–13:00 UTC)

Proofs: [DERIVATION.md](DERIVATION.md). Notes as found: [NOTES.md](NOTES.md). Scripts here, outputs in [logs/](logs/).
Notation: H = L²(0,∞), ρ_θ = {θ/x}, χ = 1_(0,1), V_N = span{ρ_{1/k}: k ≤ N}, C_λ = closed span{ρ_θ: θ ∈ [λ,1]}, λ = 1/N, L = log N,
d_N = dist(χ, V_N), D = D(1/N) = dist(χ, C_λ), P = d_N² − D² (explore38's integer-dilation penalty), C = 2 + γ − log 4π = .0461914.

## 1. Transfer
- Synonym check (rg Marcinkiewicz|Zygmund|Tchakaloff|sparsif|Batson|Spielman|Srivastava|Kadison|Weaver|discretiz|Temlyakov|Kashin over NOGO,
  TREE, Deep/*/REPORT): no N5 use (explore46's "sparse domination" is N7; MSS appears only as Ihara-zeta background).
- Tool: sampling discretization with **prescribed nodes** (Marcinkiewicz–Zygmund inequality ⇔ optimal quadrature with fixed nodes;
  Schaback's power-function bound). BSS/MSS/Weaver choose the nodes and so do not apply to the fixed set {1/k}; NOU frames likewise.
- Exact dictionary: P = dist(h*, V_N)² is the squared worst-case error of the optimal quadrature with nodes u_k = log k for the
  continuous optimizer's coefficient measure, in the RKHS of K(u,u′) = ⟨ρ_{e^{−u}}, ρ_{e^{−u′}}⟩ (h* = P_C χ). The MZ/power-function
  bound is P ≤ (∫π_N d|ν*|)², π_N(θ) = dist(ρ_θ, V_N) (vanishes at the nodes θ = 1/k).
- Node: N5, the integer-dilation penalty (TREE N5: "the remaining step is the integer-dilation penalty P(N) = o(1/log N)";
  dbarnb's THEOREM: RH ⇒ D(λ)² log(1/λ) → Σm²/|ρ|²).

## 2. Statements
- **S1 (sampling defect; elementary, certified).** Every ρ_{1/k} lies in W1 = {a_m + b/x on I_m = (1/(m+1), 1/m), b/x on (1,∞), one b}
  (Bagchi's piecewise space, math/0607733; RH ⟺ closure ∪V_N = W1). A non-integer dilation jumps inside an I_m, so
  dist(ρ_θ, V_N) ≥ ω(θ) := dist(ρ_θ, W1) > 0 for every N: the integer nodes are never an MZ/sampling set for C_λ, and the
  full-space power function never decays. Arb (s9, 3000 intervals): dist(ρ_θ, V_N)² ≥ .1325637 (θ = 2/3), .1754808 (3/4),
  .0487619 (2/5), .1133960 (3/5), .0151678 (2/9) for all N (‖ρ_{2/3}‖² = .8404).
- **S2 (exact split + Lemma B; new as far as found).** P = R + Q, R = dist(h*, W1)² (roughness of the continuous optimizer),
  Q = dist(P_W1 h*, V_N)² (its discretization inside W1); R² ≤ (D² − R)·Q, hence R·P ≤ D²·Q and, if R > 0,
  **d_N² − D² ≤ (Q/R)·D(1/N)²**. Unconditional; the same for every target in W1 and every closed C ⊇ V_N.
- **S3 (Theorem C).** d_N ≤ ‖χ − h‖ + ∫π^W_N(θ) d|ν|(θ) for every candidate h = Σc_kρ_{1/k} + ∫ρ_θ dν, with the W1-power function
  π^W_N(θ) = dist(P_W1 ρ_θ, V_N) = (π_N² − ω²)^{1/2}, which → 0 under RH (unlike π_N ≥ ω).
- **S4 (zero formula; numerical).** L·dist(P_W1 ρ_θ, V_N)² = Σ_ρ|𝓜[P_W1^⊥ρ_θ](ρ)|² (±6 %, 14 θ, N ≤ 2000), and along the continuous
  optimizer L·Q = Σ_ρ|𝓜[P_W1^⊥h*](ρ)|² (±3 %): the Burnol-type zero sum of the roughness (for the target χ it is Burnol's C) is attained.
- **Conditional (RH).** If Q(N) ≤ ε_N·R(N) with ε_N → 0 (so Q = 0 wherever R = 0) for the continuous optimizers h*_N, then
  P ≤ ε_N·D² (S2, DERIVATION 3(d)) and, with dbarnb, d_N² log N → Σm²/|ρ|² (BBLS's discrete conjecture). Sufficient: for some fixed K,
  L·Q ≤ K·Σ_ρ|𝓜g_N(ρ)|² (S4 as an upper bound) and Σ_ρ|𝓜g_N(ρ)|² ≤ Z_N·‖g_N‖² with Z_N = o(L) (an MZ upper inequality at the zeros),
  g_N = P_W1^⊥h*_N, both holding UNIFORMLY along the moving targets P_W1 h*_N (not only for each fixed target); then ε_N = K·Z_N/L, and
  limsup L²(d_N² − D²) ≤ K·Z·C if Z_N ≤ Z. Not proved: S4 is numerical; the conjunction of the two estimates is BBLS-strength and neither
  is shown easier (referee).

## 3. Tests (ζ; double precision unless certified; one process at a time, peaks ≤ 220 MB RSS, memcap ≤ 161 MB)
### 3a. Smallest case: one non-integer dilation (s3_rate, s5_zero, s10_cmp; N ≤ 2000 from one Cholesky of the Vasyunin Gram)
κ(θ) := L·(π_N² − ω²) = L·dist(P_W1 ρ_θ, V_N)²; zero sum = Σ_ρ|𝓜[P_W1^⊥ρ_θ](ρ)|² over 3000 zeros + tail (tail ≤ 2 %).

| θ | ω² | ‖ρ_θ‖² | κ, N = 160 | 640 | 2000 | zero sum | κ(2000)/zero sum |
|---|---|---|---|---|---|---|---|
| 2/3 | .13259 | .8404 | .00348 | .00332 | .00339 | .00350 | .970 |
| 3/4 | .17551 | .9455 | .00285 | .00275 | .00279 | .00283 | .986 |
| 5/6 | .19171 | 1.0506 | .00168 | .00161 | .00165 | .00165 | .997 |
| 3/5 | .11342 | .7564 | .01987 | .01903 | .01938 | .02003 | .968 |
| 5/7 | .16856 | .9005 | .00857 | .00852 | .00853 | .00839 | 1.016 |
| 2/5 | .04878 | .5043 | .02677 | .02578 | .02628 | .02714 | .968 |
| 3/8 | .04518 | .4727 | .02916 | .02826 | .02821 | .02932 | .962 |
| 2/7 | .02503 | .3602 | .02249 | .02328 | .02287 | .02228 | 1.027 |
| 2/9 | .01518 | .2801 | .01863 | .01877 | .01828 | .01729 | 1.057 |

(all 14 θ in logs/s10_cmp.log: ratios .959–1.057.) So π_N(θ)² = ω(θ)² + κ(θ)/log N: the full-space power function stalls at
the W1 floor (π_N²/ω² − 1 = 3.4·10⁻³ at θ = 2/3, N = 2000), while the in-W1 part decays exactly at the Nyman–Beurling rate,
with κ(θ) = Burnol's zero sum of the roughness (κ is flat in N to ±3 %, like L·d_N² around C: .0469, .0474, .0453, .0465 at N = 40–2000).

### 3b. At scale: head-hybrid continuous optimizer (s1_split; all θ = Q/k ∈ (1/5, 1) added to the N integer dilations)
P, R, Q for the hybrid span (a lower bound for the true penalty; explore31: 90 % of it sits in θ ∈ [1/5, 1]). S − Ω ⪰ 0 (W1 theory)
holds in every run (λ_min ≥ −4·10⁻⁹ against λ_max ≈ .1). "L·Q/zero sum" = s6_optzero (300 zeros + tail).

| N | grid Q | L²P | L²R | R/P | L³Q | c_e = LQ/R | cos | P/D² ≤ Q/R (Lemma B) | L·Q / zero sum | L·λmax(S−Ω, Ω) | Q_pw/Q | F/P | L‖e‖₁ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 40 | 8 | 0.0439 | 0.0318 | 0.723 | 0.0449 | 1.413 | 0.923 | 0.3404 ≤ 0.3830 | nan | 2.25 | 4.17 | 12.2 | 3.15 |
| 160 | 8 | 0.0480 | 0.0377 | 0.784 | 0.0527 | 1.399 | 0.939 | 0.2493 ≤ 0.2756 | 1.027 | 2.04 | 3.91 | 12.1 | 3.28 |
| 160 | 32 | 0.0576 | 0.0432 | 0.749 | 0.0734 | 1.701 | 0.960 | 0.3149 ≤ 0.3351 | nan | 2.67 | 5.24 | 19.6 | 4.89 |
| 640 | 8 | 0.0469 | 0.0387 | 0.824 | 0.0534 | 1.379 | 0.935 | 0.1908 ≤ 0.2135 | 1.000 | 2.03 | 3.94 | 12.6 | 3.35 |
| 640 | 16 | 0.0535 | 0.0431 | 0.805 | 0.0676 | 1.569 | 0.950 | 0.2235 ≤ 0.2429 | nan | 2.38 | 4.36 | 16.2 | 4.14 |
| 640 | 32 | 0.0569 | 0.0451 | 0.793 | 0.0761 | 1.688 | 0.951 | 0.2409 ≤ 0.2612 | nan | 2.67 | 5.30 | 21.5 | 5.06 |
| 2000 | 8 | 0.0501 | 0.0424 | 0.845 | 0.0589 | 1.390 | 0.943 | 0.1653 ≤ 0.1829 | 1.012 | 1.97 | 3.84 | 12.5 | 3.46 |
| 2000 | 16 | 0.0573 | 0.0474 | 0.828 | 0.0750 | 1.583 | 0.957 | 0.1935 ≤ 0.2082 | 1.015 | 2.31 | 4.27 | 16.3 | 4.30 |

- R/P = .72–.85: three quarters of the penalty is roughness (h*'s jumps inside harmonic intervals), rising with N, falling slowly
  with the grid. Lemma B holds in every run (P/D² ≤ Q/R, slack 6–12 %; cos = R/√((D² − R)Q) = .92–.96).
- c_e = L·Q/R is flat in N at each grid (1.38–1.41 at Q = 8; 1.57–1.58 at 16; 1.69–1.70 at 32) and equals the zero visibility
  Z = Σ_ρ|𝓜g(ρ)|²/‖g‖² of the roughness g (s6: 1.36–1.56): Q ≈ Z·R/log N, i.e. Q/R → 0 like 1/log N.
- Tail freedom (s13_tail: one non-integer dilation 2/(2k+1) in every gap down to 1/301, on top of the head grid Q = 8):
  N = 640: L²P .0469 → .0490, R/P .824 → .805, c_e 1.38 → 1.56, cos .935 → .893; N = 2000: L²P .0523, R/P .830, c_e 1.563 (flat in N),
  P/D² = .174 ≤ Q/R = .206. The tail carries 81–95 % of ‖e‖₁ but adds only 4 % to P.
- The profile L·e (head coefficients) is N-stable (corr .993, 640 vs 2000): the non-integer part of the optimizer is O(1/L).
- Constant: R = cos²·c_e·(D² − R)/L is an identity (definitions of cos and c_e); its content is c_e = Z (zero formula) with c_e, cos
  N-flat. Then L²P → cos²·Z·C: .067 (head, Q = 16), .071 (head, Q = 32), .059 (head Q = 8 + tail), against today's .050–.057
  (still rising with N, as D²L does) and ceilaudit6's δ = L²(d_N² − B) ≈ .068 at N = 4000.

### 3c. Theorem C against the plain MZ power-function bound (s11_thmC; h = hybrid optimizer)
| N, grid | L²P | in-W1: Q_pw/Q | bound L²(2D√Q_pw + Q_pw) | full space: F/P | bound L²(2D√F + F) |
|---|---|---|---|---|---|
| 640, 8 | .0469 | 3.94 | .212 (4.5×) | 12.6 | 1.35 (29×) |
| 2000, 8 | .0501 | 3.84 | .220 (4.4×) | 12.5 | 1.50 (30×) |
| 2000, 16 | .0573 | 4.27 | .266 (4.6×) | 16.3 | 1.99 (35×) |
| 640, 32 | .0569 | 5.30 | .305 (5.4×) | 21.5 | 2.30 (40×) |
Q_pw = (Σ|e_j|π^W_N(θ_j))², F = (Σ|e_j|π_N(θ_j))². The in-W1 bound is second order and flat in N (L³Q_pw = .21–.40); the
full-space bound is ≍ L^{−3/2} (its L² value grows like √L): the MZ/quadrature tool works only after projecting onto Bagchi's space.

### 3d. Is the zero visibility an MZ average? (s7_mz; N = 640, Q = 8, roughness Mellin transform on t ≤ 1000, step .1)
Σ_{γ≤T}|G(1/2+iγ)|² against the smooth-density average ∫_0^T|G|²θ′(t)/π dt: ratio 1.21 (T = 50), .99 (100), .92 (200), .85 (500),
.82 (1000); 65 % of ‖g‖² lies below t = 50. So Z = 1.33 vs the MZ average 1.62: Z is the mean of log(t/2π) over the roughness
spectrum, discounted ≈ 18 % because |𝓜g| dips at the zeros (𝓜h*(ρ) = 0). An MZ upper inequality at the zeros for this class is the
missing analytic step behind "Q ≤ (Z/L)·R".

## 4. Adversarial: Davenport–Heilbronn (s8_dh; natural Gram by ceilaudit6's gram.c 'per' mode, cross entries by explore38's digamma
formula; 304 on-line zeros < 420 and 9 off-line, ExploreHelio)
- The W1 structure holds verbatim (every A(1/(kx)) is constant on each I_m: V_N^DH ⊂ PC), so S1–S3 are not about ζ.
- Zero formula for DH, L(π² − ω²) at N = 640/1000 against the on-line zero sum: 2/3 .201/.193 vs .215; 3/4 .219/.217 vs .250;
  2/5 .097/.101 vs .099; 3/7 .101/.105 vs .106; 3/8 .084/.086 vs .086; 2/7 .060/.060 vs .057; 2/9 .040/.041 vs .038 (ratios .87–1.09).
- DH's off-line zeros add the H² evaluation floor (2β − 1)|𝓜g(ρ₀)|², slope 1·10⁻⁵–1.2·10⁻⁴ per unit L: ≤ .4 % at N = 1000 and 10 % only near
  N ≈ e^{100}, so the mechanism cannot tell ζ from DH at any computable N. RH-blind, like explore38's kernel part. Its only RH-sensitive term is this
  floor (Beurling's), which (heuristically) keeps Q/R ↛ 0 for DH as L → ∞.

## 5. What it gives on N5 (TREE N5: under RH the remaining step to BBLS's discrete limit is P(N) = o(1/log N))
- New exact inequality (unconditional): R·(d_N² − D²) ≤ D²·Q, i.e. d_N² ≤ (1 + Q/R)·D(1/N)² when R > 0, plus the split P = R + Q. Before: only d_N ≥ D(1/N).
- New mechanism (numerical): the penalty is 75–85 % roughness; the rest, Q, is Burnol's zero sum of the roughness divided by
  log N; hence P/D² ≈ cos²·Z/log N and the flat explore38/explore31 constant .06–.068 ≈ cos²·Z·C with Z ≈ 1.4–1.7, cos ≈ .90–.96.
  This replaces the unexplained "integer-dilation penalty" by two measured, N-flat constants of the continuous optimizer.
- Reduction (conditional on RH): BBLS's discrete conjecture ⟸ Q(N) ≤ ε_N R(N), ε_N → 0 ⟸ [S4 as an upper bound, fixed K] +
  [Z_N = o(log N), an MZ inequality at the zeros for the roughness], both uniform along the moving targets P_W1 h*_N. Measured
  Q/R = .18–.38 ≈ 1.5/L. Not a proof; the conjunction is BBLS-strength (referee).
- Barrier (certified, for the naive transfer): any MZ/frame/BSS-type discretization of the continuous optimizer onto the integers
  in the full space pays the N-independent floor ω(θ) (π_N ≥ ω > 0), so it bounds P only by (non-integer mass)² × O(1): here 12–40×
  P and ≍ L^{−3/2}. The method's full freedom (any N, any coefficients) is covered by S1, since V_N ⊂ W1 for all N.
- Side consequence of S3: d_N² log N → C under RH already follows if some continuous near-optimizers have non-integer total
  variation o(L^{−1/2}) (head grid alone: L·‖e‖₁ = 3.2–5.1; but tail midpoints carry 81–95 % of ‖e‖₁ in s13, so only the
  π-weighted form ∫π_N d|ν| = o(L^{−1/2}) is plausible; Lemma B, being scale-free in e, avoids the issue).
- Rigorous range: S1 (Arb), S2, S3 are theorems; S4, Q/R ≈ Z/L, and the constant formula are numerical (N ≤ 2000, 14 θ, double).

## 6. Novelty (searched 20 Nov 2026)
- Bagchi, math/0607733 (TeX read): Hilbert space of sequences with weight 1/(n(n+1)), γ_l = ({n/l})_n, the piecewise-constant space
  𝓜 on (1/(n+1), 1/n] and g_l = {1/(lx)} − {1/x}/l ∈ 𝓜; RH ⟺ {γ_l} total. So W1 and "RH ⟺ V_∞ = W1" are known (S1 is routine).
- arXiv API: Nyman∧Beurling∧piecewise, Bagchi∧Nyman, Nyman∧Beurling∧dilations∧continuous, Marcinkiewicz∧Zygmund∧Nyman,
  Báez∧Duarte (12 hits), Nyman∧Beurling∧Mellin (4 hits; 2510.18132 abstract: Gram decay on a 2^j3^k ladder, unrelated).
- Monthly listings Oct–Nov 2026 (math.NT, CA, CV, FA; Nyman|Beurling|Báez|Bagchi|Burnol|dilation|Marcinkiewicz|sampling discreti|
  sparsif): 17 title hits, none on NB (2611.10402 and 2610.07679 are sampling-discretization papers without ζ).
- Zenodo ("Nyman-Beurling" OR "Báez-Duarte"): 74 records; abstracts of 24382036 and 24054105 read; our explore38/explore31 lists.
- GXL evidence API (2 claims: second-order discrete–continuous gap; piecewise a + b/x structure): 0 supports, 0 refutes.
- New as far as found: S2 (Lemma B and the split), S3's W1 power function, S4 (zero formula for W1 targets and along the optimizer),
  and the decomposition of the penalty constant. S2/S3 are two-line Hilbert-space facts once W1 is in view (likely ROUTINE as math).

## 7. Caveats
- The hybrid optimizers (head grid θ = Q/k ∈ (1/5, 1), Q ≤ 32, plus tail midpoints) are not BBLS's continuous optimizer: P_hyb ≤ P, and
  c_e rises with refinement (1.38 → 1.70 with the grid, +.18 with the tail); continuum values of Z, cos and the constant formula are
  extrapolations. N ≤ 2000 (L ≤ 7.6). Double precision; Ω by Richardson in the interval cut-off (relative change ≤ 1·10⁻⁴).
- Zero sums: 3000 zeros + a smooth tail model (≤ 2 %) for single θ, 300 zeros + tail along the optimizer; DH only 304 zeros < 420.
- S2/S3 reviewed (reviewer, verdict only, claims inline): (1) W1 facts, (2) defect, (3) split, (5) Theorem C CORRECT; (4) Lemma B GAP
  (division by R; fixed: R·P ≤ D²·Q always). Reviewer's novelty view: "elementary subspace/projection consequences".
- The reduction moves the open step, it does not close it: Q = o(R) is a property of the unknown continuous optimizer, S4 is numerical,
  and an MZ upper bound at the zeros for the roughness (Z = o(L)) is not proved. ‖e‖₁ grows with tail freedom, so S3's total-variation
  corollary is fragile; Lemma B (scale-free in e) is not.
- Independent referee (20 Nov 2026, 12:55 UTC, verdict only): overall CORRECT, ROUTINE.
  (1) DERIVATION §§1, 3, 4 CORRECT; §1's sum ℂρ_1 + PC is not orthogonal (notation fixed); §4's integral needs Bochner integrability
  (stated: ν finite, θ ↦ ρ_θ continuous). (2) The reduction is CORRECT; since R = 0 can occur, it is stated as Q ≤ ε_N R, ε_N → 0
  (Q = 0 where R = 0), fixed in DERIVATION 3(d) and §2. (3) Not circular as a conditional argument: any fixed K in
  L·Q ≤ K·Σ_ρ|𝓜g_N(ρ)|² ≤ K·Z_N·R with Z_N = o(L) suffices, but the missing upper bound must hold uniformly along the moving targets
  P_W1 h*_N, not just for each fixed target (S4's 14 single-θ tests are fixed targets; the optimizer tests are N ≤ 2000 only); the
  conjunction of the two estimates is BBLS-strength and neither is shown easier. (4) ROUTINE: elementary projections relocate the
  obstruction; Q = o(R) is sufficient, not shown necessary or easier than P = o(1/L).
- RH-blind at every computable N (§4).

## Suggested NOGO line
- explore61 (Orca, exploration: sampling discretization and sparsification — Marcinkiewicz–Zygmund / prescribed-node quadrature → N5 integer-dilation penalty; BSS/MSS/Weaver/NOU choose their nodes and do not apply to {1/k}) IDEA, not a proof: integer dilations span Bagchi's piecewise space W1 (a_m + b/x on (1/(m+1), 1/m), one b; math/0607733), so they are never an MZ set for the continuous family (dist(ρ_θ, V_N)² ≥ dist(ρ_θ, W1)² > 0 for every N; Arb: ≥ .13256 at θ = 2/3) and full-space power-function bounds are 12–40× the penalty; the defect is free: exactly P = d_N² − D(1/N)² = R + Q and R·P ≤ D²·Q (R = roughness of the continuous optimizer, 75–85 % of P; Q = its in-W1 discretization), and Q obeys Burnol's zero formula L·Q = Σ_ρ|𝓜[roughness](ρ)|² (±3 % along the optimizer, ±6 % for 14 single dilations), so Q/R ≈ 1.4–1.7/L and L²P → cos²·Z·C ≈ .06–.07 (Z the roughness's zero visibility); under RH, BBLS's discrete limit ⟸ Q ≤ ε_N R with ε_N → 0 (unproved; sufficient: a fixed-K zero-sum upper bound plus Z_N = o(L), uniformly along the moving targets P_W1 h*_N, a conjunction of BBLS strength); RH-blind (Davenport–Heilbronn alike, its off-line floor ≤ .4 % at N = 10³); reviewer: CORRECT after an R > 0 fix; referee: CORRECT, ROUTINE (elementary projections relocate the obstruction; Q = o(R) sufficient, not shown necessary or easier) [ExploreSparse](Riemann/Deep/ExploreSparse/REPORT.md)

RH relevance: N5 (rate of the Báez-Duarte/Nyman–Beurling distance: the integer-dilation penalty d_N² − D(1/N)², the last step from dbarnb's continuous THEOREM to BBLS's discrete limit under RH).

## Files (run from this folder with `OMP_NUM_THREADS=1 … uv run --no-project --with numpy --with scipy python3 <script>`, under memcap)
- gram.c — copy of ceilaudit6's Gram code. Regenerate the intermediates (deleted after the run):
  `cc -O2 -o ../../../.scratch/explore61/gram gram.c -lm`; `…/gram zeta 2000 ../../../.scratch/explore61/Gz2000.bin` (3 s, 34 MB);
  `…/gram per 5 1000 ../../../.scratch/explore61/Gdh1000.bin 0 1 1.28407904384041227 1 0` (8 s, 30 MB).
- s1_split.py N Q K0 [M] — hybrid optimizer, P, R = e'Ωe, Q, S − Ω ⪰ 0 check; Vasyunin code copied from explore31's s3_hybrid.py (≤ 130 MB). run1.sh — driver.
- s3_rate.py — π_N(θ)², ω(θ)² for every N ≤ 2000 (114 MB RSS). s10_cmp.py — κ(θ) vs zero sums.
- s4_pencil.py, s12_table.py, s11_thmC.py — derived constants, tables, Theorem C bounds (from logs/s1_*.npz).
- s5_zero.py NZ M — zero sums for single θ (3000 zeros: 191 MB RSS, 206 s). s6_optzero.py — zero formula along the optimizer (≤ 112 MB).
- s7_mz.py — zero sum vs MZ average (64 MB). s8_dh.py — Davenport–Heilbronn (98 MB). s9_cert.py — Arb certificate (python-flint, 35 MB).
- s13_tail.py N Q K0 K1 — head + tail freedom (≤ 203 MB at N = 2000, 324 extra dilations).
- logs/: every output (s1_scan.log with /usr/bin/time and memcap peaks; table_split.md; *.npz with S, Ω, e for each run).
- Zeros: Riemann/Deep/ExploreHelio/results/zeta_zeros.npy (22491 ordinates < 2·10⁴) and dh_zeros.json (read only).
