# REPAIR — explore59 / ExploreDbar, answer to the referee's GAP verdict (20 Nov 2026, ~09:40 UTC)

Status: (2) FIXED · (3i) FIXED mod Lemma C · (3ii) FIXED · (3iii) FIXED mod Lemma C · (3iv) FIXED (moot) · (3v) FIXED mod Lemma C ·
(4) OPEN (Lemma R′). Lemma C is a standard contour lemma (Titchmarsh 9.6–9.7, Cartan), stated in §C. Net: RH + Lemma R′ ⇒ the limit.
The ∂̄ step is replaced everywhere by explicit Mittag-Leffler corrections with a band-limited kernel; ∂̄ remains the heuristic.
Notation as in DERIVATION: z = t − i(σ − ½); V_L = {E = ∫_[0,L] e^{−iuz}dμ(u)}, ‖E‖² = ∫∫e^{−|u−v|/2}dμdμ̄ = ∫|E|²w on the line;
K = reproducing kernel; U_Z(L) = min{‖E‖² : E ∈ V_L, E = 1 on Z with multiplicity}; L′ = L − √L; T₁ = L², T₂ = L⁶.

## (2) m² is a coalescent limit — FIXED (restatement + proof of what the limit needs)
Claim F. For a fixed finite set Z_H = {ρ : |γ| ≤ H} (true multiplicities m_ρ), L·U_{Z_H}(L) → Σ_{Z_H} m_ρ²/|ρ|² as L → ∞.
Proof. (a) One point, multiplicity m: the constraints are E(γ) = 1, E^{(j)}(γ) = 0 (1 ≤ j < m), with E^{(j)}(γ) = ∫(−iu)^j e^{−iuγ}dμ.
K(t,s) = (¼ + s²)∫_0^L e^{−iu(t−s)}du + (½ − is) + (½ + is)e^{−iL(t−s)}. The main term gives the jet Gram (¼+γ²)·Δ M Δ*,
M_jk = ∫_0^L u^{j+k}du, Δ unitary diagonal; the data (1,0,…,0) give cost [M⁻¹]₀₀/(¼+γ²) = Σ_{k<m} q_k(0)²/(¼+γ²) with
q_k = √((2k+1)/L)·P_k(2u/L − 1) orthonormal on [0,L]; q_k(0)² = (2k+1)/L (P_k(−1)² = 1), so the cost is m²/(L|ρ|²). The
endpoint terms change entry (j,k) by O(|γ|L^{j+k}) against (¼+γ²)L^{j+k+1}/(j+k+1): relative O(1/(|γ|L)) → 0.
(b) Distinct points at fixed distance d: in the jet basis scaled by L^{−j}, cross entries are O((¼+γ²)L·(1/(Ld))) against diagonal
≍ (¼+γ²)L, i.e. normalised O(1/(Ld)) → 0; so the finite normalised Gram tends to the block-diagonal limit of (a). ∎
So m² is exact in the limit for every fixed zero set; at finite L the constant depends on the rescaled configuration (referee's
point, accepted), and that dependence is a domination question, moved into Lemma R′ (4). Split-pair/cluster numerics
(dbar_pair, dbar_cluster) illustrate both regimes.

## (3i) Window, s = 1, growth for Paley–Wiener–Schwartz — FIXED modulo Lemma C
- One-sided window. Φ(z) := e^{−iz/T₁}F(z/T₁), F = h²/h(0)², h(x) = ∫φ(ξ)e^{ixξ}dξ, φ ≥ 0 an even Gevrey-2 bump on [−½,½].
  Spectrum of Φ ⊂ [0, 2/T₁] (no negative frequencies); on the line |Φ| = F(t/T₁) ≤ 1; |Φ(s)| ≤ e^{(2/T₁)max(0, ½−σ)}·C;
  in |σ − ½| ≤ 1, |Φ| ≤ C e^{−c(|t|/T₁)^{1/2}} (Gevrey decay is uniform in horizontal strips). Cost: |1 − Φ(t)| ≤ |t|/T₁ +
  C t²/T₁² (|t| ≤ T₁), ≤ 2 beyond, so ‖1 − Φ‖² ≪ 1/T₁ = L⁻².
- s = 1. D(λ) is Burnol's/BBLS's unconstrained distance in L²(0,∞) (span of the contractions of {1/t}, math/0112254 §5);
  their Mellin images −θ^sζ(s)/s lie in L²(Re s = ½) and the pole at 1 is off the line, so no condition at s = 1 is needed
  (§7's Bochner-integral argument uses none). For Nyman's constrained variant (Σc_kθ_k = 0) use G − G(1)·1 (θ = 1 term):
  G₀(1) = 0 because 1/ζ(1) = 0, so |G(1)| = |R(1)| ≪ e^{−5√L/2} by (3ii), (3iii); extra cost |G(1)|²‖ζ/s‖² = o(L⁻²).
- Growth. Need |G(s)| ≤ C_L(1+|s|)^{N}e^{(L−½)max(0,−σ)} for all s (then PWS gives supp ⊂ [0, L − ½] and §7 applies with δ < ½).
  Pieces: |Ψ(s)| ≤ 1 + ‖E_low‖·K_{L′}(z,z)^{1/2} (reproducing property; ‖E_low‖ ≤ ‖1‖ = 1) with, for z = x + iy,
  K(z,z) = (¼ + z̄²)(e^{2L′y} − 1)/(2y) + (½ − iz̄) + (½ + iz̄)e^{2L′y}: ≤ (¼+|z|²)min(L′, 1/(2|y|)) + O(1+|z|) for y < 0, and
  e^{2L′y}·poly(|z|) for y > 0; so |Ψ| ≤ C(1+|s|)√L·e^{L′max(0,½−σ)}. |1/ζ|: ≤ 2 for σ ≥ 2; under RH ≤ exp(C log t/log log t) for
  |σ − ½| ≥ 1/log log t, 0 ≤ σ ≤ 2 (Titchmarsh 14.14, and the functional equation for σ < ½); ≤ C log(|t|+2) for σ ≤ 0, |t| ≥ 1;
  circles |s + 2k| = ½ near the trivial zeros. In the thin strip |σ − ½| ≤ 1/log log t: maximum principle for the entire G on
  rectangles whose vertical sides lie in the previous region and whose horizontal sides are Lemma C's heights t_n
  (|1/ζ| ≤ t^A there). R_tr, R_high are bounded on all these contours by (3ii), (3iii). Exponents: L′ + 2/T₁ ≤ L − ½ for large L,
  kernels have type L − 1. ∎ (mod Lemma C)

## (3ii) Trivial-zero corrections — FIXED
- Terms: R_tr(s) = Σ_{k≥1} r_k e^{−(L−1)(s+2k)}/(s+2k), r_k = Ψ(−2k)Φ(−2k)/ζ′(−2k) (β ≡ 1). Referee accepted: a term with a pole
  has no compact spectrum; none is claimed. Only the entire G needs a spectrum, and it gets it from its growth ((3i), PWS).
- Prefactor without Lemma R: by the reproducing property, |E_low(−2k)| ≤ ‖E_low‖·K_{L′}(iy,iy)^{1/2}, y = 2k + ½, and
  K(iy,iy) = e^{2L′y}(y+½)²/(2y) + O(y) (the formula above with z̄ = −iy). So |Ψ(−2k)| ≤ 1 + C√k·e^{L′(2k+½)}; with
  |Φ(−2k)| ≤ Ce^{(2/T₁)(2k+½)} and |ζ′(−2k)| = (2k)!ζ(2k+1)/(2(2π)^{2k}):
  on the line |r_k e^{−(L−1)(s+2k)}/(s+2k)| ≤ C√k e^{−(√L − 1 − 2/T₁)(2k+½)} / (|ζ′(−2k)|·|2k + ½ + it|).
  Hence ‖ζR_tr‖²_w ≪ e^{−5(√L−2)} (|ζ(½+it)| ≪ t^{1/6}, w ≍ t⁻²), |R_tr(1)| ≪ e^{−5√L/2}; εL = √L beats the √k·(y+½)/√y prefactors.
- Growth: each term ≤ |r_k|e^{(L−1)max(0,−σ)}e^{−2k(L−1)}/|s+2k|, |r_k| ≤ C√k e^{(L′+2/T₁)(2k+½)}/|ζ′(−2k)|; the sum converges
  locally uniformly off {−2k} (factorial |ζ′(−2k)|) and has type L − 1 to the left, bounded for σ ≥ 0. ∎

## (3iii) High-zero corrections — FIXED mod Lemma C;  (3iv) — FIXED (moot)
- Kernel: B(x) = ((1 − e^{−ax})/(ax))², a = (L−1)/2: entire, B(0) = 1, B = Fourier–Laplace transform of the triangle on [0, L−1];
  |B(x)| ≤ min(1, 4/(a|x|)²)·e^{(L−1)max(0, −Re x)}.
- For each Lemma C curve c (around a zero cluster P_c above T₂; concentric c′ of double size): with F = ΨΦ/ζ,
  R_c(s) := Σ_{p∈P_c} Res_{w=p}[F(w)B(s−w)/(s−w)]. Since B(x)/x = 1/x + B′(0) + O(x), R_c is meromorphic with poles P_c and the
  principal parts of F (also at multiple poles). Residue theorem: R_c(s) = I_{c′}(s) := (1/2πi)∮_{c′}F(w)B(s−w)/(s−w)dw for s
  outside c′, and R_c = I_{c′} + F inside c′; so G₀ − R_c = −I_{c′} is holomorphic inside c′ and, on the line,
  ζR_c = ζI_{c′} + ΨΦ·1_{inside c′}.
- Bounds: on c′, |1/ζ| ≤ exp(C log γ log log γ) (Lemma C), |Ψ| ≤ C(1+γ)√L e^{L′r}, |Φ| ≤ Ce^{−c(γ/L²)^{1/2}} ≤ Ce^{−cL²} (γ ≥ L⁶);
  len c′ ≤ 1, r ≤ 1/log γ. So |ζR_c(½+it)| ≤ C t^{1/6}e^{−c(γ_c/L²)^{1/2} + C log γ_c log log γ_c + L}·min(1, (L|t−γ_c|)^{−2})·(1 + 1/dist);
  summing over c (≪ log γ curves per unit height): ‖ζR_high‖²_w ≪ e^{−cL²}, |R_high(1)| ≪ e^{−cL²}; growth ≤ (tiny)·e^{(L−1)max(0,½−σ+r)}.
- (3iv) No ∂̄ solution is used any more: every correction is holomorphic off its poles and is bounded on the line pointwise by the
  Cauchy integral above, so no trace estimate is needed. ∎ (mod Lemma C)

## (3v) The displayed bound — FIXED mod Lemma C
On the line 1 − ζG = (1 − Φ) + ΦE_low + ζR_tr + ζR_high; ‖ΦE_low‖ ≤ ‖E_low‖ (|Φ| ≤ 1). With ‖a + b‖² ≤ (1+η)‖b‖² + (1+1/η)‖a‖²,
η = L^{−1/2}, b = ΦE_low: D(λ)² ≤ (1 + L^{−1/2})U_{Z(L⁶)}(L − √L) + O(L^{−3/2}) (also for the constrained variant), where G is
admissible by (3i) and §7 of DERIVATION. This needs RH only for the |1/ζ| bounds in (3i) and puts no condition on ζ′(ρ).

## (4) Lemma R — OPEN; precise replacement
Accepted: average density does not give the short-interval Schur bound (under RH a window of length 1/L at height γ can hold
≍ log γ/log log γ zeros); c/L-separated pairs keep correlations ≍ 1/c; adjacent growing jet blocks can be nearly dependent;
jets are not exact for distinct zeros. What the limit needs is only tail domination:
Lemma R′ (open). Assume RH. There is ε(H) → 0 (H → ∞) such that for all large L and all H ≤ L⁶:
  L·U_{Z(L⁶)}(L − √L) ≤ L·U_{Z(H)}(L − √L) + ε(H).
(Weaker, for the order only — Lemma R″: sup_L L·U_{Z(L⁶)}(L − √L) < ∞.) Evidence, not proof: dbar_U (L·U within 1 % of C_K,
1000 zeros), dbar_riesz (Riesz bounds → 1 as density/Nyquist → 0), dbar_cluster (m = 2, 3 clusters at the m² value).

## §C  Lemma C (standard; Titchmarsh 9.6(A), 9.7; Cartan's lemma, Levin, Lectures on entire functions, Lecture 11)
For t ≥ t₀: (a) each [n, n+1] contains t_n with |ζ(σ+it_n)|⁻¹ ≤ t_n^A for −1 ≤ σ ≤ 2, chosen outside the curves in (b);
(b) the zeros with γ ≥ t₀ are enclosed by disjoint closed curves c of diameter ≤ 1/log γ, each with ≤ C log γ zeros, on whose
doubles c′ log|ζ| ≥ −C log γ · log log γ. (From log|ζ(s)| = Σ_{|γ−t|≤1}log|s−ρ| + O(log t) and Cartan's bound
|Π(s − ρ)| ≥ (h/e)^n outside discs of total radius 2h, h = 1/(100 log t).) Not re-proved here.

## Final status
RH + Lemma R′ (+ Lemma C) ⇒ lim_{λ→0} D(λ)² log(1/λ) = Σ_ρ m_ρ²/|ρ|²;  RH + Lemma R″ (+ C) ⇒ D(λ)² ≪ 1/log(1/λ).
Proved here without Lemma R′: Claim F (unconditional) and, under RH (+ Lemma C), the reduction (3v): D(λ)² ≤ (1 + L^{−1/2})·U_{Z(L⁶)}(L − √L) + O(L^{−3/2}).
