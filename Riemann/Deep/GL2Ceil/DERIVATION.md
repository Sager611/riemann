# DERIVATION — gl2ceil (ceiling of the one-piece Levinson–Conrey method, nodes N8, N12)

Sources (TeX line numbers): B = Bernard arXiv 1404.6432 (Modular_case_of_Levinson_s_theorem.tex); AT21 = Andersen–Thorner arXiv 2004.03581
(zeros-paper.tex); KRZ = arXiv 1609.03738; S = arXiv 2508.11108. Degree2Mol = Riemann/Deep/Degree2Mol/DERIVATION.md (§ refs there).

## 0. The method and its admissible class

- B l.161, AT21 l.192–194: mollifier ψ with M = T^ν, σ₀ = 1/2 − R/log T, R > 0 real, P ∈ R[x], P(0) = 0, P(1) = 1.
- B Thm th:mm2 l.770–783 (AT21 l.195–211): Q ∈ C[x], Q(0) = 1, V = Q(−(1/(2 log T)) d/ds) L(f,s); (1/T)∫|Vψ(σ₀+it)|² → c(P,Q,2R,ν/2),
  c(P,Q,r,ξ) = 1 + (1/ξ)∫₀¹∫₀¹ e^{2rv}|P′(u)Q(v) + ξP(u)Q′(v) + ξrP(u)Q(v)|² du dv (AT21 l.210: modulus squared).
- B l.785–795: zero count for Q = 1 + Σ_{n≤M} i^{n+1}λ_n[(1−2x)ⁿ − 1], λ ∈ R^M (any M); AT21 l.214–218: κ_f ≥ 1 − inf (1/R) log c(P,Q,R,ν/2).
- Lengths: AT21 Thm thm:nu l.227–231: ν < 1/4 − θ/2, θ = 7/64 → ν < 25/128 (ξ = ν/2 < 25/256); θ = 0 → ξ < 1/8. B l.796: ν < 5/27.

Definition. κ*(ξ) := sup { 1 − (1/R) log c(P,Q,R,ξ) : R > 0, P ∈ R[x] with P(0)=0, P(1)=1, Q in B's class }.
The best bound the method gives on the open range ν < ν₀ is sup_{ξ<ν₀/2} κ*(ξ) ≤ κ*(ν₀/2) (§1.7). We certify κ*(25/256) and κ*(1/8).

## 1. Reductions

1.1 Complex Q never helps (real P). B's Q = Q_r + iQ_i with Q_r = 1 + Σ_{n odd} ±λ_n[(1−2x)ⁿ − 1] (real, Q_r(0)=1, Q_r(y)+Q_r(1−y) const) and
Q_i = Σ_{n even} ±λ_n[(1−2x)ⁿ − 1] real. P real ⇒ the integrand |A_r + iA_i|² = A_r² + A_i² with A_r, A_i the brackets of Q_r, Q_i, so
c(P,Q) = c(P,Q_r) + (1/ξ)∫∫e^{2Rv}A_i² ≥ c(P,Q_r). Hence κ* = sup over the real class 𝒬 = {Q ∈ R[x]: Q(0) = 1, Q(y) + Q(1−y) ≡ β}
(Q ∈ 𝒬 ⇔ Q − β/2 odd in 1−2y ⇔ Q = 1 + Σ_{odd} b_n[(1−2y)ⁿ − 1]: exactly the real part of B's class).
Side remark: AT21 l.996 quotes inf_P c with √(AC − Im(B)²); for complex Q that is the infimum over complex P. With real P (their l.194)
the cross term is 2ξPP′Re(w w̄′) and inf_P c = (1 + e^{2R}|Q(1)|²)/2 + √(AC) coth(ξ√(C/A)) (A = ∫|w|², C = ∫|w′|²) ≥ their value;
harmless for AT21, which restricts to real Q (l.1002), where both agree.

1.2 Expansion (Degree2Mol §2): w(y) = e^{Ry}Q(y), W₀ = ∫₀¹w², W₁ = ∫₀¹w′², B_P = ∫P², C_P = ∫P′²:
   c = (1 + w(1)²)/2 + (C_P/ξ) W₀ + ξ B_P W₁,   w(0) = 1, w(1) = e^R(β − 1).                                   (E)

1.3 Closure. c is continuous on H¹(0,1) × H¹(0,1) (B_P, C_P, W₀, W₁ are continuous quadratic forms, w(1) is a continuous trace).
Polynomials with P(0)=0, P(1)=1 are H¹-dense in 𝒫 = {P ∈ H¹: P(0) = 0, P(1) = 1} (approximate P′ in L² by polynomials p_n, P_n = ∫₀^x p_n
+ x(1 − ∫₀¹p_n)). 𝒬 is H¹-dense in 𝒬̄ = {Q ∈ H¹: Q(0) = 1, Q(y) + Q(1−y) ≡ β, β ∈ R}: write Q = β/2 + g(1−2y), g odd in H¹(−1,1),
g(1) = 1 − β/2; approximate g′ (even) in L² by even polynomials, integrate from 0 (odd polynomials g_n → g in H¹), correct by
(g(1) − g_n(1))z. So sup over polynomials = sup over 𝒫 × 𝒬̄ (any degree; the H¹ closure adds nothing).

1.4 Inner minimum over P (Conrey §4; Degree2Mol §3). For a, b > 0, aC_P + bB_P over 𝒫 is a strictly convex quadratic with unique minimiser
P_λ = sinh(λx)/sinh λ, λ = √(b/a) (Euler–Lagrange P″ = λ²P). With a = W₀/ξ, b = ξW₁: inf_P c(P,Q) = c(P_{λ(Q)},Q). Hence
   κ*(ξ) = sup_{R>0, λ≥0} 1 − (1/R) log M(R,λ),   M(R,λ) := inf_{Q∈𝒬̄} J,  J = (1 + w(1)²)/2 + α∫w² + γ∫w′²,  α = C_λ/ξ, γ = ξB_λ,
C_λ = λ coth λ/2 + λ²/(2 sinh²λ), B_λ = (sinh 2λ/(4λ) − 1/2)/sinh²λ; λ = 0 is P = x (C = 1, B = 1/3).

1.5 Monotonicity in λ. For λ₁ < λ₂ minimality gives C₁ + λ₁²B₁ ≤ C₂ + λ₁²B₂ and C₂ + λ₂²B₂ ≤ C₁ + λ₂²B₁; adding, (λ₂² − λ₁²)(B₂ − B₁) ≤ 0,
so B_λ is non-increasing; then C₁ − C₂ ≤ λ₁²(B₂ − B₁) ≤ 0, so C_λ is non-decreasing. Consequences: C_λ ≥ 1, B_λ/C_λ ≤ 1/3.
Also C_PB_P > 1/4 for every P ∈ 𝒫: with a = B_P, b = C_P, 2B_PC_P ≥ inf_P' (B_P C_P' + C_P B_P') = √(B_PC_P) coth√(C_P/B_P) > √(B_PC_P).

1.6 For fixed (R, λ), J is a positive quadratic form plus linear terms on the affine space 𝒬̄; its minimiser solves the Euler–Lagrange BVP
(y-form of S §sect:gene, t = R(2y−1)): on [1/2, 1], Q″ + 2R tanh(R(2y−1)) Q′ + (R² − ω²)Q = (R² − ω²) β/(1 + e^{2R(2y−1)}), ω² = C_λ/(B_λξ²),
Q(1/2) = β/2, Q(1) = β − 1, Q(1−y) = β − Q(y); β minimises the quadratic J(β). Used only to build trial dual functions (floats).

1.7 κ* is non-decreasing in ξ: inf_P c = (1 + w(1)²)/2 + √(W₀W₁) coth(ξ√(W₁/W₀)) decreases in ξ for every Q, R (Degree2Mol §4b).

## 2. Weak duality (the certificate's engine)

Lemma D. Let R, α, γ > 0, σ ∈ H¹(0,1) (here: continuous piecewise polynomial), μ ∈ L²(0,1) with μ(y) = μ(1−y), M := ∫₀¹μ. Then
   M(R, λ) ≥ D := 1/2 + M − 2σ(0) − (2σ(1) + Me^{−R})²/2 − (1/α)∫₀¹(σ′ + μe^{−Ry})² − (1/γ)∫₀¹σ².
Proof. Let Q ∈ 𝒬̄, w = e^{Ry}Q, p := σ′ + μe^{−Ry}. Pointwise αw² ≥ 2pw − p²/α, γw′² ≥ 2σw′ − σ²/γ. Integration by parts:
∫2σw′ = 2σ(1)w(1) − 2σ(0) − ∫2σ′w (w(0) = 1). So α∫w² + γ∫w′² ≥ 2∫μe^{−Ry}w + 2σ(1)w(1) − 2σ(0) − ∫(p²/α + σ²/γ). Folding y ↔ 1−y,
2∫₀¹μe^{−Ry}w = 2∫₀^{1/2}μ(y)[Q(y) + Q(1−y)]dy = βM. With u = β − 1, w(1) = e^R u: J ≥ 1/2 + (e^{2R}/2)u² + (2σ(1)e^R + M)u + M − 2σ(0)
− ∫(p²/α + σ²/γ) ≥ (minimum over u) D. ∎
Strong duality (not needed): σ = γw′, μ = e^{Ry}(αw − γw″) at the minimiser give D = J (symmetry of μ is the Euler–Lagrange equation);
floats: gap/J ≤ 1e-13 at (R,λ) = (9.87, 0.77), (6, 0.77), (14, 1.2), (20, 2), (2, 0.3) (gapfloat.py, logs/gapfloat.log).

## 3. Cell bounds (fixed σ, μ, varying R and λ)

For fixed σ, μ: D(R; α, γ) = K₀(R) − K₁(R)/α − K₂/γ with K₀ = 1/2 + M − 2σ(0) − (2σ(1) + Me^{−R})²/2, K₁(R) = ∫(σ′ + μe^{−Ry})² ≥ 0,
K₂ = ∫σ² ≥ 0. Write X(R) = ∫σ′μe^{−Ry}, Y(R) = ∫μ²e^{−2Ry}, S₁ = ∫σ′²; K₁ = S₁ + 2X + Y.
- λ-subcell [λ_i, λ_{i+1}] (dyadic grid 2^{-HLBITS}, HLBITS = 16–21; subcells bisected only where needed): by 1.5, 1/α ≤ x := ξ/C_{λ_i}, 1/γ ≤ z := 1/(ξB_{λ_{i+1}}), so D ≥ L(R) := K₀(R) − xK₁(R) − zK₂.
- R ∈ [R_c − δ, R_c + δ], R_lo = R_c − δ, ε = R − R_c: |K₀″| ≤ F₀ := M²e^{−2R_lo} + (|2σ(1)| + |M|e^{−R_lo})|M|e^{−R_lo};
  K₁ ≤ K₁(R_c) + K₁′(R_c)ε + K₁″(R_c)ε²/2 + G₁|ε|³/6 with K₁″ = 2∫y²σ′μe^{−Ry} + 4∫y²μ²e^{−2Ry} and
  |K₁‴| ≤ 2√S₁ √Y(R_lo) + 8Y(R_lo) =: G₁ (Cauchy–Schwarz, y ≤ 1, Y decreasing in R).
  lower(ε) := K₀ + K₀′ε − x(K₁ + K₁′ε) − zK₂ + a₂ε² − F₀ε²/2 − xG₁|ε|³/6, a₂ := min(lower end of −xK₁″(R_c)/2, 0), is concave and ≤ L.
  target(R) = e^{(1−κ_up)R} is convex, so lower − target ≥ 0 on the cell iff it holds at ε = ±δ: that is what check_lam verifies (Arb).
  Then M(R,λ) ≥ e^{(1−κ_up)R}, i.e. 1 − log M/R ≤ κ_up, on the whole cell.
- Exact integrals: σ, μ are piecewise polynomials on 32 pieces of [0,1] (μ: 16 pieces on [1/2,1], reflected), local variable x ∈ [0,1];
  ∫₀¹ xʲ e^{ax} dx = ₁F₁(j+1; j+2; a)/(j+1) (Arb hypgeom_1f1). σ's piece constants are reset to the previous piece's end value
  (sums of doubles, exact at 256 bits: asserted rad = 0), so σ is exactly continuous; μ is symmetric by construction.

## 4. Tails (analytic, constants in Arb in ceil.py tails())

4.1 Closed form at R = 0: the minimiser for fixed β is Q = β/2 + (β/2 − 1) sinh(ω(y − 1/2))/sinh(ω/2), and
    J(β) = (1 + (β−1)²)/2 + (C_λ/(4ξ))β² + 2√(C_λB_λ) coth(ω/2)(1 − β/2)².  With C_λ ≥ C_Λ (Λ ≥ 0), √(CB) ≥ 1/2, coth ≥ 1:
    M(0, λ) ≥ min_β [2 − 2β + β²(3/4 + C_Λ/(4ξ))] = 2 − 1/(3/4 + C_Λ/(4ξ))  (Λ = 0: 2 − 4ξ/(3ξ + 1) = 1.6978852 at ξ = 25/256).
4.2 Slope: for fixed λ, Q: dc/dR = e^{2R}Q(1)² + 2α∫ye^{2Ry}Q² + 2γ∫ye^{2Ry}(RQ+Q′)² + 2γ∫e^{2Ry}(RQ+Q′)Q ≥ −2γ√(W₀W₁) ≥ −√(γ/α)·c,
    √(γ/α) = ξ√(B_λ/C_λ) ≤ ξ/√3. So M(R,λ) ≥ M(0,λ)e^{−sR} and κ ≤ 1 + s − log M(0)/R ≤ κ_up for R ≤ R_tail = log m₀/(1 + s − κ_up)
    (ξ = 25/256: R_tail = 0.53672 ≥ 1/2 = start of the core; for λ ≥ Λ = 3: 0.6045).
4.3 Derivative-free bound: drop ξB W₁ ≥ 0, use C ≥ C_Λ and, for each y ≤ 1/2, the pointwise minimum of e^{2Ry}q² + e^{2R(1−y)}(β − q)²:
    W₀ ≥ β² e^R gd(R)/(4R) (gd(R) = ∫₀ᴿ sech), so c ≥ 1/2 + e^{2R}K/(e^{2R} + 2K) =: crude(R, C), K = (C/ξ)e^R gd(R)/(4R).
    crude increases in C and in R (d log(e^R gd/R)/dR = 1 + sech R/gd R − 1/R > 0: trivial for R ≥ 1; for R < 1 use gd R ≤ R and
    cosh R (1 − R) < 1). Used (i) on cells [R_a, R_b] × [λ_a, λ_b]: crude(R_a, C_{λ_a}) > target(R_b); (ii) on the strip λ ≥ Λ = 3 for
    R ∈ [0.59, 28] (λ ≥ 3, R ≤ 0.6045 by 4.2 with C_3, B_3); (iii) for R ≥ 28, with C = 1 ≤ C_λ (crude increases in C, so this covers every λ ≥ 0; orchestrator review 12 Nov): log crude ≥ log K − 2Ke^{−2R}, K = e^R gd(R)/(4Rξ) ≤ e^R, so κ ≤ [log(4Rξ/gd(28)) + 2e^{−28}]/R, which decreases for
    R ≥ e gd(28)/(4ξ) = 10.9 and equals 0.0693077 < κ_up at R = 28.

## 5. Cover and bookkeeping

Core (R, λ) ∈ [R0, RMAX] × [0, Λ] (ξ = 25/256: [1/2, 28] × [0, 3]; ξ = 1/8: [1/2, 22] × [0, 3], R_tail = 0.50127, h(22) = 0.0884688;
ζ, ξ = 4/7: [1/8, 6] × [0, 6], R_tail = 0.15919, h(6) = 0.36250), initial cells 1/2 × 1/4, depth-first; a cell is closed by the crude
bound (zeroth order, 4.3) or by the dual (§3) with σ, μ built at its centre; otherwise it is bisected (in R if 13·(R-width) > λ-width,
the measured ratio of the two fixed-dual reaches, else in λ). Union = box: covcheck.py checks exact area and pairwise disjointness.
Remaining regions: R ≤ R_tail (4.2, all λ), λ ≥ Λ (4.2 with C_Λ, B_Λ for R ≤ R_tail(Λ), crude on [R_tail(Λ), RMAX]), R ≥ RMAX (4.3 iii).
So 1 − log M(R,λ)/R ≤ κ_up for all R > 0, λ ≥ 0, hence κ*(ξ) ≤ κ_up.
replay.py re-derives every cell inequality, the tails and the strip with independent code (exact-rational integrals via the finite
alternating formula of Degree2Mol §5, own bound assembly); the float proposer build() is the only shared part and is not trusted.


## Audit note (gl2audit, 12 Nov 2026)
- §1.1: the same complex-P infimum √(AC − (Im B)²) is displayed in Bernard 1404.6432 l.796–800 (as in AT l.993–996); for real P the infimum is √(AC)·coth(ξ√(C/A)). AUDIT PASSED: Riemann/Deep/GL2Audit/REPORT.md.
