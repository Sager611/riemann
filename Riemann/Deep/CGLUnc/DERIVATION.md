# CGLUnc derivation log (cglunc, 12 Nov 2026)

## D1. zeta23 = Alpöge–Furman 2608.13637 (TeX paper-v5-draft18; port Riemann/Deep/WeilEF/Z23)
Window ψ ≥ 0 even on [−1/2,1/2]; L = log(T/2π) (or λ·l, λ ≤ 1); φ ∈ C_c², supp φ = [−L/2, L/2], φ² ≈ ψ(u/L).
Grid α_k = T + 2πk/L, 0 ≤ k < d ≈ N(T,2T). Vectors v_ρ = (φ̂(γ_ρ − α_k))_k ∈ ℂ^d, γ_ρ = (ρ − 1/2)/i (complex off the line).
G̃ = (aL²)⁻¹ Σ_{Re γ_ρ ∈ I'} m_ρ v_ρ v_ρᵀ (real symmetric: v_{1−ρ̄} = conj v_ρ), a = ‖φ‖₂²/L.
Poisson–Gabor (exact, no aliasing at h = 2π/L): Σ_k φ̂(z−α_k)φ̂(z'−α_k) = L·(φ²)^(z − z') for all complex z, z'.
Inputs → three facts:
 (Z) block: P = on-line part ⪰ 0, rank P ≤ N0*, tr P ≤ N0; Q = off-line part = pull-back of ⊕ m·diag(1,−1), n₊(Q) ≤ #pairs.
 (T) tr G̃ = N + o(N) (Poisson–Gabor at z = z' + tail bound; uses N(t,t+1) ≪ log t).
 (P) ‖G̃‖²_F = Σ_{ρ,ρ'} m m' K(γ_ρ − γ_ρ') = (R(ψ) + o(1))N, K(Δ) = (ψ̂(LΔ)/∫ψ)² (Montgomery units x = ΔL/2π:
     K̂(α) = (ψ⋆ψ)(α)/(∫ψ)², supp ⊂ [−1,1]); R(ψ) = ∫K̂(α)(δ(α) + |α|)dα = (∫ψ² + ∬|u−v|ψψ)/(∫ψ)².
     Prime side: explicit formula (EF) turns G̃ into ∫φ̂φ̂ ν_X; ‖·‖² = archimedean·archimedean + prime·prime
     (diagonal n = n' by Chebyshev–Mertens Σ Λ²/n = ½log²; off-diagonal by Montgomery–Vaughan, X ≤ T) + cross.
 (L) rank–trace: P ⪰ 0, rank P ≤ r, n₊(Q) ≤ b ⇒ r ≥ 2trP + 4trQ − 4b − ‖P+Q‖² (eigenvalue form of (x−1)² ≥ 0,
     x² ≥ 0, (x−2)² ≥ 0 after von Neumann; general c: c·trP − c²r/4 + 2c·trQ − c²b ≤ ‖P+Q‖², RankTrace.lean l.160).
Chain (Thm A(i)): P₁ = simple on-line part, Q' = G̃ − P₁ (n₊ ≤ s₂ + p; doubles at 'flat charge 4'):
  N0s ≥ rank P₁ ≥ 4 tr G̃ − 2N − ‖G̃‖² = (2 − R(ψ) − o(1))N. ψ₀ = 1: R = 4/3 → 2/3; ψ_MT = cos(√2 s): R = c_MT⁻¹ = 1.3274992 → 0.6725007.
  λ < 1: R_λ(ψ₀) = 1/λ + λ/3, profile H(λ) = 2 − 1/λ − λ/3, sup_{λ≤1} = 2/3 (Lean Defs.lean Hfun).
Off-line zeros: never estimated individually; their huge entries (|φ̂(x+iy)| ≤ L e^{L|y|/2}) sit in Q, and only
  n₊(Q) ≤ p and the global tr, ‖·‖² (fixed by primes) enter. Negative eigenvalues of Q only help ((p−n)² ≥ … ).
Same chain for N0* (distinct on-line): P = all on-line points (rank ≤ N0*), doubles also on the rank side.

## D2. TEN (arena) = AF chain + Gram-defect of the simple on-line block (Solution.lean of knausgard-ten-win10-lp)
Window W: ψ(s) = cos(√2 s) + Σ_{j≤16} c_j cos(2πjs) on [−1/2,1/2] (c_j = cAMn[j]/10⁹); H_W = 2 − R(ψ) = 3/2 − cotθ₀/√2
 − Σ c_j²(1/2 − 1/(4π²j²))/(2sin²θ₀), θ₀ = 1/√2 (quadratic in c: MT is AF-optimal); TEN: H_W⁻ = 0.6720624204.
Kernel k_W(x) = ψ̂(2πx)/∫ψ (x = normalized ordinate difference ΔL/2π); retained Gram Γ_ij = k_W(x_i − x_j) on simple on-line zeros.
Exact identity: rank-side count r = #retained satisfies r ≥ 2trΓ − ‖Γ‖² + D∘, D∘ = Σ_i(λ_i−1)² = ‖Γ − I‖² ≈ Σ_{i≠j} k_W²(x_i−x_j).
AF chain + D∘: N0s ≥ r ≥ H_W·N + D∘. TEN certificate (10 points, 9 gaps, 45 pair weights a_ij, gap pressures b_r,
 reflection half-space): ∀ gap vectors g ≥ 0, c ≤ Fw(W,g) ⇒ D∘ ≥ c·r − B·span − (n−1)c (ActualDefectEstimate, l.16111).
 Then r(1 − c) ≥ (H_W − B)N: κ = (H_W⁻ − B)/(1 − c). Lean: OpethN.Counting.defect_to_dyadic (l.16116) concludes
 N0s (SIMPLE zeros on the line) ≥ (κ − ε)N; the arena statement (N0*) follows from N0s ≤ N0*.
 So the arena certificate also proves N0s/N > 0.673639 (the TREE N8 line says 'not for simple zeros': to be corrected).
Exact form of the method: with e_W(ν) := inf over configurations of density ν of (1/length)·Σ_{i≠j} k_W²(x_i − x_j),
 an exact local certificate gives at best κ*(W) = least root of ν = H_W + e_W(ν); TEN: e_W(0.6736) ≈ 0.0016 per unit length.
CCLM ceiling 0.6725007 is passed because D∘ (off-diagonal Gram energy of simple on-line zeros) is thrown away by MT/CCLM.


## D3. CGL import: what holds unconditionally, and where it breaks
Notation: zeros in a window, γ_ρ = (ρ−1/2)/i, Δ = (γ_ρ − γ_ρ')L/2π (complex if ρ or ρ' off line), S(K) := Σ_{ρ,ρ'} m m' K(Δ),
K entire, K(x) = ∫K̂(α)e^{2πiαx}dα, K̂ even, K̂(α)e^{c|α|} ∈ L¹ ∀c.
(a) F ≥ 0 survives. Reflection ρ ↦ 1−ρ̄ maps γ ↦ conj γ with equal multiplicity, so for every real-analytic weight η
 (η(z̄) = conj η(z)) and real ξ: conj(Σ_ρ m η(γ_ρ)e^{iξγ_ρ}) = Σ_ρ m η(γ_ρ)e^{−iξγ_ρ}; hence
 Σ_{ρ,ρ'} m m' η(γ_ρ)η(γ_ρ')e^{iξ(γ_ρ−γ_ρ')} = |Σ_ρ m η(γ_ρ)e^{iξγ_ρ}|² ≥ 0, i.e. F_W(α) ≥ 0 for all α, and
 S(K) ≥ 0 whenever K̂ ≥ 0. Pure symmetry: no arithmetic, no RH.
(b) |α| ≤ 1 survives: AF Thm HS (EF + MV + Chebyshev–Mertens), polarised over windows, gives
 S(K) = N(K̂(0) + ∫K̂|α| + o(1)) for K̂ ∈ span{ψ₁⋆ψ₂} ⊂ C_c([−1,1]) (complex Δ included).
(c) Hence CGL's UPPER bound is unconditional: K̂ ≤ 0 on |α| ≥ 1 − ε ⇒ S(K) ≤ N(U(K) + o(1)), U(K) = K̂(0) + ∫_{−1}^{1}K̂|α|
 (split K̂ = K̂χ + K̂(1−χ), χ ∈ C_c(−1,1), χ = 1 on [−1+ε, 1−ε]: S(K(1−χ)) ≤ 0 by (a); ε → 0).
(d) The LOWER bound is what fails. CGL use S(K) ≥ K(0)Σ_γ m_γ² (K ≥ 0 on ℝ, all Δ real). Off-line pair at depth
 y = (β−1/2)L/2π > 0, isolated: self-energy 2m²(K(0) + K(2iy)), K(2iy) = ∫K̂(α)cosh(4παy)dα.
 Lemma O1. If K̂ ≤ 0 on |α| ≥ 1 and K̂ ≠ 0 on a positive-measure subset of |α| ≥ 1 + η, then
 K(2iy) ≤ −c_η e^{4π(1+η)y} + ‖K̂‖₁ e^{4π(1+η/2)y} → −∞ (super-exponentially relative to e^{4πy}).
 Every bandwidth-one quantity (tr, HS norm of any window ψ on [−1/2,1/2], all cross terms) changes by
 O(m² e^{4πy}) per pair (|ψ̂(x+iy)| ≤ e^{π|y|}‖ψ‖₁ in these units) and N by 2m.
 Lemma O2 (deep-pair padding). Given any configuration Z, add n_p = ⌈M N e^{−4π(1+η/2)y}⌉ isolated off-line pairs at
 depth y. Then S(K) drops by ≥ (M − o(1))N while N, N0*, N0s and all bandwidth-one data change by
 O(M N e^{−2πηy}) = o(N) as y → ∞ (y ≤ L/4π allows y ≍ log T). So no inequality
 S(K) ≥ Σ_on c(m) − C·(bandwidth-one data) − o(N) holds for all configurations, and S(K) ≤ N·U(K) is vacuous for N0*/N.
(e) The HS frame cannot host it either. Any second moment the rank–trace device uses is ‖Σ_ρ m M_ρ‖²_F with M_ρ a
 translate of one matrix pattern (modulated windows on the grid); its pair kernel κ(γ−γ') = tr(M_γ M_γ') is positive
 definite (Σ c_i c̄_j κ(γ_i−γ_j) = ‖Σ c_i M_γi‖² ≥ 0), so κ̂ ≥ 0, and F_W ≥ 0 bounds the unknown |α|>1 part from BELOW only,
 while the rank–trace inequality needs an UPPER bound for ‖·‖². CGL's gain needs κ̂ < 0 beyond 1: not a pd kernel.
(f) Zero-density inputs (beyond zeta23's list) block O2 only if N(1/2+δ,T) ≪ T^{1−aδ}log T with a > 2(1+η) > 2
 (padding needs fraction e^{−4π(1+η/2)y}, density allows C e^{−2πay}); best known a = 8/7 (Conrey/Baluyot, θ < 4/7).
Freedoms tested (lesson 191): analytic weights η (|η(γ)|² multiplies the pair term; no real-analytic η is small on
 Im z = ±y and ≈ 1 on ℝ, by subharmonicity of log|η| in the strip), complex shifts K(Δ + iσ) (K̂ → K̂e^{−2πσα}: same tail),
 complex/sign-changing windows (still pd HS kernels), several windows (sums of pd kernels), more points (TEN-type
 defects are bandwidth-one and configuration-wise: covered by AF's ceiling).


## D4. TEN-class barrier (dual certificates; code in this folder)
Class C_J: AF chain for one window (or any mixture / complex combination: Gram matrix G ⪰ 0 of coefficients) in
 span{cos(√2 s), cos(2πjs), j ≤ J}, plus ANY configuration-wise bound c·r − B·span − O(1) ≤ D∘ on the simple on-line block
 (any number of points n, pair weights, gap pressures, floors, reflection rows): κ = (H_W − B)/(1 − c).
Dual lemma. If an S-periodic configuration X of density ν = p/S has e_W(X) := lim D∘/span ≤ ν − H_W, then every valid
 (c,B) has cν − B ≤ e_W(X), so κ ≤ (H_W − cν + ν − H_W)/(1 − c) = ν. (ActualDefectEstimate is proved for every ZeroConfig,
 no EF used, so it must hold for M periods of X placed at height T.)
Poisson (exact, K̂ = ψ⋆ψ/(∫ψ)² supported in [−1,1]): e_W(X) = S⁻²[K̂(0)p² + 2Σ_{1≤m<S} K̂(m/S)|Σ_x e^{2πimx/S}|²] − p/S;
 checked against a real-space sum (|d| ≤ 4000: 0.0018003 vs 0.0018005).
Arena window (knausgard-ten-win10-lp): X = (.xx)^12(.xx.xx.xxx), S = 46, p = 31, relaxed: e = 0.0017782478544798,
 ν − H_W − e = 7.1796e−5 > 0 (Arb, certify_dual.py): κ ≤ 31/46 = 0.6739130; float hull limit 0.673840 (necklaces S ≤ 16).
 Arena: 0.6736393, so the window has ≤ 2.0e−4 left (float), 2.7e−4 (certified).
All windows: 12 configurations of density 153/227 (block mixtures of .xx, .xx.xx.xxx, .xx.xx.xx.xxx, ..xxx.xx.xxx, .xx.xxx,
 relaxed at the optimum window, the arena window and MT) and λ from an SDP make
 M(λ) = [[2ν − H_MT − Q̄₀₀, −q̄ᵀ], [−q̄, W − Q̄']] ≻ 0 (e(X;c) = cᵀQ_Xc − ν, H(c) = H_MT − Σ w_j c_j²,
 w_j = (1/2 − 1/(4π²j²))/(2sin²θ₀)); Arb LDLᵀ pivots > 0 (min 3.85e−4) for J = 16 and J = 48 (verify_cover.py).
 Since tr(M G) ≥ 0 for every G ⪰ 0, this covers single windows, mixtures and complex coefficients:
 TEN class ≤ 153/227 = 0.6740088 (certified); float optimum over 16-harmonic windows 0.67394–0.67397 (winopt.py).


## D5. Review round 1 (reviewer, 15:06) and repairs
Claim B (TEN-class dual, S-procedure, Poisson normalisation): CORRECT for the stated class (J ≤ 48; admissible c < 1).
Claim A (padding): GAP as first stated — (1) smooth rapid decay does not make ∫K̂cosh(4παy) finite (needs all exponential
 moments; already the D3 hypothesis, and necessary for S(K) = N∫K̂F with complex zeros); (2) the bandwidth-one change
 needs ‖E‖²_F = o(N) for the whole padding and Cauchy–Schwarz for 2⟨G̃,E⟩, with an explicit separation; (3) scope.
Repairs (PROOF.md Prop O): base configurations of bounded depth y₀ (on-line ones suffice for the obstruction);
 I₂(y) ≍ MN^{1/2}, n_p ≍ N^{1/2} pairs at mutual distance ≍ N^{1/2}; I₁(y+y₀) ≪ I₂(y)^{1/2} (Cauchy–Schwarz) controls the
 cross terms with Z (≪ log T zeros per unit x-length); ‖E‖²_F = O(N^{1/2+1/(2A₁)}) via AF's Paley–Wiener bound (PW).
 Without the depth restriction a deep base could carry cross terms of size I₁(y + y_max): not needed for the obstruction.


## D6. Review round 2 (15:10) and repairs; single-window escape-route test
Round 2: GAP — (1) integration-by-parts constant C_k(η) ≪ (1+η)^k I₁(η) unjustified for the general exponential-moment
 class; (2) zero-density 'iff a < 2A₁' false: for Gaussian-type K̂ the padding fraction M/I₂(y) = e^{−cy²+…} is compatible
 with every slope. Repairs: Prop O restricted to C1 (band-limited C^k_c) and C2 (CGL's p(x)e^{−πx²}) with explicit
 Γ_k(η) = Σ_j (1+η)^{k−j}∫|∂^jK̂|e^{2π|α|η} bounds; zero density split: C1 needs slope > 2A to block, C2 cannot be blocked
 by any slope (padding sits at |β − 1/2| ≍ (log T)^{−1/2}); only thin-box hypotheses (BGSTB, GS: b/log T) stop it.
Escape route 'any configuration-wise lemma on one window's data' (marks.py): on-line S-periodic patterns with masses
 {0,1,2} (S ≤ 9–10), relaxed, HS/length = K̂(0) + 2S⁻²Σ_m K̂(m/S)|Σ m_j e^{2πimx_j/S}|² ≤ R_W, LP over block mixtures:
 simple fraction ≥ … caps: MT 0.676483 (.121 + .11111121), record window 0.677662 (.1121 + .111111121),
 float-optimal window 0.677965 (.11211 + .111121). So Q-side/cross-term refinements on one window have ≤ 4.0e−3 room.


## D7. Family-wide one-window cap: not obtained
cover_marks.py: 720 relaxed mark patterns (S ≤ 8, relaxed at record/MT/optimal windows); maximise t with
 Ψ − Σλ_iΦ_i ⪰ tI (Ψ = 2∫₀¹T_sym(α)α dα, Φ_i the patterns' form-factor matrices): t* = −0.0181 < 0 (SCS; the ν-capped
 runs were cut by the 1200 s watchdog). No window-independent double adversary exists in this library, so one-window
 caps are per window only (float: MT 0.676483, record 0.677662, optimal TEN window 0.677965).


## D8. Orchestrator referee (Prop D/D2: GAP repairable) — repaired in PROOF.md / REPORT.md
Class = arena ActualDefectEstimate for every ZeroConfig (Z_X on-line simple, fields checked); augmented normalised Gram
G₀₀ = 1 with certified rational λ; one common configuration X̄_B by operator-norm convergence + PD margin; complex only
Hermitian/|k|²; family J ≤ 96 (orchestrator rerun J = 96, min pivot 3.8468e−4); uncovered freedoms listed.
