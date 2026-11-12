# FULLPROOF — gl2audit (referee of discoveries/gl2-levinson-ceiling.md, nodes N8, N12)

Sources (TeX line numbers, e-prints fetched 12 Nov 2026): B = Bernard 1404.6432 (single gzip TeX, 1127 lines); AT = Andersen–Thorner
2004.03581 (zeros-paper.tex, 1107 lines). Notation: ξ = ν/2, w = e^{Ry}Q, W₀ = ∫₀¹|w|², W₁ = ∫₀¹|w′|², B_P = ∫P², C_P = ∫P′².

## 0. The input, verbatim hypotheses

- B l.159–161: ψ(s) = Σ_{n≤M} μ_f(n) n^{−s−1/2+σ₀} P(log(M/n)/log M), M = T^ν, σ₀ = 1/2 − R/log T, R > 0 real,
  "P is a real polynomial satisfying P(0)=0, P(1)=1". AT l.191–194: same, "P ∈ R[x]".
- B Thm th:mm2 l.770–783: "Q a polynomial with complex coefficients satisfying Q(0)=1", V = Q(−(1/(2 log T)) d/ds) L(f,s);
  for ν < (1−2θ)/(4+2θ): (1/T)∫₁ᵀ|Vψ(σ₀+it)|² = c(P,Q,2R,ν/2) + o(1), with B's display l.779
  c = 1 + (1/ξ)∫∫e^{2rs}[d/dx(e^{rξx}Q(s+ξx)P(x+u))|_{x=0}]² du ds. The bracket is rξPQ + ξPQ′ + P′Q; B writes [·]²,
  AT l.209 writes |·|² (the only real-valued reading for complex Q; a mean of |Vψ|² is real). Audit uses |·|².
- B l.785–795: zero count needs Q of shape (formepolynome) Q = 1 + Σ_{n≤M} i^{n+1}λ_n[(1−2x)ⁿ − 1], λ ∈ R^M, any M;
  liminf N_{f,0}/N_f ≥ limsup[1 − (1/(2R)) log((1/T)∫|Vψ|²)] ⇒ κ_f ≥ 1 − (1/(2R)) log c(P,Q,2R,ν/2).
  Renaming r = 2R (r ranges over all of (0,∞)): κ_f ≥ 1 − (1/r) log c(P,Q,r,ν/2) (B l.796, AT l.214–218). ✓ normalisation ξ = ν/2.
- AT Thm thm:nu l.227–231: the same limit for 0 < ν < 1/4 − θ/2; θ = 7/64 ⇒ 25/128 (ξ < 25/256); θ = 0 ⇒ 1/4 (ξ < 1/8).
- Freedoms the input admits: R > 0; real P, P(0)=0, P(1)=1, any degree; Q in (formepolynome), any M (complex: imaginary even
  terms); ν in the open range. NOT admitted as stated: complex P (B l.161, AT l.194 say real), shifts α, β (B l.68: the shifted moment th:mm2 is deduced from; the zero count uses th:mm2 only), other mollifier shapes, other derivative scales (the scale 1/(2 log T) and the shape are
  what make the zero count of l.790 work). P(1) = 1 is no restriction: P → tP multiplies both c and the Littlewood edge term by t².
- Caution found in the inputs: B l.796–800 (lemma "[C], part 4") and AT l.993–996 state inf_P c with √(AC − Im(B)²); that is the
  infimum over COMPLEX P (§1.6 below). Under the stated real-P hypothesis it is false when Im B ≠ 0 (real P gives √(AC)).
  Both papers then restrict to real Q (B l.801–806, AT l.1001–1004), where the two agree, so their numbers are unaffected.

## 1. Reductions (complete proofs)

1.1 Expansion. For real P and any (complex) H¹ function Q with Q(0) = 1: e^{rv}·bracket = P′(u)w(v) + ξP(u)w′(v), so
  c = 1 + (1/ξ)[C_P W₀ + 2ξ Re(∫₀¹PP′ du · ∫₀¹ w w̄′ dv) + ξ²B_P W₁].  ∫PP′ = (P(1)² − P(0)²)/2 = 1/2 and
  Re∫w w̄′ = (|w(1)|² − |w(0)|²)/2 = (|w(1)|² − 1)/2, hence   c = (1 + |w(1)|²)/2 + (C_P/ξ)W₀ + ξB_P W₁.          (E)
1.2 Complex Q never helps (real P). Q in (formepolynome): odd n give real coefficients i^{n+1} = (−1)^{(n+1)/2}, even n give
  imaginary ones. Q = Q_r + iQ_i, Q_r = 1 + Σ_odd ±λ_n[(1−2x)ⁿ−1], Q_i = Σ_even ±λ_n[(1−2x)ⁿ−1]; real P ⇒ bracket = A_r + iA_i with
  A_r, A_i real ⇒ |·|² = A_r² + A_i² ⇒ c(P,Q) = c(P,Q_r) + (1/ξ)∫∫e^{2rv}A_i² ≥ c(P,Q_r). Q_r ranges over exactly
  𝒬 = {Q ∈ R[x]: Q(0)=1, Q(y)+Q(1−y) ≡ β}: Q ∈ 𝒬 ⇔ Q − β/2 odd in (1−2y) ⇔ Q = 1 + Σ_odd b_n[(1−2y)ⁿ − 1]. ✓

1.3 Inner minimum over P. In (E) put a = W₀/ξ > 0 (w(0) = 1 ⇒ W₀ > 0), b = ξW₁ ≥ 0. F(P) = aC_P + bB_P is a strictly convex
  quadratic on the affine space 𝒫 = {P ∈ H¹: P(0)=0, P(1)=1} (C_P strictly convex in P′, B_P convex). Its Euler–Lagrange equation
  aP″ = bP with the boundary values has the solution P_λ = sinh(λx)/sinh λ, λ = √(b/a) (P_0 = x); for every h ∈ H¹₀,
  F(P_λ + h) = F(P_λ) + a∫h′² + b∫h² (the cross term 2∫(aP_λ′h′ + bP_λh) = 2[aP_λ′h]₀¹ = 0), so P_λ is the unique minimiser and
  min F = a·λcothλ = √(W₀W₁) coth(ξ√(W₁/W₀)) (for W₁ = 0: a). Hence for real P, Q ∈ 𝒬, r = R:
  c(P,Q) ≥ c(P_{λ(Q)},Q) = J_{λ(Q)}(Q) ≥ M(R,λ(Q)),  J_λ := (1 + w(1)²)/2 + α∫w² + γ∫w′², α = C_λ/ξ, γ = ξB_λ,
  M(R,λ) := inf over 𝒬̄ (H¹ closure: Q ∈ H¹ real, Q(0)=1, Q(y)+Q(1−y) ≡ β, β ∈ R free) of J_λ.
  So 1 − (1/R)log c(P,Q,R,ξ) ≤ sup_{λ≥0} [1 − (1/R)log M(R,λ)]: the barrier needs only this inequality; density of polynomials
  (DERIVATION 1.3) is needed only for equality and is not load-bearing.
  C_λ = ∫P_λ′² = λcothλ/2 + λ²/(2sinh²λ), B_λ = ∫P_λ² = (sinh2λ/(4λ) − 1/2)/sinh²λ (direct integration; C_0 = 1, B_0 = 1/3).
1.4 Monotonicity in λ. m(λ) := C_λ + λ²B_λ = min_P (C_P + λ²B_P). For λ₁ < λ₂: C₁ + λ₁²B₁ ≤ C₂ + λ₁²B₂ and C₂ + λ₂²B₂ ≤ C₁ + λ₂²B₁;
  adding: (λ₂² − λ₁²)(B₂ − B₁) ≤ 0 ⇒ B non-increasing; then C₁ − C₂ ≤ λ₁²(B₂ − B₁) ≤ 0 ⇒ C non-decreasing.
  Hence C_λ ≥ 1, B_λ ≤ 1/3, B_λ/C_λ ≤ 1/3, and B_λ > 0. Also B_PC_P ≥ 1/4 for every P ∈ 𝒫 (Cauchy–Schwarz: 1/2 = ∫PP′ ≤ √(B_PC_P)).
1.5 Monotonicity in ξ. inf_P c = (1 + w(1)²)/2 + √(W₀W₁) coth(ξ√(W₁/W₀)) (1.3) is non-increasing in ξ for each (Q, R), and the class
  of (P, Q, R) does not depend on ξ ⇒ κ*(ξ) := sup_{R,P,Q}[1 − (1/R)log c] is non-decreasing ⇒ for the open range ν < ν₀:
  sup_{ν<ν₀} κ*(ν/2) ≤ κ*(ν₀/2). Certifying at the endpoint ξ₀ = 25/256 (resp. 1/8) is therefore legitimate. ✓
1.6 Complex P (NOT in the stated class; the natural sesquilinear extension). P = a + ib, a real in 𝒫, b real in H¹₀:
  ∫P′P̄ = 1/2 + iτ, τ = ∫(ab′ − a′b) = 2∫ab′. Then c = (1 + |w(1)|²)/2 + (C_P/ξ)W₀ + ξB_P W₁ − 2τ·Im B, B := ∫w w̄′.
  (i) Real Q (Im B = 0): c(P) = c(a) + (1/ξ)∫∫e^{2rv}(b′Q + ξb(Q′+rQ))² ≥ c(a): complex P cannot help (page's claim, ✓).
  (ii) Complex Q in B's class: Q_i(0) = Q_i(1) = 0, Im B = ∫(w_i w_r′ − w_r w_i′) ≠ 0 in general, and −2τ Im B is a bilinear
  first-order coupling: jointly complex (P, Q) is NOT covered by any argument; inf over complex P gives exactly B l.796 / AT l.996:
  (1 + |w(1)|²)/2 + X coth(ξX/W₀), X = √(W₀W₁ − (Im B)²). Tested numerically in §6.

## 2. Lemma D (weak duality) — proof checked line by line

Statement: R, α, γ > 0; σ ∈ H¹(0,1); μ ∈ L²(0,1), μ(y) = μ(1−y); M := ∫μ. For every Q ∈ 𝒬̄ (any β ∈ R), with w = e^{Ry}Q:
  J(Q) ≥ D := 1/2 + M − 2σ(0) − (2σ(1) + Me^{−R})²/2 − (1/α)∫(σ′ + μe^{−Ry})² − (1/γ)∫σ².
Proof. p := σ′ + μe^{−Ry}. α(w − p/α)² ≥ 0 ⇒ αw² ≥ 2pw − p²/α; γ(w′ − σ/γ)² ≥ 0 ⇒ γw′² ≥ 2σw′ − σ²/γ. Integration by parts
(σ, w ∈ H¹, absolutely continuous): ∫2σw′ = 2σ(1)w(1) − 2σ(0)w(0) − ∫2σ′w, w(0) = Q(0) = 1. Sum: α∫w² + γ∫w′² ≥
2∫μe^{−Ry}w + 2σ(1)w(1) − 2σ(0) − ∫(p²/α + σ²/γ). e^{−Ry}w = Q and μ symmetric: 2∫₀¹μQ = 2∫₀^{1/2}μ(y)[Q(y) + Q(1−y)] = 2β∫₀^{1/2}μ
= βM. Q(1) = β − Q(0) = β − 1 =: u, w(1) = e^R u, β = u + 1: J ≥ 1/2 + e^{2R}u²/2 + (2σ(1)e^R + M)u + M − 2σ(0) − ∫(…);
min_u (e^{2R}u²/2 + bu) = −b²e^{−2R}/2 with b = 2σ(1)e^R + M gives −(2σ(1) + Me^{−R})²/2. ∎ (Sharpness, not needed: σ = γw′,
μ = e^{Ry}(αw − γw″) at the Euler–Lagrange minimiser. Float gap ≤ 1e-13 (GL2Ceil/logs/gapfloat.log) re-read, consistent.)
So M(R,λ) ≥ D for every admissible (σ, μ): any one (σ, μ) certifies a lower bound on M — no optimality of (σ, μ) is needed.

## 3. Cell bound (fixed σ, μ; R and λ vary) — derivatives recomputed

D = K₀(R) − K₁(R)/α − K₂/γ, K₀ = 1/2 + M − 2σ(0) − (2σ(1) + Me^{−R})²/2, K₁ = S₁ + 2X(R) + Y(R) ≥ 0, X = ∫σ′μe^{−Ry},
Y = ∫μ²e^{−2Ry}, S₁ = ∫σ′², K₂ = ∫σ² ≥ 0.
- λ-subcell [λ_i, λ_{i+1}]: by 1.4 1/α = ξ/C_λ ≤ ξ/C_{λ_i} =: x, 1/γ = 1/(ξB_λ) ≤ 1/(ξB_{λ_{i+1}}) =: z; K₁, K₂ ≥ 0 ⇒ D ≥ K₀ − xK₁ − zK₂. ✓
- K₀′ = Me^{−R}(2σ(1) + Me^{−R}); K₀″ = −2σ(1)Me^{−R} − 2M²e^{−2R} ⇒ |K₀″(R)| ≤ 2|σ(1)M|e^{−R_lo} + 2M²e^{−2R_lo} = F₀ for R ≥ R_lo
  (DERIVATION's F₀ = M²e^{−2R_lo} + (|2σ(1)| + |M|e^{−R_lo})|M|e^{−R_lo} is the same number). ✓
- K₁′ = −2∫yσ′μe^{−Ry} − 2∫yμ²e^{−2Ry}; K₁″ = 2∫y²σ′μe^{−Ry} + 4∫y²μ²e^{−2Ry}; K₁‴ = −2∫y³σ′μe^{−Ry} − 8∫y³μ²e^{−2Ry};
  0 ≤ y ≤ 1 and Cauchy–Schwarz: |K₁‴(R)| ≤ 2√S₁√Y(R) + 8Y(R) ≤ G₁ := 2√S₁√Y(R_lo) + 8Y(R_lo) (Y non-increasing). ✓
- Taylor with Lagrange remainder on [R_c − δ, R_c + δ] (remainder points ≥ R_lo): K₀ ≥ K₀(R_c) + K₀′ε − F₀ε²/2,
  K₁ ≤ K₁(R_c) + K₁′ε + K₁″ε²/2 + G₁|ε|³/6; x ≥ 0 and a₂ := min(inf(−xK₁″(R_c)/2), 0) give
  D ≥ lower(ε) := K₀ + K₀′ε − x(K₁ + K₁′ε) − zK₂ + (a₂ − F₀/2)ε² − xG₁|ε|³/6 (all at R_c).
- lower is concave (ε² coefficient ≤ 0, −|ε|³ concave); target t(R) = e^{(1−κ_up)R} convex; lower − t concave ⇒ its minimum on
  [−δ, δ] is at an endpoint ⇒ checking ε = ±δ (in ball arithmetic: inf(lower) > sup(t)) proves M ≥ t, i.e.
  1 − log M(R,λ)/R ≤ κ_up on the whole cell. ✓
- Crude cells: M(R,λ) ≥ crude(R, C_λ) (§4.3), crude increasing in R and C, t increasing ⇒ crude(R_a, C_{λ_a}) > t(R_b) suffices. ✓

## 4. Tails — proofs checked

4.1 R = 0 (closed form, all λ ≥ Λ). Q = β/2 + g(y − 1/2), g odd, g(±1/2) = ±(β/2 − 1); ∫Q² = β²/4 + ∫g², ∫Q′² = ∫g′² (∫g = 0).
  min_g α∫g² + γ∫g′² at g = (β/2 − 1)sinh(ωt)/sinh(ω/2), ω = √(α/γ); value γ[gg′]_{−1/2}^{1/2} = 2√(αγ)coth(ω/2)(1 − β/2)²,
  √(αγ) = √(C_λB_λ) ≥ 1/2 (1.4), coth ≥ 1. So J ≥ (1 + (β−1)²)/2 + (C_λ/(4ξ))β² + (1 − β/2)² = 2 − 2β + β²(3/4 + C_λ/(4ξ)),
  min over β: M(0,λ) ≥ m₀(Λ) := 2 − 1/(3/4 + C_Λ/(4ξ)) (C_λ ≥ C_Λ). Λ = 0, ξ = 25/256: m₀ = 2 − 4ξ/(3ξ+1) = 1.6978852. ✓
4.2 Slope. For fixed Q, λ: dJ/dR = e^{2R}Q(1)² + 2α∫ye^{2Ry}Q² + 2γ∫ye^{2Ry}(RQ+Q′)² + 2γ∫ww′ ≥ −2γ√(W₀W₁)
  ≥ −√(γ/α)(αW₀ + γW₁) ≥ −sJ, s = ξ√(B_λ/C_λ) ≤ ξ√(B_Λ/C_Λ) for λ ≥ Λ (≤ ξ/√3 for Λ = 0). Grönwall: e^{sR}J_R(Q) ≥ J_0(Q) ≥ M(0,λ);
  inf over Q (class independent of R): M(R,λ) ≥ m₀e^{−sR} ⇒ 1 − log M/R ≤ 1 + s − log m₀/R ≤ κ_up for R ≤ R_tail = log m₀/(1 + s − κ_up).
  Numbers recomputed by hand: ξ = 25/256: R_tail = 0.529384/0.986326 = 0.53672 > 1/2 = core start ✓; Λ = 3: C_3 = 1.5523,
  B_3 = 0.16251, s = 0.0316, m₀ = 1.78830, R_tail(3) = 0.6045 ✓.
4.3 Crude bound. Drop ξB_λW₁ ≥ 0; for y ≤ 1/2, min_q e^{2Ry}q² + e^{2R(1−y)}(β − q)² = β²e^R/(2cosh(R(1−2y))), so
  W₀ ≥ β²(e^R/2)∫₀^{1/2}sech(R(1−2y))dy = β²e^R gd(R)/(4R); J ≥ 1/2 + e^{2R}(β−1)²/2 + Kβ², K = (C_λ/ξ)e^R gd(R)/(4R);
  min_β = 1/2 + e^{2R}K/(e^{2R} + 2K) = crude(R, C_λ). crude = 1/2 + 1/(1/A + 1/K), A = e^{2R}/2: increasing in K, hence in C and
  in R iff e^R gd(R)/R is: d/dR log = 1 + sech R/gd R − 1/R; R ≥ 1 trivial; R < 1: gd R ≤ R ⇒ ≥ (R − 1 + sech R)/R > 0 since
  (1 − R)cosh R < 1 on (0,1) (its derivative −cosh R + (1−R)sinh R < 0). ✓
  (i) cells and (ii) the λ ≥ 3 strip on [0.59375, 28 + 1/64] (C_λ ≥ C_3): crude(R_a, C) > t(R_b) per interval ✓ (replayed).
  (iii) R ≥ RMAX, every λ ≥ 0 (C = 1 ≤ C_λ; the reviewer's point (i) — correct and needed: K ≤ e^R must hold uniformly, and with
  C_λ unbounded in λ it would not): crude > K/(1 + 2Ke^{−2R}), so log crude ≥ log K − 2Ke^{−2R}; K = e^R gd(R)/(4Rξ) ≤ e^R since
  4Rξ ≥ 10.9 > π/2 > gd(R) (not asserted in code; true for both runs: 4·28·25/256 = 10.94, 4·22/8 = 11); gd(R) ≥ gd(RMAX):
  1 − log M/R ≤ h(R) := [log(4Rξ/gd(RMAX)) + 2e^{−RMAX}]/R, decreasing for R > e·gd(RMAX)/(4ξ) (10.9 resp. 8.5);
  h(28) = 0.0693077 < 0.07005633, h(22) = 0.0884688 < 0.08971026 (hand: log(6.963)/28 = 0.06931 ✓).
4.4 Coverage of (0,∞) × [0,∞): R ≤ 0.5367 (4.2, Λ = 0, all λ); [1/2, 28] × [0, 3] (cells, §5); λ ≥ 3: R ≤ 0.6045 (4.2, Λ = 3),
  [0.59375, 28 + 1/64] (strip); R ≥ 28 (4.3 iii, all λ). Union = everything. λ → 0 is not a limit case: λ = 0 (P = x, C = 1,
  B = 1/3) lies in the box, CB(0) is exact. ✓

## 5. Code audit (replay.py, covcheck.py against §2–4; read line by line)

replay.py (python-flint, 320 bits):
- σ, μ: float proposals converted exactly to rationals (Fraction(float)); σ piece k's constant reset to piece k−1's value at x = 1
  (exact continuity ⇒ σ ∈ H¹); μ built only on [1/2,1] and reflected by p(1−x) on the mirror piece KP−1−k (y ↦ 1−y maps piece k,
  local x to piece KP−1−k, local 1−x) ⇒ μ exactly symmetric. Lemma D needs nothing else, so the untrusted build() cannot break it.
- Integrals: ∫₀¹e^{ax}G = e^aΣ(−1)^kG^{(k)}(1)/a^{k+1} − Σ(−1)^kG^{(k)}(0)/a^{k+1} (exact rationals; a = −R_c/KP etc. rational);
  piece factors E[k] = e^{−R_c k/KP} in Arb; dy = dx/KP and σ′ = KP·dσ/dx accounted for (S₁ gets ×KP, X no factor, Y gets /KP).
  Signs: X′ = −∫yσ′μe^{−Ry} (sgn −1), X″ (+1), Y′ factor −2, Y″ factor +4, Y(R_lo) with a = −2R_lo/KP. Matches §3. ✓
- K₀, K₀′ = tMe^{−R}, F₀, K₁ = S₁+2X+Y, K₁′ = 2X′+Y′, K₁″ = 2X″+Y″, G₁ = 2√S₁√Y(R_lo) + 8Y(R_lo), K₂: match §3. ✓
- check(): subcell (i,j) uses x = ξ/C_{λ_i}, z = 1/(ξB_{λ_j}) (C_λ, B_λ closed forms re-derived: C = λ²(1/2 + sinh2λ/(4λ))/sinh²λ ✓);
  q2 = min(lower(−xK₁″/2), 0); lo(∓δ) = base ∓ K₀′δ ± xK₁′δ with base = K₀ − xK₁ − zK₂ + q2δ² − F₀δ²/2 − xG₁δ³/6; accepted only
  if both arb comparisons are certain (python-flint '>' is the certain predicate); else bisect on the integer λ-grid, fail at width 1.
- crude(): gd = 2 atan tanh(R/2), K, 1/2 + eK/(e + 2K): matches §4.3; crude cells test crude(R_c − δ, C_{λ_a}) > t(R_c + δ). ✓
- tails: R_tail > R0 (core start, min over cells), h(RMAX) < κ_up, RMAX > e·gd(RMAX)/(4ξ); strip from ⌊64·max(R_tail, R_tail(Λ))⌋/64
  to RMAX + 1/64 by crude(a, C_Λ) > t(b), with the start ≤ R_tail(Λ) checked. Matches §4. Unasserted but trivially true:
  gd(R) ≤ 4Rξ for R ≥ RMAX (K ≤ e^R), and R1 = max cell edge = RMAX (covcheck prints it: 28.0 / 22.0).
- Only shared code with ceil.py: build() (floats, untrusted) and the JSON cell list (geometry, checked by covcheck).
covcheck.py: exact Fractions; every cell inside [R0, R1] × [0, LAM/HL]; Σ areas = box area and no two cells have overlapping
  interiors ⇒ the closed cells cover the box (finite union of closed rectangles of full measure in the box is the box). ✓
No error found in either file.

## 6. Falsification (own code: myceil.py, cplx.py, cplxglobal.py, spot.py; logs/)

6.1 Float ceiling, own minimiser (Chebyshev Ritz in the real class, Q = 1 + Σ_{k≤K} b_k[T_{2k−1}(1−2y) − 1], exact quadratic
  minimisation in mpmath 30 digits, 90-node Gauss–Legendre; C_λ, B_λ by quadrature of P_λ, not the closed forms; alternating
  P/Q descent to the fixed point λ = ξ√(W₁/W₀); golden section in R). ξ = 25/256 at R = 9.867773:
  K = 12, 18, 24, 30, 36, 42 → 0.06953162, 0.07004450, 0.070056052, 0.0700563108, 0.07005631656, 0.07005631669 (λ = 0.770589037);
  geometric convergence (ratio ≈ 45 per 6 K) → limit ≈ 0.0700563167 = gl2ceil's float value; κ(R) is maximal at R = 9.867773
  (κ(R ± 0.005) = 0.07005630283, symmetric; offset of the maximum ≈ 2e-9). Ramanujan ξ = 1/8, R = 7.688073, K = 36:
  0.08971025099 (gl2ceil 0.0897102510). Gaps to the certified values: 1.3e-8 and 0.9e-8. No admissible choice found above them.
6.2 Q of high degree: K = 42 is degree 83; the Ritz values increase monotonically to the limit above, which is the infimum over the
  H¹ class (density, DERIVATION 1.3) — degree cannot push past the float ceiling, which the certificate bounds for all Q ∈ 𝒬̄.
6.3 Jointly complex (P, Q) — the one freedom the reduction does not cover (requires extending B's theorem to complex P; §1.6).
  With complex P the inf over P is F_c = (1 + |w(1)|²)/2 + X coth(ξX/W₀), X² = W₀W₁ − (Im B)² (derived: Euler–Lagrange
  −(W₀/ξ)P″ + ξW₁P + 2i(Im B)P′ = 0, P = e^{iξ(ImB)(x−1)/W₀} sinh(ξXx/W₀)/sinh(ξX/W₀); this is B l.796 / AT l.996).
  (a) Local, exact: at the real optimum the imaginary directions decouple from the real ones (F_c even in Q_i), and the second
  variation in Q_i = Σ_{k≤KE} d_k[T_{2k}(1−2y) − 1] is positive definite iff t := (f_X/(2X))·ℓᵀS⁻¹ℓ < 1 (ℓ the Im B coupling,
  S the Gram cost). t = 0.23283 (ξ = 25/256, K = KE = 16 and 36) and 0.23297 (ξ = 1/8): no descent direction; the coupling
  −2τ Im B reaches only 23% of the second-order cost. Line search along the best direction: κ_c decreases (logs/cplx.log).
  (b) Global, float64: BFGS from 12 random complex starts (|d| scales 1e-3 … 1), K = KE = 8, R ∈ {4, 7, 8.6, 10}: every run
  returns to d ≈ 0 (|d| ≤ 2e-7), κ_c = κ_real to 1e-11 (logs/cplxglobal.log). So the natural sesquilinear extension does not
  beat the ceiling numerically either; it is not certified (the certificate's class is real P, which is what B and AT state).

6.4 The papers' inf_P formula, settled numerically (cplxP.py, logs/cplxP.log): R = 3, ξ = 25/256, Q = 1 − 0.9[(1−2x) − 1]
  − 0.35i[(1−2x)² − 1] + 0.2[(1−2x)³ − 1] (Bernard's shape; Im B = 160.97): direct double-quadrature minimum over real P
  (degree 9) 4722.9017460952 = √(AC)coth formula 4722.9017460951; over complex P 4720.4201831735 = the X-formula of B l.796 /
  AT l.996 to 10 digits. So that formula needs complex P; with the hypotheses as stated it overstates the method for complex Q.
6.5 Spot check of the certificate against my own minimiser (spot.py, logs/spot_25_256.log): σ, μ from ceil.build for the worst-
  margin cell (R = 9.8618 ± 0.00049, λ ∈ [0.7656, 0.7734]) and 4 random dual cells, D(R,λ) by own Gauss quadrature, M(R,λ) by
  myceil (K = 36), 5 points per cell: target ≤ D ≤ M at all 25 points (min (D − t)/t = 3.0e-7; (M − D)/M ≥ 1.3e-9 at the centres,
  ≥ 0 everywhere): weak duality respected and the certified inequality holds pointwise. An assembly or sign error common to ceil.py
  and replay.py would have shown D > M or D < t; none.
6.6 Reruns on copies in /tmp (memcap 290 MB): replay.py 25/256: REPLAY OK, 1257 dual + 19 crude, 0 failures, tails and strip OK,
  250 s, rss 88 MB, memcap peak 59 MB; replay.py 1/8: REPLAY OK, 1014 dual + 3 crude, 0 failures, 214 s, rss 87 MB, peak 58 MB;
  covcheck.py: COVER OK for 25/256 (1276 cells, R ∈ [0.5, 28]), 1/8 (1017, [0.5, 22]), ζ 4/7 (526, [0.125, 6.125]).
6.7 Reviewer (verdict-only, both items pasted inline with the numbers of 6.4): ITEM 1 (class reduction: real P ⇒ √(AC); the
  papers' X-formula is the complex-P infimum; complex Q dominated by its real part) CORRECT; ITEM 2 (concave lower bound vs convex
  target, endpoint test) CORRECT.
6.8 Arb semantics replay.py relies on (flintsem.py, logs/flintsem.log): '>' between overlapping balls is False (certain predicate),
  .lower() is an exact point, fmpq → arb encloses the rational, exp balls have radius < 1e-90 at 320 bits. ✓
6.9 Off-box scan (myscan.py, logs/myscan.log; DPS = 60): κ(R) = −1.673 (R = 0.3), −0.471 (1), −0.0736 (3), 0.0517 (6), 0.0630 (15),
  0.0533 (20), 0.0377 (28), 0.0218 (35); optimal λ* ≤ 1.99. ζ (ξ = 4/7, R = 1.291586): 0.40895932599 < 0.40896.

## 7. Verdict

Every reduction from Bernard's theorem (at AT's lengths) to sup_{R,λ}[1 − (1/R)log M(R,λ)] is proved above; Lemma D, the cell
bound, the tails and the coverage logic are proved; the two certificate checkers implement them correctly and replay cleanly;
falsification (own float ceiling, high degree, off-box R, endpoint ν, complex Q, complex P, jointly complex numerically) found no
admissible choice above the certified values. AUDIT PASSED, with scope: real P as stated by B and AT; jointly complex (P, Q) only
numerically; the papers' displayed inf_P formula is the complex-P one (finding 2 of REPORT.md).
