# DERIVATION — gl2two (two-piece Levinson mollifier for GL(2), nodes N8, N12)

Sources (TeX line numbers): KRZ = Kühn–Robles–Zeindler arXiv 1609.03738 (KRZ02_post_referee_new.tex; Math. Z. 291 (2019));
BCY = Bui–Conrey–Young arXiv 1002.4127 (Acta Arith. 150 (2011)); AT21 = Andersen–Thorner arXiv 2004.03581 (zeros-paper.tex).
Degree2Mol, GL2Ceil = Riemann/Deep/Degree2Mol, Riemann/Deep/GL2Ceil.

## 1. Mollifier, class, zero count (KRZ)

- f: holomorphic primitive cusp form, even weight, squarefree level N, trivial character (KRZ l.73–100). Narrower than AT21's class
  (holomorphic or Hecke–Maass newforms, any level, AT21 l.138).
- ψ = ψ₁ + ψ₂ (KRZ l.189–208, l.239; BCY l.225–247): ψ₁(s) = Σ_{h≤M₁} μ_f(h)h^{σ₀−1/2−s} P₁[h], ψ₂(s) = χ_f(s+1/2−σ₀) Σ_{hk≤M₂} μ_{f,2}(h)λ_f(k)
  h^{σ₀−1/2}k^{1/2−σ₀} h^{−s}k^{s−1} P₂[hk], μ_{f,2} = μ_f∗μ_f, M_i = T^{ν_i}, 0 < ν₂ ≤ ν₁ (l.208), σ₀ = 1/2 − R/log T.
- Conditions: Q(0) = 1, Q(x) + Q(1−x) = const (l.172); P₁(0) = 0, P₁(1) = 1; P₂(0) = P₂′(0) = P₂″(0) = 0 (l.208, conditiononP l.230; BCY l.243).
- Zero count (kappaineq l.260): κ_f ≥ limsup (1 − (1/(2R)) log((1/T)∫₁ᵀ|Vψ(σ₀+it)|²dt)), V = Q(−(1/(2 log T)) d/ds) L(f,s) (l.174–176).

## 2. The functional

KRZ Thm unsmoothed (l.331–343): (1/T)∫₁ᵀ|Vψ|² = c + o(1), c = Σ_{k₁+k₂=2} multinomial ∏ c_{ℓ,ℓ+k_ℓ}^{k_ℓ}, i.e. for two pieces
c = c₁₁ + 2c₁₂ + c₂₂ (= BCY Thm thm:sharp l.263–268, c = c₁ + 2c₁₂ + c₂; KRZ l.382, BCY l.390: c₁₂ real).
Term by term, KRZ's c_{ℓ,ℓ} (cellellformula l.392–401, ℓ = 2) and c_{ℓ,ℓ+1} (cellell1formula l.415–421, ℓ = 1) are BCY's eq:c2
(l.376–389) and eq:c12 (l.360–366) with θ_i = ν_i/2 (prefactors 2/3 and 4, (ν₂/ν₁)² = (θ₂/θ₁)², 2/ν₂ = 1/θ₂, P₂^{(2)} = P₂″).
KRZ §5 (l.1745–1749): the mean value "produces the constants c₁,₁, c₁,₂ and c₂,₂ with 2R and ν₁/2, ν₂/2". So, in ζ-normalised variables
R′ = 2R, θ_i = ξ_i = ν_i/2, the GL(2) functional is BCY's, and κ_f ≥ 1 − log c(P₁,P₂,Q,R′,θ₁,θ₂)/R′:

  c₁ = 1 + (1/θ₁)∫₀¹∫₀¹ e^{2Rv}(Q(v)P₁′(u) + θ₁Q′(v)P₁(u) + θ₁RQ(v)P₁(u))² du dv                         (BCY eq:c1 l.344–346)
  c₁₂ = 4(θ₂/θ₁)² e^R ∂²_{xy}[∫∫_{a,b≥0, a+b≤1}∫₀¹ u²(1−u) e^{R[θ₁(y−x)+uθ₂(a−b)]} Q(−xθ₁+auθ₂) Q(1+yθ₁−buθ₂)
        × P₁(x+y+1−(1−u)θ₂/θ₁) P₂″((1−a−b)u) du da db]_{x=y=0}                                            (BCY eq:c12 l.360–366)
  c₂ = (2/3) ∂⁴_{x²y²}[∫_{[0,1]⁴} (1−r)⁴(1/θ₂ + S) e^{−θ₂RS} Q(θ₂(−y+u(x+r)) + t(1+θ₂S)) e^{2Rt(1+θ₂S)} Q(θ₂(−x+v(y+r)) + t(1+θ₂S))
        × (x+r)(y+r) P₂″((1−u)(x+r)) P₂″((1−v)(y+r)) dt dr du dv]_{x=y=0},  S = x+y−v(y+r)−u(x+r)          (BCY eq:c2 l.376–389)

Normalisation check. KRZ's remark l.403–407 writes c₁,₁ with ν₁ where its own general formula has ν₁/2; the consistent reading is the one above:
(i) for ψ₂ = 0 it is Bernard's mean value c(P,Q,2R,ν/2) (AT21 eq. limit-cPQR l.205–212); (ii) it reproduces KRZ's printed results
0.0297607 and 0.0693872 to all digits (§5). Validation of the BCY formulas themselves: §5.

P₂ enters only through g = P₂″; P₂(0) = P₂′(0) = 0 makes g ↔ P₂ one-to-one; g(0) = 0 (P₂″(0) = 0) is imposed by the basis
φ_k(s) = s·T_{k−1}(2s−1). c₁₂ is linear and c₂ quadratic in g: c = c₁ + 2gᵀL + gᵀMg, M ⪰ 0 (limit of mean values of |Vψ₂|²).
For fixed (P₁, Q, R): min_g c = c₁ − LᵀM⁻¹L, so the first-order gain is Δκ = −log(1 − LᵀM⁻¹L/c₁)/R′.

## 3. Boundary-layer form of c₁₂ (why the second piece is weak)

With w(y) = e^{Ry}Q(y), α = auθ₂, β = 1 − buθ₂: Q(α) = e^{−Rα}w(α), RQ(α)+Q′(α) = e^{−Rα}w′(α), and e^R e^{Ruθ₂(a−b)}e^{−Rα}e^{−Rβ} = 1. Hence

  c₁₂ = 4(θ₂/θ₁)² ∫∫_{a+b≤1}∫₀¹ u²(1−u)[−θ₁²w′(α)w′(β)P₁(z) + θ₁(w(α)w′(β) − w′(α)w(β))P₁′(z) + w(α)w(β)P₁″(z)] g((1−a−b)u),

α ∈ [0, θ₂], β ∈ [1−θ₂, 1], z = 1 − (1−u)θ₂/θ₁ ∈ [1 − θ₂/θ₁, 1]. The cross term sees w only on two boundary layers of width θ₂, where
w(0) = 1 and w(1) = e^R Q(1) = O(1) (−0.46 at the unconditional optimum). At x = y = 0 the c₂ integrand is w(a₁)w(a₂)·(polynomial),
a₁, a₂ = t(1 − θ₂r(u+v)) + θ₂·(·) sweeping all of [0,1]: c₂ sees the bulk ∫w², which is ≍ θc₁ with c₁ = e^{R′(1−κ)}.
Short mollifiers (θ ≈ 0.1) force R′ ≈ 8–10, c₁ ≈ 10³–10⁴; then the leading components of L are 2·10⁻³–5·10⁻³ (per unit coefficient of g),
M₀₀ ≈ 1.4·10²–2.3·10³ (logs/resol.log), and LᵀM⁻¹L/c₁ = 1.2·10⁻¹³ (5/27, 25/149), 1.6·10⁻¹⁰ (1/4, 1/4), 2.7·10⁻¹³ (25/128, 25/149).
At ζ's θ₁ = 4/7 (R = 1.28, c₁ = 2.13) the same ratio is 2.2·10⁻³ (BCY's +0.0017). Along the R scan LᵀM⁻¹L/c₁ ∝ c₁^{−2.7…−3.5}.
Diagnostic, not a uniform bound.

## 4. Admissible lengths

Explicit theorem statements:
- KRZ Thm unsmoothed (l.331–336): ν₁ = (1−2θ)/(4+2θ) − ε, ν_ℓ = (1−2θ)/(4+6θ) − ε (ℓ ≥ 2); also l.314–322. θ = 7/64: ν₁ < 5/27, ν₂ < 25/149;
  θ = 0: ν₁, ν₂ < 1/4. KRZ's θ is defined by |λ_f(n)| ≤ τ(n)n^θ (kimsarnak l.142–144) and l.320 / l.1761 say θ = 0 is proved (Deligne);
  but the θ in (1−2θ)/(4+2θ) comes from Bernard's shifted-convolution bound, where it is the Kim–Sarnak exponent of the Maass spectrum
  (AT21 l.132, l.168–173: Bernard's 0.0693 is "under the generalized Ramanujan conjecture for Hecke–Maaß forms"). KRZ's Thm theoremofkappa
  (l.293–300) states 6.93872% "under the Ramanujan conjecture". We read θ = 0 as conditional on GRC for Maass forms (as for the one-piece record).
- AT21 Thm thm:nu (l.227–231): ν < 1/4 − θ/2 for the one-piece mean value (mollifier eqn:mollifier_def l.191).

Where each bound comes from:
- I₁₁ = ∫|Vψ₁|² (KRZ l.822–832): per pair (h₁, h₂) the off-diagonal is O((h₁h₂)^{(1+θ)/2}T^{1/2+θ+ε}) (Bernard, KRZ l.531, l.613–617); summed
  trivially: T^{1/2+θ+ν₁(θ+2)} → ν₁ < (1−2θ)/(4+2θ). AT21's Prop. prop:off-diag (l.317–326) improves the per-pair bound to (ab)^{1/2}T^{1/2+θ}(abT)^ε
  for ANY a, b; the coefficients enter only in Cor. cor:off_diag (l.328–345) via |μ_f(a)| ≪ max(1,|λ_f(a)|) and Rankin–Selberg. I₁₁ is literally
  AT21's one-piece mean value with P = P₁, so AT21's theorem gives I₁₁ for ν₁ < 25/128 (θ = 7/64).
- I₂₂ (KRZ l.838–860): the same per-pair bound with |λ_f^{∗(ℓ−1)}(k)| ≪ k^{θ+ε} → ν_ℓ < (1−2θ)/(4+6θ). With AT21's per-pair bound and Rankin–Selberg
  on average one would expect ν₂ < 25/128; nobody states it.
- I₁₂ (KRZ l.1216–1320, l.1388–1420): off-diagonal needs only ν₁ + ν₂ < 1 (Remark a l.347; l.1276–1282, l.1315–1317: "the condition
  ν_ℓ+ν_{ℓ+1}<1 is needed only for the bound of w₀"); the diagonal part uses ν₁ + ν₂ < 2 (l.1392) and ν₂ ≤ ν₁ (l.1418).

Tiers:
1. Explicit statements: (ν₁, ν₂) < (5/27, 25/149) unconditionally; (1/4, 1/4) under GRC; KRZ's class of f.
2. Consequence of published estimates (AT21 for I₁₁ + KRZ's I₁₂ and I₂₂ proofs, whose hypotheses are ν₁+ν₂ < 1, ν₂ ≤ ν₁, ν₂ < 25/149):
   (ν₁, ν₂) < (25/128, 25/149) unconditionally, KRZ's class. Not a stated theorem.
3. Plausible, unpublished: (25/128, 25/128) (AT21's per-pair bound inserted into KRZ's I₂₂).
Tier 1 is "what the published theorems justify". The gains at all three tiers are in REPORT.md.

## 5. Validation (logs/val_bcy*.log, logs/xcheck.log)

- BCY, θ₁ = 4/7, θ₂ = 1/2, R = 1.28 (l.273–284): printed Q has Q(0) = 1.002; rescaled to Q(0) = 1: κ = 0.4105121 (BCY .4105; Q(0)-shift instead
  of rescaling: 0.4104754); one-piece part 0.4087809 (Conrey .4088). The minimiser over P₂ ∈ span{x³, x⁴, x⁵} is 0.0245412x³ − 0.00635566x⁴ +
  0.0060313x⁵ = BCY's printed P₂ (l.283) to all printed digits.
- BCY simple zeros, R = 1.12, Q = 1 − 1.03x (l.285–294): the minimiser in the same span is 0.0323061x³ − 0.00553783x⁴ + 0.0076959x⁵ with
  κ* = 0.4058028 (BCY .4058). BCY print .0769594x⁵ (l.294): a factor-10 misprint (printed value gives 0.3934573 < one piece 0.4021095).
- KRZ printed parameters (l.1750–1772), R′ = 2R, θ = ν/2, Q rescaled to Q(0) = 1: 0.0693872 (θ = 0; KRZ .0693872), 0.0297607 (θ = 7/64;
  KRZ .0297607). As printed (Q(0) = 0.999999, 1.000006): 0.0693876, 0.0297588. Their P₂ (coefficients ~10⁻⁵) adds 7.8·10⁻⁹ and 1.6·10⁻⁸.
- Derivative algebra: xcheck.py codes BCY's integrands literally (complex x, y; ∂ by Cauchy's formula on |x| = |y| = 0.3, 16 × 16 points,
  same quadrature): c₁₂, c₂ agree with the jet implementation to 7·10⁻¹⁴ (BCY) and 3·10⁻¹⁴ (KRZ parameters).
- Quadrature and basis: first-order gains identical to all printed digits at (n₁₂; n_t, n_r, n_u) = (40; 120, 24, 24), (60; 160, 32, 32),
  (80; 200, 40, 40); K = 4, 8, 12 terms for g agree to 7 digits (logs/resol.log, logs/gain.log).

## 6. Complex Q

KRZ allow complex Q (l.172). V = Σ_k q_k V_k, V_k = (−1/(2 log T))^k L^{(k)}(f,s), so (1/T)∫|Vψ|² → q*Aq with A_{kl} = lim (1/T)∫V_kψ·conj(V_lψ):
Hermitian and positive semidefinite (limit of Gram matrices). A is real symmetric: with real P₁, P₂ and real coefficients (trivial character:
λ_f, μ_f real; χ_f(s̄) = conj χ_f(s)), each main term is real for real basis polynomials (KRZ l.382, BCY l.390), and the two cross terms satisfy
I₂₁(α,β) = conj I₁₂(β̄,ᾱ), so they enter A as c₁₂^{kl} + c₁₂^{lk}. Hence for q = a + ib (a₀ = 1, b₀ = 0 since Q(0) = 1): q*Aq = aᵀAa + bᵀAb ≥ aᵀAa,
i.e. c(Q) ≥ c(Re Q). Complex Q never helps (as for one piece, Degree2Mol §3, GL2Ceil §1.1).

## 7. Numerics (untrusted float; tp.py, joint.py)

- c₁ by (E) of Degree2Mol §2 (1D Gauss–Legendre). c₁₂: 3D Gauss–Legendre (simplex a = p, b = (1−p)q, Jacobian 1−p), ∂²_{xy} by expanding
  the x, y factors (A₀B₀P₁″ + (A₁B₀ + A₀B₁)P₁′ + A₁B₁P₁). c₂: 4D Gauss–Legendre; every factor is f(α + βx + γy) or (linear)·exp(linear),
  so the integrand's truncated bivariate Taylor jet (degree ≤ 2 in x and in y) is exact; ∂⁴ = 4·[x²y²]. t is integrated first:
  M_{kl} = (8/3)Σ w(1−r)⁴ Σ_{ij} Ŵ_{ij}(r,u,v) X_k[2−i](r,u) X_l[2−j](r,v), X_k the jet of (x+r)φ_k((1−u)(x+r)).
- Joint optimisation (joint.py): for fixed Q, c is a positive quadratic in (P₁, g) (P₁ = x + x(1−x)Σ_{j<12}p_jT_j(2x−1)); for fixed (P₁, g),
  c = 1/2 + vᵀHv in v = (1, h), H = H₁ + 2H₁₂ + H₂, with H₁₂, H₂ the Q-quadratic forms of c₁₂, c₂ (consistency vᵀH₁₂v = gᵀL, vᵀH₂v = gᵀMg
  checked to 10⁻¹⁷ each step). Alternating exact solves; R fixed at the one-piece R*; R scanned separately (rscan.py).
