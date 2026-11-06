# DERIVATION — primes between consecutive k-th powers from the Lean-checked half-plane and an explicit zero density (Ingham3)

Every step below is implemented in `cert.py` (Arb ball arithmetic, python-flint 0.9.0, 128 bits). Section numbers match the comments in the code.

## 0. Result

**Theorem.** (a) For every integer n ≥ 1 there is a prime p with n^7 < p < (n+1)^7.
(b) For every integer n with n^k ≥ x₀(k) there is a prime in (n^k, (n+1)^k), where

| k | log x₀(k) | n₀(k) = x₀(k)^{1/k} |
|---|---|---|
| 3 | 1881.03 | e^{627.01} ≈ 2.03·10^{272} |
| 4 | 468.84 | e^{117.21} ≈ 8.01·10^{50} |
| 5 | 288.73 | e^{57.746} ≈ 1.20·10^{25} |
| 6 | 217.45 | e^{36.242} ≈ 5.49·10^{15} |

Inputs: (I1)–(I4) of §1 only. No RH beyond T₀, no unrefereed zero count, no Chirre–Helfgott.

## 1. Inputs

- **(I1) Half-plane.** ζ(s) ≠ 0 for Re s > T* := B₂ + 2.2·10^{−254}, B₂ = (777 − √1081)/852 = 0.87338196661320751993… (J. Bebel, Lean-kernel-checked, QuasiRiemannTracker; TREE.md N7). We use **B = 0.8733819666132076**; cert.py checks B > T* in Arb. Hence every zero has β ≤ T* < B.
- **(I2) RH to height T₀ = 3·10^{12}.** Every zero with 0 < γ ≤ 3·10^{12} has β = 1/2 (Platt–Trudgian, Bull. LMS 53 (2021); verified to 3 000 175 332 800).
- **(I3) Zero count.** Hasanalizade–Shen–Wong (J. Number Theory 2022, arXiv 2107.06506, Cor. 1.2): |N(T) − (T/2π) log(T/(2πe))| ≤ 0.1038 log T + 0.2573 log log T + 9.3675 for T ≥ e. Used only through **N(t) ≤ t log t for t ≥ 14** (at t = 14: 0.1038 log 14 + 0.2573 log log 14 + 9.3675 < (1 − 1/(2π))·14 log 14, checked in Arb; both sides' difference increases in t) and N(t) = 0 for t < 14.
- **(I4) Explicit zero density.** Kadiri–Lumley–Ng, "Explicit zero density for the Riemann zeta function", J. Math. Anal. Appl. 465 (2018) 22–46 = arXiv 2101.12263, **Lemma 4.14** (TeX l. 1547–1580): for parameters α > 0, d > 0, δ > 0, η ∈ [η₀, 1/2), μ ∈ [1+η₀, 1+η], 10⁹/H₀ ≤ k ≤ 1, 1002 ≤ H ≤ H₀ (H₀ = 3.0610046·10^{10}, RH(H₀) by Platt) with U(α,k,H₀) > 1 and U(α,k,T) decreasing in T, and every σ > 1/2 + d/log H₀, T ≥ H₀:
  N(σ,T) ≤ (C₁/(2πd)) (log kT)^{2σ} (log T)^{5−4σ} T^{8(1−σ)/3} + (C₂/(2πd)) log²T,
  N(σ,T) = #{ρ: σ < β < 1, 0 < γ < T}. This is their (4.71); C₁ = C₁(α,d,δ,k,H,σ) (their (4.72)) and C₂ = C₇(η,H)(μ − σ + d/log H₀) + C₈(k,μ) (their (4.73)) are explicit (equation numbers checked in the arXiv PDF, p. 20; §4.5 lists where each formula sits in cert.py); cert.py implements every sub-definition: b₁ = 0.62, b₂ = 1.048, b₃ = 0.605, b₄ = 0.529, m₀ = (1 + (2/3)√(6/5))^{1/2}, a₁ = 0.63, a₂ = 2.851, C₁(k), C₂(k), a₃, C₃(k), C₄(k), I(A,n) = ∫₀^∞ x^A e^{−2αx²}(log x)^n dx (closed form via Γ, ψ, ψ′), J(k,T), U = 8αC₄ω₁(1/2,T,α)²J, C₅, C₆, K, V = 8αω₁(σ₂,T,α)²K, M(k,δ), b₁₂(H) = 1/(2(1−1/H)²), b₅ … b₁₁, C₇, C₈.
  - *Reproduction.* With the 20 parameter rows of KLN Table 1 (σ₀ = 0.60 … 0.99; k, μ, α, δ, d as printed), H = H₀ − 1 and η = 0.2561882 (the minimiser of C₇; η₀ = 0.2362121… from b₆(10⁹,η₀) = 1), our values of A = C₁/(2πd) and B = C₂/(2πd) agree with Table 1 within 0.002 in every row (log: logs/cert.log). We use our own Arb values, not the table.
  - *σ-dependence.* For fixed parameters log C₁ is affine in σ (each exponent in (4.72) is affine in σ) and C₂ is affine and decreasing in σ (C₇ > 0). So each row gives an explicit bound for **every** σ > 1/2 + d/log H₀, not only at σ₀. cert.py stores log A_P(σ) = c₀ + c₁σ (c₁ > 0 checked) and B_P(σ) = q₀ − q₁σ.
  - *Admissibility* (checked in Arb for every row): U(α,k,H₀) > 1; the T-dependent coefficients of J (those of 1/log T and of 1/(T^{1/6} log²T); the others are positive automatically) are positive, so U(α,k,T) decreases; V(α,k,δ,T) ≥ V_∞ := 8α(C₅I(1,0) + C₆I(2,0)/k) > 1; b₆(10⁹,η) < 1, b₆(10⁹,μ−1) < 1 (i.e. η, μ−1 > η₀), μ ≤ 1+η, k ∈ [10⁹/H₀, 1].
  - *Inputs inside KLN.* The subconvexity bound |ζ(1/2+it)| ≤ a₁t^{1/6} log t (t ≥ 3) with a₁ = 0.63 was Hiary's (2016), which rested on an error in Cheng–Graham (Johnston–Yang, arXiv 2204.01980, TeX l. 250). Hiary–Patel–Yang (J. Number Theory 256 (2024)) prove 0.618 t^{1/6} log t for t ≥ 3, and 0.618 < 0.63, while |ζ(1/2+it)| ≤ 1.461 on [0,3] gives a₂ = 2.851 as in KLN; so KLN's Lemma 3.2 and Table 1 stand. Ramaré's b₁, b₂ were rechecked numerically to X = 10⁷ (check.py: max 0.609833 ≤ 0.62 and 1.047074 ≤ 1.048).
  - *δ.* KLN Theorem 1.1 is printed with "δ ≥ 1", but the parameter conventions (TeX l. 626), Lemma 4.14 (l. 1549) and Table 1 itself use δ > 0 (δ = 0.303–0.314). The main run uses the Table 1 parameters; a second run with δ = 1 in every row (§9) keeps (a) and gives slightly larger x₀(k).

## 2. Weight and explicit formula

Fix k ≥ 3 and an integer n ≥ 2; put x = n^k, δ = k log(1 + 1/n), so xe^δ = (n+1)^k. Let m = 2 and Φ = 1_{[0,1/m]} * ⋯ * 1_{[0,1/m]} (m factors; for m = 2 the triangle Φ(v) = min(v, 1−v)₊). Φ ≥ 0 is Lipschitz, supported on [0,1], Φ(0) = Φ(1) = 0, ∫Φ = m^{−m}, ‖Φ‖_∞ ≤ m^{1−m}, so ‖Φ‖_∞/∫Φ ≤ m. Put w(u) = Φ(log(u/x)/δ) (u > 0): Lipschitz, supported on [x, xe^δ], zero at both ends.

**Explicit formula.** With W(s) = ∫₀^∞ w(u)u^{s−1}du,
  S_w := Σ_n Λ(n)w(n) = W(1) − Σ_ρ W(ρ) − ∫ w(u) du/(u³ − u).
Proof: S_w = −∫ψ(u)w′(u)du; insert ψ₀(u) = u − Σ_{|γ|<T} u^ρ/ρ − log 2π − ½ log(1 − u^{−2}) + R(u,T) (Davenport, Multiplicative Number Theory, ch. 17), whose remainder is bounded uniformly for u in the compact support and T ≥ 2 and tends to 0 off prime powers; dominated convergence and ∫(u^ρ/ρ)w′ = −W(ρ), −∫u w′ = W(1), ∫w′ = 0, ½∫log(1−u^{−2})w′ = −∫w/(u³−u). The zero sum converges absolutely (|W(ρ)| ≪ |γ|^{−m}, §3), so no truncation and no truncation error occur.
**Mellin transform.** u = xe^{δv}: W(s) = δx^s Φ̃(δs), Φ̃(z) = ∫₀¹ Φ(v)e^{zv}dv = ((e^{z/m} − 1)/z)^m. So W(1) ≥ δx∫Φ =: M₀ and |W(ρ)| = δx^β|Φ̃(δβ + iδγ)|.
**Trivial zeros.** 0 ≤ ∫w/(u³−u) ≤ W(1)/(x(x²−1)).
**Prime powers (ψ versus θ).** S_w = Σ_p log p·w(p) + PP, PP = Σ_{j≥2} Σ_p log p·w(p^j). If (n^k,(n+1)^k) has no prime then S_w = PP. Every p^j counted lies in (x, xe^δ), log p ≤ log(xe^δ)/j, w ≤ ‖Φ‖_∞, and (x^{1/j}, x^{1/j}e^{δ/j}) holds at most x^{1/j}(e^{δ/j}−1) + 1 ≤ x^{1/j}(δ/j)e^{δ/j} + 1 integers; j ≤ J := ⌊log(xe^δ)/log 2⌋. Hence
  PP/M₀ ≤ m log(xe^δ) [e^{δ/2}x^{−1/2}/4 + J e^{δ}x^{−2/3}/9 + log J/(δx)]   (Σ_{j=2}^J 1/j ≤ log J).
**Criterion.** S_w ≥ W(1)(1 − 1/(x³−x)) − Σ|W(ρ)| ≥ M₀(1 − 1/(x³−x) − Z) with
  Z := Σ_ρ x^{β−1}|Φ̃(δρ)|/∫Φ.
So **R := Z + PP/M₀ + 1/(x³−x) < 1 implies a prime in (n^k, (n+1)^k).** cert.py bounds 1/(x³−x) by 2x^{−3}.

## 3. Kernel bound

For z = s + iλ, 0 ≤ s ≤ δ ≤ δ_max ≤ m, a = 1/m:
|∫₀^a e^{zv}dv| ≤ |∫₀^a e^{iλv}dv| + ∫₀^a(e^{sv}−1)dv ≤ a|sinc(aλ/2)| + a(e^{δa}−1), and |e^{za}−1|/|z| ≤ (2|sin(aλ/2)| + e^{δa} − 1)/|λ| = a|sinc(aλ/2)| + (e^{δa}−1)/|λ|. With e^y − 1 ≤ 2y (0 ≤ y ≤ 1):
  |Φ̃(z)|/∫Φ ≤ G(|λ|),  G(λ) := (S(λ/2m) + 2δ_max min(1/m, 1/λ))^m,
S(y) := max(sin y/y, s₁) for 0 ≤ y ≤ π, min(s₁, 1/y) for y ≥ π, **s₁ = 0.2173** ≥ sup_{y≥π}|sin y|/y = 1/√(1+y*²) = 0.2172336… (tan y* = y*, y* = 4.4934); cert.py checks |sin y|/y < 0.2173 on 8000 subintervals of [π, 2π] and 1/(2π) < 0.2173. S ≥ |sinc| (sin y/y decreases on [0,π]) and S is non-increasing, so G is non-increasing; G(λ) ≤ ((2m + 2δ_max)/λ)^m for λ ≥ 2mπ. δ_max ∈ {10⁻⁹, 0.025, 0.15}; cert.py asserts δ < δ_max in every cell. (check.py: max of (|Φ̃(z)|/∫Φ)/G over 6·10⁴ random z is 0.9999999995.)
**Moments.** For 0 ≤ p < m and λ₀ ≥ 0 put Mom(p; λ₀) := ∫_{[λ₀,∞)} max(1,λ)^p d(−G(λ)) = max(1,λ₀)^p G(λ₀) + p∫_{max(1,λ₀)}^∞ G(λ)λ^{p−1}dλ. It increases in p and decreases in λ₀. cert.py: grid λ_i = r^i, r = 1 + 1/500, i ≤ N (λ_N ≥ 2m·10⁴); on [λ_i, λ_{i+1}] it uses G(λ_i) (upper bound, G non-increasing), so p∫_{λ_i}^{λ_{i+1}}Gλ^{p−1} ≤ G(λ_i)(λ_{i+1}^p − λ_i^p); beyond Λ = λ_N, p∫_Λ^∞ ≤ p(2m+2δ_max)^mΛ^{p−m}/(m−p). p runs over the grid {j/200: 1 ≤ j ≤ 300}; a needed p is rounded up to the grid. For λ₀ ≤ 1: Mom ≤ G(0) + p∫₁^∞, G(0) = (1 + 2δ_max/m)^m. For 1 ≤ λ₀ ≤ 30: Mom ≤ λ₀^pG(λ_i) + p∫_{λ_i}^∞ (λ_i ≤ λ₀ < λ_{i+1}). For λ₀ ≥ max(2mπ, λ_{i_max}): Mom ≤ (2m+2δ_max)^m λ₀^{p−m}·m/(m−p).

## 4. The zero sum

**4.1 Pairing.** Zeros are symmetric under ρ ↦ ρ̄ and ρ ↦ 1 − ρ̄, and none is real. By §3, Z ≤ 2Σ_{γ>0} x^{β−1}G(δγ). A zero with β > 1/2 pairs with 1 − β + iγ, and x^{β−1} + x^{−β} ≤ 2x^{−1/2} + (x^{β−1} − x^{−1/2}); x^{β−1} − x^{−1/2} = log x ∫_{1/2}^β x^{σ−1}dσ. Hence, with S := Σ_{γ>0}G(δγ) and F(σ) := Σ_{γ>0, β>σ}G(δγ) (Tonelli),
  Σ_{γ>0} x^{β−1}G(δγ) ≤ x^{−1/2}S + log x ∫_{1/2}^{B} x^{σ−1}F(σ)dσ,
since F(σ) = 0 for σ ≥ T* (I1). For 1/2 ≤ σ ≤ σ₀ := 0.515, F(σ) ≤ F(1/2) ≤ S/2 (off-line zeros pair at equal heights), so the σ ≤ σ₀ part is ≤ (x^{σ₀−1} − x^{−1/2})S/2, and
  **Z ≤ 2[x^{σ₀−1}S + log x ∫_{σ₀}^{B} x^{σ−1}F(σ)dσ].**
No RH is used here; RH(T₀) enters only through F (4.3).
**4.2 All zeros.** S = ∫_{(0,∞)} G(δt)dN(t) = ∫ N(t) d(−G(δt)) ≤ ∫_{14}^∞ t log t d(−G(δt)) (I3). For t₁ := 1/δ ≥ 14 and t ≥ 14, t log t ≤ t₁ log t₁·max(1, t/t₁)^{1+1/log t₁} (for t ≤ t₁ by monotonicity, for t ≥ t₁ by 1 + u/ℓ ≤ e^{u/ℓ}). With λ = δt and D := log(1/δ):
  **S ≤ δ^{−1}D·Mom(1 + 1/D; 0).**
**4.3 Zeros off the line.** For σ ≥ σ₀ > 1/2 every zero with β > σ has γ > T₀ (I2). With N_σ(t) := #{β > σ, T₀ < γ ≤ t} (N_σ(T₀) = 0), Stieltjes integration by parts gives F(σ) = ∫_{T₀}^∞ N_σ(t) d(−G(δt)) (G is continuous, so the value of N_σ at the countably many ordinates does not matter), and for t not an ordinate N_σ(t) ≤ N(σ,t) ≤ A_P(σ)(log t)^{5−2σ}t^{a} + B_P(σ)log²t (I4; log kt ≤ log t; a := 8(1−σ)/3). Put λ₀ := δT₀, ℓ₀ := log T₀, e := 5 − 2σ.
- *Case (ii), λ₀ < 1:* (log t)^e t^a = δ^{−a}(D + log λ)^e λ^a ≤ δ^{−a}D^e max(1,λ)^{a+e/D} for λ ≥ λ₀; so
  F(σ) ≤ A_P(σ)δ^{−a}D^e·Mom(a + e/D; λ₀) + B_P(σ)D²·Mom(2/D; λ₀).
- *Case (i), λ₀ ≥ 1:* t = T₀u, (ℓ₀ + log u)^e ≤ ℓ₀^e u^{e/ℓ₀}; so with p′ = a + e/ℓ₀
  F(σ) ≤ A_P(σ)T₀^aℓ₀^e λ₀^{−p′}·Mom(p′; λ₀) + B_P(σ)ℓ₀²λ₀^{−2/ℓ₀}·Mom(2/ℓ₀; λ₀).
In both cases log[δ^{−a}·(…)] = aD + e·Λ with Λ = log D (ii) or Λ = log ℓ₀ − log λ₀/ℓ₀ (i), since log T₀ − log λ₀ = D.
**4.4 σ-integration.** [σ₀, B] is cut into cells [s_j, s_{j+1}] of width 0.0025 (144 cells, the last ends at B). In each cell one KLN row P is used (the one minimising A_P at both ends; any admissible row is valid). On the cell, p(σ) ≤ p(s_j) (a and e decrease in σ) and Mom increases in p, B_P(σ) ≤ B_P(s_j), and
  log x·x^{σ−1}·A_P(σ)·exp(aD + eΛ) = log x·exp(E₀ + E₁σ),
  E₀ = −L + c₀ + (8/3)D + 5Λ, E₁ = L + c₁ − (8/3)D − 2Λ (L = log x),
an exponential in σ, integrated exactly: ∫_{s_j}^{s_{j+1}} e^{E₀+E₁σ}dσ = (e^{E₀+E₁s_{j+1}} − e^{E₀+E₁s_j})/E₁. The B-term uses log x ∫ x^{σ−1}dσ = x^{s_{j+1}−1} − x^{s_j−1}.

*Admissibility of the row in each cell.* KLN's bound holds for σ > s_min := 1/2 + d/log H₀ (log H₀ = 24.1445941). cert.py l. 214–231 (class Setup) builds the cells: for each cell [s, t] it admits a row only if the Arb comparison s > s_min is certain (l. 226; every σ in the cell is ≥ s), and among the admitted rows takes the one with the smallest A at both ends (l. 227–228). If no row is admitted, best stays None and l. 229 (best[1]) raises, so the run aborts. At the first cell, σ = 0.515, a row is admissible iff d < 0.015 log H₀ = 0.3621689; all 20 rows have d ≤ 0.341, so all are admissible there (largest s_min = 0.5141233). Rows used (H = H₀ − 1 = 30 610 045 999 and η = 0.2561882 in every row; α, d, δ, k, μ as printed in KLN Table 1):

| cells (σ-range) | row σ₀ | α | d | δ | k | μ | s_min (upper bound) | s − s_min at the lower end (lower bound) |
|---|---|---|---|---|---|---|---|---|
| [0.515, 0.6275] (first cell [0.515, 0.5175]) | 0.60 | 0.288 | 0.341 | 0.3140 | 0.5 | 1.251 | 0.5141232443 | 0.000877 |
| [0.6275, 0.675] | 0.65 | 0.256 | 0.340 | 0.3070 | 0.6 | 1.249 | 0.5140818 | 0.113 |
| [0.675, 0.7225] | 0.70 | 0.222 | 0.339 | 0.3040 | 0.8 | 1.247 | 0.5140404 | 0.161 |
| [0.7225, 0.7775] | 0.75 | 0.189 | 0.338 | 0.3030 | 1.0 | 1.245 | 0.5139990 | 0.209 |
| [0.7775, 0.8275] | 0.80 | 0.160 | 0.337 | 0.3030 | 1.0 | 1.245 | 0.5139576 | 0.264 |
| [0.8275, 0.8575] | 0.85 | 0.133 | 0.336 | 0.3030 | 1.0 | 1.245 | 0.5139162 | 0.314 |
| [0.8575, 0.8675] | 0.86 | 0.127 | 0.335 | 0.3030 | 1.0 | 1.245 | 0.5138747 | 0.344 |
| [0.8675, B] | 0.87 | 0.122 | 0.335 | 0.3030 | 1.0 | 1.245 | 0.5138747 | 0.354 |

Each row passed the admissibility checks of §1 (cert.py l. 122–125, asserted at l. 141, or l. 134 in the δ = 1 run) and reproduces its Table 1 values A, B (l. 138–141). The δ = 1 run (§9) uses the same eight rows with δ = 1 on [0.515, 0.6425], [0.6425, 0.6925], [0.6925, 0.74], [0.74, 0.7825], [0.7825, 0.8325], [0.8325, 0.86], [0.86, 0.87], [0.87, B]; its first cell again uses row 0.60.

**4.5 Where each KLN formula is implemented** (KLN equation numbers as in the arXiv PDF of 2101.12263; cert.py line numbers)

| KLN | quantity | cert.py |
|---|---|---|
| (3.5) | b₁ = 0.62, b₂ = 1.048, b₃ = 0.605, b₄ = 0.529 | l. 32 |
| (3.4), Lemma 3.2 | a₁ = 0.63, a₂ = 2.851 | l. 33 |
| (3.15) | m₀ | l. 34 |
| Lemma 4.14 | H₀ = 3.0610046·10^{10}, H = H₀ − 1 | l. 31 |
| (3.18) | ω₁ (ω₂ of (3.19) enters only through b₁₂) | l. 62 |
| (4.3), (4.4) | C₁(k), C₂(k) | l. 44, l. 45–46 |
| (4.6), (4.7), (4.8) | a₃, C₃(k), C₄(k) | l. 47, l. 49, l. 48–49 |
| (4.10), (4.11) | C₅(k,δ), C₆(k,δ) | l. 64–65, l. 66–67 |
| (4.15)–(4.16) | I(A,n), β = 2, via Γ, ψ, ψ′ | l. 38–43 (the comment there says (4.17); the formula is (4.16)) |
| (4.17) | J(k,T); its six coefficients (the positivity checks use those of 1/log T and 1/(T^{1/6} log²T)) | l. 59–61; l. 50–58 |
| (4.18) | U(α,k,T) = 4αβC₄ω₁(1/2,T,α)²J | l. 63 |
| (4.19), (4.20) | V(α,k,δ,T), K(k,δ,T); V_∞ (lower bound for V used in the check V > 1) | l. 69, l. 68; l. 70 |
| (4.23) | M(k,δ) | l. 71 |
| (4.44), (4.45) | b₅(η), b₆(X,η) | l. 79, l. 80–81 |
| (4.50), (4.51) | b₇(k,η,H₀), b₈(η,H) | l. 82–83, l. 84 |
| (4.59), (4.56) | b₉(η,H), C₇(η,H) | l. 85–88, l. 89 |
| (4.66), (4.68), (4.63) | b₁₀(k,μ), b₁₁(X,τ), C₈(k,μ) | l. 90–91, l. 92–93, l. 94–96 |
| (4.75) | b₁₂(H) = 1/(2(1 − 1/H)²) | l. 72 |
| (4.72) | C₁(α,d,δ,k,H,σ) (as log C₁) | l. 73–78 |
| (4.73) | C₂ = C₇(μ − σ + d/log H₀) + C₈, B = C₂/(2πd) = q₀ − q₁σ | l. 119–120 |
| (4.71) | A = C₁/(2πd) = exp(c₀ + c₁σ); s_min = 1/2 + d/log H₀ | l. 117–118; l. 121 |
| Table 1 | the 20 parameter rows; η | l. 99–108; l. 97 |
| Lemma 4.14 hypotheses | U(α,k,H₀) > 1, J coefficients > 0, V_∞ > 1, δ, d, α > 0, 10⁹/H₀ ≤ k ≤ 1, b₆(10⁹,η) < 1, η < 1/2, b₆(10⁹,μ−1) < 1, μ ≤ 1 + η | l. 122–125 |

## 5. Cells in x and the tail (monotonicity in x)

**5.1 Cells.** For x = n^k ∈ [e^{L_a}, e^{L_b}]: δ ∈ [δ_b, δ_a] with δ_a = k log(1 + e^{−L_a/k}), δ_b = k log(1 + e^{−L_b/k}). (i) log x·x^{σ−1} decreases in x when log x > 1/(1−σ); as 1/(1−B) < 7.9 < L_a, it is ≤ L_a e^{(σ−1)L_a}; likewise x^{σ₀−1} ≤ e^{(σ₀−1)L_a}. (ii) G(δγ) ≤ G(δ_bγ) (G non-increasing), so S, F, D, λ₀ are taken at δ_b. (iii) δ ≤ δ_a < δ_max. (iv) PP/M₀ uses log(xe^δ) ≤ L_b + δ_a, x ≥ e^{L_a}, δx ≥ δ_b e^{L_a}. The resulting R̄(L_a, L_b) bounds R for every x in the cell.
**5.2 Tail lemma.** Let L* ≥ L₀, n* = e^{L*/k}, c* = 1/(1 − 1/(2n*)). For log x = L ≥ L*: 1/δ ≤ c*n/k (log(1+y) ≥ y − y²/2), D ≤ L/k, D ≥ L/k − log k =: D_low; δT₀ ≤ kT₀c*/n* < 1 (asserted), so case (ii). With A_max = max_σ A_P(σ) (right cell ends, c₁ > 0), B_max = max B_P (left ends), M* = Mom(grid⌈8(1−σ₀)/3 + 4/D_low(L*)⌉; 0), M₁ = Mom(grid⌈1 + 1/D_low(L*)⌉; 0), Ξ(L) = L(1 − 8/(3k)) + (8/3)log(k/c*) − 2 log(L/k) > 0:
  R ≤ R_an(L) := 2[L·A_max M*(L/k)³e^{−(1−B)Ξ(L)}/Ξ(L) + B_max M*(L/k)²e^{−(1−B)L} + e^{−(1−σ₀)L}(c*/k)e^{L/k}(L/k)M₁]
     + m(L+1)[e·e^{−L/2}/4 + J·e·e^{−2L/3}/9 + log J·(c*/k)e^{L/k−L}] + 2e^{−3L},  J = (L+1)/log 2,
using e^{(σ−1)L}δ^{−a}D^{5−2σ} ≤ (L/k)³e^{−(1−σ)Ξ(L)} and ∫_{σ₀}^B e^{−(1−σ)Ξ}dσ ≤ e^{−(1−B)Ξ}/Ξ. Every term is non-increasing in L for L ≥ L* when (1−B)(1 − 8/(3k) − 2/L*) > 4/L* (main term: d log/dL ≤ 4/L − (1−B)Ξ′ with Ξ′ = 1 − 8/(3k) − 2/L > 0), 1 − σ₀ − 1/k > 2/L* and 1 − 1/k > 3/L* (checked in Arb). So R(x) ≤ R_an(log x) ≤ R_an(L*) for all x ≥ e^{L*}; cert.py checks R_an(L*) < 1.

## 6. Small n for k = 7

For 1 ≤ n ≤ 100, sympy.nextprime(n⁷) < (n+1)⁷ (all numbers < 2⁶⁴, where sympy's primality test is deterministic). The cells start at L = 7 log 100 = 32.236, where 1/δ ≥ 14 holds.

## 7. What cert.py runs (main run ≈ 100 s; with §9's δ = 1 run and the sanity check 193 s; 199 MB)

| k | cells | range of log x | cell width | max R̄ (at log x) | tail from L* | R_an(L*) |
|---|---|---|---|---|---|---|
| 7 | 592 + 3000 | [32.236, 180], [180, 240] | 0.25, 0.02 | 0.1546 (205.52) | 240 | 0.244 |
| 6 | 1651 | [217.45, 300] | 0.05 | 0.9997 (217.45) | 300 | 0.0876 |
| 5 | 2226 | [288.73, 400] | 0.05 | 0.9996 (288.73) | 400 | 0.0368 |
| 4 | 1712 | [468.84, 640] | 0.1 | 0.9997 (468.84) | 640 | 0.0169 |
| 3 | 2076 | [1881.03, 2400] | 0.25 | 0.9999 (1881.03) | 2400 | 0.0117 |

x₀(k) is the bisection point (to 0.01, rounded up) where the first cell bound drops below 1; all cells above it are then certified. Kernel order: m = 3 gives a larger bound at every probed point (k = 7, log x = 201.8: 0.170 against 0.122; k = 3, log x = 1900: 0.818 against 0.772), so m = 2 is used throughout. There is no truncation height to optimise (§2). Lower log x fail for k = 3…6 (k = 6 also holds on [27, ≈162], k = 5 on [22, ≈120], k = 4 on [17, ≈80]; not claimed).

## 8. Where the bound comes from (k = 7, check.py)

At log x = 205.5 (n = 5.62·10^{12}, δT₀ = 3.74), Z ≈ 0.152: σ ∈ [0.86, 0.87] and T ∈ [T₀, 10^{13}] give 0.062, σ ∈ [0.87, B] 0.046, σ ∈ [0.83, 0.86] 0.024, everything with σ < 0.83 below 7·10^{−4}, heights above 10^{15} below 4·10^{−5}. In this band KLN beats the trivial count (N(T) − N(T₀))/2 by far: at T = 10^{13}, N(0.85,T) ≤ 1.06·10^{11} and N(0.87,T) ≤ 2.07·10^{10} against 1.54·10^{13}; the two cross only near σ ≈ 0.78 (3.15·10^{13} vs 1.54·10^{13}). The certificate never uses the trivial count.

## 9. Checks

- check.py: an independent numpy quadrature of Z (exact KLN shape (log kt)^{2σ}(log t)^{5−4σ}t^{8(1−σ)/3}, best row at each σ, the true envelope G, no moment tables, no σ-cells) gives cert R̄ / direct Z = 1.011–1.227 at ten points (k = 7 at log x = 40 … 230; k = 6, 5, 4, 3 at log x₀(k)); logs/check.log.
- δ = 1 in every KLN row (Theorem 1.1 as printed): (a) still holds (max R̄ 0.3482, R_an(240) = 0.559); log x₀ = 1947.12, 491.86, 305.60, 231.99 for k = 3, 4, 5, 6 (n₀(3) = e^{649.04}).
- Sanity check (floats, not a proof): replacing the half-plane by Ford's region β < 1 − 1/(57.54 (log γ)^{2/3}(log log γ)^{1/3}) in the same pipeline gives log log x₀(3) = 33.99, i.e. n₀(3) = exp(e^{32.90}): Cully-Hugill's 2021 bound (x ≥ exp(e^{33.990}), n ≥ exp(e^{32.892})) and Dudek's exp(e^{33.217}) come back.
- Kernel: random test above; s₁ check in cert.py.
- Truncation: none (absolutely convergent zero sum). ψ versus θ: only prime powers inside the interval (PP). Trivial zeros: ∫w/(u³−u) ≤ W(1)/(x(x²−1)).

## 10. Constants (each appears in cert.py)

B = 0.8733819666132076; T* = B₂ + 2.2e-254; T₀ = 3e12; HSW 0.1038, 0.2573, 9.3675; KLN H₀ = 3.0610046e10, b₁ = 0.62, b₂ = 1.048, b₃ = 0.605, b₄ = 0.529, a₁ = 0.63, a₂ = 2.851, 3.006, η = 0.2561882, Table 1 rows (σ₀, k, μ, α, δ, d); σ₀ = 0.515, σ-cell 0.0025; s₁ = 0.2173; δ_max = 1e-9, 0.025, 0.15; λ-grid ratio 1 + 1/500, Λ ≥ 2m·10000, stored up to λ = 30; p-grid j/200, j ≤ 300; small-n bound 100; tail starts L* = 240, 300, 400, 640, 2400 (2600 for k = 3 in the δ = 1 run); Ford constant 57.54 (sanity check only); B₂'s 777, 1081, 852. All of these were found in cert.py with rg. The other numbers in this file are outputs of cert.py/check.py (in logs/) or quoted literature values that enter only through a₁, a₂ (0.618, 1.461).
