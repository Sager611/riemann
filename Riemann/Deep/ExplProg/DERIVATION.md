# Explicit primes in progressions from the zero-free half-plane — derivation

Every step and every numeric constant used by `cert.py` (section numbers match the comments there).
Notation: q ≥ 1 modulus, (a, q) = 1, φ Euler's function, χ a Dirichlet character, χ* the primitive character
inducing χ, q* its conductor, ρ = β + iγ nontrivial zeros with multiplicity, X = log x, L = log q,
Λ von Mangoldt, ψ(x;q,a) = Σ_{n≤x, n≡a} Λ(n), θ(x;q,a) the same over primes.

## 0. Inputs

- **(H)** Bebel (QuasiRiemannTracker, Lean-kernel-checked, 1 Nov 2026): ζ(s) ≠ 0 and L(s, χ) ≠ 0 for Re s > B_2 + 2.2·10⁻²⁵⁴, for
  every Dirichlet character χ, where B_2 = (777 − √1081)/852 = 0.87338196661320751993…. We use
  **B := 0.8733819666132076** (certified B > B_2 + 10⁻¹⁷), b := 1 − B = 0.1266180333867924, 1/b = 7.89776916646….
- **(FE)** For ζ and for primitive χ (q* ≥ 3) the functional equation makes ρ ↦ 1 − ρ̄ an involution of the multiset of
  nontrivial zeros; it fixes γ. With (H): **b ≤ β ≤ B for every nontrivial zero**. No Siegel zeros exist.
- **(BMOR)** Bennett–Martin–O'Bryant–Rechnitzer, Math. Comp. 90 (2021) 1455–1482 (arXiv 2005.02989), Thm 1.1 (checked in the
  TeX source): for χ of conductor q > 1, T ≥ 5/7, ℓ = log(q(T + 2)/2π): if ℓ ≤ 1.567 then N(T, χ) = 0; if ℓ > 1.567 then
  |N(T, χ) − ((T/π) log(qT/2πe) − χ(−1)/4)| ≤ 0.22737ℓ + 2 log(1 + ℓ) − 0.5. Here N(T, χ) = #{ρ : 0 < β < 1, |γ| ≤ T}
  with multiplicity.
- **(ROS)** Rosser, Ann. of Math. 42 (1941), Thm 19: for T ≥ 2, |N(T) − (T/2π) log(T/2πe) − 7/8| < 0.137 log T + 0.443 log log T
  + 1.588, N(T) = #{ρ of ζ : 0 < γ ≤ T}. Used only in Theorem 1.
- **(DAV)** Davenport, Multiplicative Number Theory, ch. 12: Σ_ρ Re(1/ρ) = 1 + γ_E/2 − ½ log 4π = 0.0230957089661… over the zeros of
  ζ; ch. 16–17: |L'/L(σ + iT_j)| ≪ log²(qT_j) on suitable heights T_j → ∞ for −1 ≤ σ ≤ 2 and |L'/L(s)| ≪ log(q|s|) for
  Re s ≤ −1 away from trivial zeros (qualitative, only to justify contour shifts). The first zero of ζ has γ = 14.1347… > 14.
- **(RS62)** Rosser–Schoenfeld, Illinois J. Math. 6 (1962): θ(x) < 1.01624x and ψ(x) < 1.03883x for x > 0 (Thms 9, 12);
  φ(n) > n/(e^{γ_E} log log n + 2.50637/log log n) for n ≥ 3 (Thm 15).
- Elementary: φ(n) ≥ √(n/2); ω(q) ≤ log q/log 2; Σ_{p|q} log p ≤ log q.

No GRH, no verification of GRH to a height, no unrefereed zero count is used.

## 1. Explicit formula for weights supported away from 0 (sec 1 of cert.py comments)

**Lemma EF.** Let w : (0, ∞) → [0, ∞) be continuous with support in [t₀, t₁] ⊂ (0, ∞) and w' of bounded variation, and
W(s) = ∫₀^∞ w(t) t^{s−1} dt (entire). Let χ be primitive of conductor k (k = 1: χ = 1, L = ζ). If x t₀ > 1, then

  Σ_n Λ(n) χ(n) w(n/x) = 1_{k=1} x W(1) − Σ_ρ x^ρ W(ρ) − Σ_{m∈T(χ)} x^{−m} W(−m),

T(1) = {2, 4, …}, T(χ) = {0, 2, 4, …} for even χ ≠ 1, T(χ) = {1, 3, …} for odd χ; the ρ-sum converges absolutely, and

  0 ≤ Σ_{m∈T(χ)} x^{−m} W(−m) ≤ W(0)/(1 − (x t₀)^{−2}).                                          (1.1)

*Proof.* Two integrations by parts give W(s) = (s(s+1))^{−1} ∫ t^{s+1} dw'(t), so |W(σ + it)| ≤ Var(w')·max(t₀^{σ+1}, t₁^{σ+1})/|s(s+1)|.
Mellin inversion on Re s = 2 is absolutely convergent. Shift to Re s = −U, U → ∞ through half-odd integers, along the heights
T_j of (DAV); the horizontal integrals vanish (|W| ≪ |s|^{−2}), and on Re s = −U, |x^s W(s)| ≤ (x t₀)^{−U} t₀ Var(w')/|s(s+1)| times
log(k|s|), which tends to 0 because x t₀ > 1. Residues: x W(1) at s = 1 (k = 1), −m_ρ x^ρ W(ρ) at the zeros, −x^{−m} W(−m) at the
simple trivial zeros (ζ(0) = −1/2 ≠ 0, so no residue at 0 for ζ). For (1.1): x^{−m} W(−m) = ∫ w(t) t^{−1}(x t)^{−m} dt ∈ [0, (x t₀)^{−m} W(0)],
and Σ_{m∈T} (x t₀)^{−m} ≤ 1/(1 − (x t₀)^{−2}). ∎

**Orthogonality and imprimitive characters.** For (a, q) = 1, S(x; q, a) := Σ_{n≡a} Λ(n) w(n/x) satisfies
φ(q) S = Σ_{χ mod q} χ̄(a) Σ_n Λ(n) χ(n) w(n/x). If χ is induced by χ*, the two inner sums differ by
Σ_{p|q, p∤q*} Σ_{k≥1} χ*(p^k) log p · w(p^k/x), of modulus at most

  D_q := ‖w‖_∞ Σ_{p|q} log p · #{k ≥ 1 : p^k ∈ [x t₀, x t₁]} ≤ ‖w‖_∞ (ω(q) log(t₁/t₀) + log q).                   (1.2)

The principal character is induced by 1, i.e. by ζ. With Lemma EF:

  |φ(q) S(x; q, a) − x W(1)| ≤ Ω₁ + (φ(q) − 1) Ω + φ(q) [W(0)/(1 − (x t₀)^{−2}) + D_q],                          (1.3)

Ω₁ := Σ_{ρ of ζ} |x^ρ W(ρ)|, Ω := max over primitive χ* of conductor 3 ≤ q* | q of Σ_{ρ of χ*} |x^ρ W(ρ)|.

## 2. Zero sums

**Lemma P (pairing).** Let the zeros be those of ζ or of a primitive L(s, χ). Suppose, for all real γ,
sup_{β∈[1/2,B]} x^β |W(β + iγ)| ≤ x^B G(|γ|) and sup_{β∈[b,1/2]} x^β |W(β + iγ)| ≤ x^{1/2} G_L(|γ|). Then

  Σ_ρ |x^ρ W(ρ)| ≤ ½ Σ_ρ (x^B G(|γ_ρ|) + x^{1/2} G_L(|γ_ρ|)).

*Proof.* By (FE) the multiset splits into pairs {β + iγ, 1 − β + iγ} with 1/2 < β ≤ B and singletons with β = 1/2. A pair contributes
at most x^B G + x^{1/2} G_L, which is the right side summed over its two members; a singleton contributes at most
min(x^B G, x^{1/2} G_L) ≤ ½(x^B G + x^{1/2} G_L). ∎

**Lemma Z (counting, BMOR).** Let G(t) = min(n₀, k₁/t, k₂/t²) for t > 0, G(0) = n₀. Let q ≥ 9, L = log q. For every
primitive χ of conductor 3 ≤ q* ≤ q:

  Σ_ρ G(|γ_ρ|) ≤ Z_G(L) := L z_a + 2 n₀ log(1 + L) + z_b + z_e/(1 + L),

  z_a = a(5/7) n₀ + (1/π) ∫_{5/7}^∞ G,  z_b = b(5/7) n₀ + ∫_{5/7}^∞ ρ₀ G,  z_e = 2 ∫_{max(5/7, 2π−2)}^∞ G(t)/(t + 2) dt,

a(t) = t/π + 0.22737, b(t) = (t/π) log(t/2πe) + 1/4 + 0.22737 log((t + 2)/2π) − 1/2, ρ₀(t) = (1/π) log(t/2π) + 0.22737/(t + 2);
a(5/7) = 0.4547342044, b(5/7) = −1.162576012.

*Proof.* (i) For t ≥ 5/7 put N⁺(t) := L a(t) + 2 log(1 + L) + b(t) + e(t)/(1 + L), e(t) = 2 max(0, log((t + 2)/2π)). Then
N(t, χ) ≤ N⁺(t): if ℓ* = log(q*(t + 2)/2π) > 1.567, BMOR gives N ≤ R(t; q*) := (t/π) log(q*t/2πe) + 1/4 + 0.22737ℓ* + 2 log(1 + ℓ*) − 1/2,
R is increasing in q*, and 2 log(1 + L + c) ≤ 2 log(1 + L) + 2c⁺/(1 + L) (c = log((t + 2)/2π), concavity), so R(t; q*) ≤ R(t; q) ≤ N⁺(t).
If ℓ* ≤ 1.567, N = 0 ≤ R(t; 9) ≤ N⁺(t): each term of R(t; 9) is non-decreasing on t ≥ 5/7 (5/7 > 2π/9 and ℓ(5/7; 9) > 0) and
R(5/7; 9) = 1.5520999 > 0 (certified).
(ii) For t < 5/7, N(t) ≤ N(5/7) ≤ N⁺(5/7).
(iii) Σ_ρ G(|γ_ρ|) = ∫_{[0,∞)} N(t)(−dG(t)) ≤ N⁺(5/7)(n₀ − G(5/7)) + ∫_{5/7}^∞ N⁺(−dG) = N⁺(5/7) n₀ + ∫_{5/7}^∞ (dN⁺/dt) G dt
(N⁺G → 0), and dN⁺/dt = L/π + ρ₀(t) + 2·1_{t>2π−2}/((1 + L)(t + 2)). Collecting the L-terms, the constants and the 1/(1 + L) terms gives
z_a, z_b, z_e; e(5/7) = 0.
(iv) G is piecewise c t^{−m}: n₀ on [0, t_a], k₁/t on [t_a, t_b], k₂/t² on [t_b, ∞) with t_a = k₁/n₀, t_b = k₂/k₁ if k₁² ≤ n₀ k₂
(else n₀ up to t_c = √(k₂/n₀), then k₂/t²). The integrals are exact with the antiderivatives
  ∫t^{−m}/π: t/π, (log t)/π, −1/(πt);  ∫ρ₀ t^{−m}: (t log(t/2π) − t)/π + 0.22737 log(t + 2), log²(t/2π)/(2π) + 0.113685 log(t/(t + 2)) (exactly 0.22737/2, as cert.py uses CB/2),
  −(log(t/2π) + 1)/(πt) + 0.22737(−1/(2t) + ¼ log((t + 2)/t));  ∫2t^{−m}/(t + 2): 2 log(t + 2), log(t/(t + 2)), −1/t + ½ log((t + 2)/t)
for m = 0, 1, 2 (each checked by differentiation). The conductors 3 ≤ q* < 9 are covered by using L = max(log q, log 9). ∎

**Lemma Zζ (first moment of 1/γ²).** Σ_{ρ of ζ} 1/γ² ≤ 2(1 + B²/196) Σ_ρ β/|ρ|² = 0.04637118646 (sum = (DAV) value).
*Proof.* For a pair at height γ (|γ| ≥ 14): β/|ρ|² + (1 − β)/|1 − ρ̄|² ≥ 1/(B² + γ²) because both moduli are ≤ B² + γ²; hence
2/γ² ≤ 2(1 + B²/γ²)·(pair sum). A singleton: 1/γ² = 2(1 + 1/(4γ²))·β/|ρ|² ≤ 2(1 + B²/γ²)β/|ρ|². ∎
Hence, if G ≤ k₂/t² and G_L ≤ k₂^L/t², Lemma P gives **Ω₁ ≤ 0.02318559323 (x^B k₂ + x^{1/2} k₂^L)**.

**Lemma Zζ' (Rosser count, Theorem 1 only).** N_ζ(t) := #{ρ : |γ| ≤ t} = 2N(t) ≤ (t/π) log(t/2πe) + 4.926 + 0.274 log t
+ 0.886 log log t for t ≥ 2 (ROS), and N_ζ(t) = 0 for t < 14. As in Lemma Z,
Σ_ρ G(|γ_ρ|) ≤ N_ζ⁺(14) G(14) + ∫_{14}^∞ ρ_ζ G, ρ_ζ = (1/π) log(t/2π) + 0.274/t + 0.886/(t log t), with the same antiderivatives and
∫_u^v dt/(t² log t) ≤ (1/u − 1/v)/log u, ∫_u^v dt/(t³ log t) ≤ (1/(2u²) − 1/(2v²))/log u.

## 3. Kernels

### 3.1 The Linnik / non-residue kernel
g(u) on [−ℓ, 0] is the triangle rising linearly from 0 at −ℓ to 1 at −c and falling to 0 at 0, **ℓ = 7, c = 21/10** (rise
a = 49/10). w(t) = g(log t) t^{−B}, support [e^{−7}, 1], W(s) = ∫ g(u) e^{(s−B)u} du. Closed form
∫ g e^{mu} du = (e^{−mℓ}/a − (1/a + 1/c) e^{−mc} + 1/c)/m² (two integrations by parts; dg' has masses 1/a, −(1/a + 1/c), 1/c).
  n₀ = ‖g‖₁ = 7/2, n₁ = ‖g'‖₁ = 2, n₂ = Var(g') = 2(1/a + 1/c) = 1.3605442;
  W(1) = ∫ g e^{bu} = 2.424245725, W(0) = ∫ g e^{−Bu} = 115.9708945, ‖w‖_∞ = max(e^{Bc}, e^{Bℓ−1}/(aB)) = 38.854381 (a > 1/B);
  left kernel g_L = g e^{(1/2−B)u}: m₀ = ‖g_L‖₁ = 12.707659, m₁ = 12.526061 ≥ ‖g_L'‖₁, m₂ = 12.334200 ≥ Var(g_L')
  (m₁ = ∫|g'|e^{mu} + |m| m₀, m₂ = ∫e^{mu}|dg'| + 2|m|∫|g'|e^{mu} + m² m₀, m = 1/2 − B).

**Lemma K1.** Let Y := X − ℓ = log(x e^{−7}) > 0. For β ∈ [1/2, B] and γ ≠ 0:
x^β |W(β + iγ)| ≤ x^B min(n₀, k₁/|γ|, k₂/γ²), k₁ = n₁ + n₀/(eY), k₂ = n₂ + 2n₁/(eY) + (2/(eY))² n₀; for β ∈ [b, 1/2] the same with
x^{1/2}, (m₀, m₁, m₂).
*Proof.* With ε = B − β ≥ 0: x^β W(β + iγ) = x^B ∫ h(u) e^{iγu} du, h = g e^{−ε(u+X)}. On supp g, u + X ≥ Y, so |h^{(j)}| ≤
Σ_i C(j,i)|g^{(i)}| ε^{j−i} e^{−εY} (as measures for j = 2), and sup_{ε≥0} ε^m e^{−εY} = (m/(eY))^m. Integrate by parts j = 0, 1, 2
times. Left zeros: ε' = 1/2 − β, x^β W(β + iγ) = x^{1/2} ∫ g_L e^{−ε'(u+X)} e^{iγu} du. ∎
At q = 10⁴ in Theorem 2 (Y = 83.28): k₁ = 2.01546, k₂ = 1.37849, z_a = 2.20587, z_b = −4.67563, z_e = 0.05773; at Y = ∞:
z_a = 2.1978743, so the leading coefficient of Theorem 2 tends to c_∞ = z_a/(2W(1)) = 0.453311. Of z_a, a(5/7)n₀ = 1.5916 (72%)
comes from the N(5/7) ≈ 0.455 log q zeros below height 5/7, all placed at height 0.

### 3.2 Trapezoids (Theorem 1)
α = x^{−1/2}, 0 < δ ≤ 1/2, x ≥ 64. w₊: 0 on (0, α(1 − δ)], linear to 1 at α, 1 on [α, 1], linear to 0 at 1 + δ.
w₋: 0 on (0, α], linear to 1 at α(1 + δ), 1 on [α(1 + δ), 1 − δ], linear to 0 at 1. Then 0 ≤ w₋ ≤ 1_{(α,1]} ≤ w₊,
W₊(1) = 1 − α + (1 + α)δ/2, W₋(1) = 1 − α − (1 + α)δ/2, W_±(0) ≤ log((1 + δ)/(α(1 − δ))), ‖w_±‖_∞ = 1.

**Lemma K2.** For β ∈ [1/2, B]: x^β |W_±(β + iγ)| ≤ x^B min(W_B, A₁/|γ|, A₂/γ²) with
W_B = ((1 + δ)^{B+1} − 1)/(δB(B + 1)), A₁ = (1 + δ)^B(1 + α^B), A₂ = ((1 + δ)^{B+1} + 1)(1 + α^B)/δ; for β ∈ [b, 1/2] the same with
x^{1/2} and W_h = ((1 + δ)^{3/2} − 1)/(3δ/4), A₁^L = (1 + δ)^{1/2}(1 + α^{1/2}), A₂^L = ((1 + δ)^{3/2} + 1)(1 + α^{1/2})/δ.
*Proof.* W(s) = −s^{−1} ∫ w' t^s dt = (s(s + 1))^{−1} ∫ t^{s+1} dw'. On the support xt ≥ x α(1 − δ) ≥ 4, so x^β t^{β+j} = (xt)^β t^j ≤ (xt)^B t^j.
∫|w'| t^B ≤ u₁^B + v₁^B and ∫ t^{B+1}|dw'| = (u₀^{B+1} + u₁^{B+1})/(u₁ − u₀) + (v₀^{B+1} + v₁^{B+1})/(v₁ − v₀) for ramps [u₀, u₁], [v₀, v₁];
|s|, |s + 1| ≥ |γ|. For j = 0: x^β|W| ≤ ∫ w (xt)^β dt/t ≤ x^B ∫ w t^{B−1} ≤ x^B W_δ(B), w_± ≤ w_δ := min(1, (1 + δ − t)/δ)⁺,
W_δ(s) = ((1 + δ)^{s+1} − 1)/(δs(s + 1)). Left zeros: (xt)^β ≤ (xt)^{1/2}, |s| ≥ |b + iγ| ≥ |γ|. ∎

## 4. Theorem 2 (least primes)

**Theorem 2.** For every q ≥ 2 and (a, q) = 1 there is a prime p ≡ a (mod q) with p ≤ (q log q)^{1/(1−B)} = (q log q)^{7.8978}.
For q ≥ Q, (q log q) may be replaced by (c₂ q log q) with
c₂ = 0.79 (Q = 10²), 0.75 (10³), 0.72 (10⁴), 0.67 (10⁶), 0.62 (10¹⁰), 0.56 (10²⁰), 0.53 (10³⁰), 0.51 (10⁵⁰), 0.49 (10¹⁰⁰),
0.46 (10¹⁰⁰⁰); c₂ → c_∞ ≤ 0.4534 (any c₂ > 0.453311 for q large enough).

*Proof.* Kernel 3.1, P(x) := Σ_{p≡a} log p · w(p/x). Prime powers in the class contribute at most
‖w‖_∞(ψ(x) − θ(x)) ≤ ‖w‖_∞ · 1.01624(√x + x^{1/3} log x/log 2) (RS62: ψ − θ = Σ_{k≥2} θ(x^{1/k}), k ≤ log₂ x). By (1.3), Lemmas P,
Z, Zζ, K1, (1.1), (1.2) with t₀ = e^{−7}, t₁ = 1 and φ(q) ≤ q, P(x) > 0 (hence a prime ≡ a in (e^{−7}x, x)) as soon as

  x W(1) > 0.0231856(x^B k₂ + x^{1/2} k₂^L) + q [ (x^B Z(L) + x^{1/2} Z_L(L))/2 + W(0)/(1 − e^{−2Y}) + ‖w‖_∞(7L/log 2 + L)
          + 1.01624 ‖w‖_∞(x^{1/2} + x^{1/3} X/log 2) ],                                                           (4.1)

Z = Z_G with (n₀, k₁, k₂), Z_L with (m₀, k₁^L, k₂^L), L = log q ≥ log 9.
*Tail lemma (all q ≥ e^{L_m} at once).* Put x = (c₂ q L)^{1/b}, so x^b = c₂ q L; divide (4.1) by q x^B L. It suffices that
c₂ W(1) > F(q) := Z(L)/(2L) + x^{1/2−B} Z_L(L)/(2L) + 0.0231856(k₂ + x^{1/2−B} k₂^L)/(qL) + W(0)/((1 − e^{−2Y}) x^B L)
+ ‖w‖_∞(7/log 2 + 1)/x^B + 1.01624 ‖w‖_∞(x^{1/2−B} + x^{1/3−B} X/log 2)/L.
Every term is non-increasing in q ≥ e^{L_m} when the k's are evaluated at Y_m := Y(e^{L_m}) (they decrease in Y, and Y increases with q):
Z(L)/L = z_a + (2n₀ log(1 + L) + z_b)/L + z_e/(L(1 + L)); h(L) = (2n₀ log(1 + L) + z_b)/L is non-increasing on [L_m, ∞) if
2n₀(log(1 + L_m) − L_m/(1 + L_m)) + z_b ≥ 0 (h' has numerator 2n₀(L/(1 + L) − log(1 + L)) − z_b and log(1 + L) − L/(1 + L) increases),
otherwise z_b is replaced by max(z_b, 0); z_e is replaced by max(z_e, 0); x^{1/3−B}X decreases for X > 1/(B − 1/3); the other terms
are monotone. So F(q) ≤ F̄(L_m), its value at L_m, and the certificate checks c₂ W(1) > F̄(L_m) (margins 0.26%–1.4% for the table,
21% for c₂ = 1 at Q = 100).
*Small q.* For 2 ≤ q < 100 (and as a check up to 10⁴) an exact sieve to 10⁸ gives max_a p(a, q); the worst ratio
log p(a,q)/log((q log q)^{1/b}) over 2 ≤ q ≤ 10⁴ is 0.426 (q = 2, p = 3). ∎

**Theorem 2' (φ(q) in place of q).** For every q ≥ 3, p(a, q) ≤ (φ(q) log q)^{1/(1−B)}, and for q ≥ Q the same table of c₂ holds with
φ(q) log q in place of q log q. *Proof.* Same inequality (4.1) with φ(q) instead of q in front of the bracket (no use of φ(q) ≤ q),
x = (c₂ φ(q) L)^{1/b}, divided by φ(q) x^B L. Every right-hand term is non-increasing in x and in 1/φ(q), so it may be evaluated at
x_lo = (c₂ f(q) L)^{1/b} with f(q) = q/(e^{γ_E} log log q + 2.50637/log log q) ≤ φ(q) (RS62 Thm 15); f is increasing for q ≥ 16
(f' has the sign of e^{γ_E}(log log q − 1/log q) + 2.50637(1/log log q + 1/(log q (log log q)²)) > 0), so x_lo, Y_lo and f(q)L increase
with q and the tail lemma applies at L_m = log Q (certified margins 0.096%–1.4%, 21% for c₂ = 1 at Q = 100). For 3 ≤ q < 100 the
sieve gives max_a p(a, q) ≤ (φ(q) log q)^{1/b} (see logs/sanity_output.txt). (q = 2 is excluded: φ(2) log 2 < 1.) ∎

## 5. Theorem 3 (least non-residues of characters)

**Theorem 3.** Every nonprincipal χ mod q ≥ 3 has χ(n) ∉ {0, 1} for some (prime) n ≤ (log q)^{1/(1−B)} = (log q)^{7.8978}.
For log q ≥ Λ, (log q) may be replaced by (c₃ log q): c₃ = 0.62 (Λ = 25), 0.60 (30), 0.56 (50), 0.52 (100), 0.47 (1000),
0.46 (10⁴); c₃ → c_∞ = 0.453311.

*Proof.* Let χ be nonprincipal mod q, χ* its primitive (conductor ≥ 3), D(x) = Σ_n Λ(n)(1 − Re χ(n)) w(n/x) ≥ 0 with kernel 3.1.
If χ(n) ∈ {0, 1} for all n ≤ x then only n with (n, q) > 1 contribute, so D(x) ≤ D_q. By Lemma EF for ζ and χ* and (1.2),
D(x) ≥ x W(1) − Ω₁ − Ω_{χ*} − 2W(0)/(1 − e^{−2Y}) − D_q. Contradiction once
x W(1) > Ω₁ + (x^B Z(L) + x^{1/2} Z_L(L))/2 + 2W(0)/(1 − e^{−2Y}) + 2‖w‖_∞(7L/log 2 + L). The tail lemma of §4 applies verbatim
with x = (c₃ L)^{1/b} (divide by x^B L). For 3 ≤ q ≤ e^{25.6} the trivial bound n_χ ≤ q − 1 suffices: on 7483 cells [l₀, l₁] of
log q ∈ [log 3, 25.6], (1/b) log l₀ ≥ log(e^{l₁} − 1) (certified; tightest at q = 3: 2 ≤ 2.10). If n is minimal with χ(n) ∉ {0, 1},
n is prime (χ is completely multiplicative). ∎
Sanity: max_χ n_χ over 3 ≤ q ≤ 10⁴ is 43 (q = 2706), computed as the least N with ⟨n ≤ N, (n, q) = 1⟩ = (ℤ/q)^×.

## 6. Theorem 1 (uniform prime number theorem in progressions)

**Theorem 1.** For every q ≥ 1, (a, q) = 1 and real x ≥ 2: |ψ(x; q, a) − x/φ(q)| ≤ 1.5 x^B log²(qx).
For x ≥ 64 the constant 1.5 can be replaced by 0.174, for x ≥ e^{10} by 0.0324, for x ≥ e^{20} by 0.0147, for x ≥ e^{100} by 0.0112,
for x ≥ e^{1000} by 0.0079, for x ≥ e^{10⁶} by 0.0076 (the tail of 6.4 at X_T = 100, 1000, 10⁶; its limit is
½(b/(4π) + b²/π) = 0.0076). The constant 1.5 is forced (up to 1.4854) by q = 1, x = 2.

*Proof.* Let E := ψ(x; q, a) − x/φ(q) and R := |E|/(x^B log²(qx)).
**6.1 x < q.** The class has at most one n ≤ x, so ψ(x; q, a) ≤ log x, and x/φ(q) ≤ x√(2/q) < √(2x). With log²(qx) > 4 log² x,
R ≤ max(1/(4x^B log x), √2 x^{1/2−B}/(4 log² x)), decreasing in x; at x = 2 this is 0.5681 (certified); at x = 64 it is 0.0043.
**6.2 Exact range q ≤ x < 66.** ψ(n; q, a) in ball arithmetic for all q ≤ 65, (a, q) = 1; on each [n, n + 1) split into 2000
(n = 2), 200 (n < 8) or 4 pieces [s₀, s₁]: R ≤ max(|ψ − s₀/φ|, |ψ − s₁/φ|)/(s₀^B log²(q s₀)). Maximum 1.48535 at q = 1, x ∈ [2, 2.0005];
on 64 ≤ x < 66 the maximum is 0.00849.
**6.3 Sandwich, x ≥ 64, q ≤ x.** y = √x. Since w₋ ≤ 1_{(α,1]} ≤ w₊,
ψ(y; q, a) + S₋ ≤ ψ(x; q, a) ≤ ψ(y; q, a) + S₊, S_± = Σ_{n≡a} Λ(n) w_±(n/x). With W_±(1) from 3.2 and (1.3) (t₀ ≥ α(1 − δ), so
x t₀ ≥ y/2 ≥ 4):
  |E| ≤ |ψ(y; q, a) − y/φ(q)| + (x + y)δ/(2φ(q)) + Ω₁/φ(q) + Ω + W(0)/(1 − (y(1 − δ))^{−2}) + D_q =: T₁ + … + T₆,
T₁ ≤ max(ψ(y), y) ≤ 1.03883 y; T₂ = (x + y)δ/(2φ(q)); T₃, T₄ by Lemmas P, K2 and Zζ' resp. Z (L ≥ log 9 or replaced by log 9);
T₅ ≤ (X/2 + log((1 + δ)/(1 − δ)))/(1 − (y(1 − δ))^{−2}) (≤ …/(1 − 4/x)); T₆ ≤ (L/log 2)(X/2 + log((1 + δ)/(1 − δ))) + L (T₆ = 0 for q = 1; X/2 + … + log 2 for q = 2).
**Cells.** q = 1, q = 2, and q ∈ [e^{L₀}, e^{L₁}] (L₀ ≥ log 3) with x ∈ [e^{X₀}, e^{X₁}], X₀ from 4.16 to 60, X₁ − X₀ = max(1/20, X₀/50),
L₁ − L₀ = max(1/10, L₀/20), 7307 cells. On a cell: φ(q) ≥ φ_lo = max(e^{L₀/2}/√2, e^{L₀}/(e^{γ_E} log L₁ + 2.50637/log L₀)) (RS62);
δ = min(1/2, θ φ_lo e^{−bX₁}) for θ ∈ {1/4, 1, 4, 16, 64} (best taken; δ an exact rational), so T₂/x^B ≤ δ(e^{bX₁} + e^{(1/2−B)X₀})/(2φ_lo);
T₁/x^B, T₃/x^B, T₄/x^B use x^{1/2−B} ≤ e^{(1/2−B)X₀}, α ≤ e^{−X₀/2}, L₁; T₅/x^B, T₆/x^B use X₁ above and e^{BX₀} below; the denominator is
(L₀ + X₀)². Maximum over cells: 0.1731 (at small x), 0.0323 over cells meeting X ≥ 10, 0.0146 over cells meeting X ≥ 20.
**6.4 Tail X ≥ 60, all 0 ≤ L ≤ X.** Take δ = θ e^{−bX}, θ = 1 (φ ≥ 1), δ ≤ δ_T = e^{−60b}. Then W_B ≤ N₀ := W_B(δ_T), A₁ ≤ K₁ :=
(1 + δ_T)^B(1 + e^{−BX_T/2}), t_a = A₁/W_B ∈ [1/N₀, K₁B] ⊂ [5/7, 1.25], t_b = A₂/A₁ ∈ [2/δ, (2 + δ_T)/δ], so
λ := log(t_b/2π) ∈ [0, bX + c_θ], c_θ = log((2 + δ_T)/(2πθ)) (certified: three-piece envelope, ρ₀ < 0 on [5/7, t_a]). From Lemma Z:
z_a ≤ a(5/7)N₀ + (2K₁ + K₁(λ + log 2π + log N₀))/π − (5/7)/(πB), z_b ≤ b(5/7)/B + K₁(λ²/(2π) + 0.1137 log(1 + 2N₀) + (λ + 1)/π + 0.1137),
z_e ≤ K₁(log(1 + 2/(2π − 2)) + δ_T/2); from Lemma Zζ': Σ_ζ ≤ K₁(N_ζ⁺(14)/14 + 0.0196 + 0.886/(14 log 14) + 1/π + 0.00979 + 0.886/(28 log 14))
+ K₁λ²/(2π) + K₁λ/π; the left kernels are ≤ M times the right ones (M = max(BW_h(δ_T), (1 + δ_T)^{1/2}(1 + e^{−X_T/4}),
((1 + δ_T)^{3/2} + 1)(1 + e^{−X_T/4})/2)). Dividing by (L + X)² with L ≤ X uses L/(L + X)² ≤ 1/(4X), λ/X ≤ b + c_θ⁺/X and
log(1 + max(L, log 9)) ≤ log(1 + X); every resulting term is non-increasing in X, and the value at X = 60 is 0.01420 (0.01113 at 100,
0.00788 at 1000, 0.00759 at 10⁶). W_B ≤ (1 + δ)^B/B and W_h ≤ 2(1 + δ)^{1/2} (mean value theorem) avoid cancellation for tiny δ. ∎
Sanity (q ≤ 100, 597 points x ∈ [2, 10⁸], exact ψ): max R = 1.4848 (q = 1, x = 2); for x ≥ 64 max R = 0.0093;
max |E|/(√x log²(qx)) = 0.036 for x ≥ 100 (consistent with GRH-size errors there).

## 7. Theorem 4 (logarithmic derivative at s = 1; item (iv))

**Theorem 4.** For every primitive χ mod q ≥ 3: |L'/L(1, χ)| ≤ (1/(1−B)) log log q + C with C = 10.53; C = 0.352 for log q ≥ 25,
−1.10 for log q ≥ 100, −1.99 for log q ≥ 1000, −2.14 for log q ≥ 10⁴ (certified values 10.526, 0.35129, −1.1040, −1.9962, −2.1408).
For imprimitive χ add Σ_{p|q} log p/(p − 1) (from L(s, χ) = L(s, χ*)Π_{p|q}(1 − χ*(p)p^{−s})).

**Lemma EF1.** χ primitive mod k ≥ 3, x > 1, 0 < h ≤ 1, v(t) = min(1, (1 + h − t)/h)⁺, V(s) = ((1 + h)^{s+1} − 1)/(hs(s + 1)):
  Σ_n Λ(n)χ(n) n^{−1} v(n/x) = −L'/L(1, χ) − Σ_ρ x^{ρ−1} V(ρ − 1) − Σ_{m∈T(χ)} x^{−1−m} V(−1 − m).
*Proof.* For c > 0 the left side is (2πi)^{−1}∫_{(c)} (−L'/L)(1 + s, χ) x^s V(s) ds (|V| ≪ |s|^{−2}). V is meromorphic with one simple pole,
at s = 0, residue 1 (s = −1 is removable), and |x^s V(s)| ≤ x^σ((1 + h)^{σ+1} + 1)/(h|s||s + 1|); shift to Re s = −U as in Lemma EF.
Residues: −L'/L(1, χ) at 0 (L(1, χ) ≠ 0), −x^{ρ−1}V(ρ − 1) at s = ρ − 1, −x^{−1−m}V(−1 − m) at the trivial zeros. ∎
*Bounds.* V(s) = s^{−1} h^{−1}∫₁^{1+h} t^s dt gives |V(s)| ≤ 1/|s| for Re s ≤ 0, and |V(s)| ≤ ((1 + h)^{Re s+1} + 1)/(h|s||s + 1|).
For s = ρ − 1: |s| ≥ max(1 − β, |γ|), |s + 1| ≥ |γ|. Hence, with x^{β−1} ≤ x^{−b} (β ≤ B) resp. ≤ x^{−1/2} (β ≤ 1/2), Lemma P applies with
G = min(1/b, 1/t, ((1 + h)^B + 1)/(ht²)) and G_L = min(2, 1/t, ((1 + h)^{1/2} + 1)/(ht²)) (G ≥ G_L, so the singletons β = 1/2 are covered),
and Lemma Z bounds Σ_ρ G, Σ_ρ G_L. Trivial zeros: |V(−1)| = log(1 + h)/h ≤ 1 and |V(−1 − m)| = (1 − (1 + h)^{−m})/(hm(m + 1)) ≤ 1/(m + 1)
(Bernoulli), so their sum is ≤ Σ_{m≥0} x^{−1−m} = 1/(x − 1). Prime sum: |Σ Λχv/n| ≤ Σ_{n≤y} Λ(n)/n with y = x(1 + h) ≥ 2974, and
Σ_{n≤y} Λ(n)/n ≤ log y − γ_E + 0.2/log y + 0.2/log² y, by Dusart (arXiv 1002.0442 = Ramanujan J. 45 (2018), the theorem on
Σ_{p≤x} log p/p: ≤ log x + E + 0.2/log x + 0.2/log² x for x ≥ 2974, checked in the TeX) plus Σ_p Σ_{k≥2} log p/p^k = −γ_E − E (E cancels).
So |L'/L(1, χ)| ≤ X + log(1 + h) − γ_E + 0.2/X + 0.2/X² + 1/(x − 1) + (x^{−b} Z_G(L) + x^{−1/2} Z_{G_L}(L))/2, X = log x ≥ log 2974.
*Tail.* With x = (κ log q)^{1/b} the right side minus (1/b) log log q is at most
(1/b) log κ + log(1 + h) − γ_E + 0.2/X + 0.2/X² + 1/(x − 1) + Z̄_G/(2κ) + x^{−1/2} L Z̄_{G_L}/2, every term non-increasing in q
(Z̄ = the majorant of Z/L from §4; x^{−1/2}L = (κL)^{−1/(2b)}L decreases); the certificate takes the best of 28 pairs
(κ, h) ∈ {0.2, …, 2} × {0.3, 0.5, 0.7, 1}; the optimum is κ = 0.3, h = 0.7 for log q ≥ 100. For log 3 ≤ log q ≤ 25: 240 cells
[L₀, L₁], x = max(2974, (κL₀)^{1/b}), Z at L₁, best grid point per cell; worst constant 10.526 (at small q).
Context: GRH gives 2 log log q + O(1) (Ihara–Murty–Shimura); the unconditional C₅ log log q of Hexagon 2610.00256v3 has C₅ ≈ 10¹¹
(Harvest1 evaluation); here the coefficient is 1/(1 − B) = 7.8978. Sanity: for all 4135 nonprincipal characters of prime moduli ≤ 200
(L(1, χ) = −q^{−1}Σχ(a)ψ(a/q), L'(1, χ) = −log q·L(1, χ) − q^{−1}Σχ(a)γ₁(a/q) with mpmath), max |L'/L(1, χ)| = 2.168 (q = 163, real χ),
at most 9.3% of the bound.

## 8. Certificate conventions
python-flint 0.9 Arb balls at 200 bits; every inequality is a certified sign (`pos`), otherwise the run aborts. Floats only choose δ
(converted to exact rationals) and cell boundaries. Upward roundings: all bounds are upper bounds of balls; constants of the inputs
are used as stated (0.22737, 1.567, 0.137, 0.443, 1.588, 1.01624, 1.03883, 2.50637) or rounded up (B). One command:
`uv run --no-project --with python-flint --with mpmath python3 cert.py all` (14 s, ~50 MB). The sieve part of Theorem 2 (q < 100) and the
sanity checks: `uv run --no-project --with numpy python3 sanity.py` (7 min, ~0.5 GB).

**Independent float checks (not part of the proof).** xcheck.py recomputes the zero sums Z of Lemma Z and Lemma Zζ' by fine upper
Stieltjes sums with the undecomposed BMOR and Rosser bounds: the closed forms are larger (safe) by 0.02–1.6% for L ≥ 9.2 and by 27% at
L = log 9 (logs/xcheck_output.txt). ksanity.py tests Lemmas K1, K2 and the Theorem 4 kernel bound on 20000 random (β, γ, x) each
(max ratio 0.9957, 0.9991, 0.9896 ≤ 1) and inequality (1.3) at x = 10⁷ for all q ≤ 30 and all a (max deviation/bound = 0.0044)
(logs/ksanity_output.txt).

## 9. Adversarial checklist
- x < q: §6.1 (no explicit formula needed). x between q and q^{7.9}: covered by the cells (the explicit formula is used for all x ≥ 64,
  with δ = 1/2 when x^b is small compared with φ(q); no separate "trivial range" is needed).
- q ≤ 4: q = 1, 2 are separate cell kinds (only ζ, φ = 1; for q = 2 the prime 2 enters via D₂); q = 3, 4 have conductors < 9 and use
  L = log 9 in Lemma Z; Theorem 2 for q < 100 is by sieve; Theorem 3 for q ≤ e^{25.6} is trivial.
- Prime-power and other imprimitive q: only ω(q) ≤ log q/log 2 and Σ_{p|q} log p ≤ log q enter (1.2); conductors q* | q are covered
  since Lemma Z holds for every conductor ≤ q.
- Real χ: (H) covers them, so no exceptional zero; Lemma P only uses ρ ↦ 1 − ρ̄. Even χ: the trivial zero at 0 is in T(χ) (the W(0) term).
- Principal character: induced by ζ; the main term x W(1) comes only from ζ.
- Zeros at height < 5/7 (the dominant cost): counted by N(5/7) ≤ N⁺(5/7) and placed at γ = 0 with weight G(0) = n₀ (worst case).
- Steps not used: GRH, Platt's GRH verification, Trudgian's zeta count (BMOR report a gap in it), the zero density theorems.

## 10. What limits the constants; not done
- 72% of z_a comes from the ≈ 0.455 log q zeros below height 5/7, which BMOR cannot localise better. A Hadamard bound
  Σ (1 − β)/|1 − ρ|² ≤ ½ log q + O(log log q) (Riemann/Deep/Littlewood78, or Theorem 4 here) removes zeros near β = B, γ = 0, but the
  kernel envelope is flat (= n₀) up to γ ≈ 0.57 and pairs at γ ≈ 0.5 have Hadamard weight only ≈ 1.3, so a Lagrangian combination of
  the two constraints gains little (estimated < 5%); not pursued.
- Theorem 4 uses the crude Σ_{n≤x(1+h)} Λ(n)/n for the prime part; averaging over the ramp would save about log(1 + h)/2.
