# FULLPROOF — Theorem U, line by line: Steps 5c, 6b, 6c, 6d, 6e, 6f; checks of Steps 1, 3, 7

Referee: uaudit (Orca), 12 Nov 2026. Object: Riemann/Deep/DPCUncond/PROOF.md (Theorem U: D(λ) > 0 for λ > 2/3). Notation as there: L = log T, h = 2πφ/L, φ ∈ (2/3, 1), θ ∈ (1/φ − 1, 1/2), f, η, P (deg D), c_n, R_T, R*, G_T, V = V_φ, Z_off, d_ρ = β − ½, S_a, Y_a = (a+1)/L, n_a = |S_a| ≤ σ_aTL. DPC = 2610.28527v2 Gaps_for_Dirichlet.tex, Lz = 2611.14486 SmallGapsPaper.tex.

**Constants.** A ≪ B means |A| ≤ C·B, C depending only on φ, θ, f, η, P and on the auxiliary exponents named at the place (m, p, q, r, w, δ). C never depends on ε, on T ≥ T₀, or on the zero configuration beyond the unconditional facts F1–F7. Zeros are counted with multiplicity.

Corrections found (§8) are all cosmetic. None changes a statement of PROOF.md beyond a constant.

## 0. Inputs in the exact form used

- **F1** (local count; Titchmarsh Thm 9.2). N(t + 1) − N(t) ≤ C log(t + 2) for t ≥ 0.
- **F2** (Fujii; Lz l.240–249, N counts all zeros). For integers k ≥ 1, uniformly for 0 ≤ u ≤ 100/log T': ∫_0^{T'}(N(t + u) − N(t))^{2k} ≪_k T'. Applied with T' = 3T, so it covers u ≤ 90/L.
  - (a) For real m ≥ 1 and 0 ≤ u ≤ 90/L: ∫_0^{3T}(N(t+u) − N(t))^m ≪_m T. Proof: Hölder from the even integer 2k ≥ m on an interval of length 3T.
  - (b) For 0 ≤ ℓ ≤ 1: ∫_0^{11T/4}(N(t+ℓ) − N(t))^m ≪_m T(1 + ℓL)^m. Proof: write the window as M = ⌈ℓL/90⌉ consecutive windows of length ℓ/M ≤ 90/L, each shifted by at most 1. Apply Minkowski in L^m[0, 11T/4] and (a) on [0, 3T].
  - (c) (Lz l.272–283) Cut [0, 3T] into cells I_j of length λ with h/4 ≤ λ ≤ 45/L, and let n_j = #{γ ∈ I_j}. For t ∈ I_{j−1}, I_j ⊂ (t, t + 2λ], so n_j^m ≤ λ^{−1}∫_{I_{j−1}}(N(t + 2λ) − N(t))^m. Summing, Σ_j n_j^m ≪_m T/λ ≪ TL.
  - (d) Local counts at zeros. Let Ñ_γ(W) = #{zeros γ' : |γ' − γ| ≤ W}. For W ≥ 0: Σ_{γ∈[T/2,5T/2]} Ñ_γ(W)^m ≪_m TL(1 + WL)^m. Proof: use h-cells and K = ⌈W/h⌉ + 1. For γ ∈ I_j, Ñ_γ ≤ Σ_{|j'−j|≤K} n_{j'}. Then Σ_γ Ñ^m ≤ Σ_j n_j(Σ_{|j'−j|≤K} n_{j'})^m ≤ (2K+1)^{m−1}Σ_{|j−j'|≤K}(n_j^{m+1} + n_{j'}^{m+1}) (power mean, then ab^m ≤ a^{m+1} + b^{m+1}) ≤ 2(2K+1)^m Σ_j n_j^{m+1} ≪ K^m TL by (c). For W ≥ 1 the bound follows from F1, since Ñ ≪ (W + 1)L.
- **F3** (density; Lz l.221–230: Conrey 8/7 − δ, Jutila 1 − δ, uniform in σ ≥ ½; Selberg 1/4). There is κ > 0 (κ = 1 − δ, or 8/7 − δ) with N(σ, 3T) ≪_δ T^{1−κ(σ−½)}L for all σ ≥ ½. Used in the form: #{ρ : |β − ½| ≥ δ', 0 < γ ≤ 3T} ≪ TL e^{−κδ'L} for every δ' ≥ 0. The left half is covered by the symmetry ρ ↦ 1 − ρ̄.
- **F4** (M. Riesz). Let Hg(t) = π^{−1} p.v.∫ g(u)/(t − u) du. For 1 < m < ∞, ‖Hg‖_{L^m(ℝ)} ≤ C_m‖g‖_{L^m(ℝ)}, with best constant cot(π/2m) for m ≥ 2 (Pichorides).
- **F5** (Montgomery–Vaughan). ∫_A^{A+T₀}|Σ_{n≤M} a_n n^{−it}|² dt = Σ|a_n|²(T₀ + O(n)).
- **F6** (kernels). P_y(u) = π^{−1}y/(u² + y²), ∫P_y = 1, and ∫P_y(u)e^{iuλ}du = e^{−y|λ|}. Set k_d(u) := ∫_0^d P_y(u)dy = (2π)^{−1}log(1 + d²/u²). Then ∫k_d = d and ∫k_d^q = C_q d with C_q = ∫k_1^q < ∞ (scaling u = dv). Also k_d(u) ≤ d²/(2πu²), and k_d increases with d.
- **F7** For w, J ∈ ℂ: |P(w + J) − P(w)| ≤ C_P|J|(1 + |w| + |J|)^{D−1}, since |(w + J)^k − w^k| ≤ k|J|(|w| + |J|)^{k−1}.
- **F8** (5a, checked). Let b_n := Σ_{n₁n₂=n} c_{n₁}c_{n₂}. Then |b_n| ≤ 4^k‖f‖_∞² n^{−1/2}: n has at most 2k prime factors, all exponents ≤ 2, so at most 2^{2k} ordered factorisations. Hence Σ|b_n|² ≤ 16^k‖f‖_∞⁴ ∏_{T^{η/2}<p≤T^θ}(1 + 1/p + 1/p²) ≪_η 1 (Mertens). R_T² and (R*_Y)² have length ≤ T^{2θ} = o(T). F5 gives, for 0 ≤ Y ≤ 1:
  ∫_{T/2}^{5T/2}|R_T|⁴ = (2T + O(T^{2θ}))Σ|b_n|² ≪ T,  ∫_{T/2}^{5T/2}|R*(u − iY)|⁴du = Σ|b_n|²n^{2Y}(2T + O(n)) ≪ T e^{4θYL}.

## 1. Step 5c (∫|Re V|^m ≪ T), complete

**Claim.** For real m ≥ 2 and every interval I' ⊂ [3T/4, 9T/4] (used for [T, 2T] and I := [7T/8 − 1, 17T/8 + 2]):
(i) ∫_{I'}|Re V|^m ≪_m T; (ii) ∫_{I'}|Im V|^m ≪_m T; (iii) ∫_{I'}|F̃|^m ≪_m T for the smeared count F̃ below.
The constants depend on m, φ and κ only (any κ > 0). Hence ∫_{I'}|V|^m ≪ T and ∫_{I'}|P(V)|⁴ ≪ T.

(a) **Even product.** X(z) := ξ(½ + iz) is entire, even, real on ℝ, of order 1, and X(0) ≠ 0. It is a function of z² of order ½, so Hadamard gives X(z) = X(0)∏_{γ_ρ>0}(1 − z²/z_ρ²) with no exponential factors, z_ρ := γ − i d_ρ, over all zeros with γ > 0 (both ρ and 1 − ρ̄). The product converges absolutely. For real τ off the ordinates, log|X(τ)| = log|X(0)| + Σ_ρ log|1 − τ²/z_ρ²|. Hence Re V(t) = π^{−1}Σ_ρ[log|1 − τ²/z_ρ²|]_{τ=t−h/2}^{t+h/2}.

(b) **Far and mirror parts are O(1)** for t ∈ [3T/4, 9T/4]. Write the factor as log|z_ρ − τ| + log|z_ρ + τ| − 2log|z_ρ|.
- γ > 5T/2: |d/dτ log|1 − τ²/z_ρ²|| = |Re(2τ/(z_ρ² − τ²))| ≤ 2τ/(|z_ρ|² − τ²) ≪ T/γ², since τ ≤ 9T/4 + 1 < 0.95γ. The τ-increment is h. By F1 and partial summation, Σ_{γ>5T/2} hT/γ² ≪ hL = 2πφ.
- γ < T/2: the derivative is Re(1/(τ − z_ρ) + 1/(τ + z_ρ)) ≤ 4/T + 2/T. Over N(T/2) ≪ TL zeros this gives O(hL).
- γ ∈ [T/2, 5T/2]: the mirror term has |d/dτ log|z_ρ + τ|| ≤ 1/T, which gives O(hL). The constant −2log|z_ρ| cancels.
- Hence Re V(t) = π^{−1}Σ_{γ∈[T/2,5T/2]} log|(t + h/2 − z_ρ)/(t − h/2 − z_ρ)| + E(t), |E| ≤ C_φ. Only F1 and |d_ρ| < ½ were used, not RH. This is DPC l.633–641 without RH.

(c) **Hilbert form.** For real u and d ≥ 0: log|u + id| = (P_d * log|·|)(u). Proof: log|z| = Re Log z is harmonic in Im z > 0, log-integrable on ℝ and O(log|z|) at ∞, so it equals the Poisson integral of its boundary values. Hence log|(t + h/2 − z_ρ)/(t − h/2 − z_ρ)| = (P_{|d|} * ψ)(t − γ) with ψ(u) = log|(u + h/2)/(u − h/2)| = πH[1_{[−h/2,h/2]}](u). P_{|d|}* and H are Fourier multipliers (e^{−|d||ξ|}, −i sgn ξ), so they commute on L¹ ∩ L². Therefore
  Re V = H[F̃] + E on [3T/4, 9T/4], F̃(u) := Σ_{γ∈[T/2,5T/2]}(P_{|d_ρ|} * 1_{[−h/2,h/2]})(u − γ) (P_0* = identity).
F̃ is a finite sum of bounded L¹ functions. By F4, ‖Re V‖_{L^m(I')} ≤ C_m‖F̃‖_{L^m(ℝ)} + C_φ(3T)^{1/m}.

(d) **Kernel bound.** Let D := max(d, h). Then (P_d * 1_{[−h/2,h/2]})(x) ≤ 4(h/D)·D²/(x² + D²).
- d ≤ h: P_d*1 ≤ 1. For |x| ≥ h it is ≤ h·sup_{|v|≤h/2}P_d(x − v) ≤ 4hd/(πx²) ≤ 4h²/(πx²). Both cases are ≤ 4h²/(x² + h²).
- d > h: for |v| ≤ h/2 < d/2, (x − v)² + d² ≥ (x² + d²)/2, so P_d*1 ≤ (2h/π)d/(x² + d²) = (2/π)(h/d)d²/(x² + d²).

(e) **Dyadic bins.** B'_0 := {γ ∈ [T/2, 5T/2], D_ρ < 2h}. For j ≥ 1, B'_j := {D_ρ ∈ [2^jh, 2^{j+1}h)}. These are nonempty only for 2^jh < ½. Put D_j := 2^jh. For ρ ∈ B'_j, h/D_ρ ≤ 2^{−j} and D_ρ²/(x² + D_ρ²) ≤ 4D_j²/(x² + D_j²). So the bin sum is F̃_j(u) ≤ 16·2^{−j}Σ_{B'_j} D_j²/((u − γ)² + D_j²).
With n_j(x) := #{ρ ∈ B'_j : |x − γ| ≤ D_j}: ∫1_{|x−γ|≤D_j}P_{D_j}(u − x)dx = π^{−1}[arctan(s + 1) − arctan(s − 1)], s = (u − γ)/D_j, and this is ≥ π^{−1}(1 + s²)^{−1}. Hence Σ_{B'_j} D_j²/((u−γ)² + D_j²) ≤ π(P_{D_j} * n_j)(u), and by Young (‖P_{D_j}‖₁ = 1), ‖F̃_j‖_m ≤ 16π2^{−j}‖n_j‖_m.
Cauchy–Schwarz gives ‖n_j‖_m^m ≤ (∫n_j)^{1/2}(∫n_j^{2m−1})^{1/2}.
- ∫n_j = 2D_j|B'_j|, with |B'_j| ≪ TL e^{−κ2^jhL} = TL e^{−c₀2^j}, c₀ = 2πφκ (F3; for j = 0 use |B'_0| ≪ TL). So ∫n_j ≪ 2^j e^{−c₀2^j}T.
- n_j(x) ≤ N(x + D_j) − N(x − D_j) and D_j < ½, so F2b gives ∫n_j^{2m−1} ≪ T 2^{j(2m−1)}.
Hence ‖n_j‖_m ≪ 2^j T^{1/m}e^{−c₀2^{j−1}/m}, ‖F̃_j‖_m ≪ T^{1/m}e^{−c₀2^{j−1}/m}, and ‖F̃‖_m ≤ Σ_j‖F̃_j‖_m ≪_{m,φ,κ} T^{1/m}. This proves (iii) and then (i).
(ii): Im V = N(t + h/2) − N(t − h/2) ≥ 0 (Step 1), so ∫_{I'}|Im V|^m ≤ ∫_0^{3T}(N(t + h) − N(t))^m ≪ T (F2a). ∎

## 2. Step 6b (representations), complete

(b1) R_T(x − iy) = (P_y * R_T)(x) for y ≥ 0. By F6, P_y * e^{−iλ·} = e^{−iλx}e^{−y|λ|}; with λ = log n ≥ 0 this is n^{−i(x−iy)}. The sum is finite.
(b2) R*(x − iy) = (P_{Y−y} * R*_Y)(x) for 0 ≤ y ≤ Y, R*_Y(u) := R*(u − iY) = Σc̄_n n^Y n^{iu}. Indeed P_{Y−y} * e^{iλ·} = e^{iλx}e^{−(Y−y)λ}, giving Σc̄_n n^y n^{ix} = R*(x − iy).
(b3) **ξ̂.** On 𝒰 := {Re s > ½, T/2 < Im s < 5T/2} put log ξ̂ := log ξ + Σ_{ρ∈Z_off}[Log(s − 1 + ρ̄) − Log(s − ρ)] (principal Logs, i.e. the horizontal branches).
- Each cut {Im s = γ, Re s < β} of ρ ∈ Z_off: log ξ jumps by 2πi·mult there, and −Log(s − ρ) cancels the jump (Z_off is a multiset).
- Log(s − 1 + ρ̄) is cut only in Re s < 1 − β < ½.
- So log ξ̂ is analytic on 𝒰, and ξ̂ = ξ∏(s − 1 + ρ̄)/(s − ρ) is zero-free there.
- V̂(t) := π^{−1}[log ξ̂(½ + i(t + h/2)) − log ξ̂(½ + i(t − h/2))] is analytic on (T/2 + 1, 5T/2 − 1) × (−1, 0).
- On Re s = ½, |(s − 1 + ρ̄)/(s − ρ)| = 1, so Re V̂ = Re V. The arguments of the pair (s − ρ)(s − 1 + ρ̄) = −d² − u² cancel. So Im V̂ = the on-line count + 2Σ_{Z_off}(P_d * 1_{[−h/2,h/2]})(t − γ) = F̃(t) on [3T/4, 9T/4] (zeros outside [T/2, 5T/2] do not reach these windows).
- On Im t = −1 (Re s = 3/2): |V| ≤ φ + O(1/L) (DPC l.360–380, the Euler series; no RH). Each ρ-term is ∫_{s₋}^{s₊}(2β − 1)ds/((s − ρ)(s − 1 + ρ̄)), with Re(s − ρ) ≥ ½ and Re(s − 1 + ρ̄) ≥ 1, so no cut is crossed. Since |s − ρ||s − 1 + ρ̄| ≥ (1 + (τ − γ)²)/4 and 1 + (u − γ)² ≤ 2(1 + (τ − γ)²) for |τ − u| ≤ h/2, its modulus is ≤ 16hd/(1 + (u − γ)²). With F1, Σ_{Z_off} ≤ C_φ. So |V̂| ≤ C₀(φ) on Im t = −1.
(b4) **Box representation** (point (iv) of the brief). Fix a slit point x ∈ [T − 1, 2T + 1] and 0 < y < 1.
- Sides: x₁ ∈ [x − T/8, x − T/8 + 1] and x₂ ∈ [x + T/8, x + T/8 + 1], each at distance ≥ c/L from every shifted ordinate γ' ± h/2 of every zero. This is possible: by F1 a unit interval holds ≤ C L shifted ordinates, so it has a gap ≥ 1/(CL + 1); c = c(F1).
- B := (x₁, x₂) × (−1, 0) ⊂ (7T/8 − 1, 17T/8 + 2) × (−1, 0), inside the analyticity region of V̂. [PROOF.md writes 17T/8 + 1; +2 is right; cosmetic.]
- Majorants:
  - Top edge: by the Step 1 majorant, valid for log ξ and for each finite Log (|Log(s − ζ)| ≤ |log|Im s − Im ζ|| + π + log 3 when |s − ζ| ≤ 3), |V̂(x' − iy')| ≤ T^{O(1)} + C Σ_{w}|log|x' − w|| for 0 ≤ y' ≤ 1. The sum runs over the ≪ L shifted ordinates w within 2 of x'. This majorant is uniform in y'.
  - Sides: every |x_i ± h/2 − γ'| ≥ c/L, which gives |V̂| ≤ T^{O(1)}. The Z_off sum is ≤ CL log(L/c) + Σ_{far}1/|τ − γ| ≤ CL².
  - Bottom: |V̂| ≤ C₀.
- Representation: for B_δ := (x₁, x₂) × (−1, −δ), V̂ is continuous on the closure of B_δ, so V̂(z) = ∫_{∂B_δ}V̂ dω^{B_δ}_z. As δ ↓ 0, the top-edge values converge in L¹ by dominated convergence with the uniform log majorant. The top harmonic-measure densities are ≤ P_{y−δ}(x − ·) → P_y(x − ·) uniformly. The other edges are bounded and converge. So V̂(z) = ∫_{∂B}V̂ dω^B_z, with log singularities allowed on the top edge.
- Estimate:
  - Top: for E in the top edge, ω^B(z, E) ≤ ω^{Π}(z, E), Π = {Im t < 0}, by the maximum principle. Its density is P_y(x − x'). So the top part is ≤ (P_y * v)(x), v := |V̂|_ℝ|·1_I.
  - Bottom: ≤ C₀.
  - Sides: harmonic measure ≤ that of the end of a half-strip of height 1 at distance ≥ T/8 − 1, which is ≤ (4/π)e^{−π(T/8 − 1)}. Times T^{O(1)} this is o(1).
- **Result:** |V̂(x − iy)| ≤ (P_y * v)(x) + C₁, C₁ = C₀ + 1, and ∫_ℝ v^m = ∫_I|V̂|^m ≤ 2^m(∫_I|Re V|^m + ∫_I|F̃|^m) ≪_m T by §1 with I' = I. ∎
(b5) **Δ.** On the slit domain V = V̂ + Δ, Δ = π^{−1}Σ_{ρ'∈Z_off} ℓ_{ρ'}, ℓ_{ρ'}(t) = [Log(s − ρ') − Log(s − 1 + ρ̄')]_{s₋}^{s₊}. With t' = γ' − h/2 − id': s₊ − ρ' = i(t − t'), s₊ − 1 + ρ̄' = i(t − t̄'), and for s₋ the same with t' + h, t̄' + h.
- **Far** (δ := |Re t − γ'| ≥ 2h): on the segment s₋ → s₊, |Im s − γ'| ≥ δ − h/2 ≥ 3δ/4, so no cut is crossed. ℓ_{ρ'} = ∫(2β' − 1)ds/((s − ρ')(s − 1 + ρ̄')), so |ℓ_{ρ'}| ≤ 32hd'/(9δ²). For a slit point of ρ, |γ_ρ − γ'| ≤ δ + h/2 ≤ 5δ/4, so 1 + ((γ_ρ − γ')/h)² ≤ 2(δ/h)². Hence |ℓ_{ρ'}| ≤ (8/(2πφ))(d'L)(1 + ((γ_ρ − γ')/h)²)^{−1} and Σ_far|ℓ| ≤ (2/φ)Ψ_ρ. [PROOF.md: 32Ψ_ρ; any constant C_φ suffices.]
- **Near** (δ < 2h): for A ∈ ℂ and h > 0, |log|A| − log|A − h|| ≤ log 2 + log⁺(h/|A|) + log⁺(h/|A − h|). Proof: if |A − h| ≥ h then |A| ≤ 2|A − h|; else |A| ≤ 2h. Each principal arg is bounded by π. So |ℓ_{ρ'}(t)| ≤ 2log 2 + 4π + Σ_{w∈{t', t'+h, t̄', t̄'+h}} log⁺(h/|t − w|), also for the one-sided values on slits. Pairing the two logs is necessary: separately each would cost log(1/h) ≍ log L.

## 3. Step 6c (sparse kernels), complete

Fix a and the sign + (− is identical), x_ρ := γ_ρ − h/2, Y := Y_a, d_ρ ∈ [a/L, Y). Put K_a := Σ_{S_a}k_{d_ρ}(· − x_ρ) and K'_a := Σ_{S_a}k_Y(· − x_ρ). Since k_d ≤ k_Y, K_a ≤ K'_a.

**Lemma C.** For 1 ≤ q < ∞ and 1 < w < ∞: ‖K'_a‖_{L^q(ℝ)}^q ≪_{q,w,φ} (a+2)^q σ_a^{1−1/w}T.
Proof. Split K'_a(u) = K_near + K_far, where "near" means ρ with |u − x_ρ| ≤ Y.
- **Near.** Let M(u) := #near. Power mean gives K_near(u)^q ≤ M(u)^{q−1}Σ_near k_Y(u − x_ρ)^q. If ρ is near u, every other near ρ' has |γ' − γ_ρ| ≤ 2Y, so M(u) ≤ N_ρ := #{ρ' ∈ S_a : |γ' − γ_ρ| ≤ 2Y}. Integrating, ∫K_near^q ≤ C_qY Σ_{S_a}N_ρ^{q−1} (overlaps are covered by this).
- **Far.** On 2^lY < |v| ≤ 2^{l+1}Y, k_Y(v) ≤ 4^{−l}/(2π). So K_far ≤ (2π)^{−1}Σ_{l≥0}4^{−l}M_l with M_l(u) := #{ρ : |u − x_ρ| ≤ 2^{l+1}Y}. Then ∫M_l^q = ∫M_l^{q−1}Σ_ρ1_{|u−x_ρ|≤2^{l+1}Y} ≤ 2^{l+2}YΣ_ρ(N^{(l)}_ρ)^{q−1}, with N^{(l)}_ρ the count of S_a within 2^{l+2}Y of γ_ρ.
- **Hölder on S_a**, with g ≥ 0 extended by "all zeros": Σ_{S_a}g ≤ n_a^{1−1/w}(Σ_{γ∈[T,2T]}g^w)^{1/w}. With N^{(l)}_ρ ≤ Ñ_{γ_ρ}(2^{l+2}Y) and F2d (1 + 2^{l+2}(a+1) ≤ 2^{l+3}(a+2)): Σ_{S_a}(N^{(l)})^{q−1} ≪ (σ_aTL)^{1−1/w}(TL)^{1/w}(2^{l+3}(a+2))^{q−1}. The same holds for N_ρ (l = −1).
- Hence ∫K_near^q ≪ (a+2)^qσ_a^{1−1/w}T and ‖M_l‖_q ≪ 2^l(a+2)(σ_a^{1−1/w}T)^{1/q}. By Minkowski, ‖K_far‖_q ≪ Σ_l 2^{−l}(a+2)(σ_a^{1−1/w}T)^{1/q}. ∎
[PROOF.md states the same bound. The window of Ñ in F2d may be any O((a+2)/L).]

**Poisson tails** (the integrals in 6e run over ℝ, the fourth moments over [T/2, 5T/2]). For u ∉ [T/2, 5T/2] and x_ρ ∈ [T − 1, 2T + 1], |u − x_ρ| ≥ T/3. So ∫_{ℝ∖[T/2,5T/2]}K'_a ≤ n_a·3Y²/(πT) ≤ σ_a(a+1)²/L. Also sup|R_T| ≤ 2‖f‖_∞T^{θ/2} and sup|R*_Y| ≤ 2‖f‖_∞T^{θ/2}e^{θ(a+1)} (n^Y ≤ T^{θY}). The tails are therefore ≪ T^{3θ/2}σ_a(a+1)²/L and ≪ T^{3θ/2}e^{3θ(a+1)}σ_a(a+1)²/L. v vanishes off I, so it has no tail.

**Jensen–Fubini.** For a slit point x_ρ and ψ on ℝ, write (P_y * |ψ|)^p ≤ P_y * |ψ|^p (p ≥ 1). Then ∫(P_y * |ψ|)^p(x_ρ)dμ_a = Σ_ρ∫_0^{d_ρ}(…)dy ≤ ∫|ψ|^pK_a. For R*, the kernel ∫_0^{d_ρ}P_{Y−y}dy = ∫_{Y−d_ρ}^Y P_s ds ≤ k_Y gives K'_a.

## 4. Step 6d (Δ and Ψ moments), complete

(d1) For w ∈ ℂ and real x: ∫_0^∞ log⁺(h/|x − iy − w|)^m dy ≤ ∫_ℝ log⁺(h/|y + Im w|)^m dy = 2h·Γ(m + 1).
(d2) **Correction (cosmetic).** The near ρ' of (b5) satisfy |γ' − γ_ρ| < 5h/2. PROOF.md's Ñ_ρ (window 4Y_a + h) contains them only if 4(a+1) + 2πφ ≥ 5πφ. For a = 0 and φ > 4/(3π) ≈ 0.42 this fails. **Fix:** Ñ_ρ := Ñ_{γ_ρ}(4Y_a + 3h). This window still contains the zeros at ordinate γ_ρ ± h that share a merged cut. F2d gives the same moments, Σ_{γ∈[T,2T]}Ñ_γ^p ≪_p (a+2)^pTL.
For 0 < y < d_ρ: |Δ_L(x_ρ − iy)| ≤ C[Ñ_ρ + Σ_{near ρ'}Σ_w log⁺(h/|x_ρ − iy − w|) + Ψ_ρ], with ≤ 4Ñ_ρ log terms. By the power mean and (d1):
  ∫_0^{d_ρ}|Δ_L|^m dy ≤ C_m[d_ρ(Ñ_ρ^m + Ψ_ρ^m) + (4Ñ_ρ)^{m−1}·4Ñ_ρ·2hΓ(m+1)] ≤ C'_m(d_ρ + h)(1 + Ñ_ρ + Ψ_ρ)^m.
(d3) **Ψ moments.** Use h-cells I_j on [T/2, 5T/2]. For γ ∈ I_j and γ' ∈ I_{j'}: 1 + ((γ − γ')/h)² ≥ (1 + (j − j')²)/3 (2k² − 6k + 5 > 0). Also 1 + d'L ≤ a' + 2 on B_{a'}. So Ψ_γ ≤ 3Σ_{a'}(a' + 2)(J * b_{a'})_j with J_k = (1 + k²)^{−1} and ‖J‖₁ = π coth π.
- Minkowski in ℓ^p(n_j): (Σ_γ Ψ_γ^p)^{1/p} ≤ 3Σ_{a'}(a'+2)(Σ_j n_j(J*b_{a'})_j^p)^{1/p}.
- Cauchy–Schwarz, then Young on ℤ: Σ_j n_j(J*b)_j^p ≤ (Σn_j²)^{1/2}‖J‖₁^p(Σ_j b_j^{2p})^{1/2}.
- Cauchy–Schwarz again with b_j ≤ n_j: Σ_j b_j^{2p} ≤ |B_{a'}|^{1/2}(Σ_j n_j^{4p−1})^{1/2}.
- F2c and F3 (|B_{a'}| ≪ TLe^{−κa'}) give Σ_j n_j(J*b_{a'})_j^p ≪_p TL e^{−κa'/4}. Hence Σ_γ Ψ_γ^p ≪_{p,κ} TL, using Σ_{a'}(a'+2)e^{−κa'/(4p)} < ∞.
- Hölder on S_a: Σ_{S_a}Ψ^p ≪ σ_a^{1−1/w}TL and Σ_{S_a}Ñ^p ≪ (a+2)^pσ_a^{1−1/w}TL.
(d4) ∫|Δ_L|^m dμ_a ≤ C(Y_a + h)Σ_{S_a}(1 + Ñ + Ψ)^m ≪ ((a+2)/L)·(a+2)^mσ_a^{1−1/w}TL = (a+2)^{m+1}σ_a^{1−1/w}T. Similarly ∫Ñ^m dμ_a ≤ Y_aΣ_{S_a}Ñ^m ≪ (a+2)^{m+1}σ_a^{1−1/w}T. [PROOF.md: (a+2)^C. The exact power is irrelevant.]

## 5. Steps 6a, 6e, 6f, complete

(6a) Take a merged cut x of class a, with deepest zero ρ ∈ S_a (any choice if tied). Distinct cuts have distinct (ρ, ±), since x = x_ρ^±.
- |G_T(x − iy)| = (√(2π)H)^{−1}|∫e^{−(x−iy−u)²/(2H²)}du| ≤ e^{y²/(2H²)}G_T(x) ≤ 2.
- |J_x(y)| ≤ 2(m₊ + m₋) ≤ 2Ñ_ρ.
- F7 gives |hairpin_x| ≤ 4C_P∫_0^{d_ρ}|R_T||R*|(1 + |V_L| + 2Ñ_ρ)^D dy.
- Summing over cuts: 𝓗_{T,a} ≤ 4C_P Σ_±∫|R_T||R*|(1 + |V_L| + 2Ñ)^D dμ_a^±.
- Cuts whose deepest zero has γ ∉ [T, 2T] have x ≤ T + h or x ≥ 2T − h, so |G_T(x − iy)| ≤ 2e^{−L²/33}. Their total is o(1), uniformly in the configuration, once the Step 1 majorant is read in the y-variable along the slit. With s = ½ + y + i(x ± h/2), |s − ρ''| ≥ |Re(s − ρ'')| = |y − d''| (signed d''). So |log ξ(s)| ≤ CTL² + CΣ_{|γ''−Im s|≤2}(|log|y − d''|| + C), and ∫_0^d(…)^D dy ≤ (CL)^D·C_D, whatever the ordinates. With |R_TR*| ≤ 4‖f‖²_∞T^{2θ} and ≪ TL cuts, the total is ≪ e^{−L²/33}T^{O(1)} = o(1). [PROOF.md's "integrable log singularities" is right in this reading. The x-variable reading |log|x ± h/2 − γ''|| would not be uniform.]
(6e) Hölder (3, 3, 3) on μ_a := μ_a^+ gives 𝓗^+_{T,a} ≪ (A₁A₂A₃)^{1/3}.
- **A₁** = ∫|R_T|³dμ_a ≤ ∫_ℝ|R_T|³K_a ≤ ‖R_T‖³_{L⁴[T/2,5T/2]}‖K_a‖₄ + tail ≪ T^{3/4}(a+2)(σ_a^{1−1/w}T)^{1/4} + T^{3θ/2}σ_a(a+1)² ≪ T(a+2)²σ_a^{(1−1/w)/4}.
- **A₂** = ∫|R*|³dμ_a ≤ ∫_ℝ|R*_{Y_a}|³K'_a ≪ (Te^{4θ(a+1)})^{3/4}(a+2)(σ_a^{1−1/w}T)^{1/4} + tail ≪ Te^{3θ(a+1)}(a+2)²σ_a^{(1−1/w)/4} (F8 with Y = Y_a).
- **A₃** = ∫(1 + |V_L| + 2Ñ)^{3D}dμ_a, with |V_L| ≤ (P_y * v)(x_ρ) + C₁ + |Δ_L| (b4, b5). Then A₃ ≪ μ_a(ℝ²) + ∫v^{3D}K_a + ∫|Δ_L|^{3D}dμ_a + ∫Ñ^{3D}dμ_a.
  - μ_a(ℝ²) = Σd_ρ ≤ (a+1)σ_aT.
  - ∫v^{3D}K_a ≤ ‖v^{3D}‖_r‖K_a‖_{r'} ≪ T^{1/r}(a+2)(σ_a^{1−1/w}T)^{1−1/r}.
  - The last two are ≪ (a+2)^{3D+1}σ_a^{1−1/w}T by (d4).
  - So A₃ ≪ T(a+2)^{3D+1}σ_a^{(1−1/w)(1−1/r)} (σ_a ≤ C, exponents ≤ 1).
- **Product.** 𝓗_{T,a} ≪ T(a+2)^{D+2}e^{θ(a+1)}σ_a^E, with E = (1 − 1/w)[¼ + ¼ + (1 − 1/r)]/3 = (1 − 1/w)(½ − 1/(3r)). ✓ (point (vii))
(6f) Put σ_a := min(n_tot/(TL), Ce^{−κa}).
- n_tot := #{ρ ∈ Z_off : γ ∈ [T, 2T]} ≤ εTL + O(L) by Step 2: runs of r ≥ 2 equal ordinates hold r − 1 ≥ r/2 zero gaps. So σ_a ≤ 2ε for T ≥ T₁(ε).
- n_a ≤ N(½ + a/L, 2T) ≪ TLe^{−κa} (F3).
- With integer A ≥ 0 and E' := κE − θ > 0: 𝓗_T/T ≪ (A+1)(A+2)^{D+2}e^{θ(A+1)}(2ε)^E + Σ_{a>A}(a+2)^{D+2}e^{θ(a+1)}C^Ee^{−κEa} ≪ (A+2)^{D+3}e^{θA}ε^E + (A+2)^{D+2}e^{−E'A}.
- A = ⌊c log(1/ε)⌋ with 0 < c < E/θ gives Φ(ε) ≪ (log(1/ε) + 3)^{D+3}(ε^{E−cθ} + ε^{cE'}) → 0.
- **Order of choices** (all constants fixed before ε): φ, θ → g, q, b₀, P, D, δ₀ (Step 0, DPC l.743–765) → f, η (DPC Lemma 6) → δ (κ = 1 − δ, or Conrey's 8/7 − δ), then w, r with κE > θ (possible since θ < ½ and E ↑ ½) → all C → ε₀ (Φ(ε₀) ≤ δ₀·lim‖c‖²/2 and ε₀ ≤ c₂/(4πφ)) → T₀(ε₀). ✓

## 6. Checks of Steps 1, 3, 7 (points (i)–(iii))

**Step 1.**
- With horizontal Logs, an off-line pair has constant total argument on the line: (s − ρ)(s − 1 + ρ̄) = −d² − u² < 0. Its only change is the 2π jump of Log(s − ρ) where its cut crosses the line at ½ + iγ. So arg ξ(½ + iτ) = πN(τ) + const and Im V = N(t + h/2) − N(t − h/2), with an off-line pair counted twice. ✓ [PROOF.md's "1 − ρ̄ adds 0" is right as a jump count: the pair's continuous parts cancel.]
- The majorant: |s − ρ'| ≥ |Im s − γ'|, and every principal arg is ≤ π. DPC's RH uses in Lemma 1 are the hypothesis (l.158), the zero-free continuation region, and "non-real zeros of X are at distance ≫ T". They are replaced by the horizontal convention and by |Log(ρ̃ − z)| ≤ |log|Re(ρ̃ − z)|| + π + O(1). ✓

**Step 3, identity (H)** (point (i)).
- Let Ω_δ := (a, b) × (−1, −δ) ∖ ∪slits. Ω_δ is simply connected: each slit is a vertical segment attached to the top edge.
- F = G_TR_TR*P(V) is analytic on Ω_δ, with one-sided continuous values on each slit and log singularities only at tips, intermediate tips and slit tops that meet a shifted on-line ordinate. Small circles around these contribute O(r|log r|^D) → 0.
- Positive orientation, walking the top edge from b to a: down the right side of a slit, then up the left side. This gives −i∫(F_R − F_L)dy per cut.
- δ ↓ 0: the top edge converges by dominated convergence (Step 1 majorant, uniform in δ). The hairpin integrals converge, since their integrands are log-integrable at y = 0.
- The jump is the total J_x(y) = 2i(m₊ − m₋): right of x ⇔ Im s₊ > γ for the + term (+2πi), and Im s₋ > γ for the − term (−2πi).
- **Numerical falsification (slit_identity.py, logs/slit_identity*.json).** Toy: X = e^{αs}∏(s − ρ)^m, complex cubic P, Möbius-type R of length 30, Gaussian weight. 9 configurations: single; coincident opposite cuts (γ₂ − γ₁ = h) with equal depths (J ≡ 0) and unequal depths; two off-line zeros at one ordinate; a double off-line zero (J = ±4i); an on-line zero at the off-line ordinate; an on-line shifted ordinate at a slit top (corner); a 6-zero cluster mixing all of these; depths 0.45–0.47 near the box sides.
  - Residual of (H): 2.5·10⁻¹³ to 1.5·10⁻¹¹ (scipy). The corner cases give 5.6·10⁻⁹ and 1.0·10⁻⁹ (tanh-sinh).
  - Without hairpins the residual is 3·10⁻⁴ to 50.
  - With the rejected rule (separate increments per zero), the residual is 0.16, 10.6 and 6.6 in the three coincident configurations, and agrees with (H) where no cuts coincide.
  - So the merged-jump rule is necessary and sufficient, and the orientation is right.

**Post-shift main term** (point (ii)), DPC TeX re-read.
- Section 3's standing hypothesis (l.233: zeros above Y₀ on the line) enters Prop 1 (l.262–388) only through Lemma 1 (boundary values, the side bounds at l.335–340) and through the move "By Cauchy's theorem … to Im t = −1" (l.341). Step 1 replaces the first, (H) the second. [R450's list of DPC's RH uses (l.57, 59–62, 140, 607, 610, 629) omits l.158 and l.233; these are exactly the uses Steps 1 and 3 replace.]
- After it: (1) V = G + π^{−1}Σβ(k)k^{−it} on Im t = −1 (Re s = 3/2). This is the horizontal branch: β < 1, so no cut reaches Re s ≥ 1.
- (2) Only g G^m e^{−itω} is shifted to Im t = ∓H. G is the completion-factor difference, analytic for |Im t| ≤ H (Γ-poles and s(s − 1) have |Im| ≍ T). Stirling holds since |arg s| ≤ π/2 + O(1/L).
- (3) Diagonal terms via ‖|K|‖ = O(1) and the coefficient support ≤ T^θ < T/L^{20} (l.700–704 checks exactly this).
- (4) Lemma 6 (l.686–741) uses Prop 1, Mertens and Laplace transforms. It does not use Lemma 5 (Rodgers, RH).
- No RH, any θ < 1, DPC's rough Möbius coefficients. ✓

**Uniformity in ε** (point (iii)).
- Step 5 (§1, F8) never uses Step 2. Its constants are those of F1–F5.
- In Step 6, ε enters only through σ_a ≤ 2ε (T ≥ T₁(ε)), and all implied constants are ε-free.
- In Step 7, c₂ = (δ₁/C)⁴ with C from F8 and §1, and δ₁ = δ₀ lim‖c‖²/2. So ε-free. ✓

**Step 7.**
- Re∫_ℝ G_T|R_T|²P(V) splits into: [T, 2T] ∖ E_T, where Im V ∈ {0, 1} and Re P(V) ≥ 0 (Re P(x) = q², Re P(x + i) = 0); E_T; and ℝ ∖ [T, 2T], which is O(e^{−cL²}T^{O(1)}).
- Hölder (2, 4, 4) uses only ∫_T^{2T}|R_T|⁴ ≪ T (F8) and ∫_T^{2T}|P(V)|⁴ ≪ T (§1).
- Centres covering a pair: measure h − (γ_{n+1} − γ_n) ≤ h.
- D(λ): choose 2/3 < φ < min(λ, 1) and apply U on [T/2, T]. ✓

## 7. Numerical checks of the quantitative steps (logs/)

| Test | Setup | Outcome |
|---|---|---|
| (H), 9 configurations | above | residual ≤ 1.5·10⁻¹¹ (corner cases ≤ 5.6·10⁻⁹), against 3·10⁻⁴–50 without hairpins; separate increments fail on every coincident cut |
| Lemma C on actual zeros (sparse_kernel.py) | zeros in [2000, 4000] (1960 of them); subsets of 5% and 30% (random, tightest local clusters, contiguous block); a = 0, 2, 5; depths mixed or full (K'_a); q = 2, 4 | ‖K‖_q^q/((a+2)^qσT) ≤ 0.113 in all 72 runs, even with σ in place of the claimed σ^{1−1/w}; ‖K‖_q^q ≤ 1.01× the near-part bound C_qYΣN^{q−1}; max N_ρ ≤ 5 |
| 5c on actual zeros (vmoments.py) | T = 2000, φ = 0.7, [T, 2T]; far part via mpmath log\|ξ\| at 120 points | Re V_exact − H[F̃]: mean −0.106, sd 0.114, max 0.315 (O(1)). ∫\|Re V\|^m/T = 0.295, 0.384 (m = 2, 4). M. Riesz ratio 0.50, 0.70 < cot(π/2m) = 1, 2.41 |
| 5c with planted off-line pairs (vmoments_planted.py) | J = [T, T + 400]; 5% or 30% of the zeros replaced by pairs at d = 1/L, 3/L, 0.45 | per unit length, \|Re V\|^m (m = 2/4/6) changes from 0.31/0.41/1.28 to 0.30–0.38/0.38–0.48/1.15–1.28. Smeared Im V̂ ≤ 0.82/1.46/3.73; sharp Im V ≤ 1.10/3.98/23.2. All O(1) |

## 8. Corrections (all cosmetic; none affects Theorem U)

1. 6b: I = [7T/8 − 1, 17T/8 + 2], not +1. §1 is therefore needed on I ⊄ [T, 2T]; it holds on any I' ⊂ [3T/4, 9T/4].
2. 6c/6d: the window of Ñ_ρ must be ≥ 5h/2 to cover the near ρ' of Δ. Use 4Y_a + 3h; the moment bounds are unchanged.
3. 6d: the far constant 32 becomes 2/φ, and the exponents (a+2)^C become (a+2)^{m+1}.
4. 6b near bound: valid only with the two logs of each ρ' paired. Unpaired, each would cost log(1/h) ≍ log L, which would add a power of log L to A₃.
5. Step 1 wording: "1 − ρ̄ adds 0" is a statement about jumps; the continuous arguments of the pair cancel.
6. Step 3 bookkeeping: zeros with slit in [a, b] but γ ∉ [T/2, 5T/2] (possible if a < T/2 + h/2) are not in Z_off. They contribute o(1) like all slits with γ ∉ [T, 2T]. Alternatively take a ≥ T/2 + 1.
7. Step 3, slits with γ ∉ [T, 2T]: the log majorant must be read in y (|s − ρ''| ≥ |y − d''|) for uniformity in the configuration (§5, 6a).
