# PROOF — small gaps below 2/3 unconditionally (Durkan–Pearce-Crump's two-line moments without RH)

**Theorem U.** Let 0 < γ₁ ≤ γ₂ ≤ … be the ordinates of the nontrivial zeros of ζ in the upper half-plane, counted with multiplicity (an off-line pair ρ, 1 − ρ̄ gives two equal ordinates), as in Lamzouri 2611.14486 l.177. For every φ ∈ (2/3, 1) there is c(φ) > 0 such that for all large T

  𝒩(T) := #{n : T ≤ γ_n ≤ γ_{n+1} ≤ 2T, γ_{n+1} − γ_n < 2πφ/log T} ≥ c(φ) T log T.

Hence D(λ) > 0 for every λ > 2/3 (D as in 2611.14486 l.183): choose 2/3 < φ < min(λ, 1) and apply Theorem U on [T/2, T]; for large T, 2πφ/log(T/2) < 2πλ/log T, so #{0 < γ_n ≤ T : γ_{n+1} − γ_n ≤ 2πλ/log T} ≥ 𝒩(T/2) ≥ c(φ)(T/2)log(T/2) ≫ N(T). Before: λ ≥ 0.768 (ours, explore40), 0.7851 (ours), 0.7892 (Lamzouri). Under RH: λ > 1/2 (DPC Thm B). The method's exact limit at polynomial length T^θ is 1/(1+θ) (Step 8); θ < 1/2 is forced by the unconditional fourth moment of R_T, not by off-line zeros.

TeX line numbers: DPC = 2610.28527v2 Gaps_for_Dirichlet.tex; Lz = 2611.14486 SmallGapsPaper.tex.

## Notation and fixed parameters

L = log T, h = 2πφ/L. N(t) = #{ρ : 0 < γ ≤ t} (all zeros, multiplicity). Z_off = multiset of zeros with β > 1/2, γ ∈ [T/2, 5T/2]; d_ρ = β − 1/2 ∈ (0, 1/2).

0. Fix φ ∈ (2/3, 1) and θ ∈ (1/φ − 1, 1/2) (nonempty iff φ > 2/3). Then 0 < φθ < 1 < φ(1+θ) < 2. DPC choose 1/φ − 1 < θ < 1 (l.744), and their construction (l.744–772) uses only this sign condition (l.753: p₀(0) = sin(aθ)sin(a(1+θ))/(a∫_0^θ|g|²) < 0, a = πφ, g = e^{−iax} − e^{ia}e^{iax}) and θ < 1 (Lemma 6, l.686–741), so every θ ∈ (1/φ − 1, 1/2) is admissible. It gives a real polynomial q with ∫q²p₀ < 0 (l.755) and a real polynomial b₀ with Im b₀(x + i) = Re(q(x + i)²) for all real x, solved from the highest degree down since Im(x + i)^{m+1} has leading term (m + 1)x^m (l.757–759); so deg b₀ = 2 deg q + 1. (b₀ cannot be a constant: if deg q = d ≥ 1 with leading coefficient c ≠ 0, Re q(x + i)² has leading term c²x^{2d} and does not vanish identically.) Then P = q² + i b₀ satisfies Re P(x) = q(x)² ≥ 0 (q, b₀ real on ℝ), Re P(x + i) = Re q(x + i)² − Im b₀(x + i) = 0 and Re⟨P(A_φ)g, g⟩ < 0 on L²[0, θ] (l.760–765); then smooth f with f = 0 on [θ, ∞) and η > 0 (via DPC Lemma 6, l.686–741) with limiting quotient ≤ −3δ₀ < 0. R_T(t) = Σ_n c_n n^{−it}, c_n = μ(n) n^{−1/2} f(log n/L) ∏_{p|n} ρ_η(log p/L) (DPC l.676–684): n < T^θ, prime factors > T^{η/2}, at most k = ⌊2θ/η⌋ of them. R*(t) := Σ c̄_n n^{it}. G_T as DPC l.686–689. D := deg P. All of this is operator theory, Mertens and Laplace transforms: no RH. These are fixed before ε and T.

## Step 1. Convention: Im V counts all zeros (removes DPC's RH uses at l.164/194 and l.774–785)

Define log ξ(σ+iτ) by horizontal continuation from σ = +∞ (τ not an ordinate) and V_φ(t) = π⁻¹(log ξ(½ + i(t + h/2)) − log ξ(½ + i(t − h/2))). ξ(½+iτ) is real (DPC l.156). As τ crosses an ordinate γ: an on-line zero of multiplicity m adds mπ to arg ξ; an off-line ρ adds 2π (Re(½ − ρ) < 0, the principal arg of ½+iτ−ρ jumps from −π to π) and 1 − ρ̄ adds 0. So arg ξ(½+iτ) = π N(τ) + const and

  Im V_φ(t) = N(t + h/2) − N(t − h/2)  for every t with t ± h/2 not an ordinate.

DPC's Lemma 1 majorant (l.177–194) holds unconditionally: |s − ρ'| ≥ |t − γ'| and |arg| ≤ π per linear factor under this convention, so |log ξ(½+ε+it)| ≤ C T log²T + C Σ_{|γ'−t|≤2}(|log|t−γ'|| + π), giving boundary values and dominated convergence. In the slit domain of Step 3 these boundary values from Re s > 1/2 are the horizontal ones.

## Step 2. Contradiction hypothesis

Suppose 𝒩(T) ≤ εTL. Every off-line or multiple zero with γ ∈ [T, 2T] lies in a run of r ≥ 2 equal ordinates, which contains r − 1 zero gaps; so #{off-line or multiple zeros with γ ∈ [T, 2T]} ≤ 2εTL + O(L). (Lz l.203 uses the same observation.)

## Step 3. Exact moment identity: hairpins

For ρ ∈ Z_off, s ↦ log(s − ρ) has its cut on Im s = γ, Re s < β. In t (s = ½ + i(t ± h/2)) the part in Re s > 1/2 is the vertical slit S^± = {x^± − iy : 0 ≤ y ≤ d_ρ}, x^± = γ ∓ h/2, hanging down from the real axis to the branch point x^± − i d_ρ. Cuts are merged geometrically. All shifted cuts at the same real point x are one cut, including oppositely signed ones: if γ₂ − γ₁ = h, the S⁻ of ρ₁ (at γ₁ + h/2) and the S⁺ of ρ₂ (at γ₂ − h/2) coincide. Across the merged cut at x − iy the jump of V (right minus left) is J_x(y) = 2i(m₊(x, y) − m₋(x, y)), with m_±(x, y) = #{off-line zeros with ordinate x ± h/2 and depth > y} (one zero: ±2i). P is applied to this total signed jump, never to separate increments: for P(w) = w² and equal depths the true jump is 0, while separate increments would give [(w + 2i)² − w²] + [(w − 2i)² − w²] = −8. The cut at x has depth max{d : contributing zeros}. Around each branch point (a zero's tip, where m_± changes), circles of radius r carry O(r|log r|^D) → 0 and contribute nothing. Choose a = T/2 + O(1) and b = 5T/2 + O(1) at distance ≥ c/L from every shifted ordinate γ ± h/2. This is possible because N(t + 1) − N(t) ≪ log t (unconditional), so each unit interval has a gap ≫ 1/L between shifted ordinates. On the sides Re t = a, b, then, |V| ≤ T^{O(1)}. V is analytic in the rectangle [a, b] × [−1, 0] minus the merged cuts (as DPC l.335–341). Cauchy's theorem on this simply connected domain (orientation checked numerically below), for F = G_T R_T R* P(V):

  (H) ∫_a^b F(x)dx = ∫_{a−i}^{b−i} F dt + ∫_{b−i}^{b} F dt − ∫_{a−i}^{a} F dt − i Σ_{cut points x} ∫_0^{depth(x)} (F_R − F_L)(x − iy) dy,

F_R − F_L = G_T R_T R* [P(V_L + J_x(y)) − P(V_L)] (the total signed jump). The sides are O(T^{O(1)} e^{−cL²}) (DPC l.335–341, Step 1 majorant). On Im t = −1 (Re s = 3/2) V = G + π⁻¹Σβ(k)k^{−it} whatever the zeros are, and everything in DPC l.341–388 (Prop 1) and l.698–741 (Lemma 6) after the first contour shift uses only the Euler series, the gamma factor (Stirling) and exponentials. Hence, unconditionally,

  ∫_ℝ G_T|R_T|² P(V) dx = M_T + 𝓗_T + o(T),  lim_{η↓0} lim_T M_T / ∫G_T|R_T|² = ⟨P(A_φ)f, f⟩/‖f‖²,

with ∫G_T|R_T|² ~ (T/2)Σ|c_n|² (DPC l.700–702), and 𝓗_T = −iΣ_slits∫(F_R − F_L)dy, a sum over Z_off. Zeros of Z_off with γ ∉ [T, 2T] have G_T ≤ e^{−L²/32} on their slits, and |R_T R*| ≤ T^{2θ}, |V| ≤ T^{O(1)} up to integrable log singularities, so they contribute o(1).

Toy check (toy_identity.py, log: toy_identity.json): F(s) = (planted off-line pairs)·∏(s − ½ − iγ_k)·ζ(s+1), Möbius R of length Q ∈ {1, 12}, j = 1, 2, 3, one planted zero and a cluster of three (d = 0.25, 0.1, 0.4, two of them 0.35 apart): residual of (H) 1e−7 … 1.6e−4 (real-line quadrature noise at the log singularities), against 0.03 … 38 without the hairpins.

## Step 4. Weight of one off-line zero

On the slit, R_T R*(x − iy) = M(½ + y + ix)·conj(M(½ − y + ix)) with M(s) = Σ c_n n^{½−s}: one factor at ½ + y (mean square S_−(y) = Σ|c_n|²n^{−2y} ≤ S(0)), one at the mirror point ½ − y (mean square S_+(y) = Σ|c_n|²n^{2y} ≤ T^{2θy}S(0)). So a zero at β = ½ + d carries T^{θ(β−½)} in continuous mean (Cauchy–Schwarz over x). This is a continuous-mean heuristic only, not a bound at the zero-sampled points. The rigorous aggregate factor e^{θ(a+1)} for the depth class S_a comes from 6c–6e: Poisson representations, the fourth moment of R*_{Y_a} on the line Im t = −Y_a, and sparse kernels, with no independence assumption between zeros and R_T. The weight is half of the T^{θ(2β−1)} of Lamzouri's removal (Lz l.371–377, where it comes from the pair (ρ, 1 − ρ̄) in a pair-correlation sum). Toy check (toy_weight.py, log toy_weight.json; Q = 100, 1000, 5000, 4000 random x): E|R R*(x − iy)| = (0.86–1.03)·√(S₊S₋) = S(0)·sinh(a)/a·(0.86–1.03) with a = y log Q, while E|RR*|/S₊ falls like e^{−a}. Growth in j (or deg P): |P(w ± 2im) − P(w)| ≤ C_P m(1 + |w| + 2m)^{D−1}, polynomial in V_L; that this costs no further power of T needs the weighted slit moments of V_L, proved in 6b–6e.

## Step 5. Unconditional absolute moments (replace DPC Lemma 5, l.607–671)

5a. ∫_{T/2}^{5T/2}|R_T(x + iη')|⁴ dx ≪_η T e^{4θ|η'|L} for |η'| ≤ 1, in particular ∫_T^{2T}|R_T|⁴ ≪ T and ∫|R*(x − iY)|⁴ ≪ T e^{4θYL}: R_T² = Σ b_n n^{−it} has length T^{2θ} < T, so by Montgomery–Vaughan ∫_{T/2}^{5T/2}|R_T²|² = (2T + O(T^{2θ}))Σ|b_n|²n^{±2η'}; n is T^{η/2}-rough with ≤ 2k prime factors, exponents ≤ 2, so |b_n| ≤ 3^{2k}‖f‖²_∞ n^{−1/2} and Σ|b_n|² ≤ 3^{4k}‖f‖⁴_∞ ∏_{T^{η/2}<p≤T^θ}(1 + 1/p + 1/p²) ≪_η 1 (Mertens). This is the only place θ < 1/2 is needed besides Step 6; DPC's prime-sum route (l.643–671) is not used.
5b. ∫_{T/2}^{5T/2}(N(t+u) − N(t))^m dt ≪_m T for 0 ≤ u ≤ 100/L (Fujii, as quoted Lz l.242–249, all zeros), hence ∫|Im V|^m ≪ T. Longer windows (length ℓ ≤ 1) by Minkowski over ⌈ℓL⌉ pieces: ≪ T(1 + ℓL)^m.
5c. ∫_T^{2T}|Re V|^m ≪ T for every m ≥ 2 (m = 1 by Cauchy–Schwarz; M. Riesz needs 1 < m < ∞). Pairing the Hadamard factors as DPC l.633–641 (far zeros O(1), unconditional), and using log|½ + iτ − ρ'| = (P_{d'} * log|·|)(τ − γ') (P_d the Poisson kernel, d' = |β' − ½|), πRe V(t) = π H[F̃_T](t) + O(1), H the Hilbert transform, F̃_T(u) = Σ_{γ'∈[T/2,5T/2]} (P_{d'} * 1_{[−h/2,h/2]})(u − γ'): the sharp count for on-line zeros, a Poisson-smeared count for off-line ones. M. Riesz (DPC l.641): ‖Re V‖_m ≪ ‖F̃_T‖_m + T^{1/m}. On-line part: 5b. Off-line part: P_d * 1_{[−h/2,h/2]}(x) ≤ C(h/D)D²/(x² + D²), D = max(d, h); split Z_off dyadically, D ∈ [D_j, 2D_j), D_j = 2^j h (bin j = 0 holds all d < 2h; there use the total zero count, no density). With K_j = Poisson kernel at scale D_j (‖K_j‖₁ = 1) and n_j(x) = #{bin-j zeros within D_j of x}, the bin's contribution is F_j ≤ C(h/D_j)(K_j * n_j), so by Young ‖F_j‖_m ≤ C2^{−j}‖n_j‖_m (this handles the Poisson tails). ‖n_j‖_m^m ≤ (∫n_j)^{1/2}(∫n_j^{2m−1})^{1/2} ≪ (2^j e^{−c2^j}T)^{1/2}(T 2^{j(2m−1)})^{1/2}, using 5b and the uniform density bound in the form N(½ + δ, 3T) ≪ TL·e^{−cδL} (any slope κ > 0, e.g. Selberg's; zeros left of the line by the symmetry ρ ↦ 1 − ρ̄); hence Σ_j ‖F_j‖_m ≪ T^{1/m}. So ∫|V|^m ≪ T and ∫|P(V)|⁴ ≪ T.

## Step 6. The hairpins are small: |𝓗_T| ≤ Φ(ε)T + o(T), Φ(ε) → 0 as ε → 0

Input (density, uniform in σ ≥ ½): N(σ, T) ≪ T^{1 − κ(σ−½)} log T with κ = 1 − δ for every δ > 0 (Jutila 1983, as quoted Lz l.223–225), or κ = 8/7 − δ (Conrey, proved in Baluyot's thesis; audited, Riemann/Deep/Baluyot). Needed: θ < κ/2.

6a. S_a := {ρ ∈ Z_off : γ ∈ [T,2T], d_ρL ∈ [a, a+1)}, a = 0, 1, …; Y_a = (a+1)/L. By Step 2 and density, |S_a| ≤ σ_a TL with σ_a := min(2ε + o(1), C e^{−κa}). dμ_a := Σ_{ρ∈S_a} δ_ρ ⊗ 1_{(0,d_ρ)}(y)dy (each sign ± separately, at x_ρ^± = γ_ρ ∓ h/2). Every merged cut point x of Step 3 with depth in [a/L, (a+1)/L) is the x_ρ^± of its deepest contributing zero ρ ∈ S_a. On it |J_x(y)| ≤ 2Ñ_ρ (Ñ_ρ as in 6c), and |P(w + J) − P(w)| ≤ C_P|J|(1 + |w| + |J|)^{D−1}. |G_T(x − iy)| ≤ e^{1/(2H²)}G_T(x) ≤ 2. So, keeping Ñ_ρ inside the polynomial factor, 𝓗_{T,a} ≪ ∫ |R_T||R*|·Ñ_ρ(1 + |V_L| + 2Ñ_ρ)^{D−1} dμ_a ≤ ∫ |R_T||R*|(1 + |V_L| + 2Ñ_ρ)^{D} dμ_a. Cut points whose deepest zero lies outside [T, 2T] carry G_T ≤ e^{−L²/32} (Step 3).
6b. Representations (exact): R_T(x − iy) = (P_y * R_T)(x) (frequencies −log n ≤ 0); R*(x − iy) = (P_{Y−y} * R*_Y)(x), R*_Y(u) := R*(u − iY), for 0 ≤ y ≤ Y = Y_a (frequencies log n ≥ 0); V_L = V̂ + Δ_L, where V̂ is the V-function of ξ̂ = ξ·∏_{ρ∈Z_off}(s − 1 + ρ̄)/(s − ρ): zero-free in Re s > ½ at heights [T/2, 5T/2], so V̂ is analytic in −1 < Im t < 0 for Re t ∈ [T/2 + 1, 5T/2 − 1], which contains every box [x − T/8, x + T/8] × (−1, 0) with x ∈ [T − 1, 2T + 1]; |ξ̂| = |ξ| on the line, so Re V̂ = Re V there, and Im V̂ = on-line count + smeared off-line count; on Im t = −1, |V̂| ≤ |G| + O(1/L) + Σ_{ρ'} C h d'/(1 + (u−γ')²) = O(1). Poisson representation in the box B_x = [x − T/8 + u₋, x + T/8 + u₊] × (−1, 0), with 0 ≤ u_± ≤ 1 chosen so that both vertical sides are at distance ≥ c/L from every shifted ordinate γ' ± h/2 (possible as in Step 3). Along those sides |V̂| ≤ T^{O(1)} uniformly: no zero comes within c/L of a side's top end, and the Step 1 majorant and the bound on ℓ below apply. Their harmonic measure density seen from x − iy is ≪ e^{−πT/8}, so they contribute o(1). The bottom edge contributes O(1). The top edge has kernel ≤ C P_y(x − u) and mass ≤ 1, and lies inside I := [7T/8 − 1, 17T/8 + 1] ⊂ [T/2, 5T/2]. Hence |V̂(x − iy)| ≤ C (P_y * v)(x) + C with the truncated function v := |V̂|_ℝ|·1_I, and ∫_ℝ v^m = ∫_I v^m ≪ T (Step 5). The truncation matters: without it the norm of (1 + v)^M used in 6e is infinite. Δ = π⁻¹Σ_{ρ'∈Z_off}ℓ_{ρ'}, ℓ_{ρ'}(z) = log((z − t')/(z − t' − h)) − log((z − t̄')/(z − t̄' − h)), t' = γ' − h/2 − id' (the V-function of (s − ρ')/(s − 1 + ρ̄')): |ℓ_{ρ'}(z)| ≤ 8hd'/dist² if dist := |Re z − γ'| ≥ 2h, else ≤ C(1 + Σ_{w∈{t', t'+h, t̄', t̄'+h}} log⁺(h/|z − w|)).
6c. Sparse kernels. Jensen and Fubini: ∫ (P_y*|ψ|)^p(x_ρ) dμ_a ≤ ∫|ψ|^p K_a with K_a(u) = Σ_{ρ∈S_a} k_{d_ρ}(u − x_ρ), k_d(u) = ∫_0^d P_y(u)dy = (2π)⁻¹log(1 + d²/u²) (∫k_d^q = C_q d); for R* use K'_a = Σ k_{Y_a}(· − x_ρ). Claim: ‖K_a‖_q^q, ‖K'_a‖_q^q ≪_{q,w} (a+2)^q σ_a^{1−1/w} T for all 1 < q < ∞, 1 < w < ∞. Proof: near part ∫(Σ_{|u−x_ρ|≤Y}k_ρ)^q ≤ Σ_ρ N_ρ^{q−1}∫k_ρ^q ≤ C_q Y_a Σ_{ρ∈S_a} N_ρ^{q−1}, N_ρ := #{ρ'∈S_a : |γ' − γ| ≤ 2Y_a}; far part ≤ Y_a²/u² summed over dyadic shells, same bound; Hölder Σ_{S_a}N_ρ^{q−1} ≤ |S_a|^{1−1/w}(Σ_{all γ∈[T,2T]} Ñ_ρ^{(q−1)w})^{1/w}, Ñ_ρ = #{all zeros within 4Y_a + h of γ_ρ} (this also covers the zeros at ordinate γ_ρ ± h that share a merged cut with ρ), and Σ_γ Ñ_ρ^m ≪_m TL(a+2)^m by the argument of Lz Lemma DoubleSumZeros (l.251–283: Σ_j n_j^{m+1} ≪ L∫(N(t + 4π/L) − N(t))^{m+1}) with 5b. Tails of the Poisson kernels outside [T/2, 5T/2]: for R_T they cost O(T^{3θ/2}L) = o(T); for R*_Y the growth factor stays, ∫_outside |R*_Y|³K'_a ≪ σ_a L Y_a² T^{3θ/2} e^{3θ(a+1)}, which is absorbed in the bound for A₂ below.
6d. Δ moments: ∫_0^{d_ρ}log⁺(h/|x_ρ − iy − w|)^m dy ≤ C_m h for every w (one-dimensional). The far terms sum to 8hΣ_{ρ'}d'/dist² ≤ 32Ψ_ρ with Ψ_ρ := Σ_{ρ'∈Z_off}(1 + d_{ρ'}L)(1 + ((γ − γ')/h)²)⁻¹. So ∫_0^{d_ρ}|Δ_L|^m dy ≤ C_m (d_ρ + h)(1 + Ñ_ρ + Ψ_ρ)^m.
    Moments of Ψ. The depth bins cover all of Z_off: B_{a'} := {ρ' ∈ Z_off : d_{ρ'}L ∈ [a', a' + 1)}, γ' ∈ [T/2, 5T/2], so |B_{a'}| ≪ e^{−κa'}TL by density alone. Cut [T/2, 5T/2] into cells I_j of length h. Let n_j = #{zeros with γ ∈ I_j}, b_{a',j} = #{ρ' ∈ B_{a'} : γ' ∈ I_j}, J_j = (1 + j²)⁻¹. For γ ∈ I_j, Ψ_γ ≤ C Σ_{a'}(a' + 2)(J * b_{a'})_j. Hölder and Young (‖J * b‖_{2p} ≤ ‖J‖₁‖b‖_{2p}) give
    Σ_j n_j (J * b_{a'})_j^p ≤ (Σ_j n_j²)^{1/2} ‖J‖₁^p (Σ_j b_{a',j}^{2p})^{1/2} ≪_p (TL)^{1/2} (|B_{a'}|^{1/2}(Σ_j n_j^{4p−1})^{1/2})^{1/2} ≪_p TL e^{−κa'/4},
    using Σ_j n_j^m ≪_m TL (Lz l.251–283 with 5b) and Σ_j b^{2p} ≤ (Σ_j b)^{1/2}(Σ_j b^{4p−1})^{1/2}. Minkowski in a' with the weight a' + 2: (Σ_γ Ψ_γ^p)^{1/p} ≤ C Σ_{a'}(a' + 2)(Σ_j n_j (J * b_{a'})_j^p)^{1/p} ≪_p (TL)^{1/p}, so Σ_γ Ψ_γ^p ≪_p TL. Hölder restricted to S_a: Σ_{ρ∈S_a}Ψ_ρ^p ≤ |S_a|^{1−1/w}(Σ_γ Ψ_γ^{pw})^{1/w} ≪ σ_a^{1−1/w}TL. The same Hölder gives Σ_{ρ∈S_a}Ñ_ρ^p ≪ (a+2)^p σ_a^{1−1/w}TL (6c). Hence ∫|Δ_L|^m dμ_a + ∫Ñ_ρ^m dμ_a ≪ (a+2)^C σ_a^{1−1/w} T.
6e. Hölder with exponents (3, 3, 3):
  A₁ = ∫|R_T|³dμ_a ≤ ∫|R_T|³K_a ≤ ‖R_T‖₄³‖K_a‖₄ ≪ T (a+2) σ_a^{(1−1/w)/4};
  A₂ = ∫|R*|³dμ_a ≤ ‖R*_{Y_a}‖₄³‖K'_a‖₄ ≪ T e^{3θ(a+1)}(a+2) σ_a^{(1−1/w)/4};
  A₃ = ∫(1 + |V_L| + 2Ñ_ρ)^{3D}dμ_a ≤ C[μ_a(total) + ∫(P_y * v)^{3D}dμ_a + ∫|Δ_L|^{3D}dμ_a + ∫Ñ_ρ^{3D}dμ_a] ≪ (a+1)σ_aT + ‖v^{3D}‖_{L^r(I)}‖K_a‖_{r'} + (6d) ≪ T(a+2)^C σ_a^{(1−1/w)(1−1/r)}. Here μ_a(total) = Σ_{S_a}d_ρ ≤ Y_a|S_a|, and Jensen uses mass ≤ 1 of the box kernel;
  𝓗_{T,a} ≤ (A₁A₂A₃)^{1/3} ≪ T (a+2)^C e^{θ(a+1)} σ_a^{E}, E = (1−1/w)(1/2 − 1/(3r)) ↑ 1/2.
6f. Σ_a 𝓗_{T,a} ≪ T[(A+2)^C e^{θ(A+1)}(3ε)^E + Σ_{a>A}(a+2)^C e^{θ(a+1) − κEa}] =: Φ(ε)T. Order of choices. First κ: fix δ, κ = 1 − δ (Jutila). Then E: fix w, r so that κE > θ, possible since θ < 1/2 and E ↑ 1/2. Then A = c log(1/ε) with 0 < c < E/θ. The first term is ≪ (log(1/ε) + 2)^C e^{θ}ε^{E − cθ} → 0 because cθ < E. The tail is ≪ (A + 2)^C e^{−(κE − θ)A} → 0. So Φ(ε) → 0 as ε → 0; then fix ε; then T → ∞. All implied constants here and in Step 5 depend on φ, θ, f, P, η, δ, w, r only, not on ε.

## Step 7. Conclusion (DPC l.766–787 with Steps 1–6)

Re ∫G_T|R_T|²P(V) ≤ −3δ₀(T/2)Σ|c_n|²·(1 + o(1)) + Φ(ε)T + o(T) ≤ −δ₁T once ε ≤ ε₁ (δ₁ > 0 fixed). Outside E_T (windows with ≥ 2 zeros) Im V ∈ {0, 1} (Step 1), so Re P(V) ≥ 0; Hölder with 5a and 5c: δ₁T ≤ (∫|R_T|⁴)^{1/2}(∫|P(V)|⁴)^{1/4} meas(E_T)^{1/4}, so meas(E_T) ≥ c₂T with c₂ independent of ε. Each consecutive pair with gap < h covers centres of measure ≤ h: 𝒩(T) ≥ c₂T/h − O(L) = (c₂/2πφ)TL − O(L). With ε₀ = min(ε₁, c₂/(4πφ)), 𝒩(T) ≤ ε₀TL is impossible for large T. So 𝒩(T) > ε₀ T L. ∎

## Step 8. Exactness and what blocks θ > 1/2

- On L²[0, θ], Im e^{πA_φ} = sin a·I + a cos a|1⟩⟨1| + a sin a·S after F = e^{−iax}f (DPCAudit D7), S F(x) = i∫_0^θ sgn(x−y)F(y)dy with eigenvalues 2θ/((2k+1)π); I + aS > 0 for φθ < 1/2, and ⟨1, (I + aS)⁻¹1⟩ = tan(aθ)/a (Σ_{m odd}1/(b + mπ) = −½tan(b/2)). So Im e^{πA_φ} ⋡ 0 iff −cot a·tan(aθ) > 1 iff φ(1+θ) > 1. The method with length T^θ gives exactly λ > 1/(1+θ): DPC's g is optimal (threshold.py: sign change of λ_min at 1/(1+θ) for θ = 0.45, 0.5, 4/7, 1; the identity to 1e−6).
- θ < 1/2 enters only through fourth moments (5a and 6e). Off-line zeros alone would allow θ < κ/2: 4/7 with Conrey, i.e. λ > 7/11 (the heuristic of R448), 1 under the density hypothesis. Going below 2/3 needs an unconditional bound ∫_T^{2T}|R_T|^{2+δ} ≪ T (with the matching sparse Hölder) for rough Möbius polynomials of length T^θ, θ > 1/2; Montgomery–Vaughan gives only T^{2θ}Σ|b_n|²·(…) there, and DPC's route needs RH (Rodgers).

## Inputs

| Input | Source | Used in |
|---|---|---|
| DPC Prop 1 after the shift, Lemmas 2, 3, 6, proof of Thm B | 2610.28527 l.262–388, 198–229, 494–546, 686–741, 743–787 (audited correct, Riemann/Deep/DPCAudit) | Steps 0, 3, 7 |
| Hadamard product, Stirling, N(t+1) − N(t) ≪ log t | DPC l.170–180 | Steps 1, 3, 5c, 6b |
| Montgomery–Vaughan mean value theorem | standard | 5a |
| Fujii: ∫(N(t+u) − N(t))^{2k} ≪ T, u ≤ 100/L, all zeros | Lz l.242–249 (Fujii TAMS 1974, 1981) | 5b, 5c, 6c |
| Lz Lemma DoubleSumZeros argument | Lz l.251–283 | 6c, 6d |
| M. Riesz L^m bound for H | DPC l.641 | 5c |
| Density N(σ,T) ≪ T^{1−κ(σ−½)}log T uniformly, κ = 1 − δ (Jutila) or 8/7 − δ (Conrey/Baluyot) | Lz l.223–225 | 6a, 6f (5c needs only κ > 0) |
| Poisson/strip representations, Jensen | standard | 6b, 6c |

No step is left open. Steps 5c and 6b–6d are written at the level of detail above, not line by line. Self-review (three reviewers, 12 Nov 2026, one each on Steps 1–4, 5 and 6–7): GAP-repairable, no fatal defect. Repairs made: cuts merged with their total signed jump (Step 3; extended to oppositely signed coincident cuts in the orchestrator repairs below); the Step 4 weight is only heuristic and is made rigorous at the zero-sampled points by Step 6; Poisson tails in 5c via Young's inequality, with m ≥ 2 for M. Riesz; ξ̂ localised to heights [T/2, 5T/2], plus box harmonic measure; the R*_Y tail keeps its growth factor; q, w > 1; order of choices (δ, w, r) → ε → T.

## Repairs (orchestrator 33 review)

Three reviewers: A (Steps 0–3, 7) GAP, repairable; B (Step 5) CORRECT; C (Steps 4, 6) GAP, repairable. All three: SIGNIFICANT if the repairs close. No fatal error.

1. Step 0: b₀ is a real polynomial of degree 2 deg q + 1 with Im b₀(x + i) = Re(q(x + i)²) (DPC l.757–759). It is not a constant, since Re q(x + i)² has leading term c²x^{2d}. P = q² + i b₀ has the two-line signs of DPC l.760–765. DPC take θ < 1 (l.744), and the construction uses only 0 < φθ < 1 < φ(1+θ) < 2 (l.753).
2. Step 3: all coincident shifted cuts are merged, including oppositely signed ones (γ₂ − γ₁ = h). P is applied to the total jump J_x(y) = 2i(m₊ − m₋). Tip circles contribute 0. The sides a, b are placed ≥ c/L from all shifted ordinates.
3. Step 5a: the Montgomery–Vaughan main term on [T/2, 5T/2] is 2T.
4. Step 6b: box sides moved to points ≥ c/L from all shifted ordinates, so V̂ is uniformly T^{O(1)} on them, with side kernel ≪ e^{−πT/8}. v is truncated to I = [7T/8 − 1, 17T/8 + 1]; A₃ is rewritten accordingly.
5. Steps 6a and 6d: Ñ_ρ (now counting zeros within 4Y_a + h) is kept inside the polynomial factor, |P(w + J) − P(w)| ≤ C|J|(1 + |w| + |J|)^{D−1} with |J| ≤ 2Ñ_ρ. The depth bins B_{a'} cover all of Z_off. The Ψ moments are proved via h-cells, Hölder, Young (Σ_j n_j(J * b_{a'})_j^p ≪ TL e^{−κa'/4}) and Minkowski in a', followed by Hölder on S_a.
6. Step 6f: order of choices made explicit: κ, then E with κE > θ, then A = c log(1/ε) with 0 < c < E/θ, then ε, then T.
7. Theorem U ⇒ D(λ) > 0: take 2/3 < φ < min(λ, 1) and apply U on [T/2, T].
8. Step 4: marked as a continuous-mean heuristic. The rigorous aggregate factor e^{θ(a+1)} comes from 6c–6e, with no independence assumption.

RH relevance: N8 (zero spacing).
