# LiAudit: line-by-line proofs of the LiRate claims (referee: liaudit, 12 Nov 2026)

Source audited: Riemann/Deep/LiRate/proofs.md (cited as P), discoveries/li-remainder.md. Each step: statement, proof, verdict.
Notation (as P). ρ = β + iγ zeros of ξ, multiplicity counted. Under RH ρ = ½ + iγ; list γ > 0 with multiplicity.
θ(t) := 2 arctan(1/(2t)) ∈ (0, π) for t > 0; θ'(t) = −1/(t² + ¼), θ''(t) = 2t/(t² + ¼)²; 1/t − 1/(12t³) ≤ θ(t) ≤ 1/t.
g_n(t) := 2(1 − cos nθ(t)) = 4 sin²(nθ/2) ∈ [0, min(4, n²θ²)] ⊂ [0, n²/t²].
N(T) = #{0 < γ ≤ T}; θ_RS = Riemann–Siegel theta; A(t) := θ_RS(t)/π + 1; S(t) := N(t) − A(t) (= π⁻¹ arg ζ(½ + it) by
Riemann–von Mangoldt, valid for every T > 0 not an ordinate; Titchmarsh §9.3). t(n) := (n/2)(log n − κ), κ = 1 + log 2π − γ_E.
R(n) := λ_n − t(n); y_m := R(m)/m (Keiper's y_m: λ_m/m = ½(log m + γ_E − log 2π − 1) + y_m, identical).

## F0. Li's formula in the form used (RH)
Claim. Under RH, λ_n = Σ_{γ>0} g_n(γ) for all n ∈ Z, absolutely convergent; λ_{−n} = λ_n.
Proof. 1 − 1/ρ = (ρ − 1)/ρ = (−½ + iγ)/(½ + iγ): modulus 1, argument arg(−½ + iγ) − arg(½ + iγ) = (π/2 + arctan(1/(2γ))) −
(π/2 − arctan(1/(2γ))) = θ(γ) for γ > 0. For ρ̄: e^{−iθ(γ)}. Pairing ρ, ρ̄ (the order of Li's/Bombieri–Lagarias's sum):
[1 − e^{inθ}] + [1 − e^{−inθ}] = 2 − 2cos nθ = g_n(γ). |g_n(γ)| ≤ n²/γ², Σ γ⁻² < ∞: absolute. cos even ⇒ λ_{−n} = λ_n. ✓
Verdict: correct (P "Notation").

## F1. Lemma 0 (unconditional): if RH fails, |R(n)| ≥ e^{cn} for infinitely many n.
Proof. (i) Generating function. ξ is entire of order 1, ξ(1) = ½ ≠ 0, Hadamard: ξ(s) = ξ(0)Π_ρ(1 − s/ρ) (symmetric
pairing ρ ↔ 1 − ρ). For s = 1/(1 − z): 1 − s/ρ = w_ρ (1 − z/w_ρ)/(1 − z), w_ρ := 1 − 1/ρ. Differentiating,
d/dz log ξ(1/(1 − z)) = Σ_ρ [1/(1 − z) − 1/(w_ρ − z)] = Σ_{n≥1} z^{n−1} Σ_ρ [1 − w_ρ^{−n}]  (|z| < min|w_ρ|).
w_ρ^{−1} = ρ/(ρ − 1) = 1 − 1/(1 − ρ) = w_{1−ρ}, and ρ ↦ 1 − ρ permutes the zeros, so Σ_ρ[1 − w_ρ^{−n}] = λ_n:
Σ_{n≥1} λ_n z^{n−1} = d/dz log ξ(1/(1 − z)) (Bombieri–Lagarias 1999, eq. (1.4)-type; also P's form with poles z = ρ/(ρ−1),
the same set as {w_ρ} by the symmetry).
(ii) Radius. |w_ρ|² = ((1 − β)² + γ²)/(β² + γ²) < 1 ⟺ β > ½. If RH fails there is a zero with β ≠ ½, hence (ρ ↦ 1 − ρ) one
with β > ½. |w_ρ| → 1 as |γ| → ∞ uniformly in β ∈ [0,1], so r₀ := min |w_ρ| < 1
is attained. At z = w_ρ the function has a simple pole with residue −m_ρ ≠ 0, a genuine singularity; inside |z| < r₀ the map
z ↦ 1/(1 − z) avoids s = 1 and no zero of ξ is hit, so the radius of convergence of Σ λ_n z^{n−1} is exactly r₀.
(iii) Cauchy–Hadamard: limsup |λ_n|^{1/n} = 1/r₀ > 1. Pick 1 < a < 1/r₀: |λ_n| ≥ a^n infinitely often. |t(n)| ≤ n log n, so
|R(n)| ≥ a^n − n log n ≥ e^{cn} for those n ≥ n₀, any c < log a. ✓
Verdict: correct. (Minor: P writes the series as Σ λ_n z^n "of log ξ(1/(1−z))"; the exact identity is
Σ_{n≥1} λ_n z^n/n = log ξ(1/(1 − z)) − log ξ(1) (integrate (i)); same singularities, same radius.)

## F2. Lemma A, part 1 (unconditional): λ_n^{(0)} := ∫_0^∞ g_n A' dt = t(n) − 1/4 + O(log n/n)
(A1) A' = θ_RS'/π, θ_RS'(t) = ½ Re ψ(¼ + it/2) − ½ log π: analytic on [0,∞), θ_RS(0) = 0, θ_RS'(0) = ½ψ(¼) − ½log π ≈ −2.69.
Stirling (differentiable termwise in a sector): θ_RS(t) = (t/2) log(t/2πe) − π/8 + 1/(48t) + O(t⁻³),
θ_RS'(t) = ½ log(t/2π) − 1/(48t²) + O(t⁻⁴), θ_RS''(t) = 1/(2t) + 1/(24t³) + O(t⁻⁵).
Put r(t) := A'(t) − (2π)⁻¹ log(t/2π). Then r(t) = −1/(48πt²) + O(t⁻⁴), r'(t) = 1/(24πt³) + O(t⁻⁵) (t → ∞);
r(t) = −(2π)⁻¹ log t + O(1), r'(t) = −(2πt)⁻¹ + O(1) (t → 0+). So r ∈ L¹(0, ∞).
(A2) ∫_0^∞ r = lim_T [A(T) − A(0) − (T/2π)(log(T/2π) − 1)] = lim_T [θ_RS(T)/π − (T/2π) log(T/2πe)] = −1/8. ✓
(A3) g_n ≤ min(4, n²/t²) and |A'| ≪ 1 + log(t + 2): ∫ g_n A' converges absolutely; split A' = (2π)⁻¹log(t/2π) + r:
λ_n^{(0)} = I(n) + 2∫_0^∞ r − 2J(n),  I(n) := π⁻¹∫_0^∞ (1 − cos nθ(t)) log(t/2π) dt,  J(n) := ∫_0^∞ cos(nθ(t)) r(t) dt.
(A4) Substitute t = ½ cot(v/2), v ∈ (0, π): θ(t) = v, dt = −dv/(4 sin²(v/2)), t/2π = cot(v/2)/(4π), 1 − cos nv = 2 sin²(nv/2):
I(n) = (2π)⁻¹ ∫_0^π F_n(v) ℓ(v) dv,  F_n := sin²(nv/2)/sin²(v/2) = n + 2Σ_{1≤j<n} (n − j) cos jv,  ℓ(v) := log(cot(v/2)/(4π)).
(A5) −log(2 sin(v/2)) = Σ_{k≥1} cos(kv)/k and −log(2 cos(v/2)) = Σ_{k≥1} (−1)^k cos(kv)/k in L²(0, π); subtracting,
ℓ(v) = 2Σ_{k odd} cos(kv)/k − log 4π in L². F_n is a trigonometric polynomial, so pair termwise with
∫_0^π cos jv cos kv dv = (π/2)δ_jk, ∫_0^π cos kv dv = 0 (j, k ≥ 1):
∫_0^π F_n ℓ = −nπ log 4π + 2π Σ_{k odd < n} (n − k)/k,  i.e.  I(n) = Σ_{k odd<n} (n − k)/k − (n/2) log 4π (exact). ✓
(A6) n ∈ {2m, 2m+1}: the odd k < n are 1, 3, …, 2m − 1, Σ 1/k = H_{2m} − ½H_m = ½ log m + log 2 + γ_E/2 + O(m⁻²)
(H_x = log x + γ_E + 1/(2x) + O(x⁻²); the 1/(4m) terms cancel). n = 2m: n(…) − m = (n/2)(log n + log 2 + γ_E − 1) + O(1/n).
n = 2m + 1: ½ log m = ½ log n − ½ log 2 − 1/(2n) + O(n⁻²), n(…) − m = (n/2)(log n + log 2 + γ_E) − ½ − (n − 1)/2 + O(1/n), same.
Hence I(n) = (n/2)(log n + γ_E − 1 − log 2π) + O(1/n) = t(n) + O(1/n). ✓ (Numerically I(n) − t(n): see logs/lemmaA.txt.)
(A7) J(n) = ∫_0^π cos(nv) ρ(v) dv, ρ(v) := r(½cot(v/2))/(4 sin²(v/2)). Near 0: t⁻² = 4tan²(v/2), ρ = −1/(48π cos²(v/2)) + O(v²),
ρ' = O(v) (from r'). Near π: t ~ (π − v)/4, ρ = −(8π)⁻¹ log(π − v) + O(1), ρ' = (8π(π − v))⁻¹ + O(1). Split at π − 1/n:
∫_{π−1/n}^π |ρ| ≪ log n/n; on [0, π − 1/n] integrate by parts once: ≪ (|ρ(0)| + |ρ(π − 1/n)|)/n + n⁻¹∫_0^{π−1/n} |ρ'| ≪ log n/n.
So λ_n^{(0)} = t(n) + O(1/n) − 1/4 + O(log n/n). ✓  Verdict: correct (no RH used). P's sketch says "ρ bounded near 0 and log
singularity at π: O(log n/n)" — correct, the split at π − 1/n is the missing detail.

## F3. Lemma A, part 2 (RH): the exact decomposition
Claim. For every T > 0 not an ordinate (P's hypothesis γ₁ < T is not needed):
R(n) = 7/4 + 2S(T) cos nθ(T) − 2 Re Z_T(n) + 2∫_0^T cos(nθ) A' dt − ∫_T^∞ S g_n' dt + O(log n/n),  Z_T(n) := Σ_{γ≤T} e^{inθ_γ}.
Proof. N is a right-continuous step function, so by F0 λ_n = ∫_{(0,∞)} g_n dN (absolutely). A is C¹, so
λ_n^{(0)} = ∫_{(0,∞)} g_n dA, and λ_n − λ_n^{(0)} = ∫_{(0,∞)} g_n dS. On (0, T]: ∫ g_n dN = Σ_{γ≤T} 2(1 − cos nθ_γ) = 2N(T) − 2Re Z_T(n);
∫_0^T g_n A' = 2(A(T) − A(0)) − 2∫_0^T cos(nθ)A', A(0) = 1. Difference: 2S(T) + 2 − 2Re Z_T + 2∫_0^T cos(nθ)A'.
On (T, ∞): Stieltjes integration by parts with g_n ∈ C¹, S of locally bounded variation: ∫_{(T,M]} g_n dS = g_n(M)S(M) − g_n(T)S(T)
− ∫_T^M S g_n'; S(M) = O(log M) (unconditional), g_n(M) ≤ n²/M² → 0; |g_n'| = 2n|θ'||sin nθ| ≤ 2n²θ|θ'| ≤ 2n²/t³, so ∫_T^∞ S g_n'
converges absolutely. Total: 2S(T) − g_n(T)S(T) = 2S(T) cos nθ(T); 2 − 1/4 = 7/4 with F2. ✓
Size of 2∫_0^T cos(nθ)A': in v it is 2∫_{θ(T)}^π cos(nv) a(v) dv, a(v) := A'(½cot(v/2))/(4 sin²(v/2)), C^∞ on (0, π];
a(v) = (2πv²)⁻¹ log(1/(2πv))(1 + o(1)), v³a'(v) → −∞ as v → 0, so a' < 0 on (0, v₀] (fixed v₀). One integration by parts:
≤ (|a(θ(T))| + |a(π)|)/n + n⁻¹(a(θ(T)) − a(v₀) + ∫_{v₀}^π |a'|) ≪ (1 + T² log(T + 2))/n. ✓  Verdict: correct.

## F4. Theorem 1 (RH): Σ_{n≤N} R(n)² ≥ (c₁ − ε) N^{3/2}(log N)^{1/2} for N ≥ N₀(ε), c₁ = (3/5)^{1/2}(4/5)²/(2π√2) = 0.055790
Order of limits: given ε, fix ψ (hence φ) with (∫ψ)³ ≥ 1 − ε/(2c₁); every O below depends on ψ only; then N ≥ N₀(ψ, ε).
(1) Weight. ψ ∈ C_c^∞(−½, ½) even, ψ ≥ 0, ∫ψ² = 1; φ := ψ*ψ: even, supp ⊂ (−1, 1), φ ≥ 0, φ(0) = ∫ψ² = 1, φ ≤ 1 (C–S),
∫φ = (∫ψ)² ≤ |supp ψ|·∫ψ² < 1 and → 1 along flat-top ψ. φ̂(ξ) := ∫φ e^{iuξ} du = ψ̂(ξ)² ≥ 0 (ψ̂ real: ψ even), Schwartz.
Φ(x) := Σ_n φ(n/N)e^{inx} = N Σ_k φ̂(N(x + 2πk)) (Poisson; f(y) = φ(y/N)e^{iyx} ∈ C_c^∞) ⇒ Φ real, even, Φ ≥ 0,
Φ(0) = N∫φ + O(N^{−K}), and |Φ^{(j)}(x)| ≪_{j,K} N^{1+j}(1 + N‖x‖)^{−K}, ‖x‖ := dist(x, 2πZ).
(2) R̃(n) := λ_{|n|} − t(|n|) on Z, R̃(0) = 0 (λ₀ = 0, t(0) := 0). ⟨F, G⟩ := Σ_n φ(n/N)F(n)G(n) is positive semidefinite (φ ≥ 0).
φ(±N) = 0 ⇒ ‖R̃‖² = 2Σ_{1≤n<N} φ(n/N)R(n)² ≤ 2Σ_{n≤N} R(n)². ✓
(3) For θ ∈ R: Σ_n φ(n/N)e^{−inθ}λ_{|n|} = Σ_{γ'} B(θ, θ_{γ'}), B(θ, θ') := 2Φ(θ) − Φ(θ − θ') − Φ(θ + θ'): insert F0 (valid for
negative n), finite n-sum, absolutely convergent γ'-series. B(θ, θ') = −θ'²Φ''(ξ), some |ξ − θ| ≤ θ' (second difference). ✓
(4) Trend lemma. τ(θ) := Σ_n φ(n/N)e^{−inθ} t(|n|). Claim: τ(θ) = log(2πθ)/θ² + O_ψ(1 + N^{−2}θ^{−4} log(2/θ)), N^{−1/2} ≤ θ ≤ 1.
Proof. f(x) := φ(x/N) t(|x|) is continuous, compactly supported. (a) FT of u(x) := t(|x|) away from 0: |x|^s (Re s > −1) is an
analytic family of tempered distributions with FT −2Γ(s+1) sin(πs/2)|ξ|^{−s−1} on ξ ≠ 0; ∂_s at s = 1 gives FT(|x| log|x|) =
2(log|ξ| − 1 + γ_E)/ξ², FT|x| = −2/ξ², so û = T, T(ξ) = (log|ξ| − 1 + γ_E + κ)/ξ² = log(2π|ξ|)/ξ² on ξ ≠ 0. ✓
(b) f̂ = (2π)⁻¹(Nφ̂(N·)) * û, i.e. f̂(ξ) = ⟨û, k_ξ⟩, k_ξ(η) := (2π)⁻¹Nφ̂(N(ξ − η)), ∫k_ξ = φ(0) = 1, k_ξ even about ξ,
∫(η − ξ)²k_ξ = −φ''(0)N^{−2}. Let χ be smooth, symmetric about ξ, = 1 on |η − ξ| ≤ ξ/4, supported in |η − ξ| ≤ ξ/2, |χ^{(j)}| ≪ ξ^{−j}.
On supp(1 − χ)k_ξ every seminorm of (1 − χ)k_ξ is ≪ N^{O(1)}(Nξ)^{−K} ≤ N^{−10} (ξ ≥ N^{−1/2}, K large), and û has finite order:
⟨û, (1 − χ)k_ξ⟩ = O(N^{−10}). On supp χ, û = T smooth; Taylor at ξ: ∫Tχk_ξ = T(ξ)(1 + O(N^{−10})) + 0·T'(ξ) (symmetry)
+ O(N^{−2} sup_{|η−ξ|≤ξ/2}|T''|), |T''(η)| ≪ ξ^{−4} log(2/ξ). So f̂(ξ) = log(2πξ)/ξ² + O(N^{−2}ξ^{−4}log(2/ξ) + N^{−10}) for
N^{−1/2} ≤ ξ ≤ 1; the same argument gives f̂(ξ) ≪ log(2 + |ξ|)/ξ² for |ξ| ≥ 1. (c) Poisson: Σ_n f(n)e^{−inθ} = Σ_k f̂(θ + 2πk)
(f continuous, compact support, Σ_k |f̂(θ + 2πk)| < ∞: the 2π-periodization is continuous with Fourier coefficients f(n)).
Aliases k ≠ 0: |θ + 2πk| ≥ 5|k|, Σ_{k≠0} log(2 + 5|k|)/(25k²) = O(1). ✓  [P's step 4 states this with "O(log N)"; correct, and
at θ ≥ θ(T) ≍ (log N/N)^{1/2} the error is O(1).]
(5) Dual vector W(n) := Σ_{γ≤T} cos nθ_γ, T ≤ N^{1/2}/2 (so θ_γ ≥ θ(T) ≥ 1/T − 1/(12T³) ≥ 1/(1.01T) ≥ N^{−1/2} for γ ≤ T, T ≥ 3).
φR̃ even ⇒ ⟨R̃, W⟩ = Σ_{γ≤T} [Σ_{γ'} B(θ_γ, θ_{γ'}) − τ(θ_γ)]. Bounds, X := N³:
γ' > X: θ_{γ'} < N^{−3}, |B| ≤ θ_{γ'}² N³(1 + N(θ_γ − N^{−3}))^{−K} ≪ θ_{γ'}²N^{3−K/2}, Σ_{γ'>X} θ_{γ'}² ≪ log X/X: total ≪ N^{−10}.
γ' ≤ X: 2Φ(θ_γ) − Φ(θ_γ + θ_{γ'}) with θ_γ, θ_γ + θ_{γ'} ∈ [N^{−1/2}, 2θ(γ₁)] ⊂ [N^{−1/2}, 0.15]: each ≪ N^{1−K/2}; N(X) ≪ N³ log N
terms: total ≪ N^{−10}. Remaining −Σ_{γ'≤X} Φ(θ_γ − θ_{γ'}) ≤ −Σ_{γ'≤T} Φ(θ_γ − θ_{γ'}) since Φ ≥ 0 (the only sign step). So
⟨R̃, W⟩ ≤ −Σ_{γ,γ'≤T} Φ(θ_γ − θ_{γ'}) − Σ_{γ≤T} τ(θ_γ) + O(N^{−9}).  ‖W‖² = ½Σ_{γ,γ'≤T}[Φ(θ_γ − θ_{γ'}) + Φ(θ_γ + θ_{γ'})]
(cos a cos b = ½[cos(a − b) + cos(a + b)]), the second part O(N^{−9}). With (4): ⟨R̃, W⟩ ≤ −2‖W‖² + τ',
τ' := τ_T + C_ψ N(T) log N + 1, τ_T := −Σ_{γ≤T} log(2πθ_γ)/θ_γ². ✓
(6) x := ‖W‖² ≥ D := ½Φ(0)N(T) (all terms ≥ 0; keep the N(T) diagonal pairs of the multiplicity list). If 2D > τ' (> 0):
⟨R̃, W⟩ ≤ −(2x − τ') < 0, so by Cauchy–Schwarz ‖R̃‖² ≥ (2x − τ')²/x; d/dx[(2x − τ')²/x] = (2x − τ')(2x + τ')/x² > 0 for
x > τ'/2, so ‖R̃‖² ≥ (2D − τ')²/D and Σ_{n≤N} R(n)² ≥ Φ(0)N(T)(1 − τ'/(Φ(0)N(T)))². ✓
(7) Constants. θ_γ = 1/γ − 1/(12γ³) + O(γ⁻⁵) ⇒ −log(2πθ_γ)/θ_γ² = γ² log(γ/2π) + O(log γ). Stieltjes with N = A + S, S = O(log t):
Σ_{γ≤T} γ² log(γ/2π) = (2π)⁻¹∫_0^T t² log²(t/2π)dt + O(T² log² T) = (T³/6π)(L² − 2L/3 + 2/9) + O(T² log² T), L := log(T/2π).
So τ_T = (T³L²/6π)(1 + O(1/L)); N(T) = (TL/2π)(1 + O(1/L)); C_ψN(T)log N/τ_T ≪ log N/(T²L) ≪ (log N)/N.
Take T with T²L = uN∫φ (unique for large N; T ≍ (N/log N)^{1/2} ≤ N^{1/2}/2): τ'/(Φ(0)N(T)) = u/3 + O(1/L), and the bound is
(N∫φ)^{3/2} L^{1/2} √u (1 − u/3)²(2π)⁻¹(1 + O(1/L)). max_u √u(1 − u/3)² at u = 3/5 (derivative ∝ (1 − u/3)(½ − 5u/6)) = 0.495740;
2D > τ' ⟺ u/3 < 1 ✓. L = ½ log N (1 + O(log log N/log N)) ⇒ Σ_{n≤N}R² ≥ (∫ψ)³ c₁ N^{3/2}(log N)^{1/2}(1 + O(log log N/log N)). ∎
Verdict: correct. Data at finite N (logs/thm1_chain.txt, real zeros and λ_n): every inequality of (5)–(6) holds at N = 10³–10⁵,
including the corner T = N^{1/2}/2; the trend error is 13.5 (ψ flat-top, edge 0.1) and 1384 (edge 0.02) at N = 10⁵: O_ψ, as stated.


## F5. Corollaries 1, 2 and the unconditional Ω
Cor 1 (RH): limsup |R(n)|/(n log n)^{1/4} ≥ (3c₁/2)^{1/2} = 0.28928. Proof: else there are C < (3c₁/2)^{1/2}, n₀ with
|R(n)| ≤ C(n log n)^{1/4} (n ≥ n₀); Σ_{n₀≤n≤N}(n log n)^{1/2} = (2/3)N^{3/2}(log N)^{1/2}(1 + o(1)) gives Σ_{n≤N}R² ≤ O(1) +
(2C²/3 + o(1))N^{3/2}(log N)^{1/2} < (c₁ − ε)N^{3/2}(log N)^{1/2} for ε = (c₁ − 2C²/3)/2 and N large: contradicts F4. ✓
Ω unconditional: RH ⇒ Cor 1; ¬RH ⇒ F1 (e^{cn} ≫ (n log n)^{1/4}). ✓
Cor 2 (RH): Σ_{m≤N} m^{1/2}y_m² = Σ_{m≤N} m^{−3/2}R(m)² ≥ N^{−3/2}Σ_{m≤N}R(m)² ≥ (c₁ − ε)(log N)^{1/2} → ∞. ✓
Verdict: correct.

## F6. Theorem 2 (RH): Σ_{N<n≤2N} R² ≪ N^{3/2}(log N)^{3/2}; with Montgomery's theorem ≪ N^{3/2}(log N)^{1/2} log log N
Proof. Fix φ ∈ C_c^∞((½, 4)), φ ≥ 0, φ ≥ 1 on [1, 2]; Φ_+(x) := Σ_n φ(n/N)e^{inx}, |Φ_+^{(j)}| ≪ N^{1+j}(1 + N‖x‖)^{−K}.
Take T with T² log T = N. For n ∈ (N/2, 4N), F3 gives R(n) = −2Re Z_T(n) − E(n) + O(log N), E(n) := ∫_T^∞ S g_n' (the terms 7/4,
2S(T)cos nθ(T) = O(log T), 2∫_0^T cos(nθ)A' ≪ T² log T/n ≪ 1, O(log n/n)). (a + b + c)² ≤ 3(a² + b² + c²):
Σ_{N<n≤2N} R² ≤ 3[4Σ_n φ(n/N)|Z_T(n)|² + Σ_n φ(n/N)E(n)² + O(N log² N)].
(a) Σφ|Z_T|² = Σ_{γ,γ'≤T} Φ_+(θ_γ − θ_{γ'}); mean value theorem |θ_γ − θ_{γ'}| ≥ |γ − γ'|/(T² + ¼), N/(T² + ¼) ≥ ½ log T:
≪ N Σ_{γ,γ'≤T} (1 + ½ log T|γ − γ'|)^{−K}. Unconditional count N(t + 1) − N(t) ≪ log t: inner sum ≪ log T ⇒ (a) ≪ N N(T) log T.
Under RH, Montgomery (1973; F(α) = |α| + T^{−2|α|}log T(1 + o(1)) + o(1) uniformly on |α| ≤ 1, with multiplicity) and a
majorant r(x) = C[sinc²(x/2) + sinc²((x + 1)/2)] ≥ (1 + |x|)^{−2} (zeros disjoint; supp r̂ ⊂ [−½, ½] for sinc x = sin x/x or
sin πx/(πx)): Σ_{γ,γ'≤T} r((γ − γ')log T/2π)w(γ − γ') = (T log T/2π)∫F r̂ ≪ T log T, w(u) = 4/(4 + u²). Pairs with |γ − γ'| ≤ 1
have w ≥ 4/5; pairs with |γ − γ'| > 1 contribute ≪ N(T)(log T)^{1−K} to our sum. So (a) ≪ N N(T). [P omits the w-factor
detail; harmless.] ✓
(b) E(n)² = ∫∫_{[T,∞)²} S(t)S(t')·4n²θ'(t)θ'(t') sin nθ(t) sin nθ(t'); summing, Σφ E² = ∫∫ S(t)S(t') K, K(t,t') = 2θ'(t)θ'(t') Re Σ_n φ(n/N)n²
[e^{in(θ(t)−θ(t'))} − e^{in(θ(t)+θ(t'))}]; |Σφ n²e^{inx}| ≪ N³(1 + N‖x‖)^{−K}, |θ'(t)| ≤ t^{−2}. In v' = θ(t'), dv' = |θ'(t')|dt':
∫_T^∞ |K(t,t')|dt' ≪ N³t^{−2}∫[(1 + N|θ(t) − v'|)^{−K} + (1 + Nv')^{−K}]dv' ≪ N²/t². K symmetric, |SS'| ≤ ½(S² + S'²) (Schur):
Σφ E² ≤ ∫_T^∞ S(t)²(∫|K(t,t')|dt')dt ≪ N²∫_T^∞ S² t^{−2} dt ≪ N² log log T/T (Selberg ∫_0^U S² ~ (U/2π²) log log U, unconditional;
partial integration). ✓
(c) T ≍ (N/log N)^{1/2}, N(T) ≍ (N log N)^{1/2}: (a) ≪ N^{3/2}(log N)^{3/2} [RH + Montgomery: N^{3/2}(log N)^{1/2}],
(b) ≪ N^{3/2}(log N)^{1/2} log log N, N log² N smaller. ∎  Verdict: correct (Montgomery's theorem needs RH only).
Data: Σ_{N<n≤2N}R²/(N^{3/2}(log N)^{1/2}) = 0.17–0.27 for N = 10²–5·10⁴ (logs/audit_num.txt B5).

## F7. Corollary 3 and the criterion
Claim (RH): Σ_m m^δ y_m² < ∞ ⟺ δ < 1/2. ⇐: dyadic blocks, F6: Σ_{2^j<m≤2^{j+1}} m^{δ−2}R(m)² ≪ 2^{j(δ−1/2)} j^{3/2}, summable
for δ < 1/2. ⇒ (contrapositive): δ ≥ 1/2 ⇒ m^δ ≥ m^{1/2}, divergent by F5 Cor 2. ✓
RH ⟺ (∃ or ∀) δ ∈ [0, 1/2): Σ m^δ y_m² < ∞. RH ⇒ ∀ by the above. ¬RH ⇒ |y_n| = |R(n)|/n ≥ e^{cn}/n infinitely often (F1),
so the series diverges for every real δ. ✓ δ = 1/2 diverges always (RH: Cor 2; ¬RH: F1). δ = 0 is Arias de Reyna's theorem —
verified at source (abstract, OpenAlex record of Funct. Approx. 45 (2011) 7–21: "RH is equivalent to (y_m) ∈ ℓ₂", same y_m);
his theorem numbering (Thm 1/Cor 2/Thm 6) not checked (Project Euclid PDF behind a bot wall). Exponent: Thm 1 + dyadic Thm 2 ⇒
log(Σ_{n≤N}R²)/log N → 3/2. ✓  Verdict: correct.


## F8. Theorem 3 (RH): |R(n)| ≤ (κ₃ + o(1))√n log n/(log log n)^{1/2}, κ₃ = (2/(3π))(π²/12)^{1/4} = 0.20212
Input (checked in TeX, arXiv 1309.1526 l.171, 187): S₁(t) := ∫_0^t S, −(π/24 + o(1))h ≤ S₁(t) ≤ (π/48 + o(1))h, h := log t/(log log t)²
(Carneiro–Chandee–Milinovich, RH). Proof. F3 with T = T(n) → ∞, T ≍ (n/log log n)^{1/2}: |7/4 + 2S(T)cos nθ(T)| ≪ log n,
|2∫_0^T cos(nθ)A'| ≪ T² log T/n ≪ log n/log log n, |2 Re Z_T| ≤ 2N(T) = (T/π) log T(1 + O(1/log T)).
∫_T^∞ S g_n' = −S₁(T)g_n'(T) − ∫_T^∞ S₁g_n'' (S₁g_n' → 0); |g_n'(T)| ≤ 2n/T² ⇒ first ≪ n h(T)/T² ≪ log n·log log n.
g_n'' = 2n[θ'' sin nθ + nθ'² cos nθ] ⇒ |g_n''| ≤ 4n/t³ + 2n²/t⁴; |S₁| ≤ (π/24 + o(1))h; ∫_T^∞ h t^{−4} = h(T)/(3T³)(1 + o(1)):
|∫S₁g_n''| ≤ (C₁ + o(1))(2n²/3)h(T)/T³ + O(n h(T)/T²), C₁ = π/24. f(T) := (T/π) log T + (2C₁/3)n² log T/(T³(log log T)²), log T
treated as ½ log n(1 + o(1)): minimum at T⁴ = 2πC₁n²/(log log n)², value (4/3)(T/π) log T = (2/(3π))(2πC₁)^{1/4}√n log n/√(log log n). ∎
2πC₁ = π²/12 ⇒ κ₃ = 0.2021. ✓  The parenthetical refinements check out: (i) |S| ≤ (¼ + o(1))log t/log log t alone gives
(2π)^{−1/2} = 0.399 (optimize (T/π)log T + n log T/(2T log log T)); (ii) the reviewer's 0.168: centring S₁ at −(π/96)h leaves
|S₁ − c| ≤ (π/32)h, the slowly varying centre contributes O(n h/T²) by one more partial integration, and |cos nθ(t)| averages to
2/π against the smooth weight h t^{−4} (phase n/t varies on scale t²/n ≪ t), so C₁_eff = 1/16 — correct but not claimed. Data:
max_{100≤n≤10⁵} |R|/(κ₃√n log n/√log log n) = 0.372 (n = 217; logs/audit_num.txt B3). Verdict: correct.

## F9. Theorem 4 (RH): Σ_n φ(n/N)R(n) = (7/4)Σ_n φ(n/N) + O_φ(log N), φ ∈ C_c^∞((0, ∞)), φ ≥ 0, φ ≢ 0
Proof. Φ_+(x) := Σφ(n/N)e^{inx}: Φ_+(0) = N∫φ + O(N^{−K}), |Φ_+^{(j)}(x)| ≪ N^{j+1}(1 + N‖x‖)^{−K}. (i) F2:
Σφ(λ_n^{(0)} − t(n)) = −Φ_+(0)/4 + O(Σ_{n≍N} log n/n) = −Φ_+(0)/4 + O(log N). (ii) G_N(t) := Σφ(n/N)g_n(t) = 2Φ_+(0) − 2Re Φ_+(θ(t)).
Σφ(λ_n − λ_n^{(0)}) = ∫_{(0,∞)} G_N dS (finite sum of F3's Stieltjes integrals). Parts on (ε, M]: G_N(M)S(M) → 0 (G_N ≪ N³/M²),
S(ε) = −A(ε) → −1 (S(0+) = −1: N = 0 below γ₁, A(0) = 1), G_N(ε) → G_N(0+) = 2Φ_+(0) − 2Re Φ_+(π) = 2Φ_+(0) + O(N^{−K}) (θ(0+) = π).
So ∫ G_N dS = G_N(0+) − ∫_0^∞ S G_N' = 2Φ_+(0) + ∫_0^∞ S₁G_N'' + O(N^{−K}) (S₁(0) = 0, G_N' bounded at 0, S₁G_N' → 0 at ∞).
(iii) G_N'' = −2Re[Φ_+''(θ)θ'² + Φ_+'(θ)θ'']. t ≤ N^{1/2}: θ ≥ N^{−1/2}·0.99, θ', θ'' bounded ⇒ ≪ N^{3−K/2}. N^{1/2} ≤ t ≤ N:
≪ N³t^{−4}(t/N)^K + N²t^{−3}(t/N)^K ≪ (N³/t⁴)(t/N)^K. t ≥ N: Re Φ_+'(θ) = −Σφ n sin nθ ≪ N³θ, |Φ_+''| ≪ N³ ⇒ ≪ N³/t⁴.
S₁(t) = O(log(t + 2)) (Titchmarsh 9.9, unconditional): ∫|S₁G_N''| ≪ 1 + ∫_{N^{1/2}}^N log t (N³/t⁴)(t/N)^K dt + ∫_N^∞ N³ log t/t⁴
≪ log N. Total: (2 − ¼)Φ_+(0) + O(log N). ∎  Verdict: correct; the O-constant depends on φ (through the seminorms of φ in the
Poisson bounds). Data, the theorem's own quantity [Σφ(n/N)(R(n) − 7/4)]/log N (logs/thm4_scan.txt), N = 2·10³ … 10⁵: bump (0.5,1)
in [−1.3, 1.4] (−0.10 at 10⁵); (0.9,1) in [−20, 26]; (0.98,1) in [−32, 65]; (0.01,0.02) in [−0.9, 1.4]: bounded, no growth. Earlier
normalisation (logs/audit_num.txt B4): wide windows (0.5,1), (0.1,1), x⁵·bump(0.3,1): (dev)·N/log N = −0.89, +0.05, +0.13 at
N = 10⁵; narrow windows (0.9,1), (0.98,1), (0.001,0.002), (0.01,0.011) are not yet asymptotic at N ≤ 10⁵ (dev·N/log N up to 7·10³,
decreasing where two N values are comparable) — consistent with O_φ, not a counterexample. Sharp cutoff [1, N]: mean R − 7/4 =
−0.087, −0.042, −0.061 (no smoothing, no claim).


## F10. Corollary 4 (Merlini–Sala–Sala): |λ̃_n| ≤ a log n fails for every a
λ̃_n := λ_n − λ̄_n, λ̄_n := n[z^n] log(π^{−s/2}Γ(1 + s/2))|_{s=1/(1−z)} (Maślanka's trend). Input: Voros math/0506326 eq. (TS)
(l.495–498, checked in TeX): λ̄_n ~ t(n) + 3/4 − Σ_k B_{2k}n^{1−2k}/(4k) to all orders, unconditionally; so λ̃_n = R(n) − 3/4 + O(1/n).
F5 (Ω unconditional) ⇒ λ̃_n = Ω((n log n)^{1/4}) ⇒ no a works, and O(n^ε) fails for ε < 1/4. ✓
Identity of MSS's λ_tiny with our λ̃ (falsification of a definitional mismatch): our λ̃_n equals their table at n = 5, 10, 20, 100,
500, 1000, 2000, 3000 to ≤ 2·10⁻⁶ (logs/mss_check.txt). Numerics reproduced: first n ≥ 2 with |λ̃_n| > 2 log n is n = 5083
(λ̃ = 17.158913 vs 17.067314), 7584 such n ≤ 10⁵, max |λ̃_n|/log n = 4.1814 at n = 79302, a = 3 fails from 46992, a = 4 from 79286.
Correction: P's "(Their evidence: n ≤ 30.)" is inaccurate — MSS computed n ≤ 31, quoted Maślanka's values at n ≤ 3000 and read
values off plots up to n = 8000 (their l.-table: 5080 → 16.830 vs 2 log n = 17.066). The refutation is unaffected.
(Trend residual beyond all orders: λ̄₁₀ − t − 3/4 + 1/240 = 8.4·10⁻⁵ ≫ 1/(240·10³); Voros's "mod o(n^{−N})" allows this.)
Verdict: correct (one inaccurate parenthetical).

## F11. Inputs checked at source
| Input | Where used | Check |
| --- | --- | --- |
| Li/Bombieri–Lagarias generating function, λ_{−n} = λ_n | F0, F1 | derived in F0/F1 from Hadamard + functional equation |
| Riemann–von Mangoldt N = θ_RS/π + 1 + S, S(0+) = −1 | F3, F9 | Titchmarsh §9.3; S(0+) = −1 forced by N = 0, θ_RS(0) = 0 |
| S = O(log t), S₁ = O(log t) unconditional | F3, F9 | Titchmarsh Thms 9.4, 9.9 |
| Selberg ∫_0^U S² ~ (U/2π²) log log U | F6(b) | unconditional (Selberg 1946); only ≪ U log log U used |
| Montgomery pair correlation (RH), uniform on \|α\| ≤ 1 | F6(a) PC version | statement as restated in arXiv 2108.09258 l.224–232 (Goldston–Montgomery Lemma 8); only supp r̂ ⊂ [−½, ½] needed |
| CCM: −π/24 ≤ S₁/h ≤ π/48 (+o(1)), \|S\| ≤ (¼ + o(1)) log t/log log t | F8 | arXiv 1309.1526 TeX l.171, 187, 139 |
| Lagarias: O(√n log n) under GRH | "Before" | math/0404394 TeX l.177–182 |
| Voros: RH ⟺ λ_n = t(n) + o(n); 2Z(0) = +7/4 heuristic "should not count asymptotically" | "Before", F9 | math/0506326 l.405–416, 444–451; math/0404213 l.275–283 |
| Voros (TS): λ̄_n − t(n) → 3/4 | F10 | math/0506326 l.495–498 |
| Arias de Reyna: RH ⟺ (y_m) ∈ ℓ² | F7 | abstract (OpenAlex); full text not reached |
| Merlini–Sala–Sala conjecture wording, λ_tiny = λ̃ | F10 | 1903.06233 PDF text + numeric match |
| Johansson λ_100000 | numerics | 1309.2877 TeX l.843: 4.62580782406902231409416038…·10⁵ = our FFT rerun |
Coffey is not used in any proof. Checked for contradictions: math-ph/0505052 eq. (53) and math-ph/0507042 l.468 only conjecture
|S₂(n)| = O(n^{1/2+ε}) for the arithmetic part — compatible with Thm 1 (Ω) and Thm 3 (O).


## F12. Numerical falsification (all runs under .scratch/memcap.sh 290 MB; peaks in logs)
- λ_n, n ≤ 10⁵, two independent ways: (i) rerun of LiRate fft_lambda.py (r = 0.999, 2¹⁷ points, 360 bits; peak 97 MB, 32 s;
  logs/fft_r999.log), (ii) new zeros_lambda.py: Odlyzko's first 10⁵ zeros + exact smooth tail ∫_T^∞ g_nA' − g_n(T)S(T) (RH above T), remainder
  ≤ max|S₁|(4n/T² + 2n²/(3T³)) with Trudgian's |S₁(t)| ≤ 0.059 log t + 2.067 taken at t = 10T (peak 31 MB; logs/zeros_lambda.log). Max |difference| = 4.9·10⁻⁶ (≤ 2% of the remainder bound);
  ≤ 3.5·10⁻⁹ for n ≤ 1000. λ_100000 agrees with Johansson. So the data behind every LiRate number is confirmed independently.
- Lemma A (logs/lemmaA.txt): n(I(n) − t(n)) → 1/12 (n even), −1/6 (n odd); θ_RS(T)/π − (T/2π)log(T/2πe) + 1/8 = 6.6·10⁻⁶, 6.6·10⁻⁹ at
  T = 10³, 10⁶; λ_n^{(0)} − t(n) + 1/4 + 1/(24n) = 3.7·10⁻³, 8.4·10⁻⁵, −1.3·10⁻⁶, 3.8·10⁻⁸ at n = 5, 10, 20, 50.
- Theorem 1 chain at finite N (logs/thm1_chain.txt): with real zeros and λ_n, every step holds at N = 10³, 10⁴, 10⁵, both flat-top
  ψ (edges 0.1, 0.02) and at the corner T = N^{1/2}/2: ‖W‖² − ½Σ_{γ,γ'≤T}Φ(θ_γ − θ_{γ'}) = ½Σ Φ(θ_γ + θ_{γ'}) (identity, ~10⁻⁸);
  Σ_γΣ_{γ'}B + Σ_{γ,γ'≤T}Φ(θ_γ − θ_{γ'}) = −1.4·10³ < 0 (sign step); ⟨R̃, W⟩ = −3.2133·10⁶ ≤ −(2‖W‖² − τ') = −3.2119·10⁶;
  Σ_{n≤N}R² / Theorem-1 bound = 6.1 (N = 10⁵). The trend lemma fails at N = 10³ (τ'/τ_T = 1.99: θ_γ near N^{−1/2}) and holds
  from N = 10⁴ (τ'/τ_T = 1.004, 1.00015 at 10⁵) — an asymptotic statement, as claimed.
- Constants (logs/audit_num.txt): Σ_{n≤N}R²/(c₁N^{3/2}√log N) = 26.0, 15.2, 5.2, 2.40, 1.99, 2.12, 1.70, 1.96, 2.22, 2.46 at
  N = 3, 10, 30, 10², 3·10², 10³, 3·10³, 10⁴, 3·10⁴, 10⁵; over every 100 ≤ N ≤ 10⁵ min 1.472 (N = 660), so P's data line
  "≥ 1.7 c₁ … N ≤ 10⁵" is false as stated (logs/minratio.txt), harmless (> 1); max |R|/(n log n)^{1/4} = 1.558
  (n = 79302), > 0.289 for 53% of n ∈ [5·10⁴, 10⁵]; κ₃-ratio ≤ 0.372; Thm 2 window ratio 0.17–0.27; Thm 4 above (F9).
- Discovery-file slip: "smoothed mean of R on (N/2, N) is 1.7515 (N = 10⁴)" — the value at N = 10⁴ is 1.7443 (LiRate's own
  logs/tiny.txt and ours); 1.7512 is N = 3·10⁴. No claim depends on it.

