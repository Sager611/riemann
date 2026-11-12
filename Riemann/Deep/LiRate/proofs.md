# LiRate: proofs (Theorems 1–3)

Notation. ρ = 1/2 + iγ (RH), γ > 0 listed with multiplicity; N(T) = #{0 < γ ≤ T} = A(T) + S(T),
A(t) = θ_RS(t)/π + 1, S(t) = π^{-1} arg ζ(1/2 + it) (S(0+) = −1). θ(t) := 2 arctan(1/(2t)) ∈ (0, π), so 1 − 1/ρ = e^{iθ(γ)};
θ_γ := θ(γ); θ'(t) = −1/(t² + 1/4); θ(t) = 1/t + O(t^{-3}).
g_n(t) := 2(1 − cos nθ(t)).  Li: λ_n = Σ_ρ [1 − (1 − 1/ρ)^n] = Σ_γ g_n(γ) for all n ∈ Z (absolutely: g_n(γ) ≤ n²/γ²).
t(n) := (n/2)(log n − κ), κ = 1 + log 2π − γ_E;  R(n) := λ_n − t(n);  y_m := R(m)/m.

## Lemma 0 (unconditional). If RH fails, |R(n)| ≥ e^{cn} for infinitely many n (some c > 0).
Σ_n λ_n z^{n} has the singularities z = ρ/(ρ − 1) of log ξ(1/(1−z)); for Re ρ < 1/2 these lie in |z| < 1, so the radius of
convergence is < 1 (Cauchy–Hadamard); t(n) is polynomial.

## Lemma A (RH; exact decomposition). λ_n^{(0)} := ∫_0^∞ g_n A' dt = t(n) − 1/4 + O(log n/n)
(numerically λ_n^{(0)} − t(n) = −1/4 − 1/(24n) + O(n^{−2}): residual 3.8e−8 at n = 50, 5.2e−10 at n = 200, logs/analyze.txt). For γ_1 < T:
R(n) = 7/4 + 2S(T) cos nθ(T) − 2 Re Z_T(n) + 2∫_0^T cos(nθ(t)) A'(t) dt − ∫_T^∞ S(t) g_n'(t) dt + O(log n/n),
Z_T(n) := Σ_{γ≤T} e^{inθ_γ},  |2∫_0^T cos(nθ)A'| ≪ T² log T/n,  g_n'(t) = 2n θ'(t) sin nθ(t).
Proof of the decomposition: λ_n − λ_n^{(0)} = ∫_{(0,∞)} g_n dS; split at T; Σ_{γ≤T} g_n(γ) = 2N(T) − 2Re Z_T(n);
∫_0^T g_n A' = 2(A(T) − 1) − 2∫_0^T cos(nθ)A'; ∫_{(T,∞)} g_n dS = −g_n(T)S(T) − ∫_T^∞ S g_n' (S = O(log t), g_n = O(n²/t²)).
In v = θ(t), ∫_0^T cos(nθ)A' dt = ∫_{θ(T)}^π cos(nv) a(v) dv, a(v) = A'(t(v))/(4 sin²(v/2)) ~ log(1/(2πv))/(2πv²), monotone near 0;
one integration by parts gives ≪ |a(θ(T))|/n.
Proof of c_0 = −1/4. Write A'(t) = (2π)^{−1} log(t/2π) + r(t): r(t) = −1/(48πt²) + O(t^{−4}) (t → ∞, Stirling),
r(t) = −(2π)^{−1} log t + O(1) (t → 0), so r ∈ L¹ and ∫_0^∞ r = lim_T [θ_RS(T)/π − (T/2π) log(T/2πe)] = −1/8
(θ_RS(T) = (T/2) log(T/2πe) − π/8 + O(1/T)). Then λ_n^{(0)} = I(n) + 2∫r − 2∫_0^∞ cos(nθ(t)) r(t) dt with
I(n) := π^{−1}∫_0^∞ (1 − cos nθ(t)) log(t/2π) dt. In v: I(n) = (2π)^{−1}∫_0^π F_n(v) log(cot(v/2)/(4π)) dv, F_n = sin²(nv/2)/sin²(v/2)
= n + 2Σ_{j<n}(n − j)cos jv, and log(cot(v/2)/(4π)) = 2Σ_{k odd} cos(kv)/k − log 4π (from log(2 sin(v/2)) = −Σ cos kv/k).
Hence I(n) = Σ_{k odd<n} (n − k)/k − (n/2) log 4π = (n/2)(log n + γ_E − 1 − log 2π) + O(1/n) = t(n) + O(1/n)
(Σ_{k odd ≤ 2m−1} 1/k = H_{2m} − H_m/2). The last integral is ∫_0^π cos(nv) ρ(v) dv with ρ(v) = r(t(v))/(4 sin²(v/2)) bounded
near 0 and with a log singularity at π: O(log n/n). So c_0 = 2∫r = −1/4: the constant −π/8 of θ_RS (the 7/8 = 1 − 1/8 of N(T)). ∎

## Theorem 1 (RH). For every ε > 0 and N ≥ N_0(ε): Σ_{n≤N} R(n)² ≥ (c_1 − ε) N^{3/2} (log N)^{1/2},
c_1 = (3/5)^{1/2}(4/5)²/(2π√2) = 0.05579.
Corollary 1. limsup |R(n)|/(n log n)^{1/4} ≥ (3c_1/2)^{1/2} = 0.289 under RH; with Lemma 0, R(n) = Ω((n log n)^{1/4}) unconditionally.
Corollary 2. Under RH Σ_m m^{1/2} y_m² = ∞ (Σ_{m≤N} m^{1/2} y_m² ≥ N^{−3/2} Σ_{m≤N} R(m)² ≥ c √(log N)).
Proof.
1. Weight. ψ ∈ C_c^∞(−1/2, 1/2) even, ψ ≥ 0, ∫ψ² = 1; φ := ψ * ψ: even, supp φ ⊂ [−1,1], 0 ≤ φ ≤ φ(0) = 1, ∫φ = (∫ψ)² ≤ 1
   (→ 1 as ψ → 1_{(−1/2,1/2)}), φ̂(ξ) := ∫φ(u)e^{iuξ}du = ψ̂(ξ)² ≥ 0. Φ(x) := Σ_{n∈Z} φ(n/N)e^{inx} = N Σ_{k∈Z} φ̂(N(x + 2πk)) ≥ 0
   (Poisson); Φ even, real; Φ(0) = N∫φ + O(N^{−K}); |Φ^{(j)}(x)| ≪_{j,K} N^{1+j}(1 + N‖x‖)^{−K}, ‖x‖ = dist(x, 2πZ).
2. R̃(n) := λ_{|n|} − t(|n|) on Z (t(0) := 0; λ_0 = 0). ⟨F,G⟩ := Σ_{n∈Z} φ(n/N)F(n)G(n), ‖F‖² = ⟨F,F⟩;
   ‖R̃‖² = 2 Σ_{1≤n≤N} φ(n/N) R(n)² ≤ 2 Σ_{n≤N} R(n)².
3. For 0 < θ < π (finite n-sum, interchange free): Σ_n φ(n/N)e^{−inθ}λ_{|n|} = Σ_{γ'} B(θ, θ_{γ'}),
   B(θ, θ') := 2Φ(θ) − Φ(θ − θ') − Φ(θ + θ') (= O(θ'² sup_{|x−θ|≤θ'}|Φ''|), Σ_{γ'} θ_{γ'}² < ∞).
4. Trend. τ(θ) := Σ_n φ(n/N)e^{−inθ} t(|n|) = Σ_k F(θ + 2πk), F = Fourier transform of φ(x/N)t(|x|); FT of t(|x|) is
   log(2π|ξ|)/ξ² (from FT|x|^s = −2Γ(s+1)sin(πs/2)|ξ|^{−s−1}, d/ds at s = 1), smoothing by Nφ̂(N·) changes it by
   O(N^{−2}ξ^{−4} log(1/ξ)) for ξ ≥ N^{−1/2}; aliases k ≠ 0 are O(1). So τ(θ) = log(2πθ)/θ² + O(log N) on [N^{−1/2}, 1].
   (Numerical: N = 10^5, θ = 0.003, 0.01, 0.03, 0.0707: τ − log(2πθ)/θ² = −219, −0.97, 0.34, 0.35, logs/check_tau.txt.)
5. Dual vector W(n) := Σ_{γ≤T} cos nθ_γ with T ≤ N^{1/2}/2 (so θ_γ ≥ N^{−1/2} for γ ≤ T). R̃, φ, W even ⇒
   ⟨R̃, W⟩ = Σ_{γ≤T} [Σ_{γ'} B(θ_γ, θ_{γ'}) − τ(θ_γ)].
   Split γ' ≤ X := N³ and γ' > X: the tail is ≪ N(T)·N³(N^{1/2})^{−K}·Σ_{γ'>X}θ_{γ'}² (tiny). For γ' ≤ X:
   Σ_{γ'≤X}[2Φ(θ_γ) − Φ(θ_γ + θ_{γ'})] ≪ X log X · N (N^{1/2})^{−K} (tiny, |θ_γ + θ_{γ'}| ≥ θ_γ ≥ N^{−1/2});
   −Σ_{γ'≤X} Φ(θ_γ − θ_{γ'}) ≤ −Σ_{γ'≤T} Φ(θ_γ − θ_{γ'}) since every term is ≥ 0 (Φ ≥ 0).
   Hence ⟨R̃,W⟩ ≤ −Σ_{γ,γ'≤T} Φ(θ_γ − θ_{γ'}) + τ_T + O(N(T) log N) + tiny,  τ_T := −Σ_{γ≤T} log(2πθ_γ)/θ_γ² ~ T³L²/(6π),
   L := log(T/2π).  And ‖W‖² = ½ Σ_{γ,γ'≤T}[Φ(θ_γ − θ_{γ'}) + Φ(θ_γ + θ_{γ'})] = ½ Σ_{γ,γ'≤T}Φ(θ_γ − θ_{γ'}) + tiny.
   So ⟨R̃, W⟩ ≤ −(2‖W‖² − τ'), τ' := τ_T + O(N(T) log N) + tiny.
6. ‖W‖² ≥ D := ½ Φ(0) N(T) (diagonal γ = γ' only, Φ ≥ 0; multiple zeros only increase it).
   If 2D > τ': Cauchy–Schwarz ‖R̃‖² ≥ ⟨R̃,W⟩²/‖W‖² ≥ (2x − τ')²/x at x = ‖W‖²; (2x − τ')²/x is increasing for x > τ'/2, so
   ‖R̃‖² ≥ (2D − τ')²/D, i.e. Σ_{n≤N} R(n)² ≥ Φ(0)N(T)(1 − τ'/(Φ(0)N(T)))².
7. Constants. N(T) = (T/2π)(L − 1) + O(log T), τ_T = T³L²/(6π)(1 + O(1/L)); put u := T²L/(N∫φ): τ'/(Φ(0)N(T)) = u/3 + O(1/L),
   bound = (∫φ)^{3/2} N^{3/2} L^{1/2} √u (1 − u/3)²/(2π) (1 + O(1/L)); u = 3/5 maximizes √u(1 − u/3)² = 0.49574;
   L = ½ log N + O(log log N); ∫φ → 1. ∎
   Uses RH only through λ_n = Σ_γ g_n(γ) with real θ_γ; no pair correlation, no simplicity. Data: logs/check_tau.txt
   (Σφ R̃² is 4.5–5.1 times the bound at N = 2·10^4, 10^5); logs/analyze.txt (sampled N; over every 100 ≤ N ≤ 10^5 the ratio Σ_{n≤N}R²/(c_1 N^{3/2}(log N)^{1/2}) is 1.472–2.694, minimum at N = 660, liaudit's logs/minratio.txt).

## Theorem 2 (RH). Σ_{N<n≤2N} R(n)² ≪ N^{3/2}(log N)^{3/2}; using Montgomery's pair-correlation theorem (RH),
## Σ_{N<n≤2N} R(n)² ≪ N^{3/2}(log N)^{1/2} log log N.
Corollary 3 (RH, with Corollary 2): Σ_m m^δ y_m² < ∞ ⟺ δ < 1/2. Hence RH ⟺ Σ m^δ y_m² < ∞ for some (equivalently every)
δ ∈ [0, 1/2) (⇐ by Lemma 0); δ = 0 is Arias de Reyna's ℓ² criterion (Funct. Approx. 45 (2011), Thm 1/Cor 2, Thm 6).
Mean-square exponent under RH: lim log(Σ_{n≤N}R²)/log N = 3/2 exactly.
Proof. Weight φ ∈ C_c^∞((1/2, 4)), 0 ≤ φ, φ ≥ 1 on [1,2]; Φ_+(x) := Σ_n φ(n/N)e^{inx}, |Φ_+(x)| ≪ N(1 + N‖x‖)^{−K}.
Lemma A with T² log T = N: all terms except Z_T and E(n) := ∫_T^∞ S g_n' are O(log N) for n ∈ supp.
(a) Σ_n φ(n/N)|Z_T(n)|² = Σ_{γ,γ'≤T} Φ_+(θ_γ − θ_{γ'}) ≪ N Σ_{γ,γ'≤T}(1 + ½ log T·|γ − γ'|)^{−K}
    (|θ_γ − θ_{γ'}| ≥ |γ − γ'|/(T² + 1/4)). Unconditionally (N(t + h) − N(t) ≪ (1 + h) log t) this is ≪ N·N(T) log T;
    under RH Montgomery's theorem with the majorant r(u) = C[sinc²(u/2) + sinc²((u+1)/2)] (supp r̂ ⊂ [−1/2, 1/2]) gives ≪ N·N(T).
(b) Σ_n φ(n/N) E(n)² = ∫∫_{[T,∞)²} S(t)S(t')K(t,t'), K = 2θ'(t)θ'(t') Re Σ_n φ(n/N) n²[e^{in(θ−θ')} − e^{in(θ+θ')}],
    |K| ≪ N³(tt')^{−2}[(1 + N|θ − θ'|)^{−K} + (1 + N(θ + θ'))^{−K}], ∫|K(t,t')|dt' ≪ N²/t². Schur:
    ≪ N² ∫_T^∞ S(t)² dt/t² ≪ N² log log T/T (Selberg, unconditional: ∫_0^U S² = (U/2π²) log log U (1 + o(1))).
(c) N(T) ≍ (N log N)^{1/2}, T ≍ (N/log N)^{1/2}: (a) ≪ N^{3/2}(log N)^{1/2}·{log N | 1}, (b) ≪ N^{3/2}(log N)^{1/2} log log N. ∎

## Theorem 3 (RH). |R(n)| ≤ (κ_3 + o(1)) √n log n/(log log n)^{1/2}, κ_3 = (2/(3π))(π²/12)^{1/4} = 0.202.
Proof. Lemma A: |2 Re Z_T| ≤ 2N(T) = (T/π) log T (1 + o(1)); ∫_T^∞ S g_n' = −S_1(T)g_n'(T) − ∫_T^∞ S_1 g_n'',
|g_n''| ≤ 2n²/t⁴ + 4n/t³, −(π/24 + o(1))h ≤ S_1(t) ≤ (π/48 + o(1))h, h = log t/(log log t)² (Carneiro–Chandee–Milinovich, Math. Ann. 2013), so |S_1| ≤ (π/24 + o(1))h;
total ≤ (T/π) log T + (2C_1/3) n² log T/(T³ (log log T)²) + o(·), minimized at T⁴ = 2πC_1 n²/(log log n)², C_1 = π/24. ∎ (Reviewer: centring the asymmetric CCM bounds and
averaging |cos nθ| in the leading g_n'' term gives effective C_1 = 1/16, κ_3 = (2/(3π))(π/8)^{1/4} = 0.168.)
(With |S(t)| ≤ (1/4 + o(1)) log t/log log t alone: constant (2π)^{−1/2} = 0.399.) Data: max_{100≤n≤10^5} |R(n)|/(0.202 √n log n/√(log log n)) = 0.37.

## Theorem 4 (RH). For φ ∈ C_c^∞((0,∞)), φ ≥ 0, φ ≢ 0: Σ_n φ(n/N) R(n) = (7/4) Σ_n φ(n/N) + O(log N).
So the smoothed mean of R is 7/4 = 2·(7/8), and Maślanka's oscillating part λ̃_n = λ_n − λ̄_n (λ̄_n = Li coefficients of
π^{−s/2}Γ(1+s/2) = t(n) + 3/4 − 1/(24n) + …) has smoothed mean 1, not 0.
Proof. Φ_+(x) := Σ_n φ(n/N)e^{inx}; G_N(t) := Σ_n φ(n/N) g_n(t) = 2Φ_+(0) − 2 Re Φ_+(θ(t)). By Lemma A,
Σφ(λ_n^{(0)} − t(n)) = −Φ_+(0)/4 + O(log N), and Σφ(λ_n − λ_n^{(0)}) = ∫_{(0,∞)} G_N dS = G_N(0+) − ∫_0^∞ S G_N' dt
(S(0+) = −1), G_N(0+) = 2Φ_+(0) − 2Φ_+(π) = 2Φ_+(0) + O(N^{−K}). ∫_0^∞ S G_N' = −∫_0^∞ S_1 G_N'', S_1(t) = ∫_0^t S = O(log(t+2))
(unconditional); |G_N''(t)| ≪ N^{−K} (t ≤ N^{1/2}), ≪ (N³/t⁴)(t/N)^K (N^{1/2} ≤ t ≤ N, from |Σ φ n^j e^{inθ}| ≪ N^{j+1}(Nθ)^{−K}),
≪ N³/t⁴ (t ≥ N); so ∫|S_1 G_N''| ≪ log N. ∎  Data (logs/tiny.txt, φ = bump on (1/2,1)): smoothed mean of R = 1.8874, 1.7704,
1.7443, 1.75115, 1.749898 at N = 10³, 3·10³, 10⁴, 3·10⁴, 10⁵ (|dev| = 1.0e−4 at 10⁵; log N/N = 1.2e−4); sharp mean 1.689.

## Corollary 4 (refutation). Merlini–Sala–Sala (arXiv 1903.06233, §2, "numerical conjecture"): |λ_tiny(n)| = a·log n with a ≈ 2,
"the tiny oscillations may be bounded in absolute value by a·log(n)"; λ_tiny = λ̃_n. False for every a: λ̃_n = R(n) − 3/4 + O(1/n)
= Ω((n log n)^{1/4}) by Corollary 1 (unconditional). Numerically (two DFT radii agree to 25 digits): first n with |λ̃_n| > 2 log n is
n = 5083 (λ̃_5083 = 17.158913, 2 log 5083 = 17.067314); 7584 such n ≤ 10⁵; max |λ̃_n|/log n = 4.18 at n = 79302 (λ̃ = 47.170651);
a = 3 fails from n = 46992, a = 4 from n = 79286. (Their evidence: computed values for n ≤ 31, Maślanka's values to n = 3000 and plots to n = 8000; liaudit.)

## Barrier (pointwise exponent). Theorems 1–3: 1/4 ≤ inf{c: R(n) = O(n^c)} ≤ 1/2 under RH. Lowering 1/2 needs cancellation
in Z_T(n) = Σ_{γ≤T} e^{inθ_γ} at T ≍ √n, n θ_γ ≈ n/γ ≥ √n, which RH does not supply (Lemma A is exact; the trivial bound
2N(T) is the whole loss). Same obstruction as von Koch's √x log² x for ψ(x) − x, unimproved under RH since 1901.
