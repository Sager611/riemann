# PROOF — small gaps below 7/11 + ε unconditionally (the fourth moment replaced by a zero-based pointwise bound)

**Theorem U′.** Let κ > 0 be a uniform density slope: N(σ, 3T) ≪ T^{1−κ(σ−½)} log T for all σ ≥ ½ (Conrey κ = 8/7 − δ, proved in Baluyot's thesis, audited in Riemann/Deep/Baluyot; Jutila κ = 1 − δ). For every φ ∈ (max{2/(2+κ), ½}, 1) there is c(φ) > 0 with #{n : T ≤ γ_n ≤ γ_{n+1} ≤ 2T, γ_{n+1} − γ_n < 2πφ/log T} ≥ c(φ) T log T for large T (all zeros, multiplicity, as in Theorem U). With κ = 8/7 − δ: every φ > 7/11, so D(λ) > 0 for every λ > 7/11 = 0.6364 (Theorem U: 2/3).

**What changes.** Riemann/Deep/DPCUncond/PROOF.md (Theorem U) is used verbatim except: θ ∈ (1/φ − 1, min{κ/2, 1}) instead of (1/φ − 1, ½) in Step 0 (nonempty iff φ > max{2/(2+κ), ½}; the cap θ < 1 matters only if κ > 2) (the construction needs only 0 < φθ < 1 < φ(1+θ) < 2, θ < 1); Step 5a (Montgomery–Vaughan fourth moment, the only θ < ½ input) is replaced by Lemma Z and Lemma M below; Steps 6e–6f and the Hölder step of Step 7 are replaced by 6e′ and 7′. Steps 1–4, 5b, 5c, 6a–6d are unchanged; none of them uses θ < ½ (6c's Poisson tails for R_T, which needed θ < 2/3, are no longer used). Conventions here are DPC's (2610.28527 l.614–626, 673–679): ρ_η = ψ(u)ω(u/η) smooth with supp ⊂ (η, b), b > θ; f smooth with compact support, f = 0 on [θ, ∞); c_n = μ(n) n^{−1/2} f(log n/L) ∏_{q|n} ρ_η(log q/L); n < T^θ squarefree with ≤ k = ⌊θ/η⌋ prime factors, all > T^η. X := T^θ, L := log T, d_ρ := β − ½, d⁺ := max(d, 0), κ_A(x) := min(1, x^{−A}) with A fixed large.

## Lemma Z (pointwise bound through the zeros)

Let On(t) := Σ_ρ κ_A(L|γ − t|), W(t) := Σ_ρ X^{d_ρ⁺} κ_A(L|γ − t|) (all nontrivial zeros, multiplicity), W_{−y}(t) := Σ_ρ X^{(d_ρ − y)⁺} κ_A(L|γ − t|). For t ∈ [T/2, 3T] and 0 ≤ y ≤ ½,
  |R_T(t − iy)| ≤ C (1 + W_{−y}(t)) (1 + On(t))^{k−1},   |R*(t − iy)| ≤ C X^{y} (1 + W(t)) (1 + On(t))^{k−1},
with C = C(f, ψ, ω, η, θ, b, A). (R*(t − iy) = conj R_T(t + iy) with f̄; so both are the case s = ½ ± y + it.)

*Proof.* (i) Explicit formula. φ(u) := ρ_η(u/L)/u ∈ C_c^∞((ηL, bL)), Φ(z) := ∫ φ(u) e^{zu} du. Perron and a shift to Re → −∞ (Φ decays like e^{ηL Re z} there) give, for all s,
  Q(s) := Σ_{n≥2} Λ(n) φ(log n) n^{−s} = Φ(1 − s) − Σ_ρ Φ(ρ − s) − Σ_{m≥1} Φ(−2m − s).
Λ(q)φ(log q) = ρ_η(log q/L), so Q = Σ_q ρ_η(log q/L) q^{−s} + Q_pp, Q_pp the prime powers.
(ii) Products. Fix s = ½ + y + it (R_T(t − iy)) or s = ½ − y + it (R_T(t + iy), i.e. R*(t − iy)). With x_q := ρ_η(log q/L) q^{−s}, n = q_1⋯q_j squarefree has c_n n^{1/2−s} = (−1)^j f(Σ log q_i/L) ∏ x_{q_i}, so the polynomial at s is Σ_{j≤k} (−1)^j Σ_{q_1<…<q_j} f(Σ log q_i/L) ∏ x_{q_i}. Remove distinctness by Möbius inversion on set partitions: each partition π of {1..j} gives Σ over one prime per block, block B carrying x_q^{|B|}. Remove prime powers by the binomial expansion of (Q − Q_pp). Every resulting term is (a) an absolutely convergent sum over packets: a collision block B with prime q is the packet q^h, h = |B| ≥ 2, with coefficient ρ_η(log q/L)^h q^{−hs} (T^η < q < T^b); a prime-power term from Q_pp is the packet q^h, h ≥ 2, with coefficient ρ_η(h log q/L) h^{−1} q^{−hs} (T^{η/h} < q < T^{b/h}; repair after review 2); let 𝔞 := ∏ q^h be the total packet product (repair after review 1: the shift is h log q/L, not log q/L); times (b) a free sum Σ_{q'_1..q'_m} f_𝔞(Σ log q'_i/L) ∏ x_{q'_i}, f_𝔞(v) := f(v + log 𝔞/L), of the same kind with fewer variables (f_𝔞 = 0 on [θ − log 𝔞/L, ∞), so the free variables satisfy Σu_i ≤ log(X/𝔞)). By Fourier inversion f_𝔞(Σu_i/L) = ∫ f̂_𝔞(ξ) ∏ e^{2πiξu_i/L} dξ, the free sum equals ∫ f̂_𝔞(ξ) ∏_i Q_ξ dξ (Q at s − 2πiξ/L) minus prime-power terms of type (a), and inserting (i) and Fubini (absolute convergence: Σ_ρ |Φ(ρ − s + 2πiξ/L)| ≪ T^{b} log(2 + |t| + |ξ|), polynomial in ξ, and f̂_𝔞 is Schwartz uniformly in 𝔞) it is a sum over ρ_1, …, ρ_m of
  I(ρ_1..ρ_m) := ∫_{R^m} f_𝔞(Σu_i/L) ∏ φ(u_i) e^{z_i u_i} du,  z_i := ρ_i − s,
plus pole and trivial-zero terms (z_i = 1 − s or −2m' − s), which are O(T^{−10}) since |Im s| ≍ T and IBP (iii) applies.

(iii) The decisive bound. On the support of the integrand u_i > ηL > 0 and Σu_i ≤ θL − log 𝔞 = log(X/𝔞) (f_𝔞 = 0 on [θ − log 𝔞/L, ∞)). With Re z_i = d_i − y (case s = ½ + y + it) or d_i + y (case s = ½ − y + it),
  |e^{Σ z_i u_i}| = e^{Σ (Re z_i) u_i} ≤ e^{max_i (Re z_i)⁺ · Σ u_i} ≤ (X/𝔞)^{max_i (Re z_i)⁺}.
The integrand vanishes to infinite order on the boundary of its support, and ∂_u^α of f_𝔞(Σu/L)∏φ(u_i) is ≪_α L^{−|α|−m} on a set of volume ≪ L^m; integrating by parts A times in each u_i with L|z_i| ≥ 1 gives
  |I(ρ_1..ρ_m)| ≤ C_{A,m} (X/𝔞)^{max_i (Re z_i)⁺} ∏_i κ_A(L|z_i|),  |z_i| ≥ |γ_i − t|.
This is where the length constraint is kept: an r-fold product of prime sums of length X costs X^{max d}, not X^{r·max d} as in the route |e_r| ≪ (1 + |S_1|)^r (DPC l.665–670), which is harmless under RH and fatal without it.
(iv) Summation. X^{max_i c_i} ≤ Σ_i X^{c_i} for c_i ≥ 0, so Σ_{ρ_1..ρ_m} X^{max(Re z_i)⁺} ∏ κ_A ≤ m · W_± · On^{m−1}, where W_+ := W_{−y} in the first case, and in the second case (d + y)⁺ ≤ y + d⁺ gives W_− ≤ X^{y} W. Packets carry |𝔞^{−s}| = 𝔞^{−1/2 ∓ y}. In the first case (Re z_i = d_i − y) the free part carries (X/𝔞)^{max(d_i − y)⁺} ≤ X^{max(d_i−y)⁺}. In the second case it carries (X/𝔞)^{y + max d⁺} and the product is 𝔞^{−1/2+y}(X/𝔞)^{y + max d⁺} = X^{y + max d⁺} 𝔞^{−1/2 − max d⁺} ≤ X^{y + max d⁺} 𝔞^{−1/2}; with no free variable left, 𝔞 < X gives 𝔞^{−1/2+y} ≤ X^{y}𝔞^{−1/2}. Finally Σ_𝔞 C(𝔞)𝔞^{−1/2} ≤ ∏ over the ≤ k packets of (Σ_{T^η<q<T^b} ρ_η^h q^{−h/2} + Σ_{q, h≥2} |ρ_η(h log q/L)| h^{−1} q^{−h/2}) = O_η(1), C(𝔞) the product of the |packet coefficients| times 𝔞^{Re s} (so the packet factor 𝔞^{−s} is stripped): h = 2 by Mertens (Σ_{T^η<q<T^b} 1/q and Σ_{T^{η/2}<q<T^{b/2}} 1/q are O(log(b/η))), h ≥ 3 by Σ_q Σ_{h≥3} q^{−h/2} < ∞. The number of partitions and binomial terms is ≤ C(k). This proves Lemma Z. ∎

## Lemma M (moments)

Bins B_a := {ρ : γ ∈ [T/4, 3T], d_ρ L ∈ [a, a+1)}, a ≥ 0; |B_a| ≤ C TL e^{−κa} (density). Off(t) := Σ_{ρ: d_ρ > 0} X^{d_ρ} κ_A(L|γ − t|), so W ≤ On + Off and Off ≤ Σ_a e^{θ(a+1)} F_a, F_a := Σ_{ρ∈B_a} κ_A(L|γ − t|), plus O(T^{−5}) for t ∈ [T/2, 5T/2] (zeros with γ ∉ [T/4, 3T] have L|γ − t| ≥ TL/4, so their total weight is ≤ X^{1/2} Σ_ρ κ_A(L|γ − t|) ≪ T^{−5} for A large).
(M1) ‖On‖_{L^p[T/2,5T/2]} ≪_p T^{1/p}: On ≤ Σ_m (1 + |m|)^{−A} n_m, n_m(t) := #{ρ : m ≤ L(γ − t) < m + 1}, Fujii ∫ n_m^p ≪ T (Step 5b) for |m| ≤ TL/4 (heights comparable to T), the kernel tail κ_A(|m|) ≤ |m|^{−A} times the trivial n_m ≪ log(T + |m|/L) beyond, Minkowski.
(M2) ‖F_a‖_p ≪_{p,w} T^{1/p} e^{−κa(1−1/w)/p} for every w > 1: n^{(a)}_m := #{ρ ∈ B_a : m ≤ L(γ − t) < m + 1} ≤ n_m is an integer, so (n^{(a)}_m)^p ≤ 1_{n^{(a)}_m ≥ 1} n_m^p and Hölder gives ∫ ≤ meas{n^{(a)}_m ≥ 1}^{1−1/w} ‖n_m‖_{pw}^p ≤ (|B_a|/L)^{1−1/w}(CT)^{1/w}.
(M3) Hence ‖Off‖_p ≪ T^{1/p} Σ_a e^{θ(a+1) − κa(1−1/w)/p} ≪ T^{1/p} whenever pθ < κ (choose w large), and by Lemma Z and Hölder
  ∫_T^{2T} |R_T|^q ≪ T  for every q < κ/θ.  (With κ = 8/7 − δ and θ < 4/7 there is q > 2.)
(M4) Discrete versions over the zeros (cut points x_ρ = γ_ρ ∓ h/2, γ_ρ ∈ [T, 2T]): Σ_ρ On(x_ρ)^p ≪ TL; Σ_ρ Off(x_ρ)^p ≪ TL whenever pθ < κ. Proof as M1–M3 with Lz Lemma DoubleSumZeros (Σ_j n_j^{p+1} ≪ L∫(N(t + c/L) − N(t))^{p+1}dt ≪ TL) in place of Fujii, and #{ρ : n^{(a)}_m(x_ρ) ≥ 1} ≤ Σ_{ρ'∈B_a} Ñ^{(m)}_{ρ'} ≤ |B_a|^{1−1/w}(Σ_{all ρ'} (Ñ^{(m)}_{ρ'})^w)^{1/w} ≪ (|m|+2) TL e^{−κa(1−1/w)}, Ñ^{(m)}_{ρ'} := #{ρ : |γ − γ'| ≤ (|m|+1)/L + h} (m ∈ ℤ; a B_a zero with m ≤ L(γ' − x_ρ) < m + 1 has |γ' − x_ρ| ≤ (|m|+1)/L and |x_ρ − γ_ρ| = h/2). The exponent becomes κ(1−1/w)², still → κ.

## Step 6e′ (replaces 6e–6f): hairpins, |𝓗_T| ≤ Φ(ε)T + o(T)

Notation of PROOF.md 6a–6d: S_a (zeros with γ ∈ [T, 2T], d_ρ L ∈ [a, a+1)), σ_a := min(2ε + o(1), C e^{−κa}), Y_a = (a+1)/L, cut point x_ρ, and 𝓗_{T,a} ≪ Σ_{ρ∈S_a} ∫_0^{d_ρ} |R_T R*|(x_ρ − iy) Z(x_ρ − iy) dy with Z := (1 + |V_L| + 2Ñ_ρ)^D (6a).
(a) Lemma Z with W_{−y} ≤ On + X^{−y}Off (a zero with d ≤ y contributes 1, one with d > y contributes X^{−y}X^{d}) and W ≤ On + Off:
  |R_T R*|(x − iy) ≤ C [X^{y}(1 + On) + Off](1 + On + Off)(1 + On)^{2k−2} ≤ C [X^{y} + Off](1 + Off)(1 + On)^{2k}, all at x.
(b) y-integral. Z̄_ρ := (Y_a^{−1}∫_0^{Y_a} Z(x_ρ − iy)^{r'} dy)^{1/r'} (defined for every zero, any fixed r' > 1). Hölder in y and ∫_0^{d_ρ}X^{ry}dy ≤ X^{rd_ρ}/(r log X):
  ∫_0^{d_ρ} X^{y} Z dy ≤ X^{d_ρ}(rθL)^{−1/r} Y_a^{1/r'} Z̄_ρ ≤ C (a+1) e^{θ(a+1)} L^{−1} Z̄_ρ,  ∫_0^{d_ρ} Z dy ≤ Y_a Z̄_ρ.
  With 2e^{θ(a+1)}Off ≤ e^{2θ(a+1)} + Off², for ρ ∈ S_a: ∫_0^{d_ρ}|R_T R*| Z dy ≤ C (a+2)^2 L^{−1} [e^{2θ(a+1)} + Off(x_ρ)²] Y_ρ, Y_ρ := (1 + On(x_ρ))^{2k} Z̄_ρ (also using 1 + Off ≤ 2 max(1, Off)).
(c) Discrete moments over all zeros with γ ∈ [T, 2T] (cut points x_ρ, segments 0 < y < Y_a): Σ_ρ Y_ρ^w ≪ (a+2)^C TL for every w, and Σ_ρ Off(x_ρ)^{2u} ≪ TL whenever 2uθ < κ (M4). For Y: Hölder splits On (M4) from Z̄. For Z̄ take r' ≤ m, so Z̄^m ≤ avg_y Z^m (power means), and Z ≤ C(1 + |V̂| + |Δ_L| + Ñ_ρ)^D with V_L = V̂ + Δ_L (6b). V̂ part: |V̂(x − iy)| ≤ C(P_y * v)(x) + C and Jensen give Σ_ρ avg_y (P_y * v)^{Dm}(x_ρ) ≤ ∫ v^{Dm} K, K := Σ_ρ Y_a^{−1} k_{Y_a}(· − x_ρ) (k_d as in 6c, ∫ Y^{−1}k_Y = 1), and ∫ v^{Dm}K ≤ ‖v^{Dm}‖_{L^r(I)} ‖K‖_{r''} ≪ T^{1/r} · (a+2)^C L T^{1/r''} (5c; 6c's near/far split with all zeros, keeping the logarithmic near kernel k_{Y_a} inside the L^{r''} norm, then Fujii). Δ_L part: 6d's estimate with Y_a in place of d_ρ, avg_y |Δ_L|^{Dm} ≤ C(1 + h/Y_a)(1 + Ñ_ρ + Ψ_ρ)^{Dm}, and Σ_ρ (Ñ_ρ^p + Ψ_ρ^p) ≪ (a+2)^p TL (6c, 6d).
(d) Sum over S_a, 1/v + 1/u + 1/w = 1, u > 1 with 2uθ < κ (possible iff θ < κ/2):
  Σ_{S_a} Y_ρ ≤ |S_a|^{1−1/w}(Σ Y^w)^{1/w} ≪ (a+2)^C σ_a^{1−1/w} TL,
  Σ_{S_a} Off² Y ≤ |S_a|^{1/v}(Σ Off^{2u})^{1/u}(Σ Y^w)^{1/w} ≪ (a+2)^C σ_a^{1/v} TL.
  Hence 𝓗_{T,a} ≪ T (a+2)^C [e^{2θ(a+1)} σ_a^{1−1/w} + σ_a^{1/v}] and
  Σ_a 𝓗_{T,a} ≪ T Σ_a (a+2)^C [e^{2θ(a+1)} σ_a^{1−1/w} + σ_a^{1/v}] =: Φ(ε)T.
(e) Φ(ε) → 0. Fix w with κ(1 − 1/w) > 2θ, then A = c log(1/ε) with 0 < c < (1 − 1/w)/(2θ): Σ_{a≤A} ≪ (A+2)^C e^{2θ(A+1)} ε^{1−1/w} ≪ (log 1/ε)^C ε^{1−1/w−2θc} → 0; Σ_{a>A} ≪ (A+2)^C e^{−(κ(1−1/w) − 2θ)A} → 0; the σ_a^{1/v} series → 0 by dominated convergence (majorant (a+2)^C e^{−κa/v}). Cut points of zeros outside [T, 2T]: G_T ≤ e^{−L²/32}, the trivial |R_T R*| ≤ T^{2θ} and the polynomial majorant of the integrated V-factor (Step 3). ∎

## Step 7′ (replaces the Hölder step of Step 7)

Fix q with 2 < q < κ/θ (M3) and s with 2/q + 1/s < 1. With Re P(V) ≥ 0 outside E_T (Step 1) and Step 6e′,
  δ₁T ≤ ∫_{E_T} G_T |R_T|² |P(V)| ≤ ‖R_T‖_q² ‖P(V)‖_s meas(E_T)^{1−2/q−1/s} ≪ T^{2/q + 1/s} meas(E_T)^{1−2/q−1/s}
(5c: ‖P(V)‖_s ≪ T^{1/s}), so meas(E_T) ≥ c₂T with c₂ independent of ε. The rest of Step 7 is unchanged. ∎

## Order of choices
φ ∈ (max{2/(2+κ), ½}, 1); θ ∈ (1/φ − 1, min{κ/2, 1}); P, q (Step 0), f, η (DPC Lemma 6), A; q ∈ (2, κ/θ), s; u ∈ (1, κ/(2θ)), w > u/(u − 1) with κ(1 − 1/w) > 2θ, 1/v := 1 − 1/u − 1/w, r, r'; then ε (6e′ (e), ε ≤ ε₁, ε₀ as in Step 7); then T → ∞. All constants depend on φ, θ, f, P, η, κ, δ, A, w, u, v, r, r' only, not on ε.

## Inputs (beyond those of PROOF.md for Theorem U)

| Input | Source | Used in |
|---|---|---|
| Explicit formula for Σ Λ(n)φ(log n)n^{−s}, φ ∈ C_c^∞((0,∞)), all zeros | Perron + contour shift (standard; IK §5.5, Rodgers Thm 2.1 as cited by DPC l.643) | Lemma Z (i) |
| Uniform density N(σ,3T) ≪ T^{1−κ(σ−½)} log T, κ = 8/7 − δ | Conrey, via Baluyot (audited, Riemann/Deep/Baluyot) | M2–M4, 6e′ |
| Fujii moments, Lz DoubleSumZeros | as in PROOF.md 5b, 6c | M1, M4 |
Montgomery–Vaughan is no longer used.

## Why the threshold is κ/2 (heuristic; sufficiency is proved, necessity is not)
- Heuristically an isolated off-line zero at depth a/L makes R_T ≈ X^{d} near γ (the residue of the mollifier-like polynomial; N1 shows the bound I_j ≪ X^{d} is attained up to polynomial factors by the integral itself, but no lower bound for R_T is proved). In the hairpin, R_T(x − iy) ~ X^{d−y} and R*(x − iy) ~ X^{d+y}: the product X^{2d} = e^{2θa}, against e^{−κa} zeros: 2θ < κ. In Step 7 the same zeros put ∫|R_T|^q at Σ e^{(qθ−κ)a}: q > 2 needs 2θ < κ. Both limits coincide with Step 8's heuristic 'off-line zeros alone allow θ < κ/2'.
- The gain over Theorem U is exactly that Conrey's κ = 8/7 exceeds 1: with Jutila's κ = 1 − δ the theorem gives θ < ½, i.e. Theorem U again. Conrey's κ comes from his mollifier of length T^{4/7} (Deshouillers–Iwaniec); the argument imports 'beyond T^{1/2}' through the density theorem, without any mean value theorem for R_T of length > T^{1/2}.
- Under RH (κ = ∞) Lemma Z gives DPC's Lemma 5 for all θ < 1 and all p, without Rodgers' tail estimate (Fujii's moments suffice).

## Status
Proof sketch at the level of PROOF.md for Theorem U; Lemma Z is complete; M4 and 6e′(c) reuse 6b–6d estimates for all zeros instead of S_a (same proofs with σ = 1). Review 1 (two reviewers, 12 Nov 2026, logs/review1.txt): Lemma Z GAP (packet bookkeeping: block prime q vs packet q^h), repaired in (ii)/(iv) as the reviewer proposed; Lemma M, 6e′, 7′ CORRECT (minor: |m|+2, Fujii only at comparable heights; both applied). Review 2 (logs/review2.txt): re-review of Lemma Z GAP (prime-power packets carry ρ_η(h log q/L)/h, support down to T^{η/h}), repaired in (ii)/(iv) exactly as proposed, 'no further obstruction in (i)–(iv)'; adversarial reviewer: NO FLAW FOUND (mechanism sound, no contradiction with known results; asks that the other θ < ½ uses be audited — PROOF.md Step 5a states the fourth moment is the only one besides Step 6, confirmed by text search, D5). Review 3 (logs/review3.txt): Step 6e′ with the all-zero moments (c) CORRECT (two wording repairs applied); audit of Theorem U's Steps 0, 3, 7 for other θ < ½ uses: NO OTHER θ<1/2 USE. Review 4 (logs/review4.txt): third pass on the repaired Lemma Z: CORRECT (one notation fix, applied).
Referees (orchestrator, 12 Nov 2026): A (Lemma Z) CORRECT, SIGNIFICANT; B (Lemma M, 6e′) GAP, repairable — M4's radius (m+1)/L is negative for some m ∈ ℤ; replaced by (|m|+1)/L (also in DERIVATION.md D4), no other exponent obstruction, SIGNIFICANT; C (Step 7′, order of choices, assembly) GAP, repairable, SIGNIFICANT; Step 7′ valid, order of choices consistent.
CLOSED re-reviews (orchestrator, 12 Nov 2026): B′ CLOSED (covering valid for every m ∈ ℤ, exponent κ(1 − 1/w)², kernel sum converges for A large), C′ CLOSED (every φ in the range admits θ, q, u, w; κ = 8/7 − δ gives exactly φ > 7/11; no remaining sharpness claim).
Repairs after the referees: (1) B: M4 radius (m+1)/L → (|m|+1)/L for m ∈ ℤ. (2) C: the construction needs θ < 1 and φ(1+θ) > 1, so the theorem is stated for φ ∈ (max{2/(2+κ), ½}, 1) with θ ∈ (1/φ − 1, min{κ/2, 1}) (Theorem U′, 'What changes', Order of choices); for κ = 8/7 nothing changes (7/11 > ½). (3) C: 'Lemma Z is sharp' removed; the κ/2 threshold is marked heuristic (sufficiency proved, necessity not), here and in REPORT.md.

