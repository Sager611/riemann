# FULLPROOF — line-by-line referee of Theorem U′ (λ > 7/11 unconditionally)

Referee: upaudit (Orca), 12 Nov 2026. Object: Riemann/Deep/DPCBeyond/PROOF.md (Theorem U′), on top of Riemann/Deep/DPCUncond/PROOF.md (Theorem U, "U" below) as audited in Riemann/Deep/UAudit/FULLPROOF.md. DPC = 2610.28527v2 Gaps_for_Dirichlet.tex, Lz = 2611.14486 SmallGapsPaper.tex (both re-downloaded, gzip/tar checked).

## 0. Conventions (checked against DPC l.614–626, 673–684)

- ψ, ω smooth, 0 ≤ ψ, ω ≤ 1, ψ = 1 on [0, θ], ψ = 0 on [b, ∞), b > θ; ω = 0 on (−∞, 1], ω = 1 on [2, ∞) (DPC l.673). ρ_η(u) := ψ(u)ω(u/η), so supp ρ_η ⊂ [η, b] (DPC l.679 states the weaker "⊂ (η/2, 2b)").
- f smooth, compact support, f = 0 on [θ, ∞) (DPC l.619), possibly complex. c_n := μ(n)n^{−1/2}f(log n/L)∏_{q|n}ρ_η(log q/L) (DPC l.676).
- c_n ≠ 0 ⇒ n squarefree, n < X := T^θ, every prime factor q > T^η, so j primes give T^{jη} < X: j < θ/η, j ≤ k := ⌊θ/η⌋. (U Step 0 says "> T^{η/2}, ≤ ⌊2θ/η⌋": also true, weaker. U′'s "(η, b)" is right up to closure.) **Cosmetic C1.**
- R_T(t) = Σc_n n^{−it}, R*(t) = Σc̄_n n^{it}. With M(s) := Σ_n μ(n)f(log n/L)∏_{q|n}ρ_η(log q/L)n^{−s} and s_± := ½ ± y + it:
  R_T(t − iy) = Σc_n n^{−it−y} = M(s₊);  conj R*(t − iy) = Σc_n n^{−it+y} = M(s₋).  So |R*(t − iy)| = |M(s₋)|.
- L := log T, X := T^θ, d_ρ := β − ½, d⁺ := max(d, 0), κ_A(x) := min(1, x^{−A}) (nonincreasing on [0, ∞)). On, W, W_{−y} as in PROOF.md: On(t) = Σ_ρκ_A(L|γ − t|), W(t) = Σ_ρX^{d_ρ⁺}κ_A(L|γ − t|), W_{−y}(t) = Σ_ρX^{(d_ρ−y)⁺}κ_A(L|γ − t|); all nontrivial zeros (both half-planes), multiplicity.
- φ(u) := ρ_η(u/L)/u. φ ∈ C_c^∞(ℝ), supp φ ⊂ [ηL, bL] (ω vanishes to infinite order at 1). Leibniz with u ≥ ηL: |φ^{(a)}(u)| ≤ Σ_{c≤a} C(a,c)(a−c)!‖ρ_η^{(c)}‖_∞L^{−c}(ηL)^{−1−a+c} =: C_a(ψ, ω, η)L^{−1−a}  (Z0)
- Φ(z) := ∫φ(u)e^{zu}du, entire. For Re z ≥ 0: |Φ(z)| ≤ ‖φ^{(a)}‖₁e^{bL Re z}|z|^{−a}; for Re z ≤ 0: ≤ ‖φ^{(a)}‖₁e^{ηL Re z}|z|^{−a} (a IBPs, a = 0, 1, 2, …), ‖φ^{(a)}‖₁ ≤ bL·C_aL^{−1−a}.

## 1. Lemma Z (complete proof)

**Statement.** There are C = C(f, ψ, ω, η, θ, b, A) and T₀ such that for T ≥ T₀, t ≥ T/2 (in particular t ∈ [T/2, 3T]) and 0 ≤ y ≤ ½:
  (Z+) |R_T(t − iy)| = |M(s₊)| ≤ C(1 + W_{−y}(t))(1 + On(t))^{k−1},
  (Z−) |R*(t − iy)| = |M(s₋)| ≤ C X^{y}(1 + W(t))(1 + On(t))^{k−1}.
(For k = 0 read (1 + On)^{−1} as 1; then M = f(0).)

**Z1 (explicit formula; no RH).** For every s ∈ ℂ, with Q(s) := Σ_{n≥2}Λ(n)φ(log n)n^{−s} (a finite sum: log n ∈ [ηL, bL]):
  Q(s) = Φ(1 − s) − Σ_ρΦ(ρ − s) − Σ_{m≥1}Φ(−2m − s),  all series absolutely convergent.
Proof. Fourier inversion of φ(u)e^{−cu}: φ(u) = (2πi)^{−1}∫_{(c)}Φ(w)e^{−wu}dw for every real c (|Φ(c + iτ)| ≪_{c,T} (1+|τ|)^{−2}). For c > 1 − Re s: Q(s) = (2πi)^{−1}∫_{(c)}Φ(w)(−ζ'/ζ)(s + w)dw (absolute convergence, finite n-sum anyway). Move the line to Re w = −M_j with Re(s + w) = −2M_j' − 1 (midway between trivial zeros), truncating at heights ±U_ℓ with |ζ'/ζ(σ + iU)| ≪ log²U for −2M_j−1 ≤ σ ≤ c + Re s (Davenport ch. 17; for σ ≤ −1 by the functional equation |ζ'/ζ| ≪ log(|σ| + |U|)). Horizontal pieces: ≤ (M_j + c)·sup|Φ|·log²U_ℓ ≪_j U_ℓ^{−2}log²U_ℓ → 0 (a = 2 above). Left line: |Φ(w)| ≤ ‖φ''‖₁e^{−M_jηL}|w|^{−2} and |ζ'/ζ(s + w)| ≪ log(M_j + |w|), so its integral is ≪ e^{−M_jηL}log M_j → 0 (ηL > 0; this is why supp φ ⊂ [ηL, ∞) matters). Residues of −ζ'/ζ(s + w): +1 at s + w = 1, −m_ρ at s + w = ρ, −1 at s + w = −2m. Absolute convergence: Σ_ρ|Φ(ρ − s)| ≤ C_T e^{bL}Σ_ρ min(1, |γ − Im s|^{−2}) < ∞; |Φ(−2m − s)| ≤ ‖φ‖₁e^{−(2m + Re s)ηL}, summable. ∎ (Standard: IK Thm 5.12-type; DPC l.643 uses the same with Rodgers.) Taking φ(log n) at n = q: Λ(q)φ(log q) = ρ_η(log q/L); at n = q^h: Λ(q^h)φ(h log q) = ρ_η(h log q/L)/h, nonzero only for T^{η/h} < q < T^{b/h}. ✓ (PROOF.md (i) and review 2.)

**Z2 (inversion: squarefree n → power sums; exact, finite).** Let 𝒫 := {primes q : ρ_η(log q/L) ≠ 0} ⊂ (T^η, T^b), finite. In the monoid algebra ℚ[ℕ, ·] put X_q := ρ_η(log q/L)q^{−s}[q] and let Λ_f be the linear functional [n] ↦ f(log n/L). Then M(s) = Σ_{j≥0}(−1)^jΛ_f(e_j(X)), e_j the elementary symmetric polynomial in (X_q)_{q∈𝒫}. Newton's identity in the commutative ℚ-algebra: e_j = Σ_{λ⊢j}ε_λz_λ^{−1}∏_ip_{λ_i}, ε_λ = (−1)^{j−ℓ(λ)}, z_λ = ∏_h h^{m_h}m_h!, p_h = Σ_qX_q^h = Σ_qρ_η(log q/L)^hq^{−hs}[q^h]. Applying Λ_f:
  M(s) = Σ_{j≥0}Σ_{λ⊢j}(−1)^{ℓ(λ)}z_λ^{−1} Σ_{q_1..q_ℓ∈𝒫} f(Σ_iλ_i log q_i/L) ∏_iρ_η(log q_i/L)^{λ_i}q_i^{−λ_is}.
Each inner sum vanishes unless jηL < Σλ_i log q_i < θL, so only j ≤ k occur. Σ_{λ⊢j}z_λ^{−1} = 1, so |coefficients| ≤ 1 and there are ≤ Σ_{j≤k}p(j) =: C₁(k) terms. Parts h ≥ 2 are the **collision packets** q^h with coefficient ρ_η(log q/L)^hq^{−hs} (q > T^η); parts 1 are free primes. (This is PROOF.md's set-partition Möbius inversion grouped by block type: (1/j!)Σ_{π of type λ}μ(0̂, π) = (−1)^{j−ℓ}/z_λ.) Write a term with m free primes and packet product 𝔞 = ∏_{λ_i≥2}q_i^{λ_i} as c_𝔞𝔞^{−s}·F_m(𝔞; s), |c_𝔞| ≤ 1, where
  F_m(𝔞; s) := Σ_{q'_1..q'_m∈𝒫} f_𝔞(Σ log q'_i/L) ∏ρ_η(log q'_i/L)q'_i^{−s},  f_𝔞(v) := f(v + log 𝔞/L).
Packets are summed with independent primes (two packets may share a prime): the bound in Z5 treats them as independent variables, which is what the expansion gives. ✓

**Z3 (free sums → zeros, cutoff inside).** Measures on ℝ: ν_s := Σ_nΛ(n)φ(log n)n^{−s}δ_{log n}, π_s := Σ_{q∈𝒫}ρ_η(log q/L)q^{−s}δ_{log q}, ν^pp_s := Σ_{q, h≥2}(ρ_η(h log q/L)/h)q^{−hs}δ_{h log q}; π_s = ν_s − ν^pp_s (Z1's last line). F_m(𝔞; s) = ∫f_𝔞(Σu_i/L)dπ_s^{⊗m} and (ν − ν^pp)^{⊗m} = Σ_{S⊆[m]}(−1)^{|S|}ν^{⊗([m]∖S)}⊗(ν^pp)^{⊗S}. The ν^pp factors are point masses: they become **prime-power packets** q^h (h ≥ 2, T^{η/h} < q < T^{b/h}) with coefficient (ρ_η(h log q/L)/h)q^{−hs}, absorbed into 𝔟 := 𝔞·∏_{i∈S}q_i^{h_i}. Hence F_m(𝔞; s) = Σ_S(−1)^{|S|}Σ_{pp packets}(∏coeff)·G_{m−|S|}(𝔟; s) with
  G_r(𝔟; s) := ∫f_𝔟(Σ_{i≤r}u_i/L)dν_s^{⊗r}(u),  G_0(𝔟) = f(log 𝔟/L).
Claim: G_r(𝔟; s) = Σ_{z_1..z_r∈𝒵_s}∏ε(z_i)·I_𝔟(z),  I_𝔟(z) := ∫_{ℝ^r}f_𝔟(Σu_i/L)∏φ(u_i)e^{z_iu_i}du,
𝒵_s := {1 − s} ∪ {ρ − s : ρ, multiplicity} ∪ {−2m − s : m ≥ 1}, ε = +1, −1, −1 respectively.
Proof. f_𝔟(v) = ∫f̂(ξ)e^{2πiξ(v + log𝔟/L)}dξ, f̂ Schwartz, |f̂_𝔟| = |f̂| (uniform in 𝔟). ν_s is finite, so G_r = ∫f̂_𝔟(ξ)Q(s − 2πiξ/L)^rdξ. Z1 at s' = s − 2πiξ/L: 𝒵_{s'} = 𝒵_s + 2πiξ/L, so Q(s')^r = Σ_{z∈𝒵_s^r}∏ε(z_i)Φ(z_i + 2πiξ/L), absolutely, with Σ_{z∈𝒵_s}|Φ(z + iτ)| ≤ C_T log(2 + |t| + |τ|) (Φ bound with a = 2, zero counting). This majorant to the power r is integrable against |f̂|, so sum and ξ-integral commute (dominated convergence). For each z: ∫f̂_𝔟(ξ)∏_iΦ(z_i + 2πiξ/L)dξ = ∫∏φ(u_i)e^{z_iu_i}[∫f̂_𝔟(ξ)e^{2πiξΣu_i/L}dξ]du = I_𝔟(z) (Fubini: compact u-support, f̂ ∈ L¹). ∎ ✓ The cutoff stays inside the r-fold integral: this is the whole point. (Numerical check: §6, N-A.)

**Z4 (the decisive bound: X^{max d}, not X^{r·max d}).** Let F(u) := f_𝔟(Σu_i/L)∏_{i≤r}φ(u_i) ∈ C_c^∞(ℝ^r). supp F ⊂ K_𝔟 := {u_i ∈ [ηL, bL], Σu_i ≤ log(X/𝔟)} (f_𝔟(v) = 0 for v ≥ θ − log𝔟/L). If 𝔟 ≥ X then K_𝔟 = ∅ and I_𝔟 = 0; so below 1 ≤ 𝔟 < X. On K_𝔟, since u_i > 0:
  Re Σz_iu_i = ΣRe(z_i)u_i ≤ Σ(Re z_i)⁺u_i ≤ M_z·Σu_i ≤ M_z log(X/𝔟),  M_z := max_i(Re z_i)⁺,
so |e^{Σz_iu_i}| ≤ (X/𝔟)^{M_z} on K_𝔟. Let α_i := A if L|z_i| ≥ 1, else 0. e^{z_iu_i} = z_i^{−A}∂_{u_i}^Ae^{z_iu_i}; IBP (no boundary terms, F compactly supported) gives I_𝔟(z) = ∏_{α_i=A}(−z_i)^{−A}∫(∂^αF)e^{Σz_iu_i}du. By Leibniz, (Z0) and ∂^a_v f_𝔟(v/L) = L^{−a}f^{(a)}: |∂^αF| ≤ C₂L^{−|α|−r}, C₂ = C₂(r, A, ‖f^{(≤rA)}‖_∞, C_{≤A}). vol K_𝔟 ≤ (bL)^r. Hence
  |I_𝔟(z)| ≤ C₃(X/𝔟)^{M_z}∏_iκ_A(L|z_i|),  C₃ := C₂b^r  (Z4)
independent of 𝔟, s, T. Since κ_A is nonincreasing and |z| ≥ |Im z|: κ_A(L|ρ − s|) ≤ κ_A(L|γ − t|); for z = 1 − s or −2m − s, |z| ≥ t ≥ T/2. ✓ This is the only place the length constraint enters: the r-fold integrand lives on Σu_i ≤ θL, so the exponential costs X^{M_z} once. The route |e_r| ≪ (1 + |S₁|)^r (DPC l.665–670) pays (X^{d})^r. (N1 log: I₂/e^{θa}, I₃/e^{θa} bounded, I₁^r not.)

**Z5 (summation).** κ_A-sums: (a) for c_i ≥ 0, X^{max c_i} ≤ Σ_iX^{c_i}; so for r ≥ 1, Σ_{ρ_1..ρ_r}X^{max_ic(ρ_i)}∏κ_A(L|γ_i − t|) ≤ r·(Σ_ρX^{c(ρ)}κ_A)·On^{r−1}. (b) Non-zero entries: Σ_{z∈{1−s}∪{−2m−s}}κ_A(L|z|) ≤ κ_A(LT/2) + Σ_{m≥1}κ_A(L max(T/2, 2m)) ≤ 2T(LT/2)^{−A} =: τ_T.
Case + (s = s₊, Re(ρ − s) = d_ρ − y, Re(1 − s) = ½ − y ≤ 1, Re(−2m − s) < 0). Tuples of zeros only: by (Z4) with (X/𝔟)^{c} ≤ X^{c} (c ≥ 0, 𝔟 ≥ 1) and (a): Σ|I_𝔟| ≤ C₃rW_{−y}On^{r−1}. Tuples with j ≥ 1 non-zero entries: M_z ≤ 1, so ≤ C₃X·C(r,j)τ_T^j(1 + W_{−y})(1 + On)^{r−1} (zero entries summed by (a) or by On). Since Xτ_T ≤ 2^{A+1}T^{1+θ−A} ≤ 1 for A ≥ 3: |G_r(𝔟; s₊)| ≤ C₄(1 + W_{−y})(1 + On)^{r−1} (r ≥ 1), |G₀| ≤ ‖f‖_∞. Packet moduli: |𝔟^{−s₊}| = 𝔟^{−½−y} ≤ 𝔟^{−½}. With Z2, Z3:
  |M(s₊)| ≤ C₁(k)·2^k·max_rC₄·(1 + W_{−y})(1 + On)^{k−1}·Π,  Π := Σ_{packet tuples}∏|coefficient|·𝔟^{−½} ≤ B₀^k,
where each packet (≤ k of them) contributes at most B₀ := Σ_{T^η<q<T^b}q^{−1} + Σ_{q, h≥3}q^{−h/2} + ½Σ_{T^{η/2}<q<T^{b/2}}q^{−1} ≤ 2log(b/η) + 6 (Mertens, T ≥ T₀; Σ_qΣ_{h≥3}q^{−h/2} ≤ Σ_q 2q^{−3/2} < 4). This proves (Z+).
Case − (s = s₋, Re(ρ − s) = d_ρ + y). (d + y)⁺ ≤ d⁺ + y (y ≥ 0), and X/𝔟 ≥ 1 on K_𝔟 ≠ ∅, so for zero tuples (X/𝔟)^{M_z} ≤ (X/𝔟)^{y + max d_i⁺}. Packet moduli |𝔟^{−s₋}| = 𝔟^{−½+y}. Product: 𝔟^{−½+y}(X/𝔟)^{y+max d⁺} = X^{y+max d⁺}𝔟^{−½−max d⁺} ≤ X^y𝔟^{−½}X^{max d⁺}. Then (a) with c = d⁺ gives C₃X^y𝔟^{−½}rWOn^{r−1}. r = 0: f(log𝔟/L) ≠ 0 forces 𝔟 < X, so 𝔟^{−½+y} ≤ X^y𝔟^{−½}. Non-zero entries: M_z ≤ ½ + y ≤ 1, same τ_T absorption. Packet sum as in Case +: (Z−). ∎

**Z6 (dependence, uniformity).** C = C₁(k)2^k max_{r≤k}C₄(r)·B₀^k depends on f (via ‖f^{(a)}‖_∞, a ≤ kA), ψ, ω, η, θ (via k), b, A; not on t, y, T ≥ T₀(η, b, θ). Used: t ≥ T/2 (only in τ_T), 0 ≤ y ≤ ½ (only M_z ≤ 1 and y ≥ 0). Gamma factor: in Z1 the contour goes to Re w = −∞, so the archimedean factor appears only through the trivial-zero terms, bounded in (b). Lemma Z is **correct as stated**. Remark: A ≥ 3 suffices for Z5, but M1–M4 below need A ≥ 4 + 2/p-type conditions; A = 12 fixed works for every exponent used (see §2), so "A fixed large" can be an absolute constant, chosen before p, u, w (order of choices consistent).

## 2. Lemma M (complete proof)

Inputs. (F) Fujii as quoted Lz l.242–244: ∫_0^T(N(t + ℓ) − N(t))^{2j}dt ≪_j T uniformly for 0 ≤ ℓ ≤ 100/L, N counting all zeros with multiplicity (unconditional; Lz footnote l.246: Fujii's range corrected in Fu2). Hölder (real p ≤ 2j) and Minkowski over ⌈ℓ log(4T)/100⌉ pieces give, applied at height 4T: (F′) ∫_0^{4T}(N(τ + ℓ) − N(τ))^pdτ ≪_p T(1 + ℓL)^p for 0 ≤ ℓ ≤ 1, p ≥ 1. (U) Unconditionally N(τ + H) − N(τ) ≪ (H + 1)log(τ + H + 2) for τ ≥ 0, H ≥ 0. (D) Density, Lz l.224: N(σ, T) ≪_δ T^{1−κ(σ−½)}log T uniformly in σ ≥ ½, κ = 8/7 − δ (§7). Cells: n_m(t) := #{ρ : m ≤ L(γ − t) < m + 1}, m ∈ ℤ. A zero in cell m has L|γ − t| ≥ max(|m| − 1, 0), so κ_A(L|γ − t|) ≤ 3^A(1 + |m|)^{−A} and On(t) ≤ 3^AΣ_m(1 + |m|)^{−A}n_m(t).  (C)

**M1.** ‖On‖_{L^p[T/2,5T/2]} ≪_{p,A} T^{1/p} for every p ≥ 1 (A ≥ 2). For |m| ≤ TL/4 and t ∈ [T/2, 5T/2], n_m(t) = N(t + (m+1)/L) − N(t + m/L) at heights in [T/4, 11T/4], so ∫n_m^p ≤ ∫_{T/4}^{3T}(N(τ + 1/L) − N(τ))^p ≪ T by (F′). For |m| > TL/4: n_m(t) ≤ #{ρ : ||γ| − |t + m/L|| ≤ 1/L} ≪ log(2 + 3T + |m|/L) (unconditional N(τ + 1) − N(τ) ≪ log(2+|τ|), conjugates by symmetry), so ‖n_m‖_p ≪ T^{1/p}log(T + |m|). Minkowski over m in (C). ✓ ("Fujii only at comparable heights": used only on [T/4, 3T].)

**M2.** B_a := {ρ : γ ∈ [T/4, 3T], d_ρL ∈ [a, a+1)}, a ≥ 0. By (D) at σ = ½ + a/L, height 3T: |B_a| ≤ C_δ(3T)^{1−κa/L}log 3T ≤ C₅TLe^{−κa}. F_a := Σ_{ρ∈B_a}κ_A(L|γ − t|) ≤ 3^AΣ_m(1+|m|)^{−A}n^{(a)}_m, n^{(a)}_m(t) := #{ρ ∈ B_a in cell m}. Claim ‖n^{(a)}_m‖_{L^p[T/2,5T/2]} ≤ C T^{1/p}e^{−κa(1−1/w)/p} uniformly in m ∈ ℤ, every w > 1. Proof: n^{(a)}_m is a nonnegative integer ≤ n_m, and n^{(a)}_m(t) ≥ 1 forces t + m/L ∈ [T/4 − 1/L, 3T], so (n^{(a)}_m)^p ≤ 1_{n^{(a)}_m≥1}·1_{t+m/L∈[T/4−1,3T]}n_m^p (the integer trick). Hölder (w, w′): ∫(n^{(a)}_m)^p ≤ meas{n^{(a)}_m ≥ 1}^{1−1/w}(∫_{T/4−1}^{3T}(N(τ + 1/L) − N(τ))^{pw}dτ)^{1/w}, and meas{n^{(a)}_m ≥ 1} ≤ Σ_{ρ∈B_a}meas{t : m ≤ L(γ − t) < m+1} = |B_a|/L ≤ C₅Te^{−κa}; the second factor is ≪ T^{1/w} by (F′). ∎ Minkowski in (C): ‖F_a‖_p ≪_{p,w} T^{1/p}e^{−κa(1−1/w)/p}. ✓

**M3.** Off(t) := Σ_{d_ρ>0}X^{d_ρ}κ_A(L|γ − t|). For t ∈ [T/2, 5T/2]: a zero with d_ρ > 0 and γ ∈ [T/4, 3T] lies in B_a, a = ⌊d_ρL⌋ ≤ L/2, with X^{d_ρ} = e^{θLd_ρ} < e^{θ(a+1)}; the others have |γ − t| ≥ T/4 and X^{d} ≤ X^{½}, total ≤ T^{θ/2}Σ_{|γ−t|≥T/4}(L|γ − t|)^{−A} ≪ T^{θ/2+1−A} ≤ T^{−5} (A ≥ 7). So Off ≤ Σ_ae^{θ(a+1)}F_a + T^{−5}, and by Minkowski
  ‖Off‖_{L^p[T/2,5T/2]} ≪ T^{1/p}e^θΣ_{a≥0}e^{a(θ − κ(1−1/w)/p)} ≪ T^{1/p}  if pθ < κ(1 − 1/w);
given pθ < κ choose w with κ(1 − 1/w) > pθ. ✓ Moment of R_T: (Z+) at y = 0 and W ≤ On + Off (d ≤ 0: weight 1 ≤ On's; d > 0: Off's) give |R_T| ≤ C(1 + Off)(1 + On)^k on [T/2, 3T]. For q < κ/θ pick p₁ > 1 with qp₁θ < κ, p₂ = p₁/(p₁ − 1); Hölder:
  ∫_T^{2T}|R_T|^q ≤ C^q‖1 + Off‖^q_{qp₁}‖1 + On‖^{qk}_{qkp₂} ≪ T^{1/p₁ + 1/p₂} = T. ✓ κ = 8/7 − δ, θ < 4/7: κ/θ > 2 for δ small, so some q > 2 is admissible.

**M4 (discrete, over cut points).** Points x_ρ = γ_ρ ∓ h/2, γ_ρ ∈ [T, 2T], both signs, multiplicity (≤ 2N(2T) points); Σ_ρ is over these. Ñ^{(m)}_ρ := #{ρ′ : |γ′ − γ_ρ| ≤ ℓ_m}, ℓ_m := (|m| + 1)/L + h.
(M4a) Σ_{γ_ρ∈[T/4,3T]}(Ñ^{(m)}_ρ)^p ≪_p TL(|m| + 9)^{2p+1} for every m ∈ ℤ. *Small windows* ℓ := ℓ_m ≤ ½ (Lz l.257–283 with window ℓ): cells I_j of length ℓ covering [T/4 − 1, 3T + 1], n_j := #zeros in I_j; a zero in I_j has Ñ ≤ n_{j−1} + n_j + n_{j+1}, and ab^p ≤ a^{p+1} + b^{p+1} gives Σ ≤ C_pΣ_jn_j^{p+1}. For τ ∈ I_{j−1} ⊂ [T/4 − 2, 3T + 1], I_j ⊂ (τ, τ + 2ℓ], 2ℓ ≤ 1, so Σ_jn_j^{p+1} ≤ ℓ^{−1}∫_0^{4T}(N(τ + 2ℓ) − N(τ))^{p+1}dτ ≪ L·T(1 + 2ℓL)^{p+1} by (F′) (ℓ ≥ 1/L), and 1 + 2ℓL ≤ 2(|m| + 9). *Large windows* ℓ > ½ (then |m| + 1 + 2πφ > L/2, so L < 2|m| + 16): by (U), Ñ^{(m)}_ρ ≤ #{|γ′| ≤ γ_ρ + ℓ, ||γ′| − γ_ρ| ≤ ℓ} ≪ (ℓ + 1)log(4T + ℓ) ≪ (|m| + 9)²·(1 + log(1 + |m|)) (ℓL ≤ |m| + 9), and there are ≪ TL points: Σ ≪ TL(|m| + 9)^{2p+1}. ✓ (Repair of a gap found by the reviewer in the first version of this write-up: the cell argument was applied to unbounded ℓ_m.)
(M4b) Σ_ρOn(x_ρ)^p ≪ TL: a zero in cell m of x_ρ has |γ′ − x_ρ| ≤ (|m|+1)/L, so |γ′ − γ_ρ| ≤ ℓ_m and n_m(x_ρ) ≤ Ñ^{(m)}_ρ. Minkowski in ℓ^p over points with (C): ‖On(x_·)‖_{ℓ^p} ≪ Σ_m(1+|m|)^{−A}(TL)^{1/p}(|m| + 9)^{2+1/p} ≪ (TL)^{1/p} (A ≥ 5). ✓
(M4c) Σ_ρOff(x_ρ)^p ≪ TL whenever pθ < κ. If n^{(a)}_m(x_ρ) ≥ 1, some ρ′ ∈ B_a has |γ′ − γ_ρ| ≤ ℓ_m, i.e. ρ is counted by Ñ^{(m)}_{ρ′}. Hence (factor 2: two signs; Hölder; (M4a))
  #{points : n^{(a)}_m(x_ρ) ≥ 1} ≤ 2Σ_{ρ′∈B_a}Ñ^{(m)}_{ρ′} ≤ 2|B_a|^{1−1/w}(Σ_{ρ′}(Ñ^{(m)}_{ρ′})^w)^{1/w} ≪ TL(|m| + 9)^{3}e^{−κa(1−1/w)}.
Integer trick and Hölder again: Σ_ρn^{(a)}_m(x_ρ)^p ≤ #{…}^{1−1/w}(Σ_ρ(Ñ^{(m)}_ρ)^{pw})^{1/w} ≪ TL(|m| + 9)^{2p+4}e^{−κa(1−1/w)²}. With Off(x_ρ) ≤ Σ_ae^{θ(a+1)}F_a(x_ρ) + T^{−5} (x_ρ ∈ [T/2, 5T/2]) and Minkowski in ℓ^p:
  ‖Off(x_·)‖_{ℓ^p} ≪ (TL)^{1/p}Σ_ae^{θ(a+1)−κa(1−1/w)²/p}·Σ_m(1+|m|)^{−A}(|m|+9)^{2+4/p} ≪ (TL)^{1/p}  if pθ < κ(1−1/w)² (A ≥ 8). ✓
PROOF.md's M4 bound "(|m|+2)TLe^{−κa(1−1/w)}" has a higher power of |m| + 9 with my constants; immaterial (polynomial in |m|, absorbed by κ_A). **Cosmetic C2.** Its referee-B repair (radius (|m|+1)/L for m ∈ ℤ) is right: (m+1)/L < 0 for m ≤ −2.
**Constraint on A.** Z5: A ≥ 3; M1: A ≥ 2; M3 tails: A ≥ 7; M4: A ≥ 8. Every exponent p ≥ 1 in M1–M4 enters only through (|m| + 9)^{2+4/p} ≤ (|m|+9)^6. So A := 12 works for all p, u, w at once; it can be fixed before them, as the Order of choices does. ✓ **Lemma M is correct.**

## 3. Step 6e′ (complete proof)

Setting (U 6a–6d, audited in UAudit §2–§5, unchanged): S_a, Y_a = (a+1)/L, σ_a = min(2ε + o(1), C₅e^{−κa}) (Step 2 + (D)), points x_ρ^± = γ_ρ ∓ h/2, and (U 6a, UAudit l.117–119)
  𝓗_{T,a} ≤ 4C_P Σ_{±}Σ_{ρ∈S_a}∫_0^{d_ρ}|R_TR*|(x_ρ − iy)·Z(x_ρ − iy)dy,  Z := (1 + |V_L| + 2Ñ_ρ)^D,
with |G_T(x − iy)| ≤ 2 used. Note a ≤ L/2 + 1 (d_ρ < ½), so Y_a < 1 and U 6b's box bound |V̂(x − iy)| ≤ C(P_y*v)(x) + C applies on every segment {x_ρ − iy : 0 < y < Y_a}, x_ρ ∈ [T − 1, 2T + 1].
**(a)** Zeros with d ≤ y contribute X⁰ = 1 to W_{−y}, zeros with d > y ≥ 0 contribute X^{d−y} = X^{−y}·X^d (they are in Off): W_{−y} ≤ On + X^{−y}Off, W ≤ On + Off. (Z±) at t = x_ρ:
  |R_TR*|(x_ρ − iy) ≤ C²X^y(1 + On + X^{−y}Off)(1 + On + Off)(1 + On)^{2k−2} = C²[X^y(1 + On) + Off](1 + On + Off)(1 + On)^{2k−2}
  ≤ C²(1 + Off)(1 + On)^{2k}[X^y + Off],  all of On, Off at x_ρ (independent of y). ✓
**(b)** r′ > 1, 1/r + 1/r′ = 1, Z̄_ρ := (Y_a^{−1}∫_0^{Y_a}Z^{r′}dy)^{1/r′}. Hölder and ∫_0^dX^{ry}dy ≤ X^{rd}/(rθL), d_ρ < Y_a, X^{d_ρ} < e^{θ(a+1)}:
  ∫_0^{d_ρ}X^yZ ≤ X^{d_ρ}(rθL)^{−1/r}(Y_aZ̄^{r′}_ρ)^{1/r′} ≤ (rθ)^{−1/r}(a+1)^{1/r′}e^{θ(a+1)}L^{−1}Z̄_ρ;  ∫_0^{d_ρ}Z ≤ Y_aZ̄_ρ = (a+1)L^{−1}Z̄_ρ.
With E := e^{θ(a+1)} ≥ 1: (1 + Off)(E + Off) ≤ E² + Off + E·Off + Off² ≤ 3(E² + Off²). So
  ∫_0^{d_ρ}|R_TR*|Z dy ≤ C₆(a+2)L^{−1}[e^{2θ(a+1)} + Off(x_ρ)²]·Y_ρ,  Y_ρ := (1 + On(x_ρ))^{2k}Z̄_ρ. ✓ (PROOF.md has (a+2)²; harmless.)
**(c)** All-zero moments (points x_ρ^±, γ_ρ ∈ [T, 2T], segment (0, Y_a); a fixed). Claim Σ_ρY_ρ^w ≪_w (a+2)^{C}TL. Cauchy–Schwarz: Σ_ρY^w ≤ (Σ_ρ(1 + On(x_ρ))^{4kw})^{½}(Σ_ρZ̄_ρ^{2w})^{½}; first factor ≪ (TL)^{½} (M4b). For 2w ≥ r′, power means give Z̄^{2w} ≤ Y_a^{−1}∫_0^{Y_a}Z^{2w}dy, and with V_L = V̂ + Δ_L (U 6b), Z^{2w} ≤ C[1 + |V̂|^{2wD} + |Δ_L|^{2wD} + Ñ_ρ^{2wD}].
 – V̂: |V̂|^{2wD} ≤ C[(P_y*v)^{2wD} + 1] ≤ C[P_y*(v^{2wD}) + 1] (Jensen; P_y a probability density). Fubini: Σ_ρY_a^{−1}∫_0^{Y_a}(P_y*v^{2wD})(x_ρ)dy = ∫v^{2wD}K, K := Σ_ρY_a^{−1}k_{Y_a}(· − x_ρ), k_Y(u) = ∫_0^YP_y(u)dy = (2π)^{−1}log(1 + Y²/u²), ∫k_Y^q = c_qY. Near part (|u − x_ρ| ≤ Y_a): power mean, M(u) ≤ N_ρ ≤ Ñ_ρ, so ∫(Σ_near k)^q ≤ c_qY_aΣ_ρÑ_ρ^{q−1} ≪ Y_a(a+2)^{q}TL (M4a with ℓ = 4Y_a + h). Far part: k_Y ≤ Y²/(2πu²), dyadic shells, Σ_far ≤ Σ_l4^{−l}n(u, R_l), R_l := 2^{l+1}Y_a, n(u, R) := #{points with |x_ρ − u| ≤ R}. For R_l ≤ ½: ‖n(·, R_l)‖_q ≪ T^{1/q}(1 + R_lL) (F′; support [T − 2, 2T + 2]). For R_l > ½ (large radius, by (U)): n(u, R) ≪ (R + 1)L on its support [T − R − 1, 2T + R + 1], so ‖n(·, R_l)‖_q ≪ (R_l + 1)L(T + 2R_l + 2)^{1/q}; with 2^{−l} < 4Y_a there, Σ_{R_l>½}4^{−l}2^{l+2}L(T^{1/q} + 2^{(l+3)/q}) ≪ (a+1)T^{1/q} + L. (Large-radius case written out after the review; it needs q > 1, which holds: the norm used is q = r″ = r/(r − 1) ∈ (1, ∞). Re-review: M4a CORRECT, far part CORRECT for 1 < q < ∞.) With K₀ := Y_aK: ‖K₀,near‖_q ≪ (a+2)((a+1)T)^{1/q}, ‖K₀,far‖_q ≪ Σ_l4^{−l}T^{1/q}(1 + 2^{l+1}(a+1)) ≪ (a+2)T^{1/q}; so ‖K‖_q ≪ Y_a^{−1}(a+2)²T^{1/q} ≪ (a+2)LT^{1/q}. Hölder (r, r″) on I: ∫_Iv^{2wD}K ≤ ‖v^{2wD}‖_{L^r(I)}‖K‖_{r″} ≪ T^{1/r}(a+2)^CLT^{1/r″} = (a+2)^CTL (U 5c: ∫_Iv^m ≪ T). ✓
 – Δ_L: U 6d on (0, Y_a) instead of (0, d_ρ) (its proof uses only that the segment lies in 0 < y < 1): Y_a^{−1}∫_0^{Y_a}|Δ_L|^m ≤ C_m(1 + h/Y_a)(1 + Ñ_ρ + Ψ_ρ)^m, h/Y_a ≤ 2π. Σ_ρΨ_ρ^p ≪ TL (U 6d, all zeros, before its restriction to S_a), Σ_ρÑ_ρ^p ≪ (a+2)^{p+1}TL (M4a). ✓ The near ρ′ of 6d (|γ′ − γ_ρ| < 5h/2) are counted by Ψ_ρ (each contributes ≥ (1 + 6.25)^{−1}), so the window issue of UAudit (d2) is harmless here too. **Cosmetic C4** (= UAudit d2).
 So Σ_ρZ̄^{2w} ≪ (a+2)^CTL and Σ_ρY_ρ^w ≪ (a+2)^CTL. ✓ And Σ_ρOff(x_ρ)^{2u} ≪ TL for 2uθ < κ (M4c). ✓
**(d)** Fix u ∈ (1, κ/(2θ)) (nonempty iff θ < κ/2), w > u/(u − 1), 1/v := 1 − 1/u − 1/w > 0. |S_a| ≤ 2σ_aTL points. Hölder over S_a with the all-zero sums of (c):
  Σ_{S_a}Y_ρ ≤ |S_a|^{1−1/w}(Σ_ρY^w)^{1/w} ≪ (a+2)^Cσ_a^{1−1/w}TL;  Σ_{S_a}Off²Y ≤ |S_a|^{1/v}(Σ_ρOff^{2u})^{1/u}(Σ_ρY^w)^{1/w} ≪ (a+2)^Cσ_a^{1/v}TL.
With (b): 𝓗_{T,a} ≪ L^{−1}(a+2)[e^{2θ(a+1)}Σ_{S_a}Y + Σ_{S_a}Off²Y] ≪ T(a+2)^C[e^{2θ(a+1)}σ_a^{1−1/w} + σ_a^{1/v}]. ✓
**(e)** Φ(ε) := Σ_{a≥0}(a+2)^C[e^{2θ(a+1)}σ_a^{1−1/w} + σ_a^{1/v}], σ_a := min(3ε, C₅e^{−κa}) (T ≥ T₁(ε)). Choose w also with κ(1 − 1/w) > 2θ. A_ε := c log(1/ε), 0 < c < (1 − 1/w)/(2θ): Σ_{a≤A_ε} ≤ (A_ε + 1)(A_ε + 2)^Ce^{2θ}ε^{−2θc}(3ε)^{1−1/w} → 0; Σ_{a>A_ε} ≤ C Σ_{a>A_ε}(a+2)^Ce^{2θ(a+1)−κ(1−1/w)a} → 0 (tail of a convergent series); Σ_a(a+2)^Cσ_a^{1/v} → 0 by dominated convergence (majorant (a+2)^C(C₅e^{−κa})^{1/v}). Constants: C depends on φ, θ, f, P, η, κ, δ, A, u, w, r, r′ (not on ε). Zeros with γ ∉ [T, 2T]: G_T ≤ e^{−L²/32} on their slits, |R_TR*| ≤ (Σ|c_n|)(Σ|c_n|n^{y}) ≤ 4‖f‖²_∞X^{3/2}, V-factor T^{O(1)} up to integrable log singularities (U Step 3): o(1). **6e′ is correct.**

## 4. Step 7′ and the order of choices

Re∫_ℝG_T|R_T|²P(V) ≤ −δ₁T for ε ≤ ε₁ (U Step 7 with 6e′). [T, 2T] ∖ E_T has Re P(V) ≥ 0 (Im V ∈ {0, 1}, Re P(x) = q₀(x)², Re P(x + i) = 0); ℝ ∖ [T, 2T] contributes O(e^{−cL²}T^{O(1)}) (G_T, |R_T| ≤ 2‖f‖X^{½}, ∫|V|^m polynomial). So for large T, δ₁T/2 ≤ ∫_{E_T}G_T|R_T|²|P(V)|. Fix 2 < q < κ/θ (M3) and s with 2/q + 1/s < 1. Hölder (q/2, s, (1 − 2/q − 1/s)^{−1}), G_T ≤ 1:
  δ₁T/2 ≤ ‖R_T‖²_{L^q[T,2T]}‖P(V)‖_{L^s[T,2T]}meas(E_T)^{1−2/q−1/s} ≤ C₇T^{2/q+1/s}meas(E_T)^{1−2/q−1/s}
(M3; |P(V)| ≤ C(1 + |V|)^D and U 5c ∫_T^{2T}|V|^m ≪ T), so meas(E_T) ≥ c₂T, c₂ := (δ₁/(2C₇))^{1/(1−2/q−1/s)}, independent of ε. Rest of U Step 7 unchanged (𝒩(T) ≥ c₂T/h − O(L), ε₀ = min(ε₁, c₂/(4πφ))). ✓
Order of choices (checked for circularity): φ ∈ (max{2/(2+κ), ½}, 1) → θ ∈ (1/φ − 1, min{κ/2, 1}) → Step-0 polynomial q₀, b₀, P (DPC l.744–765) → f, η, δ₀ (DPC Lemma 6) → A = 12 → q ∈ (2, κ/θ), s → u ∈ (1, κ/(2θ)), w (> u/(u−1), κ(1−1/w) > 2θ, and the M-exponents' w), v, r, r′, r″ → all constants C → ε₁ (Φ(ε₁)T ≤ δ₀·lim_T∫G_T|R_T|²/2), c₂, ε₀ → T → ∞. No constant chosen later feeds back. ✓ PROOF.md uses the letter q both for the Step-0 polynomial and the moment exponent. **Cosmetic C3.** Range of φ: θ exists iff 1/φ − 1 < min{κ/2, 1} iff φ > max{2/(2+κ), ½}; κ = 8/7 − δ gives 2/(2+κ) ↓ 7/11 as δ ↓ 0, so every φ > 7/11 is reached. ✓ D(λ) > 0 for λ > 7/11: 7/11 < φ < min(λ, 1) and U′ on [T/2, T] (as U Thm). ✓

## 5. θ-check: U's unchanged steps at θ ∈ [½, 4/7) (indeed at every θ < 1)

| Step of U | Where θ could enter | Check at θ ∈ [½, 4/7) |
|---|---|---|
| 0 | DPC construction | DPC l.744 needs 1/φ − 1 < θ < 1; l.753 sign from 0 < φθ < 1 < φ(1+θ) < 2: φ < 1, θ < 1 ✓. Lemma 6 (l.673 "0 < θ < 1"; l.702 "T^θ < T/L^{20}"; l.704–741 Mertens, the measure M_η, Laplace transforms; l.718 coprimality costs O(Σ_{q>T^η}q^{−2})) ✓ |
| 1, 2 | none (V, zeros only) | ✓ |
| 3 | sides, far slits, post-shift main term | Prop 1 setup (l.274): H = T/L, B = T/L^{20}, index sets ≤ B elements, q′/q ≤ B; integers n < T^θ satisfy this for θ < 1. sides/far slits: |R_TR*| ≤ 4‖f‖²X^{3/2} = T^{O(1)} against e^{−cL²} ✓. Main term: DPC Prop 1 l.343–388: on Im t = −1 the off-diagonal ω = log(kq/q′) ≠ 0 has |ω| ≥ 1/(2B), shifting gG^me^{−itω} to Im t = ∓H gives e^{−H|ω|} ≤ e^{−L^{19}/4}q′/(kq) (l.368–373) and total O(H‖c‖²B²(e^{−L^{19}/4} + e^{−L²/4})) (l.376–377); this needs only B = T^θ ≤ T/L^{20}, i.e. θ < 1. Diagonal ω = 0: ‖|K|‖^r‖c‖² (l.380–386), no θ. No RH after the shift (UAudit §6 point (ii), re-read here) ✓ |
| 4 | heuristic | not used ✓ |
| 5b, 5c | V and zeros only | ✓ (5c's ‖P(V)‖_s ≪ T^{1/s} uses density with any κ > 0) |
| 6a | jump bound, G_T | ✓ |
| 6b | Poisson representations of R_T, R* | exact for every θ; U′ uses only the V̂, Δ parts ✓ |
| 6c | sparse kernels; R_T-tails O(T^{3θ/2}L) (θ < 2/3) and R*_Y-tails | kernels: no θ ✓; the tails belonged to A₁, A₂, which U′ does not use ✓ |
| 6d | Δ, Ψ | ✓ |
| 7 (rest) | E_T, counting | ✓ |
Text search of U PROOF.md: θ < ½ occurs only via 5a (MV: "R_T² has length T^{2θ} < T", l.49) and its uses in 6e (A₁, A₂: ‖R_T‖₄, ‖R*_{Y_a}‖₄) and 7 (Hölder (2, 4, 4)), plus 6f's "θ < ½ and E ↑ ½". All four are replaced in U′. ✓ **No other θ < ½ use.**

## 6. Numerical falsification (code and logs in this folder; every run under memcap 290 MB, single thread)

**N-A, the identity behind Lemma Z (na_identity.py → logs/na_identity.{log,json}).** Real zeta zeros (ExploreHelio zeta_zeros.npy, 22 491 zeros to height 2·10⁴). Three settings (T, θ, η, b) = (10⁴, 0.56, 0.10, 0.66), (1.8·10⁴, 0.60, 0.10, 0.70), (1.8·10⁴, 0.52, 0.12, 0.60): 68–156 primes in 𝒫, up to 5 prime factors, prime-power packets 4, 8, 9, 16, 25, 27, 32, 49, … (support down to T^{η/h}); complex f; s = ½ + 0, ½ + 0.25, ½ − 0.25, ½ − 0.5 (+ it), t ∈ {1.02T, 1.37T, 19650.3}.
| Check | max error (|M| ∈ [0.15, 2.3]) |
|---|---|
| Z2 Newton power sums vs direct M(s) | 3.1·10⁻¹³ |
| Z3 free primes → (ν − ν^pp)^{⊗m} with packets ρ_η(h log q/L)/h, G_r direct | 4.9·10⁻¹² |
| Z3 with G_r from the zeros (cutoff inside), zeros |γ − t| ≤ 150 | 1.0·10⁻⁷ (≤ 6·10⁻⁵ at window 60: truncation, not bookkeeping) |
| control: t = 24 660 (no zeros in the table) | errors 0.07–1.8 = O(|M|): the zero side is load-bearing |
The packet bookkeeping that reviews 1–2 repaired is now verified to 10⁻¹² in both directions (s₊ and s₋, so for R_T and R*). The agent's n3 had checked only the Λ-weighted j ≤ 2 sums, not the inversion.

**N-B, the bound on real zeros (nb_ratio.py → logs/nb_ratio.*).** There W = W_{−y} = On. sup_t over [T, 2T] (step 0.05) of |M(s₊)|/(1 + On)^k and |M(s₋)|/(X^y(1 + On)^k), T = 10³, 3·10³, 9.9·10³, (θ, η) = (0.56, 0.10), (0.60, 0.15), y = 0, ¼, ½: 0.42–2.7, flat in T (2.34 → 2.51 → 2.68; 1.90 → 2.05 → 2.17). Weak test: the polynomials have only 26–76 terms and the trivial bound Σ|c_n|n^{y−½} is within a factor 1.2–1.6 of the sup, so these heights cannot separate Lemma Z from the trivial bound. No contradiction.

**N-C, planted off-line zeros in the explicit-formula toy (nc_toy.py → logs/nc_toy.{log,json}).** The toy replaces the prime measure by its zero side, ν_s(u) = −φ(u)Σ_{ρ∈𝒵}e^{(ρ−s)u}, and sets M_toy := Σ_{j≤k}((−1)^j/j!)∫f(v/L)ν_s^{*j}(v)dv (Z1–Z3 without primes, so Lemma Z's proof applies verbatim). Background: real zeta zeros near 15 000 rescaled to density L/2π; planted pairs ½ ± d + iγ_p. θ = 0.56, η = 0.1, b = 0.66, k = 5; L = 20, 30, 44; depth a = dL ∈ {2, 5, 8, 12, 0.48L} (a = 12 at L = 20 is outside the strip, kept as a stress case); configurations: single, single at 0.7/L, cluster of 3 (depths a, 0.8a, 0.6a within 0.4/L), 3 coincident pairs (same ordinate, same depth), cluster of 5, two pairs h apart (coincident cuts). Scan: 13 points t within 3/L, y ∈ {0, d/2, d, ½}, s₊ and s₋. 90 configurations, peak 33 MB.
| a | log X^d | log|M| single | coincident3 − single | cluster5 − single | multiplicative: + 2 log X^d, + 4 log X^d | max ratio A = 2 | A = 4 | A = 12 |
|---|---|---|---|---|---|---|---|---|
| 2 | 1.12 | 1.48 | 0.92 | 1.50 | 2.2, 4.5 | 0.41 | 7.8 | 19 |
| 8 | 4.48 | 2.46 | 1.27 | 1.82 | 9.0, 17.9 | 0.27 | 7.0 | 70 |
| 14.4 | 8.06 | 4.46 | 1.51 | 1.99 | 16.1, 32.3 | 0.20 | 6.5 | 607 |
| 21.1 | 11.83 | 7.07 | 1.62 | 2.02 | 23.7, 47.3 | 0.19 | 4.4 | 4.4·10³ |
(log = natural; ratios = |M_toy|/bound with the constant set to 1, max over all configurations at that a.)
- **Z4 survives:** three coincident pairs of depth d cost e^{0.9…1.6} (a bounded combinatorial factor, cf. 3^k), five clustered pairs e^{1.5…2.0}; if the zeros' contributions multiplied (the route |e_r| ≪ (1 + |S₁|)^r), 3 and 5 zeros would cost an extra X^{2d}, X^{4d} = e^{23.7}, e^{47.3} at a = 21. With A = 2 and 4 the ratio to Lemma Z's bound is ≤ 0.41 and ≤ 8.4 for every configuration, and does not grow with a.
- The toy is exactly scale-invariant in units of 1/L (identical numbers at L = 20, 30, 44 for equal a), so L only limits a ≤ L/2. A single deep zero gives |M| ≈ X^de^{−c√a}, c ≈ 0.95–1.04 for a ≥ 14 (cf. the agent's N1 for I₁), with a flat profile near the zero (a = 8, debug run at L = 20: |M| at L|γ − t| = 3 is still 90 % of its central value). Hence with A = 12 the ratio grows (19 → 4.4·10³) towards its bound C(12): the IBP constant ‖∂^AF‖ of Gevrey-class bumps is astronomically large. This is consistent with Lemma Z (C = C(A) is a constant) but shows the lemma is asymptotic: "A fixed large" carries a huge C(A). Not a counterexample.
- Exponents at the planted zero: |M(s₊)| decreases and |M(s₋)| increases with y at rates ≤ X^{y} (log ratio/(θLd) = −0.36…−0.55 and +0.43…+0.79 at y = d), below the bounds X^{(d−y)⁺} and X^{y+d} of (Z±).

**N-D, the M3 threshold in the toy (nd_moments.py → logs/nd_moments.{log,json}).** One pair of depth a/L (a = 1…21, L = 44) over the background; J_q(a) := excess ∫|M_toy(½ + it)|^q over |t − γ_p| ≤ 6/L. A density profile e^{−κa} gives the q-moment Σ_a e^{−κa}J_q(a). Fitted growth d log J_q/da over a ∈ [14, 21]: 0.78, 0.97, 1.17, 1.56 for q = 2, 2.5, 3, 4, i.e. ≈ 0.69·qθ (qθ = 1.12, 1.40, 1.68, 2.24), rising with a (≈ 0.45·qθ at a ≈ 7) as the e^{−c√a} factor predicts: growth → qθ. Partial sums with κ = 8/7: q ≤ 2.5 converged (1.52, 3.32), q = 3 borderline (7.01 → 7.19), q = 4 diverging (30.9 → 35.1 → 71.8); with κ = 1.6 all converge. So in the toy the moment series diverges exactly when the growth rate passes κ, and the growth rate tends to qθ: the threshold q < κ/θ of M3 is the asymptotic one (heuristic sharpness; at reachable depths the true moments are smaller). No contradiction with M3.

**N-E, a dense random configuration (ne_config.py → logs/ne_config.{log,json}).** Toy of N-C at L = 30. Each of the 716 background zeros within 75 of t₀ becomes an off-line pair of depth (a + U)/L with probability e^{−κa} (a = 1…14, κ = 8/7, total 0.47), giving 354 random pairs (deepest a = 7.3). Three pairs are also forced at a = 8.5, 11.5 and 14.5. Moments are taken over |t − t₀| ≤ 60 at step 1/(4L), peak 24 MB.
- Lemma Z: max |M|/((1 + W)(1 + On)^{k−1}) with A = 2 and constant 1 is 1.70 (isolated configurations in N-C: ≤ 0.41). A bounded constant, no growth mechanism visible.
- Moments per unit length (background in brackets): q = 1: 2.8 (1.1); 2: 60 (1.4); κ/θ = 2.04: 70 (1.5); 2.5: 448 (1.7); 3: 3.7·10³ (2.1). The ±6/L window of the forced a = 14.5 pair carries 10 %, 41 %, 42 %, 50 % and 57 % of these.
- Additivity: N-D's isolated-pair excess, rescaled to L = 30 (×44/30) and interpolated to a = 14.5, predicts 3.05·10³ for q = 2. The window here carries 2.93·10³, within 4 %. So clustering many shallow pairs does not amplify a deep one, and N-D's bin decomposition Σ_a|B_a|J_q(a) is the right model for M3.
- The blow-up at q ≥ 2 comes from the forced pairs, which are 1 per 716 zeros where the density profile allows e^{−κ·14.5} ≈ 6·10⁻⁸. This is the M3 mechanism (deep zeros denser than density permits), not a defect of the proof.

## 7. Inputs (brief item 3)

- **Conrey's density, slope 8/7 − δ** (Lz l.221–226: "N(σ,T) ≪_ε T^{1−(8/7−ε)(σ−1/2)} log T, which holds uniformly for σ ≥ 1/2"; footnote: announced in Conrey's BAMS note without full proof, complete proof in Baluyot's thesis §3.6). Our audit Riemann/Deep/Baluyot/REPORT.md (AUDIT PASSED): Thm 1.2.1 N(σ,T) ≪ T^{1+θ_m(1−2σ)}log T for every θ_m < 4/7, i.e. κ = 2θ_m = 8/7 − δ. Uniformity (REPORT l.47, l.50): constants depend on θ_m, δ, ε, ε₀, never on σ; σ ∈ [½ + 1/log y, ½ + ε₀] by the mollifier; ½ ≤ σ ≤ ½ + 1/log y by N(T) ≪ T log T (thesis p. 142), which is also the brief's case σ − ½ ≪ 1/log T (there T^{1−κ(σ−½)} ≥ e^{−O(1)}T, so the bound is trivial); σ ≥ ½ + ε₀ by Ingham T^{3(1−σ)/(2−σ)}log⁵T, whose exponent lies below the Conrey line by ≥ c(ε₀) > 0 (Ingham − Conrey vanishes at ½, has slope 2θ_m − 4/3 < 0 there and is concave), so the log⁵ is absorbed: **a single log T, uniformly for σ ≥ ½.** Height 3T: (3T)^{1−κ(σ−½)}log 3T ≪ T^{1−κ(σ−½)}log T. Multiplicity: Littlewood's lemma counts zeros with multiplicity. What U′ uses: |B_a| ≤ C₅TLe^{−κa} for 0 ≤ a ≤ L/2 ✓.
- **Fujii** (Lz l.240–248): ∫_0^T(N(t + h) − N(t))^{2k}dt ≪_k T uniformly for 0 ≤ h ≤ 100/log T; unconditional, all zeros; Lz's footnote: the h-range of Fujii's main theorem was corrected in Fu2. Used as (F′) in M1, M2, M4a ✓.
- **Lz Lemma DoubleSumZeros** (l.251–284): cells of length 2π/L, Hölder, ab^q ≤ a^{q+1} + b^{q+1}, Σn_j^{q+1} ≤ (L/2π)∫(N(t + 4π/L) − N(t))^{q+1}, Cauchy–Schwarz and Fujii ✓. M4a is this argument with window ℓ_m ✓.
- **Fallback.** With Jutila's κ = 1 − δ (Lz l.226) the range θ < κ/2 is θ < ½ − δ/2: Theorem U′ gives λ > 2/3, i.e. Theorem U, nothing weaker. With any κ ∈ (1, 2]: λ > 2/(2 + κ).

## 8. Corrections found

| # | Where | What | Effect |
|---|---|---|---|
| C1 | PROOF.md conventions vs U Step 0 | supp ρ_η ⊂ [η, b] (DPC l.673: ω = 0 on (−∞, 1]); U's "> T^{η/2}, ≤ ⌊2θ/η⌋" is weaker but true | none |
| C2 | M4 | power of (|m| + 2) | none (A = 12 absorbs (|m|+9)^6) |
| C3 | Order of choices | q denotes both the Step-0 polynomial and the moment exponent | notation |
| C4 | 6e′(c) via U 6d | window of Ñ_ρ (UAudit d2); the near zeros are also counted by Ψ_ρ | none |
| R1 | this write-up, M4a (first version) | cell argument applied to unbounded windows ℓ_m; reviewer GAP; repaired with (U) for ℓ > ½ | none for U′ |
| R2 | this write-up, 6e′(c) far shells | large radii need q > 1; holds (q = r″) | none |
