# DERIVATION — a prime in (x, x + x^0.856] for every real x ≥ 1, and the conditional gains from the Chourasiya–Simonič density (Ingham4)

Every step is implemented in `cert.py` (Arb ball arithmetic, python-flint 0.9.0, 128 bits). `cert.py` is a copy of `Riemann/Deep/Ingham3/cert.py`; §11 lists exactly what changed. "I3 §n" cites Ingham3's DERIVATION.

Two pipelines: **A0** = Ingham3's pipeline unchanged (triangle weight; KLN's Lemma 4.14 at their RH height H₀ = 3.0610046·10^{10}, Table 1 reproduced). **A1** = A0 with two changes: the trapezoid weight 1_{[0,1/10]} * 1_{[0,9/10]} (§3′) and KLN's Lemma 4.14 evaluated with RH(3·10^{12}) (§1, I4′).

## 0. Results

**Theorem A (main claim; refereed zero counts only).** For every real x ≥ 1 the interval (x, x + x^{0.856}] contains a prime.
- A1 certifies θ = 0.856 and not 0.855 (max bound 1.152 at log x ≈ 190.6 with cells of width 0.004). A0 certifies θ = 0.857 and not 0.856 (1.106). So 0.856 is the least multiple of 0.001 either pipeline certifies.
- Theorem A implies I3's k = 7 result: for x = n^7, (n+1)^7 − n^7 > 7n^6 > x^{0.856}.

**Theorem B (conditional).** Assume the zero-density bound of Chourasiya–Simonič (arXiv 2507.15184 v2, Cor. 1.1; unrefereed) with B₁ raised as in §7 (their b₄ = 0.54 replaced by what Ramaré–Zuniga Alterman's preprint 2603.25961 proves). Then (pipeline A1):
(i) (x, x + x^{0.842}] contains a prime for every real x ≥ 1;
(ii) (n^6, (n+1)^6) contains a prime for every integer n ≥ 1;
(iii) (n^k, (n+1)^k) contains a prime for n ≥ n₀(k): n₀(5) = e^{45.584} ≈ 6.26·10^{19}, n₀(4) = e^{91.18} ≈ 3.96·10^{39}, n₀(3) = e^{468.43} ≈ 2.71·10^{203} (I3, KLN only: e^{57.746}, e^{117.21}, e^{627.01}).

## 1. Inputs

- (I1)–(I3) as I3 §1: ζ(s) ≠ 0 for Re s > T* = B₂ + 2.2·10^{−254} (Bebel, Lean; we use B = 0.8733819666132076 > T*, checked); RH up to T₀ = 3·10^{12} (Platt–Trudgian 2021, verified to 3 000 175 332 800); N(t) ≤ t log t for t ≥ 14 (Hasanalizade–Shen–Wong 2022: 0.1038, 0.2573, 9.3675).
- (I4) Kadiri–Lumley–Ng, JMAA 465 (2018), Lemma 4.14, with I3's implementation (all of b₁ = 0.62, b₂ = 1.048, b₃ = 0.605, b₄ = 0.529, a₁ = 0.63, a₂ = 2.851, 3.006, η = 0.2561882, the 20 parameter rows of Table 1, H = H₀ − 1), evaluated at every σ (log A affine in σ). Pipeline A0: H₀ = 3.0610046·10^{10}; cert.py reproduces Table 1 (A and B within 0.002) in every run.
- (I4′) **Pipeline A1: KLN with RH(3·10^{12}).** KLN state RH(H₀) as a hypothesis (TeX l. 124–128) and fix "unless specified … H₀ = 3.0610046·10^{10}" (l. 622) for their numbers. Our reading of the TeX (every occurrence of H₀, l. 124–1098: Theorem 1.1 at l. 194–196, the parameter conventions l. 626, Lemmas 4.1–4.13, Lemma 4.14 at l. 1547–1580 as cited by I3) is that H₀ enters only through RH(H₀), T ≥ H₀ ≥ H, X = kT ≥ kH₀ ≥ 10^9, log H₀ and log log H₀ in the constants, and one monotonicity condition δ < log H₀(log log H₀ − 1)/2 (l. 1077; asserted per row: δ ≤ 0.314 against 33.9). So Lemma 4.14 holds with any H₀ for which RH(H₀) is known, here 3·10^{12} (Platt–Trudgian). Precedent: Cully-Hugill–Johnston (Int. J. Number Theory 19 (2023), proof of Lemma 7.1, refereed: "Using Platt and Trudgian's verification … up to H₀ = 3·10^{12} … we can recalculate the constants in Lemma 4.14 of [KLN]"; they also replace KLN's (3.13), which we keep) and Johnston–Yang (2204.01980). cert.py (`build_kln_H0`) re-evaluates every Table 1 row with H₀ = 3·10^{12}, H = H₀ − 1 and re-checks all admissibility conditions (as in I3 §1; η, μ − 1 > η₀ through b₆(10^9, ·) < 1, η₀ being the root of b₆(10^9, η₀) = 1); row σ₀ = 0.60 fails U(α,k,H₀) > 1 (0.976) and is dropped (19 rows remain; Setup admits a row in a σ-cell only if s > s_min = 1/2 + d/log H₀ holds in Arb, I3 §4.4). Consistency: our A(0.98) = 15.753 against Cully-Hugill–Johnston's 15.743 (their parameters and their (3.13) differ slightly). At σ = 0.86 the new A is 5.1% below the old.
- (I5) **Prime gaps.** Every gap between consecutive primes below 4·10^{18} is at most 1476 (Oliveira e Silva–Herzog–Pardi, Math. Comp. 83 (2014); maximal gap 1476 after 1425172824437699411).
- (C1, conditional) **Chourasiya–Simonič**, arXiv 2507.15184 v2 (30 Sep 2025), Cor. 1.1 and Table 1: for T ≥ H_RH = 3·10^{12} and σ₁ ≤ σ ≤ σ₂ (rows [n/32, (n+1)/32]),
  N(σ,T) ≤ B₁ T^{3(1−σ)/(2−σ)} (log T)^{(7−5σ)/(2−σ)} + B₂ log²T + B₃ log T, N(σ,T) = #{ρ: β ≥ σ, 0 < γ ≤ T}.
  Rows used (σ₁, σ₂, B₁, B₂, B₃, their δ): (0.50000, 0.53125, 5.069, 9.461, 167.8, 0.012) … (0.84375, 0.87500, 21.77, 6.240, 110.7, 0.192), (0.87500, 0.90625, 24.83, 5.745, 101.9, 0.215); all 13 in `CS_ROWS`.
- (C2, conditional) **Ramaré–Zuniga Alterman**, arXiv 2603.25961 v3 (1 Jul 2026), Thm 1.2: S_ε(X) := Σ_{d,e ≤ X} μ(d)μ(e)/[d,e]^{1+ε} ≤ 0.5501 for 10^{12} ≤ X < 10^{16}, ε ≤ 1/50; ≤ 0.5103 for 10^{16} ≤ X < 10^{33}; ≤ 0.4905 for X ≥ 10^{33} (ε ≤ 1/50). Values read from the sagetex output shipped with the source.

## 2. Weight and criterion for h = x^θ

Let x ≥ e^{42} be real, 0 < θ < 1, δ := log(1 + x^{θ−1}), so xe^δ = x + x^θ. Let Φ ≥ 0 be Lipschitz, supported on [0,1], Φ(0) = Φ(1) = 0, with ‖Φ‖_∞/∫Φ ≤ 2 (triangle: 2; trapezoid §3′: 10/9), and w(u) = Φ(log(u/x)/δ), supported on [x, x + x^θ] and zero at both ends. I3 §2 never uses x = n^k: the explicit formula Σ Λ(n)w(n) = W(1) − Σ_ρ W(ρ) − ∫w(u)du/(u³−u) (absolutely convergent: both weights are convolutions of two indicators, Φ̃(z) = ((e^{az} − 1)/z)((e^{bz} − 1)/z), so |Φ̃(z)| ≪ |z|^{−2} in the strip; no truncation), W(1) ≥ M₀ = δx∫Φ, the trivial-zero term ≤ W(1)/(x(x²−1)) and the prime-power bound
  PP/M₀ ≤ m log(xe^δ)[e^{δ/2}x^{−1/2}/4 + J e^{δ}x^{−2/3}/9 + log J/(δx)], m = 2, J = ⌊log(xe^δ)/log 2⌋,
hold for every real x ≥ 4 and 0 < δ ≤ 1 (J ≥ 2; we use x ≥ e^{42}); m enters only as ‖Φ‖_∞/∫Φ ≤ m. Hence **R := Z + PP/M₀ + 2x^{−3} < 1 implies a prime p with x < p < x + x^θ**, Z = Σ_ρ x^{β−1}|Φ̃(δρ)|/∫Φ.

## 3. Kernel A0 — unchanged (I3 §3)

Triangle Φ = 1_{[0,1/2]} * 1_{[0,1/2]}: G(λ) = (S(λ/4) + 2δ_max min(1/2, 1/λ))², s₁ = 0.2173, δ_max ∈ {10^{−9}, 0.025, 0.15}; moments Mom(p; λ₀) := ∫_{[λ₀,∞)} max(1,λ)^p d(−G) on the λ-grid r = 1 + 1/500 (tail beyond 2m·10^4), stored to λ = 30, p-grid j/200 (j ≤ 300).

## 3′. Kernel A1 — the trapezoid (new)

Φ_a = 1_{[0,a]} * 1_{[0,b]}, a = 1/10, b = 1 − a: Φ_a ≥ 0 is Lipschitz on [0,1], Φ_a(0) = Φ_a(1) = 0, ∫Φ_a = ab, ‖Φ_a‖_∞ = a, so ‖Φ_a‖_∞/∫Φ_a = 1/b ≤ 2. Φ̃_a(z) = ((e^{az} − 1)/z)((e^{bz} − 1)/z). For z = s + iλ, 0 ≤ s ≤ δ_max ≤ 1 and c ∈ {a, b}, I3 §3's two estimates with 1/m replaced by c give
  |(e^{cz} − 1)/z| ≤ c|sinc(cλ/2)| + c(e^{sc} − 1) and ≤ c|sinc(cλ/2)| + (e^{sc} − 1)/|λ|, so (e^y − 1 ≤ 2y on [0,1]) |(e^{cz} − 1)/z| ≤ c[S(c|λ|/2) + 2δ_max min(c, 1/|λ|)].
Hence |Φ̃_a(z)|/∫Φ_a ≤ G_a(|λ|) := (S(aλ/2) + 2δ_max min(a, 1/λ))(S(bλ/2) + 2δ_max min(b, 1/λ)), a product of non-negative non-increasing functions; G_a(0) = (1 + 2δ_max a)(1 + 2δ_max b); for λ ≥ 2π/a both arguments of S are ≥ π and S(y) ≤ 1/y, so G_a(λ) ≤ C₂/λ², C₂ = (2/a + 2δ_max)(2/b + 2δ_max). The moments are computed exactly as in I3 §3 (`TrapKernel`): grid to Λ ≥ (2π/a)·10^4 with tail p∫_Λ^∞ C₂λ^{p−3}dλ = pC₂Λ^{p−2}/(2−p); stored to λ = 2π/a + 1 = 63.9; for λ₀ ≥ max(2π/a, stored range) Mom ≤ C₂λ₀^{p−2}(1 + p/(2−p)); for λ₀ ≤ 1 Mom ≤ G_a(0) + p∫₁^∞G_aλ^{p−1}. In both kernels the stored partial sums are Python floats rounded up (`fup`: f ≥ the ball's upper end, exactly checked), which halves memory, and the grid index for λ₀ is checked exactly (λ_i ≤ λ₀ in Arb; a float-rounded index is moved down if needed). a = 1/2 gives the triangle. a = 1/10 was chosen by a float scan (a ∈ {0.05, 0.07, 0.1, 0.15, 0.2, …, 0.5}; the peak bound at θ = 0.856 is smallest for a ≈ 0.07–0.1). check.py: max over 6·10⁴ random z of (|Φ̃_a(z)|/∫Φ_a)/G_a(Im z) = 0.9999999994 (a = 1/10) and 0.9999999995 (a = 1/2).

## 4. Zero sum — unchanged (I3 §4)

Z ≤ 2[x^{σ₀−1}S + log x ∫_{σ₀}^{B} x^{σ−1}F(σ)dσ], σ₀ = 0.515 (needs only |Φ̃(δρ)|/∫Φ ≤ G(δ|γ|), true as Re(δρ) = δβ ∈ [0, δ_max]); S ≤ δ^{−1}D·Mom(1 + 1/D; 0) (D = log(1/δ), 1/δ ≥ 14); F(σ) from KLN in 144 σ-cells of width 0.0025 (the row with the smallest A at both ends), case (i) λ₀ = δT₀ ≥ 1 / case (ii) λ₀ < 1 as in I3 §4.3, exact σ-integrals of e^{E₀+E₁σ}. Conditional runs also put the CS bound (§7) in every σ-cell and use the smaller of the two (both bound the same F on that cell).

## 5. Cells in log x and the tail (θ)

**5.1 Cells.** δ(L) = log(1 + e^{−(1−θ)L}) is decreasing in L = log x. For x ∈ [e^{L_a}, e^{L_b}], δ ∈ [δ_b, δ_a] with δ_a = δ(L_a), δ_b = δ(L_b). I3 §5.1 (i)–(iv) use only this and x ≥ e^{L_a}: x-powers at L_a (log x·x^{σ−1} and x^{t−1} − x^{s−1} = ∫_s^t log x·x^{σ−1}dσ decrease in x for log x > 1/(1−σ); 1/(1−B) < 7.9 and cert.py asserts L_a(1 − B) > 1 in every cell), kernel, D and λ₀ at δ_b (G non-increasing), δ ≤ δ_a < δ_max (asserted), log(xe^δ) ≤ L_b + δ_a, δx ≥ δ_b e^{L_a}. So R̄(L_a, L_b) bounds R for every real x in the cell. cert.py: `cell_bound(S, dfun, La, Lb, m)`, `dfun = delta_theta(θ)` (θ as the exact rational j/1000). Widths: 0.25 on [42, 150), 0.05 on [150, 165) and [215, 300), 0.02 on [165, 215).

**5.2 Tail lemma.** Let L ≥ L* and y = e^{−(1−θ)L} ≤ y* := e^{−(1−θ)L*}, c* := 1/(1 − y*/2). From log(1+y) ≥ y(1 − y/2): 1/δ ≤ c*/y; so D ∈ [D_low, D_max], D_low = (1−θ)L, D_max = (1−θ)L + log c*, and λ₀ = δT₀ ≤ T₀y* < 1 (asserted): case (ii). With A_max, B_max the largest KLN A (right cell ends, c₁ > 0) and B (left ends) (exact maxima of exact numbers, `exact_max`), p_max = grid⌈(8/3)(1−σ₀) + 4/D_low(L*)⌉, M* := G(0) + p_max∫₁^∞G(λ)λ^{p_max−1}dλ (Kernel.mom's λ₀ ≤ 1 branch; ≥ Mom(p; λ₀) for every λ₀ ∈ [0, 1] and p ≤ p_max, since there Mom(p; λ₀) = G(λ₀) + p∫₁^∞Gλ^{p−1} and G(λ₀) ≤ G(0)), M₁ = Mom(grid⌈1 + 1/D_low(L*)⌉; 0) and Ξ(L) := L − (8/3)D_max − 2 log D_max:
  R ≤ R_an(L) := 2[L·A_max M* D_max³ e^{−(1−B)Ξ}/Ξ + B_max M* D_max² e^{−(1−B)L} + e^{−(1−σ₀)L} c* e^{(1−θ)L} D_max M₁]
              + m(L+1)[e·e^{−L/2}/4 + J·e·e^{−2L/3}/9 + log J·c* e^{−θL}] + 2e^{−3L},  J = (L+1)/log 2.
Steps: e^{(σ−1)L}δ^{−a}D^{5−2σ} ≤ D_max³ e^{−(1−σ)Ξ(L)} (δ^{−a} = e^{(8/3)(1−σ)D}, D^{5−2σ} = D³e^{2(1−σ)log D}, D ≤ D_max); ∫_{σ₀}^{B} e^{−(1−σ)Ξ}dσ ≤ e^{−(1−B)Ξ}/Ξ; every moment exponent a + e/D ≤ (8/3)(1−σ₀) + 4/D and 2/D is ≤ p_max; S ≤ (c*/y)D_max M₁; δ ≤ y ≤ 1 in PP (e^{δ/2}, e^δ ≤ e; log(xe^δ) ≤ L + 1) and 1/(δx) ≤ c*e^{−θL}.
Monotonicity for L ≥ L*: main term d/dL log ≤ 1/L + 3(1−θ)/D_max − (1−B)Ξ′ ≤ 4/L − (1−B)(κ − 2/L) with κ = 1 − 8(1−θ)/3 and Ξ′ = κ − 2(1−θ)/D_max ≥ κ − 2/L > 0; so it is non-increasing when (1−B)(κ − 2/L*) > 4/L* (checked in Arb, with κ − 2/L* > 0, θ − σ₀ > 2/L*, θ > 3/L*); the B-term needs 2/L < 1 − B, the S-term θ − σ₀ > 1/L, the PP terms L ≥ 3 (all far below L* = 300). Hence R(x) ≤ R_an(log x) ≤ R_an(L*) for all x ≥ e^{L*}; cert.py checks R_an(300) < 1 (0.0142 for θ = 0.856 in A1).

**5.3 k-th powers** (conditional runs): I3 §5 unchanged (`delta_pow(k)`: δ = k log(1 + e^{−L/k}); `tail_pow` = I3's tail lemma; cell widths 0.05 (k = 6, 5), 0.1 (k = 4), 0.25 (k = 3)).

## 6. Small x (θ = 0.856, 1 ≤ x < e^{42})

(i) 1 ≤ x < 2: x^θ ≥ 1, so 2 ∈ (x, x + x^θ].
(ii) 2 ≤ x < p_last, p_last = 9999991 the largest prime below 10^7: with p_n ≤ x < p_{n+1} ≤ p_last, x + x^θ ≥ p_n + p_n^θ (t + t^θ increases), so g_n := p_{n+1} − p_n ≤ p_n^θ suffices. A sieve to 10^7 gives the maximal-gap records (g, p) = (1, 2), (2, 3), (4, 7), (6, 23), (8, 89), (14, 113), (18, 523), (20, 887), (22, 1129), (34, 1327), (36, 9551), (44, 15683), (52, 19609), (72, 31397), (86, 155921), (96, 360653), (112, 370261), (114, 492113), (118, 1349533), (132, 1357201), (148, 2010733), (154, 4652353) (the sieve finds exactly 22 records, logged, and cert.py asserts that they begin with this 22-term list, hence equal it), and g^{1000} ≤ p^{856} in exact integers for each; a gap starting in [p_r, p_{r+1}) (consecutive records) is ≤ g_r ≤ p_r^θ ≤ p_n^θ. Worst case (4, 7): ratio log 4/log 7 = 0.71241…; any θ < log 4/log 7 fails at x = 7.
(iii) p_last ≤ x < e^{42} ≈ 1.739·10^{18}: p_n ≤ x and p_{n+1} ≤ 2p_n < 4·10^{18} (Bertrand), so (p_n, p_{n+1}) is a gap between consecutive primes below 4·10^{18}, hence p_{n+1} ≤ p_n + 1476 ≤ x + 1476 (I5), and 1476^{1000} ≤ p_last^{856} (exact), so 1476 ≤ x^θ.
Joint: the cells start at log x = 42 exactly. Every real x is covered (no integrality used).

## 7. The Chourasiya–Simonič density (conditional runs only; full audit: AUDIT-2507.15184.md)

**7.1 Audit of 2507.15184 v2 (TeX).** Not refereed: arXiv v2 of 30 Sep 2025, "Comments are most welcome", no journal reference or DOI (arXiv API, 6 Nov 2026). Cor. 1.1 (TeX l. 91–99) follows from Thm 2.1 (l. 228–247) with T₀ = H₀ = H = 3·10^{12}, h = 1 (l. 249–275): B₁ = (𝔄₁ + 𝔄₂(log H_RH)^{−6(1−σ₂)/(2−σ₂)})/(2πd), B₂ = 𝔄₃ + 𝔄₄, B₃ = 𝔄₃ log H_RH + 𝔄₅. Thm 2.1 uses Littlewood's lemma, Gabriel's convexity theorem (proved in their §3), KLN's Lemmas 4.1, 4.11–4.13 (refereed), their Lemma 4.2 and their own explicit fourth moment (Thm 2.2, Cor. 2.3–2.4: M₂(10^{10}, 2·10^{10}) ≤ 9.5079·10^{15}, 𝔐₁(10^{10}) = 34.842; Ramachandra's method, parameters by a randomized Mathematica search; not rechecked here).
**7.2 The unsupported constant.** Their Lemma 4.2 (l. 644–670) needs Σ_n λ_X(n)²/n^τ = ζ(τ)S_{τ−1}(X) ≤ b₄e^{γ(τ−1)}/(τ−1) with **b₄ = 0.54** for X ≥ 3·10^{10}, τ ∈ (1, 1.1], cited as a private communication (footnote, l. 659). This removes a factor (log kT)^{2σ−1} from KLN (KLN's (3.8) has b₄τ²e^{γ(τ−1)}log X/(τ−1)). The result appeared as (C2), which gives 0.689 for 10^8 ≤ X < 10^{12} and **0.5501 for 10^{12} ≤ X < 10^{16}** (ε ≤ 1/50): above 0.54. CS use X = T ≥ 3·10^{12} and ε ∈ {δ/log X, 2δ/log X} ≤ ε_max := 2δ/log(3·10^{12}) ≤ 0.0150 < 1/50 (asserted per row). With ζ(1+ε) ≤ 1 + 1/ε and e^{γε} ≥ 1 + γε, (C2) gives the inequality with
  b₄″ = 0.5501(1 + (1−γ)ε_max)  (0.55321 for the row [0.84375, 0.875]).
**7.3 Raised B₁.** K₃ = C₃ + C₄(Γ(3/2)/√(2κ) + πm₀/T₀) (their eq:K3, TeX l. 736–739, with h = 1) is affine in b₄ with non-negative coefficients (C₃ = πm₀b₄e^{2δγ/log X₀}/(2δ), C₄ = b₄/(5δe^δ) + b₃e^{−2δ}/log X₀), so K₃(b₄″) ≤ (b₄″/0.54)K₃; K₃ ≥ 5.96 > 1 in every row (cs_repro.py), so the factors max{1, K₃^{−3d/log T₀}}, max{1, K₃^{−3δ/log(hT₀)}} stay 1 and max{K₃^{(2σ₁−1)/(2−σ₁)}, K₃^{(2σ₂−1)/(2−σ₂)}} = K₃^{(2σ₂−1)/(2−σ₂)}. K₃ enters nothing else. Hence
  **B₁′ := B₁(b₄″/0.54)^{(2σ₂−1)/(2−σ₂)}** (22.1236 for [0.84375, 0.875]; `CSRow.B1`, rounded up),
and B₂′ := B₂ + B₃/log T₀ (log t ≥ log T₀ for t ≥ T₀).
**7.4 Reproduction (cs_repro.py, mpmath, 40 digits; a consistency check, not a proof).** From CS's formulas (Lemmas 4.1–4.5, Cor. 2.4, Thm 2.1) and their parameters (d, κ, δ, log₁₀A₀), our B₁, B₂, B₃ agree with Table 1 within 0.1% in all 13 rows (table/ours = 1.0000–1.0009); with b₄″ the recomputed B₁ is ≤ B₁′ in every row (22.1126 ≤ 22.1236).
**7.5 Use in cert.py.** The 13 rows are [n/32, (n+1)/32], n = 16 … 28, contiguous (asserted: σ₂ of each row equals σ₁ of the next, exact dyadic decimals). For a σ-cell [s, t]: the rows meeting it (they cover [s, t]: s ≥ their smallest σ₁ and t ≤ their largest σ₂, asserted), max B₁′ and max B₂′ (both exact numbers, rounded up). The exponents a(σ) = 3 − 3/(2−σ) and e(σ) = 5 − 3/(2−σ) are concave and decreasing; the tangent lines at the cell midpoint m (an exact double) are affine majorants for all σ < 2: a ≤ α₀ + α₁σ, e ≤ ε₀ + ε₁σ, α₁ = ε₁ = −3/(2−m)². I3 §4.3 cases (i)/(ii) hold verbatim for any exponents a, e ≥ 0 ((log t)^e t^a ≤ δ^{−a}D^e max(1,λ)^{a+e/D} in case (ii), resp. ≤ T₀^aℓ₀^e(λ/λ₀)^{p′} with p′ = a + e/ℓ₀ in case (i)), giving log x∫x^{σ−1}B₁′δ^{−a}D^e(…)dσ ≤ log x·Mom(a(s) + e(s)/pden; λ₀)∫_s^t e^{E₀+E₁σ}dσ with E₀ = −L_a + log B₁′ + α₀D + ε₀Λ, E₁ = L_a + α₁D + ε₁Λ (D, Λ as in I3 §4.3), and the B₂′-term as KLN's B-term. The tail lemmas use KLN only.

## 8. What cert.py runs (`logs/cert.log`; four parts, each in its own process)

| run | pipeline | range of log x (cells) | max R̄ (at log x) | tail | result |
|---|---|---|---|---|---|
| θ = 0.856 | A1 (KLN with RH(3·10^{12}), trapezoid) | [42, 300], 4932 cells | 0.9711 (191.82) | R_an(300) = 0.0142 | CERTIFIED (+ small x §6) |
| θ = 0.855 | A1 | stops at 187.70 (1.0005); width 0.004 on [184, 198]: max 1.152 (190.61) | — | — | fails |
| θ = 0.857 | A0 (I3 unchanged) | [42, 300], 4932 cells | 0.9319 (191.68) | 0.0084 | CERTIFIED |
| θ = 0.856 | A0 | stops at 187.84 (1.0027); width 0.004 on [183, 197]: max 1.106 (190.39) | — | — | fails |
| θ = 0.861, δ = 1 in every KLN row | A1 | [42, 300] | 0.8497 (198.62) | 0.0163 | CERTIFIED; 0.860 fails |
| θ = 0.842 | A1 + CS | [42, 300], 4932 cells | 0.8666 (175.20) | 0.0815 | CERTIFIED (conditional); 0.841 fails (1.0005 at 173.96) |
| k = 6, every n | A1 + CS | [6 log 100, 300], 5448 cells (width 0.05); n ≤ 100 direct | 0.7032 (177.18) | 0.129 | CERTIFIED (conditional) |
| k = 5 | A1 + CS | [227.92, 400], 3442 cells | 0.9997 (227.92) | 0.0528 | n₀ = e^{45.584} |
| k = 4 | A1 + CS | [364.71, 640], 2753 cells | 0.9997 (364.71) | 0.0238 | n₀ = e^{91.18} (364.71/4 = 91.1775) |
| k = 3 | A1 + CS | [1405.27, 2400], 3979 cells | 0.9999 (1405.27) | 0.0162 | n₀ = e^{468.43} (1405.27/3 = 468.423) |

x₀(k) = e^{L₀} is I3's bisection point (to 0.01, rounded up) where the first cell bound drops below 1; all cells above it are certified. The one command runs the four parts (main, a0, delta1, cs) one after another, each in its own process (peak 104, 92, 104, 164 MB; parent 33 MB); 632 s in total, one thread.

## 9. Checks

- **check.py** (independent float quadrature: exact KLN shape (log kt)^{2σ}(log t)^{5−4σ}t^{8(1−σ)/3} with the best row at each σ, exact CS shape B₁′t^a(log t)^e + B₂log²t + B₃log t, the true kernel envelope, no moment tables, no σ-cells), pipeline A1: cert/direct = 1.012–1.207 at 19 points (θ = 0.856 at log x = 42 … 280; θ = 0.842, k = 6 (CS); k = 5, 4, 3 at log x₀); pipeline A0 (PIPE=A0): `logs/check_A0.log`. Kernel random test: §3′.
- **cs_repro.py**: §7.4.
- **Reviewer subagents** (each given the step's text and code inline): tail lemma, small x, B₁ correction, generalization to real x, CS in σ-cells, trapezoid kernel, KLN at RH(3·10^{12}). Fixed after review: Bertrand for the gap below 4·10^{18}; the assertion log x (1 − B) > 1 per cell; CS rows asserted contiguous; exact maxima instead of float-keyed ones; exact grid index in the moment lookup (both inherited from I3; no result changed); the KLN condition of l. 1077 asserted; wording (the ratio log 4/log 7, case (i) of §7.5, the |z|^{−2} decay).
- **Why θ stops at 0.856** (check.py band breakdown, A1, θ = 0.856, log x = 191.82, x ≈ 2·10^{83}, δT₀ = 3.03; Z ≈ 0.949): σ ∈ [0.86, 0.87] gives 0.232 (T ∈ [3·10^{12}, 10^{13}]) + 0.205 (T ∈ [10^{13}, 10^{15}]), σ ∈ [0.87, B] 0.159 + 0.136, σ ∈ [0.83, 0.86] 0.107 + 0.101; σ < 0.83 below 0.008; T > 10^{15} below 10^{−3}. The binding zeros are those KLN allows at β ∈ [0.83, B] at heights 3·10^{12}–10^{15}; the trapezoid trades a lower weight just above the RH height for a longer tail. One step of 0.001 in θ changes R by about 19%.
- **δ = 1** (KLN Thm 1.1 as printed, δ ≥ 1): θ = 0.861 (A1); I3 explains why Lemma 4.14 holds for δ > 0.
- **Small x**: the sieve records equal the published initial maximal-gap records; worst ratio 0.7124 at p = 7.
- Table 1 reproduction (A0 rows) is re-run and asserted in every part that uses KLN at H₀ = 3.0610046·10^{10}, and in the main part as an implementation check.

## 10. Constants (each appears in cert.py)

B = 0.8733819666132076; T* = B₂ + 2.2e-254 (777, 1081, 852); T₀ = 3e12; HSW 0.1038, 0.2573, 9.3675; KLN H₀ = 3.0610046e10 (A0) and 3e12 (A1, `H0_RH`), b₁ = 0.62, b₂ = 1.048, b₃ = 0.605, b₄ = 0.529, a₁ = 0.63, a₂ = 2.851, 3.006, η = 0.2561882, Table 1 rows; σ₀ = 0.515, σ-cell 0.0025; s₁ = 0.2173; δ_max 1e-9, 0.025, 0.15; λ-grid 1 + 1/500, 10000, stored to 30 (triangle) and 2π/a + 1 (trapezoid); p-grid j/200, j ≤ 300; trapezoid a = 0.1 (`A_TRAP`); θ = 856/1000 (main), 0.855; 0.857, 0.856 (A0); 0.860–0.863 (δ = 1); 0.842, 0.841 (conditional); L₁ = 42; L* = 300 (θ), 300, 400, 640, 2400 (k = 6, 5, 4, 3); widths 0.25, 0.05, 0.02 and 150, 165, 215, fine 0.004 on [184, 198] and [183, 197]; sieve 10^7; 1476; 4e18; the 22-term records list; CS rows (13 × 6 numbers), b₄ = 0.54, 0.5501, ε ≤ 1/50; k = 6: n ≤ 100 direct; brackets (200, 240), (330, 400), (1300, 1500); memory limit 300 MB per part.

## 11. What changed in cert.py relative to Ingham3's

Unchanged (copied): inputs, Nplus, all KLN functions and `KP`/`build_kln` (Table 1 reproduction and admissibility asserts), `check_S1`, `S_up`, `Kernel` (except that its stored partial sums are floats rounded up and its grid index is checked exactly), `p_up`, the σ-cell construction, the arithmetic of `cell_bound` (the exact σ-integral moved into `expint`), `tail_pow` (= I3 `tail_bound`), `scan`, `find_L0`, `run_k`.
Changed: `cell_bound` takes a δ-function of log x instead of k (`delta_pow(k)` reproduces I3; `delta_theta(θ)` is new), asserts L_a(1−B) > 1, and, when CS is enabled, takes the smaller of the KLN and CS contributions per σ-cell; `scan` can continue past a failure (fine-cell peaks); `direct_small_k` generalizes I3's k = 7 check; the driver runs four parts in separate processes.
New: `TrapKernel` (§3′), `build_kln_H0` (I4′), `fup`, `exact_max` (A_max, B_max and the CS maxima were float-keyed maxima in I3), `delta_theta`, `tail_theta`, `gap_records`, `small_x_theta`, `run_theta`, `fine_peak`, `CS_ROWS`, `CSRow`, `cs_for_cell`, `run_k_alln`, `run_part`; exploration modes `probe`, `probek`, `probe3`, `fine`, `probe2`, `probe2b` (not part of the certificate).
Removed: I3's Ford sanity check and its δ = 1 run for k (the δ = 1 run is now for θ).
