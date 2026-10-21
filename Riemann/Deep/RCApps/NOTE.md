# RCApps: primes in ring classes for almost all discriminants (short intervals, sectors): proofs

Date 2026-10-21. Brief: `.scratch/orca/rcapps-brief.md`. Files: `.scratch/rcapps/` (`params.py` → `params.out`: thresholds and
parameter checks; `axq.py`, `axq2.py` → `axq.out`, `axq2.out`: arXiv sweeps; `gxl_*.json`: evidence API; `oa.py`: citations of Ditchen).
Project references: [C] = [Combos NOTE](../Combos/NOTE.md) (Theorem 1 = ring class zero density, Lemma 4 = second moment),
[D2] = [Density2 NOTE](../Density2/NOTE.md) (Theorem A, Proposition 6), [LP] = [LeastPrimeForm NOTE](../LeastPrimeForm/NOTE.md)
(Theorem B). Nothing here uses the 2026 half-planes (L), (O), (Hu).

## Verdicts

| Target | Result | Status |
|---|---|---|
| 3. almost-all-D Linnik | for every A > A₀ = 1.653595… (real root of A³−5A²+14A−14), all but O(X^{1−c}) fundamental D < 0 with \|D\| ≤ X: every form of discriminant D represents a prime ≤ 2\|D\|^A (Theorem 2) | **new**. Before: \|D\|^{20/3+ε} outside a density-zero set (Ditchen 2013); [LP] Thm B (this project, unpublished): A > 1.7242 outside a density-zero set, superseded |
| 1. short intervals | same D: every ring class has ≫ H/(\|D\|^{1/2+δ}log x) prime ideals of norm in (x, x+H], all x ≥ \|D\|^A, H ≥ x\|D\|^{−κ} (Theorem 2) | **new** |
| 2. sectors | same, with the angle of the generator in any arc of length ≥ \|D\|^{−κ} (Theorem 3); p = a²+nb² for almost all n (Corollary 4) | **new**; for all n: NOGO (log-free density near σ = 1) |
| tool | Σ_{X/2<\|D\|≤X} Σ_{χ²≠1} N(σ,T;χ) ≪ (XT)^{10(1−σ)/(3−2σ)+ε} (Proposition 1) | new for this family; the general GL₂ bound (Pascadi–Thorner) gives X^{11} |
| uniform in D (\|D\| ≤ x^δ, all D) | needs an explicit log-free density near σ = 1 for each D | no (NOGO unchanged) |

Explicit constants: for every A ≥ 1.7 one may take c = κ = 10^{−4}, δ = 1/2600 (checked in `params.py`).

## 0. Notation and inputs

- 𝔇: negative fundamental discriminants; 𝔇(X) := {D ∈ 𝔇 : X/2 < |D| ≤ X}. For D ∈ 𝔇: K = Q(√D), G = Cl(K), h = |G|,
  w = |O_K^×|, Ĝ the class group characters (all primitive, D fundamental), L(s,χ) their Hecke L-functions. Ĝ₂ := {χ : χ² = 1}
  (genus characters), Ĝ_c := Ĝ ∖ Ĝ₂. N(σ,T;χ) := #{ρ = β+iγ : L(ρ,χ) = 0, β ≥ σ, |γ| ≤ T} with multiplicity;
  N^c_D(σ,T) := Σ_{χ∈Ĝ_c} N(σ,T;χ).
- 𝔉(D): the 2^{ω(D)} integers d | D with d and D/d each equal to 1 or a fundamental discriminant; χ_d := (d/·), χ_1 := 1, L(s,χ_1) = ζ(s).
- (G) Genus theory (Cox, *Primes of the form x²+ny²*, §6; [LP] (Gen)): Ĝ₂ = {χ_{d,D/d}}, one character per unordered pair
  {d, D/d} ⊂ 𝔉(D), with L(s,χ_{d,D/d}) = L(s,χ_d)L(s,χ_{D/d}). Taking −L'/L: for every weight w,
  Σ_{χ∈Ĝ₂} χ̄(C) S_χ(w) = Σ_{d∈𝔉(D)} ε_d(C) ψ_w(χ_d), where S_χ(w) := Σ_𝔞 Λ_K(𝔞)χ(𝔞)w(N𝔞), ψ_w(χ) := Σ_n Λ(n)χ(n)w(n),
  ε_d(C) := χ_{d,D/d}(C) ∈ {±1}, ε_1 = ε_D = 1.

Inputs.
- (I1) [C, Thm 1]: for 0.6 ≤ σ ≤ 1 and t ∈ R, Σ_{χ∈Ĝ} #{ρ : β ≥ σ, |γ−t| ≤ 1} ≤ C(ε)|D|^{2σ(1−σ)+ε}(|t|+3)³, and
  ≤ C|D|^{1/2}log³(3|D|)log(|t|+3) for every σ.
- (I2) [C, Lemma 4] with mollifier length 1: Σ_{χ∈Ĝ, χ≠1} |L(½+iu,χ)|² ≤ C(ε)|D|^{1/2+ε}(|u|+3)^{1+ε}.
- (I3) Humphries–Thorner, Compos. Math. 160 (2024) (arXiv 2103.05634), §5 and Cor 4.4, with n = 2, π′ = 1, F = Q, θ₂ = 0.
  θ₂ = 0 is allowed when every π in the family satisfies the Ramanujan conjecture (their Thm 4.1 uses "the best exponent
  towards GRC for π ∈ 𝒮"); Pascadi–Thorner (arXiv 2508.14888, corollary `cor:MVT_mu` in their proof of Thm 1.3) remove that condition. Statement
  used: let 𝒮 be a finite set of cuspidal automorphic representations of GL₂(A_Q) satisfying GRC, Q ≥ max_π C(π) (analytic
  conductor), T ≥ 3, ε > 0 and V := Q^{2+ε}T^{3+ε}|𝒮|. For every 0 ≤ σ ≤ 1 and every Y ≥ 3 with log Y ≍ log V,
  Σ_{π∈𝒮} N_π(σ,T) ≪_ε V^ε [Y^{2(1−σ)} + Y^{1/2−σ}(V·I₂)^{1/2} + Y^{1−σ}],  I₂ := Σ_{π∈𝒮}∫_{−T}^{T}|L^{ur}(½+iv,π)|²dv.
  This is HT's inequality (bound_int) in §5 with both μ-integrals bounded by their Cor 4.4 (tail ≪ log V, mollifier mean
  square ≪ V log V); HT then bound I₂ by convexity, and we do not.
- (I4) Thorner–Zaman, arXiv 2208.11123, Thm 1.2 and §2. ℒ(s,Q) := Π_{q≤Q}Π*_{χ mod q} L(s,χ) (q = 1 gives ζ), β₁(Q) its
  largest real zero, ν(Q) := min{1, (1−β₁(Q))log Q}; N(σ,Q), N*(σ,Q) count zeros with β > σ, |γ| ≤ Q (N* omits β₁(Q)). For
  Q ≥ 3, σ ≥ 39/40: N ≤ 10^{88}(10^{421}Q^{99})^{1−σ}, N* ≤ 10^{93}ν(Q)(10^{466}Q^{170})^{1−σ}. Zero-free region (their
  McCurley, Platt and Landau lemmas; Mossinghoff–Trudgian for ζ): every zero ρ ≠ β₁(Q) of ℒ(s,Q) with |γ| ≤ Q has
  1−β ≥ c₀/log Q, c₀ := 1/19.3; if β₁(Q) > 1−c₀/log Q it is simple and belongs to a unique real primitive χ₁ of conductor ≤ Q.
- (I5) If χ is real primitive of conductor q ≥ 3 and β is a real zero of L(s,χ), then 1−β ≥ c_S q^{−1/2}log^{−2}q
  (class number formula and |L′(σ,χ)| ≪ log²q near 1; Davenport, Ch. 14).
- (I6) Explicit formula (IK Thm 5.11; [C] (F6)): for w ∈ C_c^∞([x, x+H]), x ≥ 2, and χ a Dirichlet character or χ ∈ Ĝ of
  conductor q: Σ Λχw = δ_{χ=1}W(1) − Σ_ρ W(ρ) + R_χ, W(s) := ∫w(u)u^{s−1}du, |R_χ| ≤ C(H/x + x^{−1/2})log(qx).

## 1. A density estimate for the whole family

**Proposition 1.** For X, T ≥ 3, 0 ≤ σ ≤ 1 and ε > 0: Σ_{D∈𝔇(X)} N^c_D(σ,T) ≤ C(ε)(XT)^{10(1−σ)/(3−2σ)+ε}.

*Proof.* For χ ∈ Ĝ_c let π_χ be the cuspidal representation of GL₂(A_Q) with L(s,π_χ) = L(s,χ) (Hecke's weight-one theta
series of level |D|, nebentypus χ_D; cuspidal since χ ≠ χ̄). π_χ = π_{χ′} iff χ′ ∈ {χ, χ̄}; GRC holds (Satake parameters
χ(𝔭), χ(𝔭̄), or ±1 at inert p, all unimodular); C(π_χ) ≤ 12|D|. Let 𝒮 := {π_χ : D ∈ 𝔇(X), χ ∈ Ĝ_c}; then
|𝒮| ≤ Σ_{|D|≤X} h(D) ≤ CX^{3/2}log X, Q := 12X, and the left side is ≤ 2Σ_{π∈𝒮}N_π(σ,T). The only new input is I₂:
L^{ur}(s,π_χ) = L(s,χ)Π_{p|D}(1 − χ(𝔭_p)p^{−s}) with 𝔭_p² = (p), so |L^{ur}(½+iv,π_χ)| ≤ 2^{ω(D)}|L(½+iv,χ)|, and (I2) gives
I₂ ≤ C(ε)Σ_{D∈𝔇(X)} 4^{ω(D)}|D|^{1/2+ε}T^{2+ε} ≤ C(ε)X^{3/2+2ε}T^{2+ε}. (HT's convexity bound is |𝒮|Q^{1/2}T² ≍ X²T².)
Now V ≤ C(ε)X^{7/2+3ε}T^{3+ε}. Choose Y := (V·X^{3/2}T²)^{1/(3−2σ)}; log Y ≍ log V since 1 ≤ 3−2σ ≤ 3. Then
Y^{1/2−σ}(V I₂)^{1/2} ≤ C(ε)(XT)^{O(ε)}Y^{2(1−σ)}, Y^{1−σ} ≤ Y^{2(1−σ)}, and Y^{2(1−σ)} ≤ C(ε)(XT)^{(10+O(ε))(1−σ)/(3−2σ)}.
Rename ε. □

## 2. Main theorem

**Theorem 2.** Let A₀ = 1.653595… be the real root of A³−5A²+14A−14. For every A > A₀ there are effective c, κ, δ > 0 (for
A ≥ 1.7: c = κ = 10^{−4}, δ = 1/2600) and an effective X₀ with the following property. For X ≥ X₀ let E_A(X) be the set of
D ∈ 𝔇 with |D| ≤ X for which some C ∈ G, x ≥ |D|^A and H ∈ [x|D|^{−κ}, x] violate

  #{𝔭 ∈ C prime of degree one : x < N𝔭 ≤ x+H} ≥ ν_D H / (16 h(D) log 2x).   (2.1)

Then #E_A(X) ≪_A X^{1−c}. Here ν_D := min{1, (1−β₁(|D|^δ))·δlog|D|} if the exceptional character χ₁ of (I4) with
Q = |D|^δ exists and equals χ_d for some d ∈ 𝔉(D), and ν_D := 1 otherwise; by (I5), ν_D ≥ c|D|^{−δ/2}log^{−1}|D|.

The statement for A implies it for every A′ > A (same exceptional set), so take A₀ < A ≤ 1.7.

### 2.1 Parameters
φ(σ) := (10/(3−2σ) − A)(1−σ) is decreasing on [0.8, 1] (φ′ ≤ −(10/1.4 − 2) + 0.2·20/1.96 < 0) and φ(A/2) < 1 ⟺
(2−A)(A²−3A+10) < 2(3−A) ⟺ A³−5A²+14A−14 > 0 ⟺ A > A₀. Fix σ₀ ∈ [0.8, A/2) ⊂ [0.8, 0.85) with φ(σ₀) < 1 (σ₀ = 0.84 at A = 1.7) and set
ε′ := (1−φ(σ₀))/4, s₁ := (1−σ₀)(A−2σ₀) > 0, δ := min{A/4400, s₁/2, ε′}, κ := δ/4 (κ := 10^{−4} allowed for A ≥ 1.7),
ε₀ := ε₁ := δ/100, σ_g := 0.995, σ_j := σ₀ + jε′/A (0 ≤ j ≤ J, J minimal with σ_J ≥ 1−ε′/A), T_X := 2X^{κ+ε₁}, Q₀(X) := (X/2)^δ.
`params.py` checks every inequality used below at A = 1.7 with these values (and for Theorem 3).

### 2.2 Exceptional sets (X a power of 2, D ∈ 𝔇(X))
- E₁(X): some j has N^c_D(σ_j, T_X) > X^{A(1−σ_j)−2ε′}. By Proposition 1 (with ε = ε₀) and Markov, since φ(σ_j) ≤ φ(σ₀) and
  T_X^{10(1−σ_j)/(3−2σ_j)} ≤ CX^{(10/7)(κ+ε₁)} (σ_j ≥ 0.8): #E₁(X) ≤ Σ_j X^{2ε′−A(1−σ_j)}C(XT_X)^{10(1−σ_j)/(3−2σ_j)+ε₀}
  ≤ CJ X^{1−4ε′+2ε′+3ε₀+(10/7)(κ+ε₁)} ≤ CJ X^{1−ε′}.
- E₂(X): some d ∈ 𝔉(D) with |d| > Q₀(X) has L(ρ,χ_d) = 0 with β > σ_g, |γ| ≤ T_X. If |d| ∈ (Z, 2Z] with Z ≥ Q₀(X) ≥ T_X, such
  a zero is counted by N(σ_g, 2Z) of (I4), so there are ≤ 10^{88}(10^{421}(2Z)^{99})^{0.005} ≤ 10^{91}Z^{0.495} such d; each
  divides ≤ X/Z+1 discriminants of 𝔇(X). Summing over dyadic Z: #E₂(X) ≤ C X Q₀(X)^{−0.505} ≤ C X^{1−0.505δ}.
- Union over dyadic X: #E_A(X) ≤ C X^{1−min{ε′,0.505δ}} + O(X₀) (D with |D| < D₀ go into E_A). c := δ/2 ≤ min{ε′, 0.505δ}.

### 2.3 Weights and decomposition
Fix D ∈ 𝔇(X) ∖ (E₁ ∪ E₂) with |D| ≥ D₀, C ∈ G, x ≥ |D|^A, H ∈ [x|D|^{−κ}, x]; Q₀ := |D|^δ, T₁ := 2|D|^{κ+ε₁} ≤ T_X.
Let w₀ ∈ C_c^∞((0,1)), 0 ≤ w₀ ≤ 1, w₀ = 1 on [¼,¾], w(u) := w₀((u−x)/H). Then W(1) ≥ H/2, |W(ρ)| ≤ Hx^{β−1} (as u ≥ x,
β ≤ 1), and k-fold integration by parts gives |W(ρ)| ≤ C_k Hx^{β−1}(2x/(H|ρ|))^k. Let S_C := Σ_{𝔞∈C}Λ_K(𝔞)w(N𝔞). By
orthogonality and (G),

  h S_C = Σ_{d∈𝔉(D)} ε_d(C)ψ_w(χ_d) + Σ_{χ∈Ĝ_c} χ̄(C)S_χ(w).   (2.2)

Split 𝔉(D) = 𝔉_s ⊔ 𝔉_l, 𝔉_s := {d : |d| ≤ Q₀} (so 1 ∈ 𝔉_s, D ∈ 𝔉_l). Put β₁ := β₁(Q₀), and M(C) := W(1) − ε_{d₀}(C)W(β₁) if
χ₁ exists and χ₁ = χ_{d₀} with d₀ ∈ 𝔉_s, M(C) := W(1) otherwise. Apply (I6) to every term of (2.2) and move W(1) (from d = 1)
and W(β₁) (from d₀) into M(C). Write E_a, …, E_f for the moduli of the following parts of h S_C − M(C):
(a) zeros of χ ∈ Ĝ_c with β ≥ σ₀, |γ| ≤ T₁; (b) zeros of χ ∈ Ĝ_c with β < σ₀, |γ| ≤ T₁; (c) zeros of any character with
|γ| > T₁; (d) zeros of χ_d, d ∈ 𝔉_l, |γ| ≤ T₁; (e) zeros of χ_d, d ∈ 𝔉_s, |γ| ≤ T₁, other than the zero β₁ moved into M(C);
(f) all remainders R.

### 2.4 Bounds (|D| ≥ D₀)
- (a) Since D ∉ E₁ and x^{σ−1} ≤ |D|^{−A(1−σ)} ≤ 2^A X^{−A(1−σ)}: E_a ≤ HΣ_{j<J}N^c_D(σ_j,T_X)x^{σ_{j+1}−1}
  ≤ H J 2^A X^{Aε′/A−2ε′} = J2^A H X^{−ε′}, and N^c_D(σ_J,T_X) < 1 because A(1−σ_J) ≤ ε′.
- (b) For a unit window |γ−k| ≤ 1, with N(σ) := #{zeros of all χ ∈ Ĝ : β ≥ σ, |γ−k| ≤ 1}: the zeros with β < σ₀ give
  Σ x^{β−1} ≤ x^{−1/2}N(0) + log x∫_{1/2}^{σ₀}x^{σ−1}N(σ)dσ. Write x = |D|^A y, y ≥ 1; for σ ≤ σ₀,
  x^{σ−1}log x = |D|^{−A(1−σ)}y^{−(1−σ)}(A log|D| + log y) ≤ |D|^{−A(1−σ)}(A log|D| + 10) (1−σ ≥ 0.1). By (I1), the part
  σ ∈ [½, 0.6] is ≤ |D|^{1/2−0.4A+o(1)} ≤ |D|^{−0.16}, and on [0.6, σ₀], x^{σ−1}N(σ) ≤ C(ε₀)|D|^{−(1−σ)(A−2σ)+ε₀}(|k|+3)³ with
  (1−σ)(A−2σ) ≥ s₁ there (decreasing, as σ₀ < A/2 < (A+2)/4). Summing over |k| ≤ T₁: E_b ≤ H|D|^{−s₁+ε₀+4(κ+ε₁)+o(1)}.
- (c) Each window holds ≤ C h log(|D|(|k|+3)) zeros (all of Ĝ, or of the ≤ 2^{ω(D)} χ_d). With k ≥ 2/ε₁ in the decay bound
  and x/H ≤ |D|^κ: E_c ≤ H|D|^{−1}.
- (d) Since D ∉ E₂, these zeros have β ≤ σ_g, and there are ≤ 2^{ω(D)}·CT₁log|D| of them: E_d ≤ H|D|^{−0.005A+κ+ε₁+o(1)}.
- (e) Every such zero is a zero of ℒ(s,Q₀) with |γ| ≤ Q₀ (T₁ ≤ Q₀; distinct d give distinct primitive characters).
  Zeros with β ≤ 39/40: ≤ CQ₀³log Q₀ of them, total ≤ HQ₀³x^{−1/40}log Q₀ ≤ H|D|^{−0.04}. Zeros with β > 39/40:
  *Case 1* (χ₁ exists, i.e. β₁ > 1−c₀/L₀): the zeros counted in (e) are zeros ≠ β₁, counted by F := N*(·,Q₀). Put
  ℓ := log x ≥ (A/δ)L₀, L₀ := log Q₀ ≥ 100, Λ := ℓ − 1073.1 − 170L₀ ≥ 0.9ℓ (as 170δ/A ≤ 0.04, |D| ≥ D₀); F(σ) = 0 for
  σ > 1−c₀/L₀ by (I4). Partial summation and (I4):
  Σ x^{β−1} = x^{−1/40}F(39/40) + ℓ∫_{39/40}^{1−c₀/L₀}x^{σ−1}F(σ)dσ ≤ 10^{93}ν(Q₀)[e^{−0.9ℓ/40} + 1.12·e^{−c₀Λ/L₀}].
  Since Λ/L₀ ≥ A/δ − 181 and A/δ ≥ 4400 > 4383 (`params.py`), the bracket is ≤ 1/(16·10^{93}): E_e ≤ H(ν(Q₀)/16 + |D|^{−0.04}).
  *Case 2* (no χ₁): every zero has 1−β ≥ c₀/L₀; the same computation with N(·,Q₀) (constants 10^{88}, 10^{421}, Q^{99})
  gives E_e ≤ H(1/16 + |D|^{−0.04}).
- (f) Σ|R| ≤ C(h + 2^{ω(D)})(H/x + x^{−1/2})log(|D|x) ≤ H|D|^{−1}.

### 2.5 Conclusion
If χ₁ = χ_{d₀} with d₀ ∈ 𝔉(D): M(C) ≥ W(1) − W(β₁) ≥ (1−x^{β₁−1})W(1) ≥ min{½, (1−β₁)ℓ/2}W(1) ≥ ν_D W(1)/2 ≥ ν_D H/4
(W(β₁) ≤ x^{β₁−1}W(1), ℓ ≥ L₀, ν_D = ν(Q₀)), and E_e ≤ Hν_D/16 + H|D|^{−0.04}. Otherwise M(C) = W(1) ≥ H/2, ν_D = 1 and
E_e ≤ H/16 + H|D|^{−0.04}. The other errors are ≤ H|D|^{−δ/2−η} for some η > 0, because ε′ > δ/2,
s₁ − ε₀ − 4(κ+ε₁) > δ/2 and 0.005A − κ − ε₁ > δ/2 (§2.1, `params.py`), while ν_D ≥ c|D|^{−δ/2}log^{−1}|D|. Hence
h S_C ≥ ν_D H/4 − ν_D H/16 − ν_D H/16 = ν_D H/8 for |D| ≥ D₀. Prime powers 𝔭^k (k ≥ 2), inert and ramified primes contribute
≤ C(Hx^{−1/2} + log x)log x ≤ ν_D H/(16h) to S_C. So Σ log N𝔭 ≥ ν_D H/(16h) over degree-one 𝔭 ∈ C with N𝔭 ∈ (x, x+H],
which gives (2.1). □

## 3. Corollaries
- Least prime: for D ∉ E_A, x = |D|^A, H = x: every class has a degree-one prime of norm ≤ 2|D|^A, so every primitive positive
  definite form of discriminant D represents a prime p ∤ D with p ≤ 2|D|^A ([C] (F1), (F2)). Exceptional set O(X^{1−c}).
- Short intervals: every such form represents ≫ H/(|D|^{1/2+δ}log x) primes in (x, x+H] for all x ≥ |D|^A,
  H ≥ x|D|^{−κ}; ≫ H/(h(D)log x) unless a genus factor χ_d, |d| ≤ |D|^δ, has an exceptional zero.
- **Corollary 4.** For A ≥ 1.7, all but O(N^{1−10^{−4}}) squarefree n ≤ N with n ≡ 1, 2 (mod 4) have: for all x ≥ (4n)^A,
  H ≥ x(4n)^{−10^{−4}} and every arc I ⊂ R/πZ of length ≥ (4n)^{−10^{−4}}, there is a prime p = a²+nb² in (x, x+H] with
  arg(a+ib√n) ∈ I. (D = −4n, O_K = Z[√−n], C = principal class, 𝔠 = O_K, β_𝔭 = a+b√−n; Theorem 3.) Without the arc: the
  least prime a²+nb² is ≤ 2(4n)^{A}. Previously: ≪ n^{L₂}, L₂ unspecified, for all n (Kowalski–Michel 2002, Cor 7).

## 4. Sectors

**Theorem 3.** Same A, c, κ, δ, ν_D (for A ≥ 1.7 the same numbers). For all D ∈ 𝔇, |D| ≤ X, outside a set of size ≪ X^{1−c}:
for every C ∈ G, every integral 𝔠 ∈ C^{−1}, every arc I ⊂ R/(2π/w)Z of length θ ≥ |D|^{−κ}, every x ≥ |D|^A and
H ∈ [x|D|^{−κ}, x], the number of degree-one 𝔭 ∈ C with x < N𝔭 ≤ x+H and arg β_𝔭 ∈ I, where 𝔭𝔠 = (β_𝔭), is
≥ wθν_D H/(128π h log 2x). (Changing 𝔠 rotates all arg β_𝔭 by one fixed angle, so the statement does not depend on 𝔠.)

*Proof (changes to §2).* (i) g ∈ C^∞(R/(2π/w)Z), 0 ≤ g ≤ 1, supported in I, = 1 on its middle half, g = Σ_k ĝ(k)e^{iwkφ};
ĝ(0) ≥ wθ/4π, |ĝ(k)| ≤ ĝ(0), |ĝ(k)| ≤ C_m θ(θ|k|)^{−m}. (ii) With 𝒜(D,M) as in [D2] (f = 1; χ_ℓĜ = unramified characters of
∞-type ℓ ∈ wZ): for 𝔞 ∈ C, (β_𝔞/|β_𝔞|)^ℓ = ψ(𝔞)ψ(𝔠) for ψ ∈ χ_ℓĜ, hence
Σ_{𝔞∈C}Λ_K w(N𝔞)g(arg β_𝔞) = ĝ(0)S_C + Σ_{k≠0}ĝ(k)h^{−1}Σ_{ψ∈χ_{wk}Ĝ}ψ(𝔠)S_ψ(w). (iii) |k| > K := θ^{−1}|D|^{ε₂}
(ε₂ := δ/100): |Σ_{𝔞∈C}Λ_K w e^{iwk arg β_𝔞}| ≤ 2Σ_{x<n≤x+H}Λ(n) ≤ 5H (Brun–Titchmarsh, H ≥ x^{0.99}); with m ≥ 1 + 12/ε₂ the total
is ≤ 5H·C_m(θK)^{1−m} ≤ ĝ(0)H|D|^{−10} (ĝ(0) ≥ |D|^{−κ}/2π). (iv) 0 < |k| ≤ K: L(s,ψ) is entire (CM form
of weight |ℓ|+1, ψ ≠ ψ^c cuspidal, GRC, multiplicity ≤ 2 in the family since ψ^c = ψ̄); no genus characters occur. The total
is ≤ ĝ(0)h^{−1}Σ_{ψ∈𝒜(D,M)∖Ĝ, M = wK ≤ 6|D|^{κ+ε₂}}Σ_ρ|W(ρ)|, bounded as (a)–(c) with: [D2, Thm A] (𝔔 = |D|(M+3)²,
per window T^{1+ε}) in place of (I1), so s₁ becomes (1−σ₀)(A − 2σ₀(1+2μ)), μ := κ+ε₂+o(1); and Proposition 1′ in place
of Proposition 1: Σ_{D∈𝔇(X)}Σ_{ψ∈𝒜(D,M)∖Ĝ₂}N(σ,T;ψ) ≤ C(ε)(X⁵M⁶T⁵)^{2(1−σ)/(3−2σ)+ε}, proved like Proposition 1 with
|𝒮′| ≤ CX^{3/2}M log X, C(π_ψ) ≤ 12X(M+4)², and I₂′ ≤ C(ε)X^{3/2}(M+3)T²(XMT)^ε from [D2, (2.1)] with mollifier length 1.
The exceptional set E₁′ uses M_X := 6X^{κ+ε₂}; φ becomes φ′(σ) := ((10+12μ)/(3−2σ) − A)(1−σ) and ε′ := (1−φ′(σ₀))/4;
at A = 1.7 (σ₀ = 0.84, μ = 1.01·10^{−4}) φ′(σ₀) = 0.94027, ε′ = 0.0149, s₁′ = 0.00315, and every check of §2.1 holds (`params.py`). (v) D ∉ E₁ ∪ E₂ ∪ E₁′ gives h S_C ≥ ν_D H/8 (§2.5), so Σ_{𝔞∈C}Λ_K w g(arg β_𝔞) ≥ ĝ(0)ν_D H/(8h) − ĝ(0)H|D|^{−δ/2−η}/h
− ĝ(0)H|D|^{−10} ≥ ĝ(0)ν_D H/(16h). Remove prime powers as in §2.5 and use 0 ≤ w, g ≤ 1, supp g ⊂ I. □

## 5. What fixes A₀
The binding point is σ = A/2, where (I1) (per D) hands over to Markov (average over D). With family exponent F(σ):

| Family density for Σ_{\|D\|~X}N^c_D | A₀ |
|---|---|
| Kowalski–Michel 2002, Thm 2 (c₀ > 13/2): X^{c₀(1−σ)/(2σ−1)} | 1.72419 (= [LP] Thm B) |
| Pascadi–Thorner 2508.14888, Thm 1.3 for GL₂ (Q⁵\|𝒮\|⁴ = X^{11}) | 1.70555 |
| **Proposition 1**: X^{10(1−σ)/(3−2σ)} | **1.65360** |
| if the GL₂ family large sieve had threshold \|𝒮\| = X^{3/2} instead of Q²\|𝒮\| = X^{7/2}: X^{6(1−σ)/(3−2σ)} | 1 (condition (A−1)(A²−4A+6) > 0) |
| GRH | 1 + ε for every D |

So the method reaches the GRH exponent for almost all D once the threshold of the GL₂ large sieve over this family drops to
its size; the remaining loss X² is the Rankin–Selberg conductor of π_χ × π̃_{χ′} for D ≠ D′.

## 6. Novelty audit (2026-10-21)
- arXiv API, 41 queries (`axq.out`, `axq2.out`): binary quadratic forms × almost all / short interval(s) / Linnik / least prime /
  Bombieri–Vinogradov / BDH; ring class × primes; class group × short intervals / large sieve / zero density / Deuring–
  Heilbronn; x²+ny²; Hecke angles / sectors; Chebotarev × almost all × families; au:Pascadi; au:Thorner AND quadratic.
  Relevant hits only: Ditchen 1312.1502 (read: §1, Lemmas; density-one 20/3, BV with (log X)^{−A} savings, large sieve
  N + N^{1/2}Q^{5/2}); Zaman 1710.08914 (upper bounds); Pierce–Turnage-Butterbaugh–Wood 1709.09637, Chen An 2004.14410,
  Lemke Oliver–Thorner–Zaman 2012.14422, Lemke Oliver–Smith 2405.08383 (almost-all-fields Chebotarev for fixed Galois group;
  ring class fields have degree 2h → ∞, not covered); Panidapu–Thorner 2105.11093 (sectors, fixed field); Gaudet 2408.01610.
- OpenAlex: Ditchen's paper has 3 citing works (TZ17, Zaman 2016, Zaman 2018); none improves it.
- GXL evidence API, 2 claims: only Ditchen supports (3+ε for most forms); nothing below 20/3 for all forms of almost all D,
  nothing on power-saving exceptional sets or short intervals.
- Pascadi–Thorner 2508.14888 (tex) and Humphries–Thorner 2103.05634 (tex): no quadratic forms, class groups or dihedral
  families. Hexagon 2610.00178: §54 (all D, uses (O), A unspecified), no almost-all statements. Kowalski–Michel 2002 (pdf):
  Cor 7 only. Duke (Acta Arith. 52): angles in the T-aspect for a fixed field.
- Read for inputs: Thorner–Zaman 2208.11123 (Thm 1.2, §2 lemmas), Thorner–Zaman 2108.10878 (statements).

## 7. Dead ends (NOGO candidates) and leads
- Uniform-in-D versions (all D, |D| ≤ x^δ, or least prime in a sector for every n): need an explicit log-free density near
  σ = 1 for each D; (I1) is |D|^ε-lossy there (unchanged NOGO, [LP] §1).
- Per-D large sieve inside the HT detector (to shorten the mollifier below X^{7/2}): the tail term then costs a factor X
  (each D separately), so the family bound exceeds X. Mollifier length X^v < X^{7/2} gives Σ N ≥ X^{1+…}.
- Only the D-aspect improves: (I1) handles σ < A/2 for every D; Markov needs the family sum below X^{1+A(1−σ)}.
- Leads: (1) a GL₂ large sieve for {π_χ : |D| ≤ X} with threshold below X^{7/2} lowers A₀ toward 1 (§5); (2) non-fundamental
  D (orders): Proposition 1 needs the multiplicity of each π_{χ*} over D = f²d_K; (3) explicit X₀ and larger c (here c is set
  by δ ≤ A/4400, i.e. by the constants 10^{93}, Q^{170} of (I4)).
