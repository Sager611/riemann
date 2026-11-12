# Natural-number Nyman–Beurling (Báez-Duarte) criterion in L^p(0,1), 1 < p < 2

Node N5. Author: Orca agent `nblp`, 2026-11-05. Paper proof (no Lean). Every analytic input is either proved
here or quoted with its exact statement in §2.

## 0. Notation

* `s = σ + it`; `{x} = x − ⌊x⌋`; `μ` = Möbius; `ρ = β + iγ` runs over the non-trivial zeros of ζ (with
  multiplicity); `Θ := sup β`.
* For `k ≥ 1`, `x > 0`: `ρ_k(x) := {1/(kx)}` (Báez-Duarte's functions on (0,∞)).
* For `0 < x < 1`: `f_k(x) := {1/(kx)} − k^{-1}{1/x}`. This is the brief's `ρ_{1/k}` (`ρ_a(x) = {a/x} − a{1/x}`,
  `a = 1/k`). Note `f_1 ≡ 0`.
* `χ := 1_{(0,1)}` (a function on (0,1)); `χ₊ := 1_{(0,1]}` (a function on (0,∞)).
* `I_n := (1/(n+1), 1/n]`, `n ≥ 1`, a partition of (0,1]; `|I_n| = 1/(n(n+1))`.
* `V := span{f_k : k ≥ 2}`; `W_p :=` closure of `V` in `L^p(0,1)`;
  `S_p := {F ∈ L^p(0,1) : F is a.e. constant on every I_n}`;
  `d_{N,p} := inf_{c_2,…,c_N ∈ ℂ} ‖χ − Σ_{k=2}^N c_k f_k‖_{L^p(0,1)}`.
* `1 < p < 2`, `r := 2p/(2−p)` (so `1/p = 1/2 + 1/r`).
* `L²_σ := L²((0,∞), x^{2σ−1} dx)`.
* `γ(s) := ½ s(s−1) π^{−s/2} Γ(s/2)`, `ξ(s) := γ(s) ζ(s)`; for `ε > 0`: `Γr_ε(s) := γ(s+ε)/γ(s)`.
* `P_N(s) := Σ_{k ≤ N} μ(k)(1 − k/N) k^{−s}` (Riesz-smoothed Möbius polynomial).

## 1. Results

**Proposition C (structure).**
(a) For `k ≥ 1`, `n ≥ 1`, `x ∈ I_n`: `f_k(x) = {n/k}`.
(b) `V ⊂ S_p`; `S_p` is a closed subspace of `L^p(0,1)` of infinite codimension. Hence **no** set of the
functions `f_k` spans a dense subspace of `L^p(0,1)`, for any `p`; the literal "dense in L^p(0,1)" of the brief is false.
(c) For `m ≥ 1` let `(D_m F)(x) := F(mx)` for `0 < x < 1/m`, `:= 0` for `1/m ≤ x < 1`. Then
`‖D_m F‖_p = m^{−1/p}‖F‖_p` and `D_m f_k = f_{km} − k^{−1} f_m` a.e.; so `D_m W_p ⊂ W_p`.
(d) `χ ∈ W_p ⟺ W_p = S_p ⟺ d_{N,p} → 0`.

**Theorem A (exact L^p analogue of Báez-Duarte).** Let `1 < p < 2`. The following are equivalent:
(i) `ζ(s) ≠ 0` for `Re s > 1/p` (i.e. `Θ ≤ 1/p`);
(ii) `χ ∈ W_p`, i.e. `d_{N,p} → 0`;
(iii) `W_p = S_p`, i.e. `span{ρ_{1/k} : k ≥ 2}` is dense in the step-function space `S_p`.
Moreover, under (i), with `ε_0 := min{1/2, (2/p − 1)/2, 1 − 1/p}`, for every `ε ∈ (0, ε_0]` there is
`N_0(p, ε)` such that for all `N ≥ N_0(p,ε)`
`‖ χ + Σ_{k=2}^{N} μ(k) k^{−ε} (1 − k/N) f_k ‖_{L^p(0,1)} ≤ A_p ε^{1 − 1/p} + 2ε`,
with `A_p := K_p^{1/2} (2/r)^{1/r}`, `K_p := 4p/(2−p) + C_Γ/(2π)`,
`C_Γ := e^4 ∫_ℝ (8 + log(2+|t|)) (2+|t|)^{1/2} (1/4 + t²)^{−1} dt = 5993.05…`.

**Theorem B (explicit coefficients, power rate, open range).** Let `1/2 ≤ θ < 1`, assume `ζ(s) ≠ 0` for
`Re s > θ`, let `1 < p < 1/θ` and `0 < η < 1/p − θ`. Then there is `C = C(p, θ, η)` with
`‖ χ + Σ_{k=2}^{N} μ(k)(1 − k/N) f_k ‖_{L^p(0,1)} ≤ C N^{−η}` for all `N ≥ 2`.

**Corollaries.**
1. (Unconditional.) By [C9] (`ζ ≠ 0` for `Re s > 7/8`): `χ ∈ W_p` for every `1 < p ≤ 8/7`, and for
   `1 < p < 8/7` the explicit bound `‖χ + Σ_{k≤N} μ(k)(1−k/N) f_k‖_p ≪_{p,η} N^{−η}` holds for every `η < 1/p − 7/8`.
2. For `1 < p < 2`: `χ ∈ W_p ⟺ p ≤ 1/Θ`. For `p > 2`: `χ ∉ W_p` (unconditionally). (`p = 2`: Báez-Duarte [C10].)
3. RH `⟺` `χ ∈ W_p` for every `p ∈ (1,2)`. Under RH, for `1 < p < 2` and `η < 1/p − 1/2`:
   `‖χ + Σ_{k≤N} μ(k)(1−k/N) f_k‖_p ≪ N^{−η}` (a power rate, unlike `d_{N,2} ≍ (log N)^{−1/2}`).

## 2. Cited results (exact statements)

* **[C1]** (Titchmarsh, *The Theory of the Riemann Zeta-Function*, 2nd ed., eq. (2.1.5).) For `0 < σ < 1`:
  `ζ(s) = −s ∫_0^∞ {u} u^{−s−1} du`.
* **[C2]** (Davenport, *Multiplicative Number Theory*, 3rd ed., Ch. 12.) `ξ` is entire of order 1, its zeros are
  exactly the non-trivial zeros `ρ` of ζ (`0 < β < 1`), `Σ_ρ |ρ|^{−2} < ∞`,
  `ξ(s) = e^{A + Bs} Π_ρ (1 − s/ρ) e^{s/ρ}` (absolutely convergent), and `B = −Σ_ρ Re(1/ρ)`, the series being
  absolutely convergent (`Re(1/ρ) = β/|ρ|² ≤ |ρ|^{−2}`).
* **[C3]** (Binet; DLMF §5.9(ii), Whittaker–Watson §12.31.) For `Re z > 0`:
  `log Γ(z) = (z − ½) log z − z + ½ log(2π) + ∫_0^∞ (½ − 1/t + 1/(e^t − 1)) e^{−tz} t^{−1} dt`.
  (DLMF 4.36.3.) `coth z = 1/z + Σ_{k≥1} 2z/(z² + k²π²)`.
* **[C4]** (Borel–Carathéodory; Titchmarsh, *The Theory of Functions*, §5.5.) If `f` is analytic on
  `|z − z_0| ≤ R` and `Re f(z) ≤ U` on `|z − z_0| = R`, then for `|z − z_0| ≤ r < R`:
  `|f(z)| ≤ (2r/(R − r)) U + ((R + r)/(R − r)) |f(z_0)|`.
* **[C5]** (Hadamard three circles; ibid. §5.3.) If `f` is analytic on `|z − z_1| ≤ r_3` and
  `M(ϱ) := max_{|z − z_1| = ϱ} |f(z)|`, then for `0 < r_1 < r_2 < r_3`:
  `log M(r_2) ≤ [log(r_3/r_2) log M(r_1) + log(r_2/r_1) log M(r_3)] / log(r_3/r_1)`.
* **[C6]** (Titchmarsh, *Zeta-Function*, Thm 4.11.) For `σ ≥ σ_0 > 0`, `|t| < 2πx/C` (`C > 1` fixed):
  `ζ(s) = Σ_{n ≤ x} n^{−s} − x^{1−s}/(1−s) + O(x^{−σ})`. Consequence used: with `σ_0 = 1/2`, `C = 2`,
  `x = |t| + 1`: there is an absolute `C_6` with `|ζ(σ+it)| ≤ C_6 |t|^{1/2}` for `σ ≥ 1/2`, `|t| ≥ 1`
  (`Σ_{n≤x} n^{−1/2} ≤ 2√x`, `|x^{1−s}/(1−s)| ≤ x^{1/2}/|t|`).
* **[C7]** (Plancherel.) `h ↦ ĥ(t) = ∫_ℝ h(u) e^{−iut} du` extends to `L²(ℝ) → L²(ℝ)` with
  `‖ĥ‖² = 2π‖h‖²`; `(1/2π)∫_{−T}^{T} ĥ(t) e^{iut} dt → h` in `L²(du)` as `T → ∞`; for `h ∈ L¹ ∩ L²` the
  integral formula holds.
* **[C8]** (Elementary, `σ > 1`.) `ζ(s) = Π_p (1 − p^{−s})^{−1}`, `log ζ(s) = Σ_p Σ_{m≥1} p^{−ms}/m`,
  `|ζ(s)| ≤ ζ(σ)`, `|1/ζ(s)| = |Σ μ(n) n^{−s}| ≤ ζ(σ)`, `|log ζ(s)| ≤ log ζ(σ)`,
  `1/(σ−1) ≤ ζ(σ) ≤ σ/(σ−1)`; `ζ` is decreasing on `(1, ∞)`; `ζ(s̄) = conj ζ(s)`.
* **[C9]** (OpenAI 2026; Lean-checked in `Riemann/Deep/QRH`.) `ζ(s) ≠ 0` for `Re s > 7/8`.
* **[C10]** (Báez-Duarte, Rend. Mat. Acc. Lincei 2003, arXiv math/0202141.) RH `⟺` `χ ∈ W_2`.

## 3. Structure: proof of Proposition C

(a) Let `x ∈ I_n`, `y := 1/x ∈ [n, n+1)`. For an integer `k ≥ 1`, `⌊y/k⌋ = ⌊⌊y⌋/k⌋ = ⌊n/k⌋`. Hence
`{y/k} − k^{−1}{y} = (y/k − ⌊n/k⌋) − (y − n)/k = n/k − ⌊n/k⌋ = {n/k}`.

(b) By (a) each `f_k` is constant on every `I_n`, so `V ⊂ S_p`. The map `F ↦ (F|_{I_n})_n` is an isometry of
`S_p` onto the weighted sequence space `ℓ^p(w)`, `w_n = 1/(n(n+1))`; `S_p` is closed (an `L^p` limit of functions
constant on `I_n` is a.e. constant on `I_n`). `S_p ≠ L^p(0,1)`: e.g. `g := 1_{(1/2, 3/4)} − 1_{(3/4, 1)}` is in
`L^q` and annihilates `S_p` (`∫_{I_1} g = 0` and `g = 0` off `I_1`), but `∫ g · 1_{(1/2,3/4)} = 1/4 ≠ 0`.
Infinitely many such `g` (one per `I_n`) are independent, so the codimension is infinite.

(c) `‖D_m F‖_p^p = ∫_0^{1/m} |F(mx)|^p dx = m^{−1} ‖F‖_p^p`. For `x < 1/m`:
`f_{km}(x) − k^{−1} f_m(x) = {1/(kmx)} − (km)^{−1}{1/x} − k^{−1}{1/(mx)} + (km)^{−1}{1/x} = f_k(mx)`; for
`x > 1/m` both sides vanish (`y := mx > 1`: `{1/(ky)} − k^{−1}{1/y} = 1/(ky) − 1/(ky) = 0`); `x = 1/m` is
null. So `D_m V ⊂ V`, and by continuity `D_m W_p ⊂ W_p`.

(d) If `W_p = S_p` then `χ ∈ S_p = W_p`. Conversely if `χ ∈ W_p`, then by (c) `D_m χ = 1_{(0,1/m)} ∈ W_p`,
so `1_{I_m} = 1_{(0,1/m)} − 1_{(0,1/(m+1))}` a.e. lies in `W_p` for all `m`; finite combinations of the `1_{I_m}`
are dense in `S_p ≅ ℓ^p(w)` (`p < ∞`), so `W_p = S_p`. Finally `χ ∈ W_p ⟺ d_{N,p} → 0` by definition
(`d_{N,p}` is non-increasing in `N`, and `f_1 = 0`). ∎

**Proposition N (necessity; = (ii) ⟹ (i)).** Let `1 < p < ∞`. If `ζ(w_0) = 0` with `β_0 := Re w_0 > 1/p`,
then `χ ∉ W_p`.

*Proof.* `g(x) := x^{w_0 − 1}` lies in `L^q(0,1)` (`q = p/(p−1)`), because `q(β_0 − 1) > −1 ⟺ β_0 > 1/p`.
By Lemma T below (with `s = w_0`, `0 < β_0 < 1`), `∫_0^1 f_k(x) x^{w_0−1} dx = ∫_0^∞ (ρ_k − k^{−1}ρ_1) x^{w_0 − 1} dx
= −ζ(w_0) w_0^{−1} (k^{−w_0} − k^{−1}) = 0` (`ρ_k − k^{−1}ρ_1 = 0` on `(1,∞)`; here `ρ_k` are the functions of §0),
so `g` annihilates `V`, hence `W_p`; but `∫_0^1 χ g = 1/w_0 ≠ 0`. ∎

(For `p > 2` take `w_0 = 1/2 + 14.1347…i`: `χ ∉ W_p` unconditionally. In sequence form Proposition N is
Gallardo-Gutiérrez–Seco, Theorem A with `α = −2`.)

## 4. Mellin–Plancherel toolkit

For `F ∈ L²_σ` put `φ(u) := F(e^{−u})` and `𝓜_σ F := (φ(u) e^{−σu})^∧` (Fourier transform [C7]). Since
`∫_0^∞ |F|² x^{2σ−1} dx = ∫_ℝ |φ(u) e^{−σu}|² du`, `𝓜_σ : L²_σ → L²(ℝ)` is a bijection with
`‖F‖²_{L²_σ} = (2π)^{−1} ∫_ℝ |𝓜_σ F(t)|² dt`, and if `x^{σ−1}F ∈ L¹(0,∞)` then
`𝓜_σ F(t) = ∫_0^∞ F(x) x^{s−1} dx`, `s = σ + it` (substitute `x = e^{−u}`).

**Lemma T.** (a) For `0 < σ < 1`, `k ≥ 1`: `x^{σ−1}ρ_k ∈ L¹(0,∞)`, `ρ_k ∈ L²_σ`,
`‖ρ_k‖_{L²_σ} = k^{−σ}‖ρ_1‖_{L²_σ}`, `‖ρ_1‖²_{L²_σ} ≤ 1/(2σ) + 1/(2 − 2σ)`, and
`𝓜_σ ρ_k(t) = −k^{−s} ζ(s)/s`. (b) For `σ > 0`: `χ₊ ∈ L²_σ`, `𝓜_σ χ₊(t) = 1/s`.

*Proof.* (a) `ρ_k(x) = ρ_1(kx)`; `|ρ_1| ≤ 1` on `(0,1]` and `ρ_1(x) = 1/x` for `x > 1`, giving integrability and
`‖ρ_1‖² ≤ ∫_0^1 x^{2σ−1} dx + ∫_1^∞ x^{2σ−3} dx`; scaling `x ↦ kx` gives `k^{−2σ}`. With `u = 1/(kx)`:
`∫_0^∞ {1/(kx)} x^{s−1} dx = k^{−s} ∫_0^∞ {u} u^{−s−1} du = −k^{−s} ζ(s)/s` by [C1]. (b) `∫_0^1 x^{s−1} dx = 1/s`. ∎

**Lemma M (Hölder).** Let `1 ≤ p < 2`, `a < 1/p`. For measurable `F` on `(0,1)`:
`‖F‖_{L^p(0,1)} ≤ (1 − r(a − 1/2))^{−1/r} ‖x^{a − 1/2} F‖_{L²(0,1)}`. In particular for `a = 1/p − ε/2`
(`ε > 0`) the constant is `(2/(rε))^{1/r}`.

*Proof.* `|F| = |x^{a−1/2}F| · x^{1/2−a}`; Hölder with `1/p = 1/2 + 1/r`:
`‖F‖_p ≤ ‖x^{a−1/2}F‖_2 ‖x^{1/2−a}‖_r` and `‖x^{1/2−a}‖_r^r = ∫_0^1 x^{r(1/2 − a)} dx = 1/(1 − r(a − 1/2))`,
finite because `r(a − 1/2) < r(1/p − 1/2) = 1`. For `a = 1/p − ε/2`: `1 − r(a − 1/2) = rε/2`. ∎

**Lemma PW (shifting the line).** Let `a < b`, let `Ψ` be analytic on a neighbourhood of
`{a ≤ Re s ≤ b}` with `M_Ψ := sup_{a ≤ σ ≤ b} ∫_ℝ |Ψ(σ+it)|² dt < ∞`. If `F ∈ L²_b` and
`𝓜_b F(t) = Ψ(b + it)` a.e., then `F ∈ L²_a` and `𝓜_a F(t) = Ψ(a + it)` a.e.

*Proof.* For `σ ∈ {a, b}` let `ψ_σ ∈ L²(ℝ)` be the inverse transform of `t ↦ Ψ(σ+it)`. By Fubini,
`∫_ℝ ∫_a^b (|Ψ(σ+iT)|² + |Ψ(σ−iT)|²) dσ dT ≤ 2(b−a) M_Ψ < ∞`, so there are `T_n → ∞` with
`h_n := ∫_a^b (|Ψ(σ+iT_n)|² + |Ψ(σ−iT_n)|²) dσ → 0`. Fix `u ∈ ℝ`; Cauchy's theorem for `Ψ(s)e^{us}` on the
rectangle `[a,b] × [−T_n, T_n]` gives
`|∫_{−T_n}^{T_n} Ψ(a+it) e^{u(a+it)} dt − ∫_{−T_n}^{T_n} Ψ(b+it) e^{u(b+it)} dt| ≤ e^{|u|(|a|+|b|)} (b−a)^{1/2} · 2 h_n^{1/2}`
(Cauchy–Schwarz on each horizontal side). By [C7] the truncated inverse transforms
`(2π)^{−1}∫_{−T_n}^{T_n} Ψ(σ+it)e^{iut}dt` converge to `ψ_σ` in `L²(du)`; along a subsequence a.e. Hence
`e^{ua}ψ_a(u) = e^{ub}ψ_b(u)` for a.e. `u`. By hypothesis `ψ_b(u) = F(e^{−u})e^{−bu}`, so
`F(e^{−u})e^{−au} = ψ_a(u) ∈ L²(du)`, i.e. `F ∈ L²_a` with `𝓜_a F = Ψ(a + i·)`. ∎

## 5. Analytic lemmas

**Lemma Γ.** (a) If `Re w ≥ 1/2`, `w ≠ 1`, then `Re (γ'/γ)(w) ≤ 4 + ½ log(1+|w|)` when `Re w < 1`, and
`Re (γ'/γ)(w) ≤ 5 + ½ log(1+|w|)` when `|Im w| ≥ 1`.
(b) If `0 < ε ≤ 1`, `Re s ≥ 1/2` and every point of the segment `[s, s+ε]` satisfies one of the two conditions in
(a), then `|Γr_ε(s)| ≤ e^{5ε} (2 + |s|)^{ε/2}`; if moreover `Re(s + ε) < 1`, then
`|Γr_ε(s)|² ≤ e^{8ε} (2 + |t|)^{ε}`.

*Proof.* (a) `γ'/γ(w) = 1/w + 1/(w−1) − ½ log π + ½ ψ(w/2)`. Differentiating [C3] under the integral
(dominated by `t e^{−t Re z}/12`, see below) gives, for `Re z > 0`,
`ψ(z) = log z − 1/(2z) − ∫_0^∞ q(t) e^{−tz} dt`, `q(t) := ½ − 1/t + 1/(e^t − 1)`. From [C3] (coth with
`z = t/2`, and `½ coth(t/2) = ½ + 1/(e^t−1)`): `q(t) = Σ_{k≥1} 2t/(t² + 4k²π²) ∈ [0, t/12]`. For `Re z ≥ 1/4`:
`|∫_0^∞ q e^{−tz}| ≤ ∫_0^∞ (t/12) e^{−t/4} dt = 4/3` and `|1/(2z)| ≤ 2`, so `Re ψ(z) ≤ log|z| + 4`.
With `z = w/2`: `½ Re ψ(w/2) ≤ ½ log|w| − ½ log 2 + 2 ≤ ½ log(1+|w|) + 2`. Further `Re(1/w) ≤ 1/Re w ≤ 2`,
`−½ log π < 0`, and `Re 1/(w−1) = (Re w − 1)/|w−1|² ≤ 0` if `Re w < 1`, `≤ 1/|Im w| ≤ 1` if `|Im w| ≥ 1`.
(b) `γ` is analytic and zero-free on `{Re w ≥ 1/2} \ {1}` (its only zero there is `w = 1`; `Γ` has no zeros),
so `log|Γr_ε(s)| = ∫_0^ε Re(γ'/γ)(s+u) du ≤ ε(5 + ½ log(2 + |s|))` (`|s + u| ≤ |s| + ε ≤ |s| + 1`). If
`Re(s+ε) < 1`, use the first bound and `|s + u| ≤ 1 + |t|`: `log|Γr_ε(s)|² ≤ 2ε(4 + ½ log(2+|t|))`. ∎

**Lemma H (Hadamard midline inequality).** Suppose every zero `ρ` of ζ has `β ≤ θ*`. Let `ε > 0`,
`Re s ≥ θ* − ε/2`, `s ≠ 1`, `Re s > 0`. Then `|ξ(s)| ≤ |ξ(s+ε)|`, hence
`|ζ(s)| ≤ |ζ(s+ε)| |Γr_ε(s)|` (right side read as `|ξ(s+ε)|/|γ(s)|` if `s + ε = 1`).

*Proof.* If `ξ(s) = 0` there is nothing to prove. Otherwise `ξ(s+ε) ≠ 0` (`Re(s+ε) > θ*`). By [C2], for `w`
not a zero, `log|ξ(w)| = Re A + Re(Bw) + Σ_ρ (log|1 − w/ρ| + Re(w/ρ))`. Subtracting at `w = s+ε` and `w = s`
and using `Σ_ρ Re(1/ρ) = −B` (absolutely convergent; `ε` real):
`log|ξ(s+ε)| − log|ξ(s)| = Σ_ρ log|(ρ − s − ε)/(ρ − s)|`. Each term is `≥ 0`, because
`|ρ − s − ε|² − |ρ − s|² = (β − σ − ε)² − (β − σ)² = ε(ε − 2(β − σ)) ≥ 0` as `β ≤ θ* ≤ σ + ε/2`.
Finally `γ(s) ≠ 0` (`Re s > 0`, `s ≠ 1`) and `ζ = ξ/γ`. ∎

**Lemma A (Lindelöf bound in a zero-free half-plane).** Let `1/2 ≤ θ < 1`, assume `ζ(s) ≠ 0` for `Re s > θ`.
For `δ ∈ (0, (1−θ)/2]` and `η' > 0` there is `C_A = C_A(θ, δ, η')` such that for `σ ≥ θ + 2δ`:
`|ζ(s)|^{−1} ≤ C_A (2+|t|)^{η'}`, and, if also `|t| ≥ 1`, `|ζ(s)| ≤ C_A (2+|t|)^{η'}`.

*Proof.* By [C8] (`ζ(s̄) = conj ζ(s)`) take `t ≥ 0`. *Small t.* On the compact set
`{θ + 2δ ≤ σ ≤ 3, 0 ≤ t ≤ 6}` the function `1/ζ` is continuous (no zeros; `1/ζ(1) = 0`), and `ζ` is continuous
on its part with `t ≥ 1`; for `σ ≥ 3`, `|ζ|^{±1} ≤ ζ(3)` [C8]. *Large t* (`t ≥ 6`). Let
`Ω := {Re s > θ} \ (θ, 1]`, simply connected, free of zeros and poles; let `log ζ` be the branch on `Ω` equal to
`Σ_p Σ_m p^{−ms}/m` for `Re s > 1`.
(1) Borel–Carathéodory [C4], for any height `t ≥ 4`, with `z_0 = 2 + it`, `R = 2 − θ − δ/2`, `r = 2 − θ − δ`
(`R − r = δ/2`, `r ≤ 3/2`, `R + r ≤ 3`): the disc `|z − z_0| ≤ R` lies in `{Re z ≥ θ + δ/2 ≥ 1/2, Im z ≥ 2} ⊂ Ω`
(take `C_6 ≥ 1`, so `U > 0`), where
`Re log ζ = log|ζ| ≤ log C_6 + ½ log(t + 2) =: U` by [C6]; `|log ζ(z_0)| ≤ log ζ(2) < 1` [C8]. Hence for
`|z − z_0| ≤ r`: `|log ζ(z)| ≤ (6/δ) U + 6/δ`. Since every `σ + it` with `θ + δ ≤ σ ≤ 2` satisfies
`|σ + it − z_0| ≤ r`, and `|log ζ| ≤ log ζ(2) < 1` for `σ ≥ 2`:
(★) `|log ζ(z)| ≤ (C_7/δ) log(Im z + 2)` for all `Re z ≥ θ + δ`, `Im z ≥ 4`, with `C_7` absolute.
(2) Three circles [C5] for `log ζ` with centre `z_1 = 1 + 2δ + it`, `r_1 = δ`, `r_3 = 1 + δ − θ`, and for a target
`σ ∈ [θ + 2δ, 1 + δ]`, `r_2 = 1 + 2δ − σ ∈ [δ, 1 − θ]`. As `r_3 ≤ 3/4`, the disc `|z − z_1| ≤ r_3` lies in
`Ω ∩ {Re z ≥ θ + δ, 4 ≤ Im z ≤ t + 1}`. On `|z − z_1| = r_1`: `Re z ≥ 1 + δ`, so `M(r_1) ≤ log ζ(1+δ) ≤ log(1 + 1/δ)` [C8]. On
`|z − z_1| = r_3`: by (★) (and `|log ζ| < 1` where `Re z ≥ 2`) `M(r_3) ≤ (C_7/δ) log(t+4)`. With
`a(σ) := log(r_2/r_1)/log(r_3/r_1) ≤ a_0 := log((1−θ)/δ)/log((1+δ−θ)/δ) < 1`:
`|log ζ(σ + it)| ≤ M(r_2) ≤ max(1, log(1+1/δ))^{1−a} ((C_7/δ) log(t+4))^{a} ≤ C_9 (log(t+4))^{a_0}` (`C_9 = C_9(θ,δ)`).
For `σ ≥ 1 + δ`, `|log ζ| ≤ log ζ(1+δ)`. By Young's inequality `C_9 y^{a_0} ≤ η' y + C_10` (`C_10 = C_10(C_9, a_0, η')`)
with `y = log(t+4)`, so `|ζ(s)|^{±1} ≤ e^{|log ζ(s)|} ≤ e^{C_10 + log ζ(1+δ)} (t+4)^{η'}`. ∎

**Lemma P (Riesz–Perron).** Assume `ζ(s) ≠ 0` for `Re s > θ` (`1/2 ≤ θ < 1`). Let `η ∈ (0,1)`,
`δ ∈ (0,(1−θ)/2]`, `η' ∈ (0, 1/2]`. If `σ − η ≥ θ + 2δ` then for all `N ≥ 1`
`|P_N(s) − 1/ζ(s)| ≤ C_P N^{−η} (2 + |t|)^{η'}`, `C_P = C_P(θ, δ, η, η')`. (`1/ζ(1) := 0`.) If `σ − η > 1` the
hypothesis on zeros is not needed.

*Proof.* For `y > 0`, `c > 0`: `(2πi)^{−1} ∫_{(c)} y^w dw/(w(w+1)) = (1 − 1/y)_+` (for `y > 1` close left: residues
`1` at `0`, `−1/y` at `−1`; for `y ≤ 1` close right: no poles). Take `c > max(0, 1 − σ)`, `y = N/k`, multiply by
`μ(k)k^{−s}`, sum over `k` (absolutely convergent: `Σ_k k^{−σ−c} N^c ∫ |dw|/|w(w+1)| < ∞`):
`P_N(s) = (2πi)^{−1} ∫_{(c)} ζ(s+w)^{−1} N^w dw/(w(w+1))`. Shift to `Re w = −η`. The integrand is analytic
in `−η ≤ Re w ≤ c` except the simple pole `w = 0` (residue `1/ζ(s)`): `1/ζ(s+w)` is analytic because
`Re(s+w) ≥ σ − η > θ`, and `w = −1` is outside. On `Im w = ±T`, `|1/ζ(s+w)| ≤ C_A(2 + |t| + T)^{η'}` (Lemma A
for `Re(s+w) ≤ 3`, [C8] beyond) and `|w(w+1)| ≥ T²`, so the horizontal integrals vanish as `T → ∞`. Thus
`P_N(s) − 1/ζ(s) = (2πi)^{−1} ∫_{(−η)} ζ(s+w)^{−1} N^w dw/(w(w+1))` and, with `w = −η + iv`,
`|P_N(s) − 1/ζ(s)| ≤ (N^{−η}/2π) ∫_ℝ C_A (2+|t|)^{η'} (1+|v|)^{η'} dv/(|w||w+1|) =: C_P N^{−η}(2+|t|)^{η'}`,
using `2 + |t + v| ≤ (2+|t|)(1+|v|)`; the `v`-integral is finite (`|w||w+1| ≥ max(η,|v|) max(1−η,|v|)`,
`η' < 1`). If `σ − η > 1`, directly: `P_N(s) − 1/ζ(s) = −N^{−1}Σ_{k≤N} μ(k)k^{1−s} − Σ_{k>N} μ(k)k^{−s}`, so
`|P_N(s) − 1/ζ(s)| ≤ N^{−1}Σ_{k≤N} k^{−η} + Σ_{k>N} k^{−1−η} ≤ (1/(1−η) + 1/η) N^{−η}`, with no hypothesis on zeros. ∎

## 6. Proof of Theorem B

Fix `η ∈ (0, 1/p − θ)`; put `σ_0 := (θ + η + 1/p)/2 ∈ (θ + η, 1/p)` and
`δ := min{(σ_0 − η − θ)/2, (1 − θ − η)/2, (1−θ)/2} > 0`, so `σ_0 − η ≥ θ + 2δ`, `1 − η ≥ θ + 2δ`, `σ_0 ≥ θ + 2δ`.
Let `F_N := χ + Σ_{k=2}^N μ(k)(1 − k/N) f_k` and `Φ_N := χ₊ + Σ_{k=1}^N μ(k)(1 − k/N) ρ_k` on `(0,∞)`. Since
`f_1 = 0` and `Σ_k c_k f_k = Σ_k c_k ρ_k − (Σ_k c_k/k) ρ_1` on `(0,1)`:
`F_N = Φ_N|_{(0,1)} − P_N(1) ρ_1|_{(0,1)}`, hence (Lemma M with `a = σ_0`, `‖ρ_1‖_{L^p(0,1)} ≤ 1`)
`‖F_N‖_p ≤ (1 − r(σ_0 − 1/2))^{−1/r} ‖Φ_N‖_{L²_{σ_0}} + |P_N(1)|`.
By Lemma T, `𝓜_{σ_0} Φ_N(t) = (1 − ζ(s)P_N(s))/s = ζ(s)(1/ζ(s) − P_N(s))/s` on `Re s = σ_0`. Lemma A (`η' = 1/8`)
gives `|ζ(σ_0 + it)| ≤ C_A (2+|t|)^{1/8}` for `|t| ≥ 1` and `|ζ(σ_0+it)| ≤ C_11` for `|t| ≤ 1` (continuity,
`σ_0 < 1`); Lemma P (`η' = 1/8`) gives `|1/ζ(s) − P_N(s)| ≤ C_P N^{−η}(2+|t|)^{1/8}`. Therefore
`‖Φ_N‖²_{L²_{σ_0}} ≤ (2π)^{−1} max(C_A, C_11)² C_P² N^{−2η} ∫_ℝ (2+|t|)^{1/2} (σ_0² + t²)^{−1} dt = C_12 N^{−2η}`.
Lemma P at `s = 1` (`1 − η ≥ θ + 2δ`): `|P_N(1)| = |P_N(1) − 1/ζ(1)| ≤ 2^{1/8} C_P N^{−η}`. Hence
`‖F_N‖_p ≤ C N^{−η}`. ∎

## 7. Proof of Theorem A

(ii) ⟹ (i) is Proposition N; (ii) ⟺ (iii) is Proposition C(d). It remains to prove (i) ⟹ (ii) together
with the quantitative bound. Assume (i), fix `ε ∈ (0, ε_0]` and put `σ' := 1/p − ε/2`. The choice of `ε_0`
gives: `σ' ≥ 1/2`; `σ' + ε = 1/p + ε/2 < 1`; `2σ' − 1 ≥ (2/p − 1)/2 > 0`; `1 − ε > σ'`; `ε ≤ 1/2`.

**A1 (the functions).** For `Re ε > 0` (here `ε` real) define `φ_ε(j) := Σ_{d|j} μ(d) d^{−ε}`,
`E_ε(y) := Σ_{j≤y} φ_ε(j) − y/ζ(1+ε)` (`y > 0`), `Φ_ε(x) := χ₊(x) − E_ε(1/x)` (`x > 0`) and
`H_ε(x) := −E_ε(⌊1/x⌋)` (`0 < x < 1`), a function in `S`. For each `x > 0`, `y = 1/x`:
`Σ_{k≥1} μ(k)k^{−ε}{y/k} = Σ_{k≤y} μ(k)k^{−ε}{y/k} + y Σ_{k>y} μ(k)k^{−1−ε}` converges absolutely, and
equals `y/ζ(1+ε) − Σ_{k≤y} μ(k)k^{−ε}⌊y/k⌋ = −E_ε(y)` (`Σ_{k≤y} μ(k)k^{−ε}⌊y/k⌋ = Σ_{j≤y} φ_ε(j)`). Hence
pointwise `Φ_ε = χ₊ + Σ_k μ(k)k^{−ε}ρ_k`, and on `(0,1)`
`H_ε(x) = Φ_ε(x) − 1 − {1/x}/ζ(1+ε)` (check: `−E_ε(y) − (y − ⌊y⌋)/ζ(1+ε) = −E_ε(⌊y⌋)`).
For `N ≥ 1` let `H_{ε,N} := Σ_{k=2}^N μ(k)k^{−ε}(1 − k/N) f_k ∈ V` and
`Φ_{ε,N} := χ₊ + Σ_{k=1}^N μ(k)k^{−ε}(1−k/N)ρ_k`; as in §6, on `(0,1)`:
`H_{ε,N} = Φ_{ε,N} − 1 − P_N(1+ε){1/x}`.

**A2 (transform on a right line).** Let `b ∈ (1 − ε, 1)`. By Lemma T, `Σ_k |μ(k)| k^{−ε}‖ρ_k‖_{L²_b} ≤
‖ρ_1‖_{L²_b} Σ_k k^{−ε−b} < ∞`, so `Σ_k μ(k)k^{−ε}ρ_k` converges in `L²_b`; its sum equals the pointwise sum
a.e., so `Φ_ε ∈ L²_b`, and (unitarity + absolute convergence of `Σ μ(k)k^{−s−ε}`, `Re(s+ε) > 1`)
`𝓜_b Φ_ε(t) = Ψ_ε(b + it)`, where `Ψ_ε(s) := (1 − R_ε(s))/s`, `R_ε(s) := ζ(s)/ζ(s+ε)`.

**A3 (shift to the left).** `Ψ_ε` is analytic on `{1/p − ε < Re s < 1}` (there `Re(s+ε) > 1/p`, so `1/ζ(s+ε)`
is analytic by (i); `ζ(s)` is analytic for `s ≠ 1`), a neighbourhood of the strip `σ' ≤ Re s ≤ b`. In the
strip: Lemma H (`θ* = 1/p`, `Re s ≥ 1/p − ε/2`) gives `|R_ε(s)| ≤ |Γr_ε(s)|`; for `|t| ≥ 1`, Lemma Γ(b) gives
`|Γr_ε(s)| ≤ e^{5/2}(3 + |t|)^{1/4}`; for `|t| ≤ 1`, `Γr_ε` is continuous on the compact set
`{σ' ≤ σ ≤ b, |t| ≤ 1}` (`γ(s) ≠ 0` as `s ≠ 1`), hence bounded by some `C_13(ε)`. So
`sup_{σ' ≤ σ ≤ b} ∫ |Ψ_ε(σ+it)|² dt ≤ ∫_ℝ (1 + max(C_13, e^{5/2})(3+|t|)^{1/4})² (1/4 + t²)^{−1} dt < ∞`.
Lemma PW: `Φ_ε ∈ L²_a` with `𝓜_a Φ_ε = Ψ_ε(a + i·)` for every `a ∈ [σ', b]`.

**A4 (the key estimate `B(ε) := ‖Φ_ε‖²_{L²_{σ'}} ≤ K_p ε`).** On `Re s = σ'`, write `R = R_ε(s)`. Since
`|1 − R|² = 1 − 2 Re R + |R|²` and `|R| ≤ |Γr_ε(s)|` (Lemma H),
`2π B(ε) = ∫_ℝ |1 − R|² |s|^{−2} dt ≤ ∫ |s|^{−2} dt − 2 Re J + ∫ |Γr_ε(s)|² |s|^{−2} dt`, `J := ∫_ℝ R(σ'+it) |s|^{−2} dt`.
(a) `∫_ℝ dt/(σ'² + t²) = π/σ'`.
(b) *Exact evaluation of `J`.* On the line, `|s|² = s(2σ' − s)`, so `J = (1/i) ∫_{(σ')} f(s) ds`,
`f(s) := R(s)/(s(2σ' − s))`. In `Re s ≥ σ'`, `f` is meromorphic with only two simple poles, `s = 1` and
`s = 2σ'` (`2σ' > 1`): `R` is analytic there except at `s = 1` (`ζ(s+ε) ≠ 0` since `Re(s+ε) > 1/p`; at
`s = 1 − ε`, `R = 0`). Integrate over the rectangle `[σ', X] × [−T, T]`, `X > 2σ' + 1`. On the horizontal sides,
`|f| ≤ |Γr_ε(s)|/T² ≤ e^{5ε}(2 + X + T)^{ε/2}/T²` (Lemma H; Lemma Γ(b), `|Im| = T ≥ 1`) `→ 0` uniformly, times
length `X`. On `Re s = X`: `|R(s)| ≤ ζ(X) ζ(X+ε) ≤ ζ(X)²` [C8] and
`∫_ℝ |ds|/(|s||2σ'−s|) ≤ ∫ dt/((X − 2σ')² + t²) = π/(X − 2σ')`, which `→ 0` as `X → ∞`. The integrals over
`Re s = σ'` converge absolutely (`|f| ≤ max(C_13, e^{5/2})(3+|t|)^{1/4}/|s|²` by A3). By the residue theorem,
`∫_{(σ')} f ds = −2πi (Res_{s=1} f + Res_{s=2σ'} f)`, with
`Res_{s=1} f = 1/(ζ(1+ε)(2σ' − 1))` and `Res_{s=2σ'} f = −R(2σ')/(2σ')`. Hence
`J = 2π [ R(2σ')/(2σ') − 1/((2σ' − 1) ζ(1+ε)) ]` (real).
(c) By Lemma Γ(b) (`Re(s+ε) < 1`): `|Γr_ε|² − 1 ≤ e^{8ε}(2+|t|)^ε − 1 ≤ (8ε + ε log(2+|t|)) e^{8ε}(2+|t|)^{ε}`
(`e^x − 1 ≤ x e^x`, `x ≥ 0`), so, with `ε ≤ 1/2`, `σ' ≥ 1/2`:
`∫ |Γr_ε|² |s|^{−2} dt ≤ π/σ' + ε e^4 ∫_ℝ (8 + log(2+|t|))(2+|t|)^{1/2}(1/4 + t²)^{−1} dt = π/σ' + C_Γ ε`.
Adding (a), (b), (c):
`2π B(ε) ≤ (2π/σ')(1 − R(2σ')) + 4π/((2σ' − 1) ζ(1+ε)) + C_Γ ε`.
Now `R(2σ') = ζ(2σ')/ζ(2σ' + ε) ≥ 1` (`ζ` decreasing on `(1,∞)`), so the first term is `≤ 0`;
`1/ζ(1+ε) ≤ ε` [C8] and `2σ' − 1 ≥ (2 − p)/(2p)` (choice of `ε_0`). Therefore
`B(ε) ≤ (4p/(2 − p)) ε + (C_Γ/2π) ε = K_p ε`.

**A5 (`H_ε` approximates `−χ`).** By A1, on `(0,1)`: `χ + H_ε = Φ_ε − {1/x}/ζ(1+ε)`. Lemma M with
`a = σ'` and A4:
`‖χ + H_ε‖_{L^p(0,1)} ≤ (2/(rε))^{1/r} ‖x^{σ'−1/2}Φ_ε‖_{L²(0,1)} + ε ≤ (2/(rε))^{1/r} (K_p ε)^{1/2} + ε
= A_p ε^{1/2 − 1/r} + ε = A_p ε^{1 − 1/p} + ε`.

**A6 (`H_{ε,N} → H_ε`).** Put `σ'' := 1/p − ε/4 ∈ (σ', 1/p)`. By A3 and Lemma T,
`Φ_ε − Φ_{ε,N} ∈ L²_{σ''}` with transform `(ζ(s)/s)(P_N(s+ε) − 1/ζ(s+ε))` on `Re s = σ''` (the transform of
the finite sum `Φ_{ε,N}` is `(1 − ζ(s)P_N(s+ε))/s`). On this line:
* `|ζ(s)| ≤ |ζ(s + ε/2)| |Γr_{ε/2}(s)|` (Lemma H with `ε/2`: `σ'' = 1/p − (ε/2)/2`);
  `Re(s + ε/2) = 1/p + ε/4 = θ + 2δ` with `θ := 1/p`, `δ := ε/8 ≤ (1−θ)/2`; Lemma A (`η' = 1/8`):
  `|ζ(s + ε/2)| ≤ C_A (2+|t|)^{1/8}` for `|t| ≥ 1`, and `≤ C_14(ε)` for `|t| ≤ 1` (continuity; `1/p + ε/4 < 1`);
  Lemma Γ(b): `|Γr_{ε/2}(s)| ≤ e^{2}(3 + |t|)^{1/8}`.
* Lemma P (`σ = Re(s+ε) = 1/p + 3ε/4`, shift `η = ε/2`, `σ − η = θ + 2δ`, `η' = 1/8`):
  `|P_N(s+ε) − 1/ζ(s+ε)| ≤ C_P N^{−ε/2} (2+|t|)^{1/8}`.
Hence `‖Φ_ε − Φ_{ε,N}‖²_{L²_{σ''}} ≤ C_15(p,ε) N^{−ε} ∫_ℝ (3+|t|)^{3/4} (1/4 + t²)^{−1} dt`. Lemma P at
`s = 1 + ε` (`η = ε/2`, `σ − η > 1`): `|1/ζ(1+ε) − P_N(1+ε)| ≤ (1/(1−ε/2) + 2/ε) N^{−ε/2}`. With A1 and Lemma M
(`a = σ''`, constant `(4/(rε))^{1/r}`):
`‖H_ε − H_{ε,N}‖_{L^p(0,1)} ≤ (4/(rε))^{1/r} ‖Φ_ε − Φ_{ε,N}‖_{L²_{σ''}} + |1/ζ(1+ε) − P_N(1+ε)| ≤ C_17(p,ε) N^{−ε/2}`.

**A7 (conclusion).** Choose `N_0(p,ε)` with `C_17 N_0^{−ε/2} ≤ ε`. For `N ≥ N_0`, A5 and A6 give
`‖χ + H_{ε,N}‖_p ≤ A_p ε^{1−1/p} + 2ε`. Since `H_{ε,N} ∈ V`, `d_{N,p} ≤ A_p ε^{1−1/p} + 2ε` for `N ≥ N_0(p,ε)`;
letting `ε → 0` gives `d_{N,p} → 0`, i.e. (ii). ∎

## 8. Proof of the corollaries

1. [C9] gives (i) for every `p ≤ 8/7`, so Theorem A applies for `1 < p ≤ 8/7`; Theorem B with `θ = 7/8`
   applies for `1 < p < 8/7`.
2. (i) `⟺ Θ ≤ 1/p ⟺ p ≤ 1/Θ`. For `p > 2` use Proposition N with a zero on `Re s = 1/2 > 1/p`.
3. RH `⟺ Θ = 1/2 ⟺` (i) for all `p ∈ (1,2)`; the rate is Theorem B with `θ = 1/2`. ∎

## 9. Remarks

* *What is new and what is routine.* Dimitrov–Oliveira (arXiv 1608.07887, 2016) prove, for Dirichlet
  `L`-functions and `p ∈ (1,2]`: `L ≠ 0` on `Re s > 1/p` `⟺ inf_{A_n} ∫_{Re s = 1/p} |1 − L A_n|² |s|^{−2}|ds| → 0`
  (a Hilbert-space, i.e. weighted-`L²`, natural-number criterion). Lemma M turns their sufficiency, applied on the
  line `Re s = θ`, into Theorem A(i ⟹ ii) for `p < 1/θ` in one line, so the open range is routine (Theorem B
  adds explicit coefficients and the power rate). Lemma M cannot reach the endpoint `p = 1/θ` from their
  endpoint statement (weighted `L²` on the line `1/p` only controls the Lorentz norm `L^{p,2} ⊋ L^p`). The
  endpoint needs A4: on the line `1/p − ε/2`, symmetric to `1/p + ε/2` about the line `Re s = 1/p` that bounds
  the zeros, the Hadamard product gives `|ζ(s)| ≤ |ζ(s+ε)Γr_ε(s)|`, and the exact residue evaluation of
  `∫ Re R / |s|²` turns this into `B(ε) = O(ε)`, which beats the Hölder loss `ε^{−1/r}` exactly because `p > 1`.
* *Constants.* `C_Γ = 5993.05`, `K_{8/7} = 959.16`, `A_{8/7} = 27.80`, so the bound in Theorem A is
  qualitative. The illustration in REPORT.md shows `‖χ + H_ε‖_{1.1} ≈ (0.9–1.15) ε` numerically (the true zeros
  sit at `β = 1/2 < 1/p`).
* *p = 2.* The same argument (with `2σ' = 1 − ε < 1`, `R(1−ε) = 0`, and Lemma Γ constants adjusted for
  `Re w ≥ 1/4`) gives, under RH, `‖χ + H_ε‖_2² ≤ (4ε/(1−ε²)) + O(ε)`: a short proof of Báez-Duarte's theorem
  with rate `ε^{1/2}`. Not needed here; not checked in detail.
* *Not covered.* The analogous converse in Noor's Hardy-space model `H^p(𝔻)` (Ghosh–Kremnizer–Noor–Santos,
  Adv. Math. 455 (2024), "It is unclear whether the converse … holds") is a different space and is not
  addressed. Weighted `ℓ^p_α` spaces (Gallardo-Gutiérrez–Seco) with `α ≠ −2` are likely accessible by the
  same method (weight `x^{−α−2}` in Lemma M); not done.
