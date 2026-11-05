# Natural-number Nyman–Beurling–Báez-Duarte criterion in L^p(0,1), 1 < p < 2, for Dedekind zeta functions

Node N5. Author: Orca agent `nflp`, 2026-11-05. Paper proof (no Lean). It adapts [NBLp] = `Riemann/Deep/NBLp/PROOF.md`
(the case K = ℚ) and [NBLp2] = `Riemann/Deep/NBLp2/PROOF.md` lemma by lemma; each item is marked **unchanged**,
**changed (how)** or **new**, with its K-dependence. Inputs not proved here are quoted in §2.

## 0. Notation

* K is a number field of degree `n = r_1 + 2r_2`, discriminant `d_K`, ring of integers `O_K`. `𝔞, 𝔟, 𝔠` run over non-zero
  integral ideals, `N𝔞 = |O_K/𝔞|`, `μ_K` is the Möbius function of the ideal monoid. `K = ℚ` is exactly [NBLp].
* `a_K(m) := #{𝔞 : N𝔞 = m}`, `A_K(y) := Σ_{m≤y} a_K(m) = #{𝔞 : N𝔞 ≤ y}` (so `A_K(y) = 0` for `y < 1`);
  `b_K(m) := Σ_{N𝔞=m} μ_K(𝔞)`, so `|b_K(m)| ≤ a_K(m)` and `b_K(m) = 0` unless `m ∈ 𝒩_K := {N𝔞}` (a multiplicative monoid).
* `ζ_K(s) = Σ_𝔞 N𝔞^{−s} = Σ_m a_K(m) m^{−s}`, `1/ζ_K(s) = Σ_m b_K(m) m^{−s}` (`σ > 1`); `c_K := Res_{s=1} ζ_K(s) > 0`.
* **The functions.** `κ_K(y) := c_K y − A_K(y)` (`y ≥ 0`); this is de Roton's complementary function
  `Ψ_F(y) = Res_{s=1}(y^s F(s)/s) − Σ_{n≤y} a_F(n)` for `F = ζ_K` [dR09, (1)]. For `K = ℚ`: `κ(y) = {y}`.
  For integers `k ≥ 1` (real dilations `θ ≥ 1` occur only in Remark 9.1), `x > 0`: `r_k(x) := κ_K(1/(kx))`, `λ_k := r_k − k^{−1} r_1`; for an ideal 𝔟: `λ_𝔟 := λ_{N𝔟}`.
  Since `A_K = 0` on `[0,1)`: `r_k(x) = c_K/(kx)` for `x > 1/k`, hence **`λ_k = 0` on `(1,∞)`**; `λ_1 = 0`.
  For `K = ℚ`, `λ_k = f_k = ρ_{1/k}` of [NBLp].
* `α_K := 1 − 2/(n+1) ∈ [0,1)`; by [K2], `|κ_K(y)| ≤ C_K y^{α_K}` for `y ≥ 1` (only `α_K < 1` is used).
* `𝟙 := 1_{(0,1)}` on `(0,1)`, `𝟙₊ := 1_{(0,1]}` on `(0,∞)`. `I_m := (1/(m+1), 1/m]`; `S_p`, `r := 2p/(2−p)`, `L²_σ := L²((0,∞), x^{2σ−1}dx)`,
  `𝓜_σ` as in [NBLp].
* **Membership condition** `(M_p)`: `κ_K(1/x) ∈ L^p(0,1)`, i.e. `∫_1^∞ |κ_K(u)|^p u^{−2} du < ∞`. Under `(M_p)` every `λ_k|_{(0,1)}`
  lies in `L^p(0,1)` (`r_k(x) = r_1(kx)`, and `r_1 = c_K/x` is bounded on `[1,k]`). For `K = ℚ` it is trivial (`0 ≤ κ < 1`).
* `V_K := span{λ_k|_{(0,1)} : k ∈ ℕ, k ≥ 2}`; `V_K^id := span{λ_𝔟 : 𝔟 ≠ O_K} = span{λ_k : k ∈ 𝒩_K, k ≥ 2}`. Under `(M_p)`:
  `W_{K,p}, W^id_{K,p}` := their closures in `L^p(0,1)`; `d_{N,p}(K) := inf_{c_k} ‖𝟙 − Σ_{2≤k≤N} c_k λ_k‖_p`, and `d^id_{N,p}(K)`
  the same with `k ∈ 𝒩_K`.
* `Γ_ℝ(s) := π^{−s/2}Γ(s/2)`, `Γ_ℂ(s) := 2(2π)^{−s}Γ(s)`;
  `γ_K(s) := s(s−1)|d_K|^{s/2} Γ_ℝ(s)^{r_1} Γ_ℂ(s)^{r_2}`, `ξ_K := γ_K ζ_K`; for `ε > 0`: `Γr_ε(s) := γ_K(s+ε)/γ_K(s)`.
  `c_1(K) := 3 + ½ log|d_K| + 2n` (`c_1(ℚ) = 5`).
* `R_ε(s) := ζ_K(s)/ζ_K(s+ε)`; `P_N(s) := Σ_{k≤N} b_K(k)(1 − k/N) k^{−s} = Σ_{N𝔟≤N} μ_K(𝔟)(1 − N𝔟/N) N𝔟^{−s}`; `1/ζ_K(1) := 0`.
* `C(K) := sup_{0<ε≤1/2} c_K/(ε ζ_K(1+ε)) < ∞` (continuous and positive on `(0, 1/2]`, limit 1 at 0 since
  `ζ_K(1+ε) ~ c_K/ε`; `C(ℚ) ≤ 1`).
* `Θ_K := sup{Re ρ : ζ_K(ρ) = 0, Re ρ > 0} ∈ [1/2, 1]` (`≤ 31/32` by [CT]; `[K1]` alone gives only `≤ 1`).

## 1. Results

**Proposition C_K (structure)** — changed (values on `I_m`; `(M_p)` needed for `V_K ⊂ L^p`).
(a) For `k, m ≥ 1`, `x ∈ I_m`: `λ_k(x) = A_K(m)/k − A_K(⌊m/k⌋)`. For an ideal 𝔟:
`λ_𝔟(x) = A_K(m)/N𝔟 − #{𝔠 ⊆ 𝔟 : N𝔠 ≤ m}` (expected minus actual number of multiples of 𝔟 of norm ≤ m; for `K = ℚ`: `{m/k}`).
(b) Each `λ_k` is constant on every `I_m`; under `(M_p)`, `V^id_K ⊂ V_K ⊂ S_p`, a closed subspace of infinite codimension.
(c) `D_m λ_k = λ_{km} − k^{−1} λ_m` a.e. on `(0,1)` (`D_m` as in [NBLp]); so `D_m V_K ⊂ V_K` for all `m ≥ 1` and
`D_m V^id_K ⊂ V^id_K` for `m ∈ 𝒩_K`.
(d) Under `(M_p)`: `𝟙 ∈ W_{K,p} ⟺ W_{K,p} = S_p ⟺ d_{N,p}(K) → 0`.

**Proposition M_K (membership)** — new. (a) `(M_p)` holds for `1 ≤ p < 1/α_K = (n+1)/(n−1)` (every `p` if `n = 1`).
(b) If `ζ_K ≠ 0` on `Re s > θ` (`1/2 ≤ θ < 1`), then `(M_p)` holds for every `1 ≤ p < 1/θ`.
(c) If `1 < p < 2` and `ζ_K ≠ 0` on `Re s > 1/p`, then `(M_p)` holds (the endpoint).
(d) (de Roton's claim [dR07, §5], via [dR1]: `Ψ^{(1)}_F ∈ L²(0,∞)` for `F` in the Selberg class of degree `< 4`; not re-checked
here.) For `n ≤ 3`, `(M_p)` holds for every `p ≤ 2`.

**Theorem A_K.** Let `1 < p < 2`. The following are equivalent:
(i) `ζ_K(s) ≠ 0` for `Re s > 1/p` (i.e. `Θ_K ≤ 1/p`);
(ii) `(M_p)` holds and `𝟙 ∈ W^id_{K,p}`: the indicator of `(0,1)` is an `L^p`-limit of combinations of the ideal dilations `λ_𝔟`;
(iii) `(M_p)` holds and `𝟙 ∈ W_{K,p}` (all natural-number dilations), i.e. `d_{N,p}(K) → 0`;
(iv) `(M_p)` holds and `W_{K,p} = S_p`.
Moreover, under (i), with `ε_0(p,K) := min{1/2, (2/p − 1)/2, 1 − 1/p, 1/(2n)}`, for every `ε ∈ (0, ε_0]` there is `N_0(p,K,ε)`
such that for `N ≥ N_0`
`‖ 𝟙 + Σ_{N𝔟≤N} μ_K(𝔟) N𝔟^{−ε} (1 − N𝔟/N) λ_𝔟 ‖_{L^p(0,1)} ≤ A_{p,K} ε^{1−1/p} + (1 + τ_{p,K}) ε`,
`A_{p,K} := K_{p,K}^{1/2}(2/r)^{1/r}`, `K_{p,K} := 4pC(K)/(2−p) + C_{Γ,K}/(2π)` (A4), `τ_{p,K} := C(K) ‖κ_K(1/x)‖_{L^p(0,1)}/c_K`.
For `K = ℚ` this is [NBLp] Theorem A.

**Theorem B_K (open range, explicit Riesz–Möbius approximants, power rate).** Let `1/2 ≤ θ < 1`, `ζ_K(s) ≠ 0` for `Re s > θ`,
`1 < p < 1/θ`, `0 < η < 1/p − θ`. Then `(M_p)` holds and
`‖ 𝟙 + Σ_{N𝔟≤N} μ_K(𝔟)(1 − N𝔟/N) λ_𝔟 ‖_{L^p(0,1)} ≤ C(K,p,θ,η) N^{−η}` for all `N ≥ 2`.
On `I_m` the approximant equals `1 + P_N(1)A_K(m) − Σ_{N𝔟≤N} μ_K(𝔟)(1 − N𝔟/N) A_K(m/N𝔟)`.

**Proposition N_K (necessity; every `1 < p < ∞`)** — unchanged argument, plus `(M_p)`. If `(M_p)` holds and `ζ_K(w_0) = 0` with
`Re w_0 > 1/p`, then `𝟙 ∉ W_{K,p}` (a fortiori `𝟙 ∉ W^id_{K,p}`).

**Corollaries.**
1. (Cureton–Thorner, Lean-checked, unrefereed [CT].) For **every** number field K and every `1 < p ≤ 32/31`: `(M_p)` holds and
   `𝟙 ∈ W^id_{K,p}`. For `1 < p < 32/31` and `η < 1/p − 31/32`: `‖𝟙 + Σ_{N𝔟≤N} μ_K(𝔟)(1 − N𝔟/N) λ_𝔟‖_p ≪_{K,p,η} N^{−η}`.
2. (Conditional.) If Liu's half-plane holds [L2] (`ζ_K = L_K(s, 1) ≠ 0` on `Re s > B_2 = (777 − √1081)/852 = 0.8733819666`, every K;
   unrefereed, audited sound given 2610.26835v2): the same for every K and every `1 < p ≤ 1/B_2 = 1.1449744…`, with the rate for
   `p < 1/B_2`, `η < 1/p − B_2`. If Hu's 7/8 [Hu] holds (unrefereed, not Lean-checked): every `p ≤ 8/7`.
3. For `1 < p < 2`: (iii) `⟺ p ≤ 1/Θ_K`. For `p > 2`, (iii) fails for every K (unconditionally).
4. GRH for `ζ_K` ⟺ (ii) holds for every `p ∈ (1,2)`; under GRH the rate `N^{−η}`, `η < 1/p − 1/2`, holds for every `p < 2`.
5. If `K/ℚ` is abelian, `ζ_K = Π_χ L(s,χ)` (`χ` over the primitive characters attached to K), so (i) for `ζ_K` is the conjunction
   of [NBLp2]'s condition (i) for these `χ`; for non-abelian K, `ζ_K` involves non-abelian Artin L-functions and no reduction to
   Dirichlet L-functions is available.

## 2. Cited results

[C1]–[C10] as in [NBLp §2] (in particular [C3] Binet/coth, [C4] Borel–Carathéodory, [C5] three circles, [C7] Plancherel).
K-specific inputs:
* **[K1]** (Hecke; Neukirch, *Algebraic Number Theory*, VII Cor. 5.11 and (5.10); Iwaniec–Kowalski [IK] §5.10 and Thm 5.6.)
  `ζ_K` continues meromorphically to ℂ with a single, simple pole at `s = 1`, residue `c_K = 2^{r_1}(2π)^{r_2} h_K R_K/(w_K |d_K|^{1/2}) > 0`.
  `Λ_K(s) := |d_K|^{s/2} Γ_ℝ(s)^{r_1} Γ_ℂ(s)^{r_2} ζ_K(s)` satisfies `Λ_K(1 − s) = Λ_K(s)` and has simple poles at `s = 0, 1` only, so
  `ξ_K(s) = s(s−1)Λ_K(s) = γ_K(s)ζ_K(s)` is entire of order 1 and real on ℝ. Its zeros are exactly the non-trivial zeros `ρ = β + iγ`
  of `ζ_K`, all with `0 < β < 1` (`ζ_K ≠ 0` on `Re s ≥ 1`); they are infinitely many and symmetric under `ρ ↦ 1 − ρ̄` and `ρ ↦ ρ̄`;
  `Σ_ρ |ρ|^{−2} < ∞`; `ξ_K(s) = e^{A+Bs} Π_ρ (1 − s/ρ) e^{s/ρ}` with `Re B = −Σ_ρ Re(1/ρ)` ([IK] Thm 5.6; Stark, Invent. Math. 23
  (1974) 135–152, for `ζ_K`). On `Re s > 0`, `γ_K` is analytic with the single (simple) zero `s = 1` (`Γ` has no zeros, no poles there).
* **[K2]** (Ideal theorem with remainder; Landau, *Einführung in die elementare und analytische Theorie der algebraischen Zahlen und
  der Ideale* (1918); the weaker exponent `1 − 1/n` in Lang, *Algebraic Number Theory*, VI §3 Thm 3.) For `y ≥ 1`:
  `A_K(y) = c_K y + O_K(y^{1−2/(n+1)})`. Only "some exponent `α_K < 1`" is used.
* **[K3]** (Elementary, `σ > 1`.) `ζ_K(s) = Π_𝔭 (1 − N𝔭^{−s})^{−1}`, `log ζ_K(s) = Σ_𝔭 Σ_{m≥1} N𝔭^{−ms}/m`; `|ζ_K(s)|, |1/ζ_K(s)| ≤ ζ_K(σ)`,
  `|log ζ_K(s)| ≤ log ζ_K(σ) ≤ n log ζ(σ)` (each rational `p` has at most `n` prime factors 𝔭, each with `N𝔭 ≥ p`); `ζ_K` is positive and
  decreasing on `(1,∞)`; `ζ_K(s̄) = conj ζ_K(s)`. Abel summation: `ζ_K(s) = s ∫_1^∞ A_K(u) u^{−s−1} du`.
* **[K4]** (Convexity; Rademacher, Math. Z. 72 (1959) 192–204, Thm 4: for `0 < η ≤ 1/2`, `−η ≤ σ ≤ 1 + η`,
  `|ζ_K(s)| ≤ 3 |(1+s)/(1−s)| (|d_K| (|1+s|/2π)^n)^{(1+η−σ)/2} ζ(1+η)^n`.) Consequence used (take `η = 1/2`; `|ζ_K| ≤ ζ(3/2)^n` for
  `σ ≥ 3/2`): there is `C_6(K) ≥ 1` with `|ζ_K(σ+it)| ≤ C_6(K)(|t|+2)^{n/2}` for `σ ≥ 1/2`, `|t| ≥ 1`. Any polynomial bound would do.
* **[CT]** (Cureton–Thorner, Hexagon 2610.00200v2, formalized in Lean as AxiomMath/DedekindZetaZeroFree; see
  `Riemann/Deep/LiuHecke/NOTE.md` §5 and `COROLLARIES.md`.) For every number field K, `ζ_K(s) ≠ 0` for `Re s > 31/32` (also every
  finite-order Hecke L-function). Unrefereed.
* **[L2]** (Liu, Hexagon 2611.00082, on top of arXiv 2610.26835v2; `COROLLARIES.md`, Artin line.) Every finite-order Hecke L-function
  of every number field is zero-free on `Re s > B_2 = (777 − √1081)/852 = 0.8733819666`. `ζ_K` is the Hecke L-function of the trivial
  character of K, so no Brauer induction is needed for `ζ_K`. Unrefereed; our audit (LiuHecke) found the field step sound given
  2610.26835v2.
* **[Hu]** (Hu, Hexagon 2610.00216v1.) 7/8 for every unitary Hecke character of every number field. Unrefereed, not Lean-checked.
* **[dR07]** (A. de Roton, *Généralisation du critère de Beurling–Nyman pour l'hypothèse de Riemann*, Trans. AMS 359 (2007) 6079–6110,
  HAL hal-00091952, Thm III.) For `F` in the Selberg class: GRH for `F` ⟺ `Ψ^{(1)}_F ∈ L²(0,∞)` and `χ` lies in the `L²(0,1)`-closure of
  `B_F = span{Ψ_F(θ_k/t)}`, `0 < θ_k ≤ 1`. §5: "pour une fonction F de la classe de Selberg de degré inférieur à 4, la fonction
  `Ψ^{(1)}_F` est de carré intégrable" (proved in [dR1] = de Roton, Acta Arith. 126 (2007) 27–55).
* **[dR09]** (A. de Roton, *Une approche séquentielle de l'hypothèse de Riemann généralisée*, J. Number Theory 129 (2009) 2647–2658,
  HAL hal-00091966, Thm I and Cor. 1.) For `F ∈ S` satisfying GRH, `−Σ_{n≤N} a^{−1}(n) e^{−a log n/log log N} Ψ_F(1/(nx))` → `χ` in
  `L²(0,∞)` at rate `O((log log N)^{−1/3})`; hence GRH for `F` ⟺ `Ψ^{(1)}_F ∈ L²(0,∞)` and `χ ∈ \overline{B_F^{nat}}`,
  `B^{nat}_F = span{Ψ_F(1/(kt)) : k ∈ ℕ*}`. This is the `p = 2` case of the present question (in `L²(0,∞)`, without side condition).

## 3. Structure

**Proof of Proposition C_K** — changed in (a); (b)–(d) as in [NBLp].
(a) `y := 1/x ∈ [m, m+1)`. Since every `N𝔞` is an integer, `A_K(y/k) = #{𝔞 : kN𝔞 ≤ y} = #{𝔞 : kN𝔞 ≤ m} = A_K(⌊m/k⌋)` and
`A_K(y) = A_K(m)`, so `λ_k(x) = c_K y/k − A_K(⌊m/k⌋) − k^{−1}(c_K y − A_K(m)) = A_K(m)/k − A_K(⌊m/k⌋)`. For `k = N𝔟`:
`#{𝔞 : N(𝔞𝔟) ≤ m} = #{𝔠 ⊆ 𝔟 : N𝔠 ≤ m}` (unique factorization: `𝔠 ⊆ 𝔟 ⟺ 𝔠 = 𝔞𝔟`).
(b) Constancy is (a). Under `(M_p)` each `λ_k ∈ L^p(0,1)`, so `V_K ⊂ S_p`; closedness and infinite codimension of `S_p`: [NBLp] C(b).
(c) `r_{km}(x) = r_k(mx)` for all `x > 0`, so `λ_{km}(x) − k^{−1}λ_m(x) = r_k(mx) − (km)^{−1}r_1(x) − k^{−1}r_1(mx) + (km)^{−1}r_1(x) = λ_k(mx)`.
For `x < 1/m` this is `(D_mλ_k)(x)`; for `x > 1/m` it is `λ_k(mx) = 0` (`mx > 1`). If `m ∈ 𝒩_K` and `k ∈ 𝒩_K` then `km ∈ 𝒩_K`.
(d) As in [NBLp] C(d) (it uses `D_m` for every `m ≥ 1`, hence the full family `V_K`). ∎

(For `V^id_K` only `D_m`, `m ∈ 𝒩_K`, is available, and `1_{I_m}` need not be reachable when `m + 1 ∉ 𝒩_K`; Theorem A_K shows
nevertheless that `𝟙 ∈ W^id_{K,p} ⟺ 𝟙 ∈ W_{K,p}`.)

**Proof of Proposition N_K** — unchanged argument; the Mellin identity now needs analytic continuation.
Let `w_0` be a zero with `Re w_0 > 1/p`; by [K1], `Re w_0 < 1`. Put `g(x) := x^{w_0−1} ∈ L^{p'}(0,1)` (`p' = p/(p−1)`;
`p'(Re w_0 − 1) > −1 ⟺ Re w_0 > 1/p`). Since `λ_k = 0` on `(1,∞)`, for `α_K < Re w < 1` Lemma T_K gives (absolutely convergent)
`∫_0^1 λ_k x^{w−1} dx = −ζ_K(w)(k^{−w} − k^{−1})/w`. The left side is analytic on `Re w > 1/p` because `λ_k ∈ L^p(0,1)` by `(M_p)`
(locally uniform domination by `‖λ_k‖_p ‖x^{w−1}‖_{p'}`); the right side is analytic on `Re w > 0` (at `w = 1` the pole of `ζ_K` is cancelled
by `k^{−w} − k^{−1}`). They agree on the non-empty strip `max(α_K, 1/p) < Re w < 1`, hence on the connected half-plane `Re w > 1/p`.
At `w_0`: `∫_0^1 λ_k g = 0` for all `k`, so `g` annihilates `W_{K,p}`, while `∫_0^1 𝟙 g = 1/w_0 ≠ 0`. ∎

## 4. Mellin–Plancherel toolkit

**Lemma T_K** — changed (strip `α_K < σ < 1` instead of `0 < σ < 1`).
(a) For `α_K < σ < 1` and `k ≥ 1`: `x^{σ−1}r_k ∈ L¹(0,∞)`, `r_k ∈ L²_σ`, `‖r_k‖_{L²_σ} = k^{−σ}‖r_1‖_{L²_σ}`,
`‖r_1‖²_{L²_σ} ≤ C_K²/(2(σ − α_K)) + c_K²/(2 − 2σ)`, and `𝓜_σ r_k(t) = −k^{−s} ζ_K(s)/s`. (b) `𝓜_σ 𝟙₊ = 1/s` (`σ > 0`) — unchanged.
(c) If `r_1 ∈ L²_a` with `𝓜_a r_1 = −ζ_K(s)/s` for some `0 < a < 1`, then the same holds for every `r_k` with `k^{−s}` and `k^{−a}` (scaling).

*Proof.* (a) `r_1(x) = c_K/x` for `x > 1` and `|r_1(x)| = |κ_K(1/x)| ≤ C_K x^{−α_K}` for `0 < x ≤ 1` [K2]: this gives integrability and the
bound; `r_k(x) = r_1(kx)` gives the scaling. With `u = 1/(kx)`: `∫_0^∞ r_k x^{s−1} dx = k^{−s} ∫_0^∞ κ_K(u) u^{−s−1} du`. Here
`∫_0^1 c_K u^{−s} du = c_K/(1 − s)` (`σ < 1`), and `∫_1^∞ κ_K(u) u^{−s−1} du` is analytic on `σ > α_K` and equals
`c_K/(s−1) − ζ_K(s)/s` for `σ > 1` ([K3], Abel summation). The right side is analytic on `σ > α_K` (at `s = 1`, `ζ_K(s)/s` has residue
`c_K`), so the identity holds on `σ > α_K`; adding, `∫_0^∞ κ_K(u) u^{−s−1} du = −ζ_K(s)/s` for `α_K < σ < 1`. (c) Scaling. ∎

**Lemma M (Hölder)** and **Lemma PW (shifting the line)** — unchanged ([NBLp §4]).

**Lemma MEM_K (lines left of `α_K`)** — new.
(a) If `ζ_K ≠ 0` on `Re s > θ` (`1/2 ≤ θ < 1`), then for every `a ∈ (θ, 1)`: `r_1 ∈ L²_a` and `𝓜_a r_1 = −ζ_K(s)/s`.
(b) If `1 < p < 2`, `ζ_K ≠ 0` on `Re s > 1/p` and `0 < ε ≤ ε_0(p,K)`, the same holds for every `a ∈ [1/p − ε/2, 1)`.
In both cases Lemma T_K(c) extends to all `r_k` on these lines.

*Proof.* Fix `b ∈ (max(α_K, a), 1)`. By Lemma T_K, `r_1 ∈ L²_b` with `𝓜_b r_1 = Ψ(b + i·)`, `Ψ(s) := −ζ_K(s)/s`, which is analytic on a
neighbourhood of the strip `a ≤ Re s ≤ b` (`b < 1`). By Lemma PW it suffices that `sup_{a≤σ≤b} ∫_ℝ |ζ_K(σ+it)|² |s|^{−2} dt < ∞`. On
`{a ≤ σ ≤ b, |t| ≤ 1}`, `ζ_K` is continuous (no pole). For `|t| ≥ 1`:
(a) with `δ := min{(a − θ)/2, (1 − θ)/2}` every `σ ≥ a` satisfies `σ ≥ θ + 2δ`, and Lemma A_K (`η' = 1/8`) gives `|ζ_K| ≤ C_A(2 + |t|)^{1/8}`.
(b) For `σ ≥ 1/p + ε/2` use Lemma A_K with `θ = 1/p`, `δ = ε/4 ≤ (1 − 1/p)/2`. For `1/p − ε/2 ≤ σ < 1/p + ε/2` Lemma H_K (`θ* = 1/p`)
gives `|ζ_K(s)| ≤ |ζ_K(s+ε)| |Γr_ε(s)|` with `Re(s+ε) ≥ 1/p + ε/2`, so `|ζ_K(s+ε)| ≤ C_A(2+|t|)^{1/8}`, and Lemma Γ_K(b) (every point
of `[s, s+ε]` has real part `≥ 1/2` and imaginary part `t`, `|t| ≥ 1`) gives `|Γr_ε(s)| ≤ e^{c_1 ε}(2 + |s|)^{nε/2} ≤ e^{c_1/2}(3 + |t|)^{1/4}`
(`nε ≤ 1/2`, `|s| ≤ 1 + |t|`). In both cases `|ζ_K(σ+it)|² ≤ C (3 + |t|)^{3/4}`, integrable against `|s|^{−2}`. ∎

**Proof of Proposition M_K.** (a) By [K2], `∫_1^∞ |κ_K(u)|^p u^{−2} du ≤ C_K^p ∫_1^∞ u^{pα_K − 2} du < ∞` iff `pα_K < 1`.
(b) For `p < 1/θ` pick `a ∈ (θ, 1/p)`; Lemma MEM_K(a) and Lemma M give `‖r_1‖_{L^p(0,1)} ≤ (1 − r(a − ½))^{−1/r} ‖r_1‖_{L²_a} < ∞`
(for `p ≤ 1` use `L^{p_1} ⊂ L^p` on `(0,1)` with some `p_1 ∈ (1, 1/θ)`). (c) The same with `a = 1/p − ε_0/2` from Lemma MEM_K(b). (d) is
quoted. ∎

## 5. Analytic lemmas

**Lemma Γ_K** — changed: degree `n` and `|d_K|` enter. (a) If `Re w ≥ 1/2` and (`Re w < 1` or `|Im w| ≥ 1`), then
`Re(γ_K'/γ_K)(w) ≤ c_1(K) + (n/2) log(1 + |w|)`, and `≤ c_1(K) − 1 + (n/2) log(1 + |w|)` if `Re w < 1`.
(b) If `0 < ε ≤ 1`, `Re s ≥ 1/2` and every point of `[s, s+ε]` satisfies the condition in (a), then
`|Γr_ε(s)| ≤ e^{c_1 ε}(2 + |s|)^{nε/2}`; if moreover `Re(s + ε) < 1`, then `|Γr_ε(s)|² ≤ e^{2(c_1−1)ε}(2 + |t|)^{nε}`.
(For `K = ℚ`: `e^{5ε}(2+|s|)^{ε/2}` and `e^{8ε}(2+|t|)^ε`, as in [NBLp].)

*Proof.* `(γ_K'/γ_K)(w) = 1/w + 1/(w−1) + ½ log|d_K| + (r_1/2)(ψ(w/2) − log π) + r_2(ψ(w) − log 2π)`. By [NBLp] Lemma Γ (proof),
`Re ψ(z) ≤ log|z| + 4` for `Re z ≥ 1/4`; apply it with `z = w/2` and `z = w`:
`(r_1/2)(log|w| − log 2π + 4) + r_2(log|w| − log 2π + 4) ≤ (r_1/2 + r_2)(log(1+|w|) + 4) = (n/2)log(1+|w|) + 2n`.
Further `Re(1/w) ≤ 1/Re w ≤ 2`, and `Re 1/(w−1) ≤ 0` if `Re w < 1`, `≤ 1/|Im w| ≤ 1` if `|Im w| ≥ 1`. Summing gives (a).
(b) `γ_K` is analytic and zero-free on `{Re w ≥ 1/2} ∖ {1}` [K1], so `log|Γr_ε(s)| = ∫_0^ε Re(γ_K'/γ_K)(s+u) du`; use (a) with
`|s + u| ≤ |s| + 1`, resp. (if `Re(s+ε) < 1`) the second bound with `|s + u| ≤ 1 + |t|`. ∎

**Lemma H_K (Hadamard midline inequality)** — changed only in the input [K1].
Suppose every non-trivial zero of `ζ_K` has `β ≤ θ*`. Let `ε > 0`, `Re s ≥ θ* − ε/2`, `Re s > 0`, `s ≠ 1`. Then `|ξ_K(s)| ≤ |ξ_K(s+ε)|`,
hence `|ζ_K(s)| ≤ |ζ_K(s+ε)| |Γr_ε(s)|` (read as `|ξ_K(1)|/|γ_K(s)|` if `s + ε = 1`).

*Proof.* [NBLp] Lemma H verbatim with the product of [K1]: if `ξ_K(s) ≠ 0`, then `ξ_K(s+ε) ≠ 0` and
`log|ξ_K(s+ε)| − log|ξ_K(s)| = ε(Re B + Σ_ρ Re(1/ρ)) + Σ_ρ log|(ρ − s − ε)/(ρ − s)| = Σ_ρ log|(ρ − s − ε)/(ρ − s)|` (only `Re B`
enters since `ε` is real; `Σ Re(1/ρ)` converges absolutely), and each term is `≥ 0` because
`|ρ − s − ε|² − |ρ − s|² = ε(ε − 2(β − σ)) ≥ 0` for `β ≤ θ* ≤ σ + ε/2`. Divide by `|γ_K(s)| ≠ 0` [K1]. ∎

**Lemma A_K (Lindelöf bound in a zero-free half-plane)** — changed: polynomial input [K4] of degree `n/2`; `K`-constants.
Let `1/2 ≤ θ < 1` and `ζ_K ≠ 0` on `Re s > θ`. For `δ ∈ (0, (1−θ)/2]` and `η' > 0` there is `C_A = C_A(K, θ, δ, η')` such that for
`σ ≥ θ + 2δ`: `|ζ_K(s)|^{−1} ≤ C_A (2 + |t|)^{η'}`, and, if also `|t| ≥ 1`, `|ζ_K(s)| ≤ C_A (2 + |t|)^{η'}`.

*Proof.* By [K3] (`ζ_K(s̄) = conj ζ_K(s)`) take `t ≥ 0`. *Small t* and `σ ≥ 3`: as in [NBLp] Lemma A (compactness; `|ζ_K|^{±1} ≤ ζ_K(3)`).
*Large t* (`t ≥ 6`): `Ω := {Re s > θ} ∖ (θ, 1]` is simply connected and free of zeros and poles; `log ζ_K` is the branch equal to the Euler
logarithm on `Re s > 1`. (1) Borel–Carathéodory [C4] with `z_0 = 2 + it'` (any `t' ≥ 4`), `R = 2 − θ − δ/2`, `r = 2 − θ − δ`: the disc
`|z − z_0| ≤ R` lies in `{Re z ≥ θ + δ/2 ≥ 1/2, Im z ≥ 2} ⊂ Ω`, where `Re log ζ_K = log|ζ_K| ≤ log C_6(K) + (n/2) log(t' + 4) =: U` by
[K4], and `|log ζ_K(z_0)| ≤ log ζ_K(2) < n` [K3]. Hence `|log ζ_K(z)| ≤ (6/δ)(U + n)` for `|z − z_0| ≤ r`, and so
(★) `|log ζ_K(z)| ≤ (C_7(K)/δ) log(Im z + 4)` for `Re z ≥ θ + δ`, `Im z ≥ 4`.
(2) Hadamard three circles [C5] for `log ζ_K`, centre `z_1 = 1 + 2δ + it`, `r_1 = δ`, `r_3 = 1 + δ − θ ≤ 3/4`, target `σ ∈ [θ + 2δ, 1 + δ]`,
`r_2 = 1 + 2δ − σ`: `M(r_1) ≤ log ζ_K(1 + δ) ≤ n log(1 + 1/δ)` [K3], `M(r_3) ≤ (C_7/δ) log(t + 5)` by (★), so with
`a_0 := log((1−θ)/δ)/log((1+δ−θ)/δ) < 1`: `|log ζ_K(σ+it)| ≤ C_9(K,θ,δ)(log(t + 5))^{a_0}`; for `σ ≥ 1 + δ`, `|log ζ_K| ≤ n log(1 + 1/δ)`.
Young's inequality `C_9 y^{a_0} ≤ η' y + C_10` (`y = log(t+5)`) gives `|ζ_K(s)|^{±1} ≤ e^{|log ζ_K(s)|} ≤ C_A (2 + t)^{η'}`. ∎

**Lemma P_K (Riesz–Perron)** — changed: `μ → b_K` (equivalently `μ_K` on ideals); constants.
Assume `ζ_K ≠ 0` on `Re s > θ` (`1/2 ≤ θ < 1`). Let `η ∈ (0,1)`, `δ ∈ (0, (1−θ)/2]`, `η' ∈ (0, 1/2]`. If `σ − η ≥ θ + 2δ`, then for all
`N ≥ 1`: `|P_N(s) − 1/ζ_K(s)| ≤ C_P N^{−η}(2 + |t|)^{η'}`, `C_P = C_P(K, θ, δ, η, η')`. If `σ − η > 1`, no hypothesis is needed and
`|P_N(s) − 1/ζ_K(s)| ≤ 2 ζ_K(σ − η) N^{−η}`.

*Proof.* [NBLp] Lemma P verbatim: `P_N(s) = (2πi)^{−1} ∫_{(c)} ζ_K(s+w)^{−1} N^w dw/(w(w+1))` (absolutely convergent, `|b_K| ≤ a_K`,
`c > max(0, 1−σ)`), shift to `Re w = −η`; the only pole crossed is `w = 0` with residue `1/ζ_K(s)` (none if `s = 1`, where
`1/ζ_K(1+w)` vanishes at `w = 0`); horizontal integrals vanish by Lemma A_K and [K3]; the remaining integral is bounded by Lemma A_K and
`2 + |t + v| ≤ (2 + |t|)(1 + |v|)`. If `σ − η > 1`: `P_N(s) − 1/ζ_K(s) = −N^{−1}Σ_{k≤N} b_K(k)k^{1−s} − Σ_{k>N} b_K(k)k^{−s}`, and
`N^{−1}k^{1−σ} = N^{−η}(k/N)^{1−η}k^{−(σ−η)} ≤ N^{−η}k^{−(σ−η)}` (`k ≤ N`), `k^{−σ} ≤ N^{−η}k^{−(σ−η)}` (`k > N`); sum with `|b_K| ≤ a_K`. ∎

## 6. Proof of Theorem B_K — changed (the `r_1` correction needs `(M_p)`; lines may lie left of `α_K`)

Fix `η ∈ (0, 1/p − θ)`, `σ_0 := (θ + η + 1/p)/2 ∈ (θ + η, 1/p)`, `δ := min{(σ_0 − η − θ)/2, (1 − θ − η)/2, (1 − θ)/2}`. `(M_p)` is
Proposition M_K(b). Let `c_k := b_K(k)(1 − k/N)`, `F_N := 𝟙 + Σ_{k≤N} c_k λ_k` (this is the ideal sum of the statement: group the ideals by
norm; `λ_1 = 0`) and `Φ_N := 𝟙₊ + Σ_{k≤N} c_k r_k`. Since `Σ_k c_k/k = P_N(1)`: `F_N = Φ_N|_{(0,1)} − P_N(1) r_1|_{(0,1)}` and
`‖F_N‖_p ≤ (1 − r(σ_0 − ½))^{−1/r} ‖Φ_N‖_{L²_{σ_0}} + |P_N(1)| ‖r_1‖_{L^p(0,1)}` (Lemma M). By Lemma MEM_K(a) (`σ_0 > θ`) and Lemma T_K(b,c),
`𝓜_{σ_0}Φ_N = (1 − ζ_K(s)P_N(s))/s = ζ_K(s)(1/ζ_K(s) − P_N(s))/s` on `Re s = σ_0`. Lemma A_K (`η' = 1/8`; continuity for `|t| ≤ 1`, `σ_0 < 1`)
and Lemma P_K (`σ_0 − η ≥ θ + 2δ`) give `‖Φ_N‖²_{L²_{σ_0}} ≤ C N^{−2η} ∫ (2+|t|)^{1/2}(σ_0² + t²)^{−1} dt`, and Lemma P_K at `s = 1`
(`1 − η ≥ θ + 2δ`) gives `|P_N(1)| ≤ 2^{1/8} C_P N^{−η}`. Hence `‖F_N‖_p ≤ C N^{−η}`. On `I_m`:
`F_N = 1 + Σ_k c_k (A_K(m)/k − A_K(⌊m/k⌋)) = 1 + P_N(1)A_K(m) − Σ_{N𝔟≤N} μ_K(𝔟)(1 − N𝔟/N) A_K(m/N𝔟)`. ∎

## 7. Proof of Theorem A_K

(ii) ⟹ (iii): `V^id_K ⊂ V_K`. (iii) ⟹ (i): Proposition N_K. (iii) ⟺ (iv): Proposition C_K(d). It remains to prove (i) ⟹ (ii) with the
bound. Assume (i). Then `(M_p)` holds (Proposition M_K(c)). Fix `ε ∈ (0, ε_0(p,K)]`, `σ' := 1/p − ε/2`. As in [NBLp]: `σ' ≥ 1/2`;
`σ' + ε < 1`; `2σ' − 1 ≥ (2 − p)/(2p) > 0`; `1 − ε > σ'`; `ε ≤ 1/2`; new: `nε ≤ 1/2`.

**A1 (the functions)** — changed. `φ_{ε,K}(j) := Σ_{kl=j} b_K(k) k^{−ε} a_K(l)`, the coefficients of
`ζ_K(s)/ζ_K(s+ε) = Σ_𝔠 Π_{𝔭|𝔠}(1 − N𝔭^{−ε}) N𝔠^{−s}`. For `y > 0`: `κ_K(y/k) = c_K y/k` for `k > y`, so
`Σ_{k≥1} b_K(k)k^{−ε}κ_K(y/k) = Σ_{k≤y} b_K(k)k^{−ε}κ_K(y/k) + c_K y Σ_{k>y} b_K(k)k^{−1−ε}` converges absolutely, and equals
`c_K y/ζ_K(1+ε) − Σ_{j≤y} φ_{ε,K}(j)` (`Σ_{k≤y} b_K(k)k^{−ε} A_K(y/k) = Σ_{j≤y} φ_{ε,K}(j)`). Put
`Φ_ε(x) := 𝟙₊(x) + Σ_k b_K(k)k^{−ε} r_k(x)` (pointwise, `x > 0`) and `H_ε := Σ_k b_K(k)k^{−ε} λ_k` on `(0,1)`. Then on `(0,1)`
`H_ε = Φ_ε − 1 − r_1/ζ_K(1+ε)`, and on `I_m`: `H_ε = A_K(m)/ζ_K(1+ε) − Σ_{j≤m} φ_{ε,K}(j)` (a step function; for `K = ℚ` it is
[NBLp]'s `−E_ε(m)`). For `N ≥ 1`: `H_{ε,N} := Σ_{k≤N} b_K(k)k^{−ε}(1 − k/N) λ_k = Σ_{N𝔟≤N} μ_K(𝔟)N𝔟^{−ε}(1 − N𝔟/N) λ_𝔟 ∈ V^id_K` (since
`b_K(k) = 0` for `k ∉ 𝒩_K`), `Φ_{ε,N} := 𝟙₊ + Σ_{k≤N} b_K(k)k^{−ε}(1 − k/N) r_k`, and on `(0,1)`: `H_{ε,N} = Φ_{ε,N} − 1 − P_N(1+ε) r_1`.

**A2 (transform on a right line)** — changed only through Lemma MEM_K. Let `b ∈ (1 − ε, 1)` (`⊂ [σ', 1)`). By Lemma MEM_K(b),
`Σ_k |b_K(k)| k^{−ε} ‖r_k‖_{L²_b} ≤ ‖r_1‖_{L²_b} Σ_k a_K(k) k^{−ε−b} = ‖r_1‖_{L²_b} ζ_K(b + ε) < ∞`, so `Σ_k b_K(k)k^{−ε} r_k` converges in
`L²_b` (to the pointwise sum, a.e.), `Φ_ε ∈ L²_b` and `𝓜_bΦ_ε(t) = Ψ_ε(b + it)`, `Ψ_ε(s) := (1 − R_ε(s))/s`
(`Σ_k b_K(k) k^{−s−ε} = 1/ζ_K(s+ε)`, `Re(s+ε) > 1`).

**A3 (shift to the left)** — unchanged argument; bounds via Lemmas H_K, Γ_K. `Ψ_ε` is analytic on `{1/p − ε < Re s < 1}` (there
`Re(s+ε) > 1/p`, so `ζ_K(s+ε) ≠ 0` by (i); `ζ_K(s)` is analytic for `s ≠ 1`). On `{σ' ≤ σ ≤ b, |t| ≥ 1}`: `|R_ε(s)| ≤ |Γr_ε(s)| ≤
e^{c_1/2}(3 + |t|)^{1/4}` (Lemma H_K with `θ* = 1/p`; Lemma Γ_K(b); `nε ≤ 1/2`); on `{σ' ≤ σ ≤ b, |t| ≤ 1}`, `R_ε` is continuous (`b < 1`), hence
bounded by `C_13(ε,K)`. So `sup_{σ'≤σ≤b} ∫|Ψ_ε(σ+it)|² dt < ∞` and Lemma PW gives `Φ_ε ∈ L²_a`, `𝓜_aΦ_ε = Ψ_ε(a + i·)` for `a ∈ [σ', b]`.

**A4 (key estimate `B(ε) := ‖Φ_ε‖²_{L²_{σ'}} ≤ K_{p,K} ε`)** — changed: pole term with `c_K`; `n` and `|d_K|` in `C_{Γ,K}`. On `Re s = σ'`,
`R := R_ε(s)`, `|R| ≤ |Γr_ε(s)|` (Lemma H_K), so
`2πB(ε) = ∫_ℝ |1 − R|²|s|^{−2} dt ≤ ∫|s|^{−2}dt − 2 Re J + ∫|Γr_ε(s)|²|s|^{−2} dt`, `J := ∫_ℝ R(σ' + it)|s|^{−2} dt`.
(a) `∫ dt/(σ'² + t²) = π/σ'` — unchanged.
(b) *Exact evaluation of J.* `|s|² = s(2σ' − s)` on the line, so `J = (1/i)∫_{(σ')} f(s) ds`, `f(s) := R(s)/(s(2σ' − s))`. In `Re s ≥ σ'`, `f` is
meromorphic with exactly two simple poles, `s = 1` and `s = 2σ'` (`2σ' > 1`): `ζ_K(s+ε) ≠ 0` for `Re(s+ε) ≥ 1/p + ε/2`, and at `s = 1 − ε`
the pole of `ζ_K(s+ε)` makes `R = 0`. Rectangle `[σ', X] × [−T, T]`, `X > 2σ' + 1`: on the horizontal sides
`|f| ≤ e^{c_1ε}(2 + X + T)^{nε/2}/T² → 0` (Lemma H_K on all of `Re s ≥ σ'`; Lemma Γ_K(b), `|Im| = T ≥ 1`); on `Re s = X`:
`|R| ≤ ζ_K(X)ζ_K(X+ε) ≤ ζ_K(X)²` [K3] and `∫|ds|/(|s||2σ' − s|) ≤ π/(X − 2σ') → 0`; the integral over `Re s = σ'` converges absolutely (A3).
Residues: `Res_{s=1} f = c_K/(ζ_K(1+ε)(2σ' − 1))` (`Res_{s=1} ζ_K = c_K`), `Res_{s=2σ'} f = −R(2σ')/(2σ')`; the rectangle is clockwise
for the region to its right, so `∫_{(σ')} f ds = −2πi(Res_{s=1} + Res_{s=2σ'})` and
`J = 2π [ R(2σ')/(2σ') − c_K/((2σ' − 1) ζ_K(1+ε)) ]` (real).
(c) Lemma Γ_K(b) (`Re(s+ε) < 1`): `|Γr_ε|² ≤ e^{x}`, `x := ε(2(c_1 − 1) + n log(2 + |t|))`; `e^x − 1 ≤ x e^x` and
`e^x ≤ e^{(c_1−1)/n}(2+|t|)^{1/2}` (`ε ≤ 1/(2n)`). With `σ' ≥ 1/2`:
`∫|Γr_ε|²|s|^{−2}dt ≤ π/σ' + C_{Γ,K} ε`, `C_{Γ,K} := e^{(c_1−1)/n} ∫_ℝ (2(c_1−1) + n log(2+|t|))(2+|t|)^{1/2}(1/4 + t²)^{−1} dt`
(`C_{Γ,ℚ} = C_Γ = 5993.05…` of [NBLp]).
Adding: `2πB(ε) ≤ (2π/σ')(1 − R(2σ')) + 4πc_K/((2σ'−1)ζ_K(1+ε)) + C_{Γ,K}ε`. Here `R(2σ') = ζ_K(2σ')/ζ_K(2σ'+ε) ≥ 1` (`ζ_K` decreasing on
`(1,∞)`, [K3]), `c_K/ζ_K(1+ε) ≤ C(K)ε` and `2σ' − 1 ≥ (2−p)/(2p)`. Therefore `B(ε) ≤ (4pC(K)/(2−p) + C_{Γ,K}/(2π)) ε = K_{p,K} ε`.

**A5 (`H_ε` approximates `−𝟙`)** — changed (correction term). On `(0,1)`: `𝟙 + H_ε = Φ_ε − r_1/ζ_K(1+ε)`. Lemma M (`a = σ'`), A4, `(M_p)`:
`‖𝟙 + H_ε‖_p ≤ (2/(rε))^{1/r}(K_{p,K}ε)^{1/2} + ‖r_1‖_{L^p(0,1)} C(K)ε/c_K = A_{p,K} ε^{1−1/p} + τ_{p,K} ε`.

**A6 (`H_{ε,N} → H_ε`)** — unchanged argument; constants depend on K. `σ'' := 1/p − ε/4 ∈ (σ', 1/p)`. By A3 and Lemma MEM_K(b),
`Φ_ε − Φ_{ε,N} ∈ L²_{σ''}` with transform `(ζ_K(s)/s)(P_N(s+ε) − 1/ζ_K(s+ε))` (the transform of `Φ_{ε,N}` is `(1 − ζ_K(s)P_N(s+ε))/s`). On this
line, for `|t| ≥ 1`: `|ζ_K(s)| ≤ |ζ_K(s+ε/2)| |Γr_{ε/2}(s)|` (Lemma H_K with `ε/2`); `Re(s+ε/2) = 1/p + ε/4 = θ + 2δ` with `θ = 1/p`, `δ = ε/8`,
so `|ζ_K(s+ε/2)| ≤ C_A(2 + |t|)^{1/8}` (Lemma A_K); `|Γr_{ε/2}(s)| ≤ e^{c_1/4}(3 + |t|)^{1/8}` (Lemma Γ_K(b), `nε/4 ≤ 1/8`); and Lemma P_K
(`σ = 1/p + 3ε/4`, `η = ε/2`, `σ − η = θ + 2δ`, `η' = 1/8`) gives `|P_N(s+ε) − 1/ζ_K(s+ε)| ≤ C_P N^{−ε/2}(2 + |t|)^{1/8}`. For `|t| ≤ 1` use
continuity. Hence `‖Φ_ε − Φ_{ε,N}‖²_{L²_{σ''}} ≤ C_15(p,K,ε) N^{−ε}`. Lemma P_K at `s = 1 + ε` (`η = ε/2`, `σ − η > 1`):
`|1/ζ_K(1+ε) − P_N(1+ε)| ≤ 2ζ_K(1 + ε/2) N^{−ε/2}`. With A1 and Lemma M (`a = σ''`):
`‖H_ε − H_{ε,N}‖_p ≤ (4/(rε))^{1/r}‖Φ_ε − Φ_{ε,N}‖_{L²_{σ''}} + ‖r_1‖_{L^p(0,1)} |1/ζ_K(1+ε) − P_N(1+ε)| ≤ C_17(p,K,ε) N^{−ε/2}`.

**A7 (conclusion; order of limits)** — unchanged. Fix `ε`, choose `N_0(p,K,ε)` with `C_17 N_0^{−ε/2} ≤ ε`. For `N ≥ N_0`:
`d^id_{N,p}(K) ≤ ‖𝟙 + H_{ε,N}‖_p ≤ A_{p,K}ε^{1−1/p} + (1 + τ_{p,K})ε`. `d^id_{N,p}` is non-increasing in `N`, so its limit is at most this for
every `ε ∈ (0, ε_0]`; `N → ∞` is taken at fixed `ε`, then `ε → 0` (`K_{p,K}`, `τ_{p,K}` do not depend on `ε`). The hypothesis (i) is used
only on `Re s ≥ 1/p + ε/4` (Lemmas A_K, P_K) and through Lemma H_K; nothing is used on the line `Re s = 1/p`. ∎

## 8. Proof of the corollaries

1. [CT] gives (i) for every K and every `p ≤ 32/31` (`1/p ≥ 31/32`); Theorem A_K applies for `1 < p ≤ 32/31`, Theorem B_K with `θ = 31/32`
   for `p < 32/31`. (`(M_p)` comes from Proposition M_K(c), so no lattice-point input beyond [K2] is needed, for any degree.)
2. The same with [L2] (`θ = B_2`, `1/B_2 = 852/(777 − √1081) = 1.14497440…`) or [Hu] (`θ = 7/8`), conditionally.
3. (i) ⟺ `Θ_K ≤ 1/p`. For `p > 2`: [K1] gives a non-trivial zero, and by the symmetry `ρ ↦ 1 − ρ̄` one with `β ≥ 1/2 > 1/p`; if `(M_p)`
   fails, (iii) fails by definition, otherwise Proposition N_K applies.
4. GRH ⟺ `Θ_K = 1/2` ⟺ (i) for all `p < 2`; the rate is Theorem B_K with `θ = 1/2`.
5. For abelian `K`, `ζ_K(s) = Π_{χ} L(s,χ)` over the primitive Dirichlet characters of the character group of K (conductor–discriminant
   formula); the non-trivial zeros of `ζ_K` are the union of those of the factors (no factor vanishes on `Re s ≥ 1`). ∎

## 9. Remarks

**9.1 Why this formulation.**
* *The function.* `κ_K = Ψ_{ζ_K}` is de Roton's complementary function [dR07, dR09]: the unique "error term" whose Mellin transform is
  `−ζ_K(s)/s` (Lemma T_K), exactly as `{y}` for ζ. Its dilates are counting errors of ideal multiples (Proposition C_K(a)):
  `λ_𝔟 = A_K(m)/N𝔟 − #{𝔠 ⊆ 𝔟 : N𝔠 ≤ m}` on `I_m`, the ideal analogue of `{m/k} = m/k − ⌊m/k⌋`.
* *Ideals versus norms.* Dilations indexed by integral ideals depend only on the norm, and the Möbius coefficients `μ_K(𝔟)` group into
  `b_K(k)`, which vanish off `𝒩_K`. So the ideal family `{λ_𝔟}` is the natural one (smallest, and it carries the explicit approximants), while
  the full natural-number family `{λ_k}` (de Roton's `B_F^{nat}`) gives the same criterion and, in addition, the structure statement (iv).
* *The side condition.* As in [NBLp], `λ_k = r_k − r_1/k` (Báez-Duarte's/Beurling's `Σ c_k/θ_k = 0`) makes every function vanish on
  `(1,∞)`, so the problem lives in `L^p(0,1)` with target `𝟙`, and the span sits in the step-function space `S_p`.
* *The membership condition `(M_p)`.* For `n ≥ 2`, `κ_K` is unbounded (for `K = ℚ(i)` it is ¼ of the circle-problem error,
  `Ω(y^{1/4})` by Hardy), and whether `κ_K(1/x) ∈ L^p(0,1)` is a mean-value question about the ideal-counting error. It holds unconditionally
  for `p < (n+1)/(n−1)` (Landau), and whenever (i) holds (Proposition M_K(c)); so it never restricts the corollaries. Without it the closure
  is not defined. De Roton includes the same condition (`Ψ^{(1)}_F ∈ L²(0,∞)`) in her `p = 2` theorems and removes it only by smoothing
  `Ψ_F` (Riesz means of order `k`, [dR07, §5]); a smoothed `L^p` version would change both the functions and the target and is not
  pursued. In Theorem A_K, `(M_p)` is a genuine part of (ii)–(iv): (i) implies it, and it is needed for (ii) to be meaningful.
* *Why not Dimitrov–Oliveira's class.* Their `L^p` theorem (arXiv 1608.07887, Thm GNB, continuous dilations, all `p > 1`) assumes
  `κ_η` **bounded** (their class `𝒟`, TeX l. 317–323), which excludes every `ζ_K` with `n ≥ 2`. Theorem A_K also yields their continuous
  version for `ζ_K` (functions `Σ_k b_k κ_K(1/(θ_k x))`, real `θ_k ≥ 1`, `Σ b_k/θ_k = 0`): sufficiency follows from `V_K ⊂` that span, and necessity from the
  same functional `x^{w_0−1}` (transform `−ζ_K(w)(Σ b_k θ_k^{−w})/w` by Lemma T_K).

**9.2 Where K enters (lemma by lemma).**

| Item | K-dependence |
|---|---|
| Prop. C_K, N_K | values `A_K(m)/k − A_K(⌊m/k⌋)`; `(M_p)` needed for `V_K ⊂ L^p`; Mellin identity continued from `α_K < Re w < 1` |
| Lemma T_K | strip `α_K < σ < 1` (Landau's exponent `α_K = 1 − 2/(n+1)`; only `α_K < 1` used) |
| Lemma MEM_K (new) | moves the Mellin identity to lines left of `α_K` under a zero-free half-plane (Lindelöf + Hadamard) |
| Lemma Γ_K | `c_1(K) = 3 + ½ log|d_K| + 2n`; growth `(2+|t|)^{nε/2}` forces `ε ≤ 1/(2n)` |
| Lemma H_K | Hadamard product of `ξ_K` ([K1]); unchanged otherwise |
| Lemma A_K | convexity input `(|t|+2)^{n/2}` ([K4]); `log ζ_K(2) < n` |
| Lemma P_K | coefficients `b_K` (`|b_K| ≤ a_K`) |
| A4 | pole term `c_K/ζ_K(1+ε) ≤ C(K)ε` (`ζ_K(1+ε) ≥ c_K/ε` may fail when the Euler–Kronecker constant is negative, hence `C(K)`); `C_{Γ,K}` |
| A5, Thm B_K | correction `‖κ_K(1/x)‖_{L^p(0,1)} · |P_N(1)|` (bounded only via `(M_p)`) |
| `ε_0(p,K)` | `min{1/2, (2/p − 1)/2, 1 − 1/p, 1/(2n)}` |

No uniformity in K is claimed. The degree never restricts `p`: it only shrinks `ε_0` and enlarges constants.

**9.3 What is routine and what is not.** Given [NBLp], the proof is a change of inputs, with two non-cosmetic points: (1) for
`n ≥ 2` the functions are unbounded and, for large `n`, not known to lie in `L²_σ` on the lines `σ ≈ 1/p` that the argument uses; Lemma
MEM_K shows that the zero-free hypothesis itself puts them there, so the endpoint argument and the corollaries hold for every degree; (2) the
gamma factor has degree `n`, so the Hadamard-midline loss is `|t|^{nε/2}` and the endpoint shift must satisfy `ε ≤ 1/(2n)`, which is harmless
because `ε → 0` anyway.

**9.4 `p = 2`.** Not covered (A4 uses `2σ' > 1`). De Roton [dR09] proves the `p = 2` natural-number criterion for the Selberg class in
`L²(0,∞)` (no side condition) with rate `(log log N)^{−1/3}`.
