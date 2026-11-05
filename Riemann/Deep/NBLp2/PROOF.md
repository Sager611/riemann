# Natural-number Nyman–Beurling (Báez-Duarte) criterion in L^p(0,1) for Dirichlet L-functions, 1 < p < 2

Node N5. Author: Orca agent `nblp2`, 2026-11-05. Paper proof (no Lean). It adapts [NBLp] = `Riemann/Deep/NBLp/PROOF.md`
(the case q = 1) lemma by lemma; every item is marked **unchanged**, **changed (how)** or **new**. Inputs that are
not proved here are quoted in §2.

## 0. Notation

* `χ` is a Dirichlet character mod `q ≥ 1`, induced by the primitive character `χ*` mod `q* | q`; `a ∈ {0,1}` with
  `χ(−1) = (−1)^a`; `δ = δ_χ := 1` if `χ` is principal (then `χ* = 1`, `q* = 1`), else `0`; `c_q := φ(q)/q`.
  `q = 1` (`χ ≡ 1`, `L = ζ`) is exactly [NBLp].
* `S_χ(y) := Σ_{1≤n≤y} χ(n)` (`y ≥ 0`; `S_χ(y) = 0` for `y < 1`).
* `κ_χ(y) := δ c_q y − S_χ(y)` (`y ≥ 0`). For `q = 1`: `κ(y) = {y}`. `C_κ := sup_{y≥0} |κ_χ(y)| ≤ φ(q)`: if `δ = 0`, `S_χ` is
  `q`-periodic (the sum over a period vanishes) and `|S_χ| ≤ φ(q)` on `[0,q]`; if `δ = 1`, `κ_χ(y+q) = κ_χ(y)` and on `[0,q]` both
  `c_q y` and `S_χ(y)` lie in `[0, φ(q)]`.
* For `k ≥ 1`, `x > 0`: `r_k(x) := κ_χ(1/(kx))` and `λ_k := r_k − δ k^{−1} r_1` (so `λ_1 ≡ 0` when `δ = 1`). Since `S_χ(y) = 0`
  for `y < 1`: `r_k(x) = δ c_q/(kx)` for `x > 1/k`, hence **`λ_k = 0` on `(1, ∞)`** in both cases.
  For `q = 1`, `λ_k = f_k = ρ_{1/k}` of [NBLp].
* `𝟙 := 1_{(0,1)}` on `(0,1)`; `𝟙₊ := 1_{(0,1]}` on `(0,∞)` ([NBLp] writes `χ`, `χ₊`; renamed because `χ` is the character).
* `V_χ := span{λ_k|_{(0,1)} : k ≥ 1}`; `W_{χ,p}` := its closure in `L^p(0,1)`;
  `d_{N,p}(χ) := inf_{c_k ∈ ℂ} ‖𝟙 − Σ_{k=1}^N c_k λ_k‖_{L^p(0,1)}`. `I_n`, `S_p`, `r = 2p/(2−p)`, `L²_σ`, `𝓜_σ` as in [NBLp].
* In the notation of Dimitrov–Oliveira (arXiv 1608.07887, lines 317–330), with `η(s) = L(s,χ)` itself (their
  `α = δ`, `β = δ c_q`): `κ_η = κ_χ` and `V_χ = C_{η,ℕ}` (for `δ = 1` their side condition
  `β Σ b_k/k^α = 0` is built into `λ_k = r_k − r_1/k`), so `W_{χ,p} = C^p_{η,ℕ}`.
* `L(s,χ) = Σ χ(n) n^{−s}`; `E_χ(s) := Π_{p | q, p ∤ q*} (1 − χ*(p) p^{−s})`, so `L(s,χ) = L(s,χ*) E_χ(s)`
  (for `δ = 1`: `L(s,χ) = ζ(s) Π_{p|q}(1 − p^{−s})`).
* Gamma factor of `χ*`: if `χ* ≠ 1`, `γ(s) := (q*/π)^{(s+a)/2} Γ((s+a)/2)`; if `χ* = 1`,
  `γ(s) := ½ s(s−1) π^{−s/2} Γ(s/2)` ([NBLp]). `ξ(s,χ*) := γ(s) L(s,χ*)`; `Γr_ε(s) := γ(s+ε)/γ(s)`.
* `R_ε(s) := L(s,χ)/L(s+ε,χ)`; `P_N(s) := Σ_{k≤N} μ(k) χ(k) (1 − k/N) k^{−s}` (Riesz mean of `1/L(s,χ)`); `1/L(1,χ) := 0`
  if `δ = 1`.
* `Θ_χ := sup{Re ρ : L(ρ,χ) = 0, Re ρ > 0}` (= sup of the real parts of the non-trivial zeros of `L(s,χ*)`; `Θ_χ ≥ 1/2`).

## 1. Results

**Proposition C_χ (structure)** — changed (values on `I_n`), argument unchanged.
(a) For `k, n ≥ 1`, `x ∈ I_n`: `λ_k(x) = δ k^{−1} S_χ(n) − S_χ(⌊n/k⌋)`.
(b) `V_χ ⊂ S_p`; literal density of `V_χ` in `L^p(0,1)` is false for every `p`.
(c) `D_m λ_k = λ_{km} − δ k^{−1} λ_m` a.e. on `(0,1)` (`D_m` as in [NBLp]); so `D_m W_{χ,p} ⊂ W_{χ,p}`.
(d) `𝟙 ∈ W_{χ,p} ⟺ W_{χ,p} = S_p ⟺ d_{N,p}(χ) → 0`.

**Theorem A_χ.** Let `1 < p < 2` and `χ` any Dirichlet character. Equivalent:
(i) `L(s,χ) ≠ 0` for `Re s > 1/p` (i.e. `Θ_χ ≤ 1/p`);
(ii) `𝟙 ∈ W_{χ,p}`, i.e. `d_{N,p}(χ) → 0`;
(iii) `W_{χ,p} = S_p`.
Moreover, under (i), with `ε_0 := min{1/2, (2/p − 1)/2, 1 − 1/p}`, for every `ε ∈ (0, ε_0]` there is `N_0(p,χ,ε)` with
`‖𝟙 + Σ_{k≤N} μ(k) χ(k) k^{−ε} (1 − k/N) λ_k‖_{L^p(0,1)} ≤ A_{p,χ} ε^{1−1/p} + (1 + qδ) ε` for `N ≥ N_0`,
`A_{p,χ} := K_{p,χ}^{1/2} (2/r)^{1/r}`, `K_{p,χ}` from A4 (for `q = 1` it is [NBLp]'s `K_p`).

**Theorem B_χ (open range, explicit coefficients, power rate).** Let `1/2 ≤ θ < 1`, `L(s,χ) ≠ 0` for `Re s > θ`,
`1 < p < 1/θ`, `0 < η < 1/p − θ`. Then `‖𝟙 + Σ_{k≤N} μ(k) χ(k) (1 − k/N) λ_k‖_{L^p(0,1)} ≤ C(p,θ,η,χ) N^{−η}` (`N ≥ 2`).

**Proposition N_χ (necessity, all `1 < p < ∞`)** — unchanged argument. If `L(w_0,χ) = 0`, `Re w_0 > 1/p`, then `𝟙 ∉ W_{χ,p}`.

**Corollaries.**
1. (Unconditional.) By [C9_χ] (`L(s,χ) ≠ 0` for `Re s > 7/8`, every `χ`; Lean port in our repo): (ii) holds for every `χ` and
   every `1 < p ≤ 8/7`; with [C9'] (`Re s > θ_B = 0.87338196661320752`) for every `1 < p ≤ 1/θ_B = 1.1449744…`.
   Theorem B_χ with `θ = 7/8` gives the power rate for `p < 8/7`.
2. For `1 < p < 2`: `𝟙 ∈ W_{χ,p} ⟺ p ≤ 1/Θ_χ`. For `p > 2`: `𝟙 ∉ W_{χ,p}` for every `χ` (unconditionally).
3. GRH for `L(s,χ)` ⟺ `𝟙 ∈ W_{χ,p}` for all `p ∈ (1,2)`; under it the rate `N^{−η}`, `η < 1/p − 1/2` holds.
4. For principal `χ_0` mod `q`, (i) is `ζ(s) ≠ 0` on `Re s > 1/p`: every modulus gives an `L^p` criterion for ζ with the step
   functions `S_{χ_0}(n)/k − S_{χ_0}(⌊n/k⌋)`. The side condition can be dropped (Remark 9.1).

## 2. Cited results

[C1]–[C8], [C10] as in [NBLp §2]. In addition:
* **[C1_χ]** (Abel summation.) For `χ` non-principal and `σ > 0`: `L(s,χ) = s ∫_1^∞ S_χ(u) u^{−s−1} du` (absolutely
  convergent, `|S_χ| ≤ C_κ`); in Lean: `LFunction_eq_mul_Ichi` (`Riemann/Deep/QRHCor`). For every `χ` and `σ > 1` the
  same formula holds. Also, for `σ > 0`, `x ≥ 1`: `L(s,χ) = Σ_{n≤x} χ(n) n^{−s} − S_χ(x) x^{−s} + s∫_x^∞ S_χ(u) u^{−s−1} du`
  (non-principal).
* **[C2_χ]** (Davenport, *Multiplicative Number Theory*, 3rd ed., Ch. 9, Ch. 12 (17)–(18), Ch. 16.) For primitive `χ* ≠ 1`,
  `ξ(s,χ*)` is entire of order 1; its zeros are exactly the non-trivial zeros `ρ` of `L(s,χ*)`, all with `0 < β < 1`;
  `Σ_ρ |ρ|^{−2} < ∞`; `ξ(s,χ*) = e^{A + Bs} Π_ρ (1 − s/ρ) e^{s/ρ}`; `Re B = −Σ_ρ Re(1/ρ)`. Functional equation:
  `ξ(1 − s, χ̄*) = w(χ*) ξ(s,χ*)`, `|w| = 1`, so the zeros are symmetric under `ρ ↦ 1 − ρ̄`, and there are infinitely
  many. `γ` is analytic and zero-free on `Re s > 0` (poles of `Γ((s+a)/2)` at `s = −a − 2m ≤ 0`). `L(s,χ) ≠ 0` for
  `Re s ≥ 1` (all `χ`; for `δ = 1`, `s ≠ 1`); `L(s,χ_0)` has a simple pole at `1` with residue `c_q`.
* **[C8_χ]** (Elementary, `σ > 1`.) `L(s,χ) = Π_p (1 − χ(p)p^{−s})^{−1}`, `log L(s,χ) = Σ_p Σ_m χ(p)^m p^{−ms}/m`,
  `|log L(s,χ)| ≤ log ζ(σ)`, `|L(s,χ)|^{±1} ≤ ζ(σ)` (`1/L = Σ μ(n)χ(n) n^{−s}`); `conj L(s̄, χ) = L(s, χ̄)`;
  `L(σ,χ_0) = Σ_{(n,q)=1} n^{−σ}` is positive and decreasing on `(1,∞)`.
* **[C9_χ]** (OpenAI 2026; Lean port `OAI.DirichletCharacter.LFunction_ne_zero_of_seven_eighths_lt_re` and
  `OAI.riemannZeta_ne_zero_of_seven_eighths_lt_re`, `Riemann/Deep/QRH`, sorry-free.) For every `q ≥ 1`, every `χ` mod `q`
  and `Re s > 7/8` with `¬(χ = 1 ∧ s = 1)` (Mathlib's `χ = 1` is the principal character): `L(s,χ) ≠ 0`.
  **[C9']** (J. Bebel, QuasiRiemannTracker, github.com/GTimG/QuasiRiemannTracker, 31 Oct 2026; Lean-kernel-checked per the
  tracker, axiom list not inspected by us; see `TREE.md` N7.) `ζ(s) ≠ 0` and `L(s,χ) ≠ 0` for `Re s > 0.87338196661320752`.
* **[DO]** (Dimitrov–Oliveira, arXiv 1608.07887, Thm "GNB" (TeX line 654) and Thm "GBNew".) For `η ∈ 𝒟` (in particular
  `η = L(·,χ)`) and every `p > 1`: `η ≠ 0` on `Re s > 1/p ⟺ 𝟙 ∈ C^p_η` (continuous dilations `θ_k ≥ 1`); and
  `𝟙 ∈ C^2_{η,ℕ}` ⟺ `η ≠ 0` on `Re s > 1/2` (the case `p = 2` of Theorem A_χ; also de Roton for the Selberg class).

## 3. Structure

**Proof of Proposition C_χ** — changed in (a), (c); unchanged in (b), (d).
(a) `y := 1/x ∈ [n, n+1)`; `⌊y/k⌋ = ⌊n/k⌋`, so `S_χ(y/k) = S_χ(⌊n/k⌋)`, `S_χ(y) = S_χ(n)`, and
`λ_k(x) = δc_q y/k − S_χ(⌊n/k⌋) − δk^{−1}(c_q y − S_χ(n)) = δ k^{−1} S_χ(n) − S_χ(⌊n/k⌋)`. (For `q = 1`: `{n/k}`.)
(b) As in [NBLp]: each `λ_k` is constant on every `I_n`; `S_p` is closed of infinite codimension.
(c) `r_{km}(x) = κ_χ(1/(kmx)) = r_k(mx)` for all `x > 0`, so for all `x > 0`:
`λ_{km}(x) − δk^{−1}λ_m(x) = r_k(mx) − δ(km)^{−1}r_1(x) − δk^{−1}r_1(mx) + δ(km)^{−1}r_1(x) = λ_k(mx)`.
For `x < 1/m` this is `(D_mλ_k)(x)`; for `x > 1/m` it is `λ_k(mx) = 0` (`mx > 1`, `λ_k = 0` on `(1,∞)`) `= (D_mλ_k)(x)`.
(d) As in [NBLp]: `D_m 𝟙 = 1_{(0,1/m)} ∈ W`, so `1_{I_m} ∈ W` for all `m`, and these span a dense subspace of `S_p`. ∎

**Proof of Proposition N_χ** — unchanged. Zeros with `Re w_0 > 0` have `Re w_0 < 1` [C2_χ]. `g(x) := x^{w_0−1} ∈ L^{p'}(0,1)`
(`p' = p/(p−1)`; `p'(Re w_0 − 1) > −1 ⟺ Re w_0 > 1/p`). As `λ_k = 0` on `(1,∞)`, Lemma T_χ gives
`∫_0^1 λ_k g = 𝓜[λ_k](w_0) = −L(w_0,χ)(k^{−w_0} − δk^{−1})/w_0 = 0`, while `∫_0^1 𝟙 g = 1/w_0 ≠ 0`. ∎

## 4. Mellin–Plancherel toolkit

**Lemma T_χ** — changed (the functions; for `δ = 0` valid for all `σ > 0`).
(a) Let `0 < σ < 1` (any `σ > 0` if `δ = 0`), `k ≥ 1`. Then `x^{σ−1} r_k ∈ L¹(0,∞)`, `r_k ∈ L²_σ`,
`‖r_k‖_{L²_σ} = k^{−σ}‖r_1‖_{L²_σ}`, `‖r_1‖²_{L²_σ} ≤ C_κ²/(2σ)` (`δ = 0`), `≤ C_κ²/(2σ) + c_q²/(2 − 2σ)` (`δ = 1`), and `𝓜_σ r_k(t) = −k^{−s} L(s,χ)/s`.
(b) `𝓜_σ 𝟙₊ = 1/s` (`σ > 0`) — unchanged.

*Proof.* `r_k(x) = r_1(kx)`; `|r_1| ≤ C_κ` on `(0,1]` and `r_1(x) = δc_q/x` on `(1,∞)`: integrability, the bound, scaling.
With `u = 1/(kx)`: `∫_0^∞ r_k x^{s−1} dx = k^{−s} ∫_0^∞ κ_χ(u) u^{−s−1} du`. If `δ = 0` this is
`−k^{−s}∫_1^∞ S_χ u^{−s−1} du = −k^{−s}L(s,χ)/s` [C1_χ]. If `δ = 1`: `∫_0^1 c_q u · u^{−s−1} du = c_q/(1−s)`
(`σ < 1`); for `σ > 1`, [C1_χ] gives `∫_1^∞ κ_χ u^{−s−1} du = c_q/(s−1) − L(s,χ)/s`, and both sides are analytic on
`{σ > 0} ∖ {1}` (`κ_χ` bounded; [C2_χ]), so it holds for `0 < σ < 1` too. Adding: `−L(s,χ)/s`. ∎

**Lemma M (Hölder)** — unchanged. **Lemma PW (shifting the line)** — unchanged.

## 5. Analytic lemmas

**Lemma Γ_χ** — unchanged if `χ* = 1` ([NBLp] Lemma Γ); changed (q-factor, no exceptional point) if `χ* ≠ 1`.
Let `χ* ≠ 1`. (a) For `Re w ≥ 1/2`: `Re(γ'/γ)(w) ≤ 2 + ½ log q* + ½ log(1 + |w|)`.
(b) For `0 < ε ≤ 1`, `Re s ≥ 1/2`: `|Γr_ε(s)| ≤ e^{2ε} (q*)^{ε/2} (2 + |s|)^{ε/2}`; if also `Re(s+ε) < 1`,
`|Γr_ε(s)|² ≤ e^{4ε} (q*)^{ε} (2 + |t|)^{ε}`.
Hence in all cases (with [NBLp]'s side conditions when `χ* = 1`): `|Γr_ε(s)| ≤ e^{5ε}(q*)^{ε/2}(2+|s|)^{ε/2}` and, if
`Re(s+ε) < 1`, `|Γr_ε(s)|² ≤ e^{8ε}(q*)^{ε}(2+|t|)^{ε}`.

*Proof.* `(γ'/γ)(w) = ½ log(q*/π) + ½ ψ((w+a)/2)`. `z := (w+a)/2` has `Re z ≥ 1/4`, `|z| ≤ (1 + |w|)/2`, and
`Re ψ(z) ≤ log|z| + 4` for `Re z ≥ 1/4` ([NBLp] Lemma Γ proof). So
`Re(γ'/γ)(w) ≤ ½ log q* − ½ log π + ½ log(1+|w|) − ½ log 2 + 2`. (b) `γ` is analytic and zero-free on `Re w > 0`
[C2_χ], so `log|Γr_ε(s)| = ∫_0^ε Re(γ'/γ)(s+u) du` with `|s + u| ≤ |s| + 1`; if `Re(s+ε) < 1`, `|s+u| ≤ 1 + |t|`. ∎

**Lemma E (Euler factors)** — new. `E_χ` is entire with zeros only on `Re s = 0` (`|χ*(p)| = 1` for `p ∤ q*`). For
`Re s ≥ 1/2` and `0 ≤ ε ≤ 1`: `|E_χ(s)/E_χ(s+ε)| ≤ q^{3ε}`, `|E_χ(s)| ≤ 2^{ω(q)} ≤ q`.

*Proof.* Put `z := χ*(p)p^{−s}`, `|z| = p^{−σ} ≤ p^{−1/2}`. Then `|1 − z| ≤ |1 − zp^{−ε}| + |z|(1 − p^{−ε})`,
`1 − p^{−ε} ≤ ε log p` and `|1 − zp^{−ε}| ≥ 1 − p^{−1/2}`, so the factor ratio is `≤ 1 + ε log p/(√p − 1) ≤ p^{2.42ε}`
(`1/(√2 − 1) < 2.42`). The product over `p | q` is `≤ (rad q)^{2.42ε} ≤ q^{3ε}`. `|1 − z| ≤ 2` gives the last bound. ∎

**Lemma H_χ (Hadamard midline inequality)** — changed: [C2_χ] for `ξ(·,χ*)`; Euler factors via Lemma E.
Suppose `L(s,χ) ≠ 0` for `Re s > θ*` (some `θ* ≥ 1/2`). Let `0 < ε ≤ 1` and `Re s ≥ max(θ* − ε/2, 1/2)`; if `δ = 1` let
`s ∉ {1, 1 − ε}`. Then `|L(s,χ)| ≤ q^{3ε} |Γr_ε(s)| |L(s+ε,χ)|`, i.e. `|R_ε(s)| ≤ q^{3ε}|Γr_ε(s)|`. (If `δ = 1`:
`R_ε(1−ε) = 0`.)

*Proof.* Non-trivial zeros of `L(·,χ*)` with `Re > 0` are zeros of `L(·,χ)` (Lemma E), so all have `β ≤ θ*`. If
`χ* ≠ 1`, the argument of [NBLp] Lemma H applies verbatim to the product of [C2_χ]: for `w` not a zero,
`log|ξ(w,χ*)| = Re A + Re(Bw) + Σ_ρ (log|1 − w/ρ| + Re(w/ρ))`; subtracting at `w = s+ε` and `w = s` (`ε` real), the
terms `ε Re B + ε Σ_ρ Re(1/ρ)` cancel (only `Re B` enters; `Σ Re(1/ρ)` converges absolutely), leaving
`Σ_ρ log|(ρ − s − ε)/(ρ − s)| ≥ 0` as `(β−σ−ε)² − (β−σ)² = ε(ε − 2(β − σ)) ≥ 0`. So `|ξ(s,χ*)| ≤ |ξ(s+ε,χ*)|`; divide by
`|γ(s)| ≠ 0`. If `χ* = 1`, this is [NBLp] Lemma H. Multiply by `|E_χ(s)| ≤ q^{3ε}|E_χ(s+ε)|` (Lemma E). ∎

**Lemma A_χ (Lindelöf bound in a zero-free half-plane)** — changed: convexity bound with a `q`-constant; no slit if `δ = 0`;
negative `t` via `χ̄`. Let `1/2 ≤ θ < 1`, `L(s,χ) ≠ 0` for `Re s > θ`. For `δ' ∈ (0,(1−θ)/2]`, `η' > 0` there is
`C_A(θ,δ',η',χ)` such that for `σ ≥ θ + 2δ'`: `|L(s,χ)|^{−1} ≤ C_A(2+|t|)^{η'}`, and `|L(s,χ)| ≤ C_A(2+|t|)^{η'}` if also
`δ = 0` or `|t| ≥ 1`.

*Proof.* (0) *Convexity* (replaces [C6]): for `σ ≥ 1/2`, `|t| ≥ 1`: `|L(s,χ)| ≤ C_6(χ)(|t|+2)^{1/2}`, with
`C_6(χ) := 2 + 5C_κ` if `δ = 0` and `q C_6` if `δ = 1`. Indeed for `δ = 0`, `σ ≤ 3`, `x := |t| + 2` in [C1_χ]:
`|Σ_{n≤x}| ≤ 2√x`, `|S_χ(x)x^{−s}| ≤ C_κ`, `|s∫_x^∞| ≤ |s| C_κ x^{−σ}/σ ≤ 2(2x) C_κ x^{−1/2}` (`|s| ≤ 3 + |t| ≤ 2x`);
for `σ > 3`, `|L| ≤ ζ(3) < 2`. For `δ = 1`, `|L| = |ζ||E_χ| ≤ q C_6 |t|^{1/2}` ([C6], Lemma E).
(1) *Sign of t.* The zeros of `L(·,χ̄)` are the conjugates of those of `L(·,χ)` [C8_χ], so `χ̄` satisfies the hypothesis,
and `|L(σ − it, χ)| = |L(σ + it, χ̄)|`; prove `t ≥ 0` for both `χ` and `χ̄` and take the larger constant.
(2) *Small t*: compactness as in [NBLp] (for `δ = 0` both `L` and `1/L` are continuous on `{θ+2δ' ≤ σ ≤ 3, |t| ≤ 6}`).
(3) *Large t* (`t ≥ 6`): [NBLp] Lemma A steps (1)–(2) verbatim with `ζ → L(·,χ)`, `Ω := {Re s > θ}` if `δ = 0` (no
pole), `Ω := {Re s > θ} ∖ (θ,1]` if `δ = 1`, `log L` the branch equal to the Euler logarithm on `Re s > 1`,
`U := log C_6(χ) + ½ log(t + 4)` (heights in the Borel–Carathéodory disc are `≤ t + 3/2`), `|log L(2+it,χ)| ≤ log ζ(2) < 1`,
and `|log L| ≤ log ζ(1+δ')` on `Re z ≥ 1 + δ'` [C8_χ]. Only `C_7, C_9, C_10` change (they depend on `C_6(χ)`). ∎

**Lemma P_χ (Riesz–Perron)** — changed: `μ → μχ`, Lemma A_χ. Assume `L(s,χ) ≠ 0` for `Re s > θ`. For `η ∈ (0,1)`,
`δ' ∈ (0,(1−θ)/2]`, `η' ∈ (0,1/2]`: if `σ − η ≥ θ + 2δ'` then `|P_N(s) − 1/L(s,χ)| ≤ C_P N^{−η}(2+|t|)^{η'}` (`N ≥ 1`). If
`σ − η > 1`, no hypothesis is needed and `C_P = 1/(1−η) + 1/η`.

*Proof.* [NBLp] Lemma P verbatim: `P_N(s) = (2πi)^{−1}∫_{(c)} L(s+w,χ)^{−1} N^w dw/(w(w+1))` (the Dirichlet series of
`1/L` has coefficients `μ(k)χ(k)`, `|μχ| ≤ |μ|`), shift to `Re w = −η`; the only pole is `w = 0` with residue `1/L(s,χ)`
(removable, value 0, at `s = 1` if `δ = 1`); horizontal integrals vanish by Lemma A_χ and [C8_χ]. ∎

## 6. Proof of Theorem B_χ — changed (the `r_1` correction only if `δ = 1`; constants depend on `χ`)

Choose `σ_0, δ'` as in [NBLp §6]. `F_N := 𝟙 + Σ_{k≤N} c_k λ_k`, `c_k := μ(k)χ(k)(1 − k/N)`;
`Φ_N := 𝟙₊ + Σ_{k≤N} c_k r_k` on `(0,∞)`. Since `Σ c_k λ_k = Σ c_k r_k − δ(Σ c_k/k) r_1` and `Σ_k c_k/k = P_N(1)`:
`F_N = Φ_N|_{(0,1)} − δ P_N(1) r_1|_{(0,1)}`, so (Lemma M, `‖r_1‖_{L^p(0,1)} ≤ C_κ`)
`‖F_N‖_p ≤ (1 − r(σ_0 − ½))^{−1/r} ‖Φ_N‖_{L²_{σ_0}} + δ C_κ |P_N(1)|`. By Lemma T_χ,
`𝓜_{σ_0}Φ_N = (1 − L(s,χ)P_N(s))/s = L(s,χ)(1/L(s,χ) − P_N(s))/s`; Lemma A_χ (`η' = 1/8`; for `|t| ≤ 1` continuity,
`σ_0 < 1`) and Lemma P_χ give `‖Φ_N‖²_{L²_{σ_0}} ≤ C N^{−2η}`; if `δ = 1`, Lemma P_χ at `s = 1` gives
`|P_N(1)| ≤ 2^{1/8} C_P N^{−η}`. ∎

## 7. Proof of Theorem A_χ

(ii) ⟹ (i) is Proposition N_χ; (ii) ⟺ (iii) is Proposition C_χ(d). Assume (i); fix `ε ∈ (0, ε_0]`,
`σ' := 1/p − ε/2`. As in [NBLp]: `σ' ≥ 1/2`; `σ' + ε < 1`; `2σ' − 1 ≥ (2 − p)/(2p)`; `1 − ε > σ'`; `ε ≤ 1/2`.

**A1 (the functions)** — changed. `φ_ε(j) := Σ_{d|j} μ(d)d^{−ε}` (unchanged);
`E_{ε,χ}(y) := Σ_{j≤y} χ(j)φ_ε(j) − δ c_q y / L(1+ε,χ)`; `Φ_ε(x) := 𝟙₊(x) − E_{ε,χ}(1/x)` (`x > 0`). For `y = 1/x`:
`Σ_{k≥1} μ(k)χ(k)k^{−ε} κ_χ(y/k) = Σ_{k≤y} μ(k)χ(k)k^{−ε}κ_χ(y/k) + δc_q y Σ_{k>y} μ(k)χ(k)k^{−1−ε}` (for `k > y`,
`S_χ(y/k) = 0`) converges absolutely, and equals `δc_q y/L(1+ε,χ) − Σ_{k≤y} μ(k)χ(k)k^{−ε}S_χ(y/k) = −E_{ε,χ}(y)`,
because `Σ_{k≤y} μ(k)χ(k)k^{−ε} Σ_{n≤y/k} χ(n) = Σ_{j≤y} χ(j)φ_ε(j)` (complete multiplicativity). So pointwise
`Φ_ε = 𝟙₊ + Σ_k μ(k)χ(k)k^{−ε} r_k`. Put `H_ε := Σ_{k≥1} μ(k)χ(k)k^{−ε} λ_k` (pointwise on `(0,1)`); then on `(0,1)`
`H_ε = Φ_ε − 1 − δ r_1/L(1+ε,χ)`, and `H_ε(x) = −Σ_{j≤m} χ(j)φ_ε(j) + δ S_χ(m)/L(1+ε,χ)`, `m = ⌊1/x⌋` (a step function;
for `q = 1` it is [NBLp]'s `−E_ε(m)`). If `δ = 0`, `Φ_ε = 𝟙₊ + H_ε` vanishes on `(1,∞)`. For `N ≥ 1`:
`H_{ε,N} := Σ_{k≤N} μ(k)χ(k)k^{−ε}(1 − k/N) λ_k ∈ V_χ`, `Φ_{ε,N} := 𝟙₊ + Σ_{k≤N} μ(k)χ(k)k^{−ε}(1−k/N) r_k`, and on `(0,1)`
`H_{ε,N} = Φ_{ε,N} − 1 − δ P_N(1+ε) r_1`.

**A2 (transform on a right line)** — changed only through Lemma T_χ. For `b ∈ (1 − ε, 1)`:
`Σ_k |μ(k)χ(k)| k^{−ε}‖r_k‖_{L²_b} ≤ ‖r_1‖_{L²_b} Σ_k k^{−ε−b} < ∞`, so `Φ_ε ∈ L²_b` and
`𝓜_bΦ_ε(t) = Ψ_ε(b+it)`, `Ψ_ε(s) := (1 − R_ε(s))/s` (`Σ_k μ(k)χ(k)k^{−s−ε} = 1/L(s+ε,χ)`, `Re(s+ε) > 1`).

**A3 (shift to the left)** — unchanged argument; bounds via Lemmas H_χ, Γ_χ. `Ψ_ε` is analytic on
`{1/p − ε < Re s < 1}` (`L(s+ε,χ) ≠ 0` there by (i); `L(s,χ)` analytic for `s ≠ 1`). On `{σ' ≤ σ ≤ b, |t| ≥ 1}`:
`|R_ε(s)| ≤ q^{3ε}|Γr_ε(s)| ≤ q^{3/2} e^{5/2} q^{1/4} (3 + |t|)^{1/4}`; on `{σ' ≤ σ ≤ b, |t| ≤ 1}`, `R_ε` is continuous
(`b < 1`), hence `≤ C_13(ε,χ)`. So `sup_{σ'≤σ≤b} ∫|Ψ_ε|² dt < ∞` and Lemma PW gives `Φ_ε ∈ L²_a`,
`𝓜_aΦ_ε = Ψ_ε(a + i·)` for every `a ∈ [σ', b]`.

**A4 (key estimate `B(ε) := ‖Φ_ε‖²_{L²_{σ'}} ≤ K_{p,χ} ε`)** — changed (pole term only if `δ = 1`; new step for
`R_ε(2σ')` if `δ = 0`; `q`-constants). On `Re s = σ'`, `R := R_ε(s)`:
`2πB(ε) = ∫|1 − R|²|s|^{−2} dt ≤ ∫|s|^{−2}dt − 2 Re J + ∫|R|²|s|^{−2}dt`, `J := ∫_ℝ R(σ'+it)|s|^{−2} dt` (`J` may be complex).
(a) `∫ dt/(σ'² + t²) = π/σ'` — unchanged.
(b) `J = (1/i)∫_{(σ')} f(s) ds`, `f(s) := R(s)/(s(2σ' − s))`. In `Re s ≥ σ'`, `f` is meromorphic with poles only at
`s = 2σ'` and, if `δ = 1`, `s = 1` (`L(s+ε,χ) ≠ 0` as `Re(s+ε) > 1/p`; `R(1−ε) = 0` if `δ = 1`). Rectangle
`[σ', X] × [−T, T]`: on the horizontal sides `|f| ≤ q^{3ε}e^{5ε}q^{ε/2}(2 + X + T)^{ε/2}/T² → 0` (Lemma H_χ on all of
`Re s ≥ σ'`, Lemma Γ_χ); on `Re s = X`: `|R| ≤ ζ(X)ζ(X+ε)` [C8_χ] and the integral is `≤ ζ(X)²π/(X − 2σ') → 0`;
the left integral converges absolutely (A3). Residues: `Res_{s=2σ'} f = −R(2σ')/(2σ')` and, if `δ = 1`,
`Res_{s=1} f = c_q/(L(1+ε,χ)(2σ' − 1))` (`Res_{s=1}L(s,χ_0) = c_q`). The rectangle is clockwise relative to the
right-hand region, so `∫_{(σ')} f ds = −2πi Σ Res` and
`J = 2π [ R(2σ')/(2σ') − δ c_q/((2σ' − 1) L(1+ε,χ)) ]`.
(c) Lemmas H_χ and Γ_χ (`Re(s+ε) < 1`, `q* ≤ q`): `|R|² ≤ q^{6ε}e^{8ε}q^{ε}(2+|t|)^{ε} = e^{x}`,
`x := ε(8 + 7 log q + log(2+|t|))`, so `|R|² − 1 ≤ x e^{x} ≤ ε e^4 q^{7/2}(8 + 7 log q + log(2+|t|))(2+|t|)^{1/2}`
(`ε ≤ 1/2`) and `∫|R|²|s|^{−2}dt ≤ π/σ' + C_{Γ,q} ε`,
`C_{Γ,q} := e^4 q^{7/2} ∫_ℝ (8 + 7 log q + log(2+|t|))(2+|t|)^{1/2}(1/4 + t²)^{−1} dt` (`C_{Γ,1} = C_Γ`).
Adding: `2πB(ε) ≤ (2π/σ')(1 − Re R(2σ')) + δ · 4π c_q/((2σ' − 1)L(1+ε,χ)) + C_{Γ,q}ε`.
*δ = 1* (unchanged): `R(2σ') ≥ 1` ([C8_χ]: `L(·,χ_0)` positive decreasing on `(1,∞)`); `L(1+ε,χ_0) = ζ(1+ε)Π_{p|q}(1 − p^{−1−ε}) ≥ c_q/ε`;
so `B(ε) ≤ (4p/(2−p) + C_{Γ,q}/(2π)) ε`.
*δ = 0* (new step): no pole term. With `u := 2σ' ≥ u_0 := 1/2 + 1/p > 1` (as `ε ≤ (2/p − 1)/2`):
`|L(u,χ) − L(u+ε,χ)| ≤ Σ_n n^{−u}(1 − n^{−ε}) ≤ ε|ζ'(u)| ≤ ε|ζ'(u_0)|` and `|L(u+ε,χ)|^{−1} ≤ ζ(u_0)`, so
`1 − Re R(2σ') ≤ |1 − R(2σ')| ≤ ε|ζ'(u_0)|ζ(u_0)` and, with `σ' ≥ 1/2`, `B(ε) ≤ (2|ζ'(u_0)|ζ(u_0) + C_{Γ,q}/(2π)) ε`.
Let `K_{p,χ}` be the bracket of the relevant case. (`|ζ'(u)| ≤ (u−1)^{−2} + 0.35`, `ζ(u) ≤ u/(u−1)`.)

**A5 (`H_ε` approximates `−𝟙`)** — changed (correction only if `δ = 1`). On `(0,1)`: `𝟙 + H_ε = Φ_ε − δ r_1/L(1+ε,χ)`; Lemma M
(`a = σ'`), A4, `‖r_1‖_{L^p(0,1)} ≤ C_κ ≤ φ(q)`, `1/L(1+ε,χ_0) ≤ ε/c_q`:
`‖𝟙 + H_ε‖_p ≤ (2/(rε))^{1/r}(K_{p,χ}ε)^{1/2} + δ q ε = A_{p,χ} ε^{1−1/p} + δqε`.

**A6 (`H_{ε,N} → H_ε`)** — unchanged argument, constants depend on `χ`. `σ'' := 1/p − ε/4`. By A3 and Lemma T_χ,
`Φ_ε − Φ_{ε,N} ∈ L²_{σ''}` with transform `(L(s,χ)/s)(P_N(s+ε) − 1/L(s+ε,χ))`. On this line:
`|L(s,χ)| ≤ q^{3ε/2}|Γr_{ε/2}(s)||L(s + ε/2,χ)|` (Lemma H_χ with `ε/2`); `Re(s + ε/2) = 1/p + ε/4 = θ + 2δ'` with `θ = 1/p`,
`δ' = ε/8`, so Lemma A_χ (`η' = 1/8`) bounds `|L(s+ε/2,χ)| ≤ C_A(2+|t|)^{1/8}` (`|t| ≥ 1`; continuity for `|t| ≤ 1`,
`1/p + ε/4 < 1`); Lemma Γ_χ: `|Γr_{ε/2}(s)| ≤ e^{5/4}q^{1/8}(3+|t|)^{1/8}`; Lemma P_χ (`σ = 1/p + 3ε/4`, `η = ε/2`):
`|P_N(s+ε) − 1/L(s+ε,χ)| ≤ C_P N^{−ε/2}(2+|t|)^{1/8}`. Hence `‖Φ_ε − Φ_{ε,N}‖²_{L²_{σ''}} ≤ C_15 N^{−ε}`. If `δ = 1`, Lemma
P_χ at `s = 1 + ε` (no hypothesis): `|1/L(1+ε,χ) − P_N(1+ε)| ≤ (1/(1−ε/2) + 2/ε)N^{−ε/2}`. With A1 and Lemma M:
`‖H_ε − H_{ε,N}‖_p ≤ (4/(rε))^{1/r}‖Φ_ε − Φ_{ε,N}‖_{L²_{σ''}} + δC_κ|1/L(1+ε,χ) − P_N(1+ε)| ≤ C_17(p,χ,ε) N^{−ε/2}`.

**A7 (conclusion; order of limits)** — unchanged. Fix `ε`; choose `N_0(p,χ,ε)` with `C_17 N_0^{−ε/2} ≤ ε`. For `N ≥ N_0`:
`d_{N,p}(χ) ≤ ‖𝟙 + H_{ε,N}‖_p ≤ A_{p,χ}ε^{1−1/p} + (1 + δq)ε`. `d_{N,p}` is non-increasing in `N`, so
`lim_N d_{N,p}(χ) ≤ A_{p,χ}ε^{1−1/p} + (1+δq)ε` for every `ε ∈ (0,ε_0]`: `N → ∞` is taken at fixed `ε`, and only then
`ε → 0`. Every constant depends on `ε` only through fixed-`ε` strips strictly right of `Re s = 1/p − ε`; nothing is
used on the line `Re s = 1/p` itself. ∎

## 8. Proof of the corollaries

1. [C9_χ] gives (i) for every `χ` and `p ≤ 8/7`; [C9'] for `p ≤ 1/θ_B` (`1/0.87338196661320752 = 1.14497440…`).
2. (i) ⟺ `Θ_χ ≤ 1/p`. For `p > 2`: `L(s,χ*)` has infinitely many non-trivial zeros, symmetric under `ρ ↦ 1 − ρ̄`
   [C2_χ] (for `χ* = 1`: `1/2 + 14.1347…i`), so one has `β ≥ 1/2 > 1/p`; it is a zero of `L(s,χ)`; apply Prop. N_χ.
3. GRH for `χ` ⟺ `Θ_χ = 1/2` ⟺ (i) for all `p < 2`; the rate is Theorem B_χ with `θ = 1/2`.
4. `L(s,χ_0) = ζ(s)Π_{p|q}(1 − p^{−s})` and the product has no zeros in `Re s > 0`. ∎

## 9. Remarks

**9.1 (principal χ: the side condition can be dropped).** For `δ = 1` let `V'_χ := span{r_k|_{(0,1)} : k ≥ 1} = V_χ + ℂ r_1`.
Then `𝟙 ∈ closure(V'_χ) ⟺ (i)` too. ⟸: `V_χ ⊂ V'_χ`. ⟹: let `L(w_0,χ_0) = 0`, `Re w_0 > 1/p`. `w_0 ∉ ℝ` (on
`(0,1)`, `ζ < 0` and `Π(1 − p^{−σ}) > 0`) and `w̄_0` is a zero. For `ℓ_w(F) := ∫_0^1 F x^{w−1}dx`:
`ℓ_w(λ_k) = 0`, `ℓ_w(r_1) = −L(w,χ_0)/w − c_q/(1 − w) = −c_q/(1−w)` at `w ∈ {w_0, w̄_0}`. The functional
`ℓ := (1 − w_0)ℓ_{w_0} − (1 − w̄_0)ℓ_{w̄_0}` is bounded on `L^p(0,1)`, kills `V'_χ`, and
`ℓ(𝟙) = 1/w_0 − 1/w̄_0 ≠ 0`. (For `q = 1` this is the `L^p(0,1)` analogue of Báez-Duarte's removal of
`Σ b_k/θ_k = 0`.) For `δ = 0` no side condition arises.

**9.2 (what 1608.07887 actually asks).** Problem 1 there (TeX lines 247–251) reads: "For any fixed `n ∈ ℕ` and `p ∈ (1,2]`,
determine the best approximation of `1_{(−∞,∞)}` by products of the form `L(1/p+it)A_n(1/p+it)` in `L²(ℝ, ω)`, where
`ω(t) = 1/(1/p² + t²)` is the weight function." It is a fixed-`n`, weighted-`L²` extremal problem (its answer is an
orthogonal projection, a Gram-determinant quotient), **not** an `L^p` natural-number problem; Theorem A_χ does not solve it.
Their "p" criterion (Thm GBNew/GBLp) is the `p = 2` Báez-Duarte criterion for the shifted series
`η(s) = L(s + 1/p − 1/2, χ)` (TeX lines 787–796), i.e. a Hilbert-space statement; the `L^p(0,1)` natural-number closure
`C^p_{η,ℕ}` is defined there (line 330) but used only for `p = 2`. Theorem A_χ is the natural-number (Báez-Duarte)
strengthening of their Theorem GNB in `L^p(0,1)`, `1 < p < 2`, for `η = L(·,χ)`.

**9.3 (their p-weighted family is the wrong family in L^p).** With D–O's `κ_{χ,p}(x) = βx^α − Σ_{k≤x} χ(k)k^{1/2−1/p}`
the Mellin transforms carry `η(s) = L(s + 1/p − 1/2, χ)`; in `L^p(0,1)` the annihilators `x^{w−1}` of Prop. N_χ
(`Re w > 1/p`, `η(w) = 0`) detect only zeros of `L` with `Re > 2/p − 1/2 > 1/p`, so that closure is governed by a different
half-plane (its sufficiency would need Theorem A for the shifted series; not attempted). The weights are a device to move the line `1/p` to `1/2` for `L²(0,∞)`; in
`L^p(0,1)` the unweighted `κ_χ` is the right family (Prop. N_χ).

**9.4 (uniformity in q).** All constants (`C_{Γ,q}`, `K_{p,χ}`, `N_0`) depend on `q`; `C_{Γ,q}` grows like `q^{7/2} log q`
(crude: Lemma E and `(q*)^ε` were bounded with `ε ≤ 1/2`). No uniformity in `χ` is claimed or needed.
