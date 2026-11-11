FORMALIZED: Theorem 1 (Σ_{n≤N} R(n)² ≥ (c₁ − ε)N^{3/2}(log N)^{1/2}, c₁ = (3/5)^{1/2}(4/5)²/(2π√2)) and Corollaries 1–4 (limsup |R(n)|/(n log n)^{1/4} ≥ 0.289; R(n) = Ω((n log n)^{1/4}) unconditionally; Σ m^{1/2}y_m² = ∞; Σ m^δ y_m² < ∞ ⟺ δ < 1/2; RH ⟺ Σ m^δ y_m² < ∞ for some/every δ ∈ [0, 1/2); Merlini–Sala–Sala's a·log n bound and O(n^ε), ε < 1/4, both fail), sorry-free, standard axioms only, from named hypotheses: Li's formula with real phases (RH), Riemann–von Mangoldt, Theorem 2's upper bound, ¬RH ⇒ radius of Σλ_n zⁿ < 1, and four Fourier estimates of proofs.md for the ψ∗ψ weights (Φ ≥ 0 with Φ(0) ≥ (1−η)N; rapid decay of Φ; the trend transform τ(θ) = log(2πθ)/θ² + O(log N); τ_T ≤ (1+η)T³L²/(6π)). Proved in Lean, not assumed: the Fourier-positivity/Cauchy–Schwarz core (steps 2–6), the favourable-sign step, the optimisation u = 3/5 (step 7), all corollaries, Lemma 0's elementary part, and the Fejér weights as a witness that the weight hypothesis is satisfiable.

# LiRate in Lean 4 + Mathlib (lilean, 2026-11-11)

Lean v4.35.0-rc2, Mathlib v4.35.0-rc2 (repo `lakefile.toml`). Files in `Riemann/Deep/LiRate/Lean/`: Core (253 lines), Asym (336), Reduce (56), Trend (175), Cor (295), Cor3 (245), Lemma0 (100), Fejer (112), Main (128, with the axiom prints); 1,700 lines. No `sorry`, `admit`, `axiom`, `native_decide`.

## Main statements (namespace `LiRate`)

`R(n) = lam n − t n`. Zeros: any family `θ : κ → ℝ` of real phases (RH: 1 − 1/ρ = e^{iθ_γ}); `Z T : Finset κ` = the zeros with γ ≤ T; `w N n = φ(n/N)`; `P` = RH (any `Prop`; instantiate with Mathlib's `RiemannHypothesis`).

```lean
-- Core.lean (proofs.md Theorem 1, steps 1–6)
def Phi N w x := w 0 + 2 * ∑ n ∈ Icc 1 N, w n * cos (n * x)            -- Φ(x) = Σ_{|n|≤N} φ(n/N)e^{inx}
def ip N w F G := w 0 * F 0 * G 0 + 2 * ∑ n ∈ Icc 1 N, w n * F n * G n -- ⟨F,G⟩, F, G even
def B N w a b := 2 * Phi N w a - Phi N w (a - b) - Phi N w (a + b)
def Wfun θ J n := ∑ j ∈ J, cos (n * θ j)                                -- W(n)
def trend N w t a := ip N w (fun n => cos (n * a)) t                    -- τ(θ)
def tail N w θ lam J a := ip N w (fun n => cos (n * a)) lam - ∑ k ∈ J, B N w a (θ k)
def tauPrime N w θ lam t J := 2 * J.card * ∑ j ∈ J, Phi N w (θ j)
    + ∑ j ∈ J, tail N w θ lam J (θ j) - ∑ j ∈ J, trend N w t (θ j)       -- exact residual τ'
lemma tail_eq_tsum_compl (hlam) : tail N w θ lam J a = ∑' k : ↑((J : Set κ)ᶜ), B N w a (θ k)
lemma ip_R_W : ip N w (fun n => lam n - t n) (Wfun θ J) = tauPrime N w θ lam t J - 2 * ip N w (Wfun θ J) (Wfun θ J)
theorem core_bound (hw0 : ∀ n, 0 ≤ w n) (hw1 : ∀ n, w n ≤ 1) (hΦ : ∀ x, 0 ≤ Phi N w x)
    (hlam : ∀ n : ℕ, HasSum (fun k => 2 * (1 - cos (n * θ k))) (lam n)) (ht0 : t 0 = 0)
    (J : Finset κ) (τ : ℝ) (hτ : tauPrime N w θ lam t J ≤ τ) (hτ0 : 0 ≤ τ) (hτM : τ < J.card * Phi N w 0) :
    (J.card * Phi N w 0 - τ) ^ 2 / (J.card * Phi N w 0) ≤ ∑ n ∈ Icc 1 N, (lam n - t n) ^ 2

-- Asym.lean (step 7)
def RvM Z := ∃ C T₀, ∀ T ≥ T₀, |#(Z T) - T/(2π) * log (T/(2π·e))| ≤ C * log T
def GoodWeights w η := (∀ N n, 0 ≤ w N n) ∧ (∀ N n, w N n ≤ 1) ∧ (∀ N x, 0 ≤ Phi N (w N) x)
    ∧ ∀ᶠ N in atTop, (1 - η) * N ≤ Phi N (w N) 0
def TrendTail w θ lam t Z η := ∀ᶠ N in atTop, ∀ T, √N / log N ≤ T → T ≤ √N / 2 →
    tauPrime N (w N) θ lam t (Z T) ≤ (1 + η) * (T^3 * log (T/(2π))^2 / (6π))
theorem thm1 (hlam) (ht0) (hN : RvM Z) (hW : ∀ η > 0, ∃ w, GoodWeights w η ∧ TrendTail w θ lam t Z η) :
    Thm1Bound (fun n => lam n - t n)

-- Reduce.lean, Trend.lean (TrendTail from its ingredients)
lemma tail_le (hΦ) (hJX : J ⊆ X) : tail N w θ lam J a ≤ ∑ k ∈ X \ J, (2 * Phi N w a - Phi N w (a + θ k)) + tail N w θ lam X a
def TinyTerms w θ lam Z Xs η :=   -- Z T ⊆ Xs N (γ ≤ X = N³), and for √N/log N ≤ T ≤ √N/2:
    2 #Z(T) Σ_{j∈Z T} Φ(θ_j) + Σ_{j∈Z T} (Σ_{k∈Xs N \ Z T} (2Φ(θ_j) - Φ(θ_j+θ_k)) + tail_{Xs N}(θ_j)) ≤ η/2 · T³L²/(6π)
def TrendPt w θ t Z := ∃ C, ∀ᶠ N, ∀ T ∈ [√N/log N, √N/2], ∀ j ∈ Z T, log (2π θ j)/θ j^2 - C log N ≤ trend N (w N) t (θ j)
def TauT θ Z := ∀ η > 0, ∀ᶠ T in atTop, -∑ j ∈ Z T, log (2π θ j)/θ j^2 ≤ (1 + η) T³L²/(6π)
theorem trendTail_of … : TrendTail w θ lam t Z η          -- from TinyTerms-form and the trend sum
theorem htrend_of (hN : RvM Z) (h4 : TrendPt w θ t Z) (h7 : TauT θ Z) (hη : 0 < η) :
    ∀ᶠ N, ∀ T ∈ [√N/log N, √N/2], -∑ j ∈ Z T, trend N (w N) t (θ j) ≤ (1 + η/2) T³L²/(6π)
theorem thm1_of_lemmas [DecidableEq κ] (hlam) (ht0) (hN : RvM Z) (h7 : TauT θ Z)
    (hW : ∀ η > 0, ∃ w Xs, GoodWeights w η ∧ TinyTerms w θ lam Z Xs η ∧ TrendPt w θ t Z) :
    Thm1Bound (fun n => lam n - t n)

-- Cor.lean, Cor3.lean, Lemma0.lean, Fejer.lean
def c₁ := √(3/5) * (4/5)^2 / (2 * π * √2)
def Thm1Bound R := ∀ ε > 0, ∀ᶠ N in atTop, (c₁ - ε) * (N * √N * √(log N)) ≤ ∑ n ∈ Icc 1 N, R n ^ 2
def Lemma0 R := ∃ c > 0, ∃ᶠ n in atTop, exp (c * n) ≤ |R n|
def Thm2Upper R := ∃ C N₀, ∀ N ≥ N₀, ∑ n ∈ Ioc N (2N), R n ^ 2 ≤ C * (N * √N) * (log N * √(log N))
def liTrend n := n/2 * (log n - (1 + log (2π) - eulerMascheroniConstant))
theorem cor1 (h : Thm1Bound R) (hc : c < √(3 * c₁ / 2)) : ∃ᶠ n in atTop, c * √(√(n * log n)) ≤ |R n|
theorem cor1_num (h : Thm1Bound R) : ∃ᶠ n in atTop, 0.289 * √(√(n * log n)) ≤ |R n|
theorem omega_uncond (h1 : P → Thm1Bound R) (h0 : ¬P → Lemma0 R) : ∃ c > 0, ∃ᶠ n in atTop, c * √(√(n * log n)) ≤ |R n|
theorem cor2 (h : Thm1Bound R) : ¬ Summable (fun m => √m * (R m / m) ^ 2)
theorem cor3 (h1 : Thm1Bound R) (h2 : Thm2Upper R) (δ : ℝ) : Summable (fun m => m ^ δ * (R m / m) ^ 2) ↔ δ < 1/2
theorem rh_iff_exists (h2 : P → Thm2Upper R) (h0 : ¬P → Lemma0 R) : P ↔ ∃ δ ∈ Set.Ico 0 (1/2), Summable (fun m => m ^ δ * (R m / m) ^ 2)
theorem rh_iff_forall (h2) (h0) : P ↔ ∀ δ ∈ Set.Ico 0 (1/2), Summable (fun m => m ^ δ * (R m / m) ^ 2)
theorem cor4 (h1 : P → Thm1Bound R) (h0 : ¬P → Lemma0 R) (hlt : ∃ C, ∀ n, |lt n - R n| ≤ C) (a : ℝ) :
    ¬ ∀ᶠ n in atTop, |lt n| ≤ a * log n
theorem cor4_rpow (…) (hε : ε < 1/4) (C' : ℝ) : ¬ ∀ᶠ n in atTop, |lt n| ≤ C' * n ^ ε
theorem lemma0_of_unbounded (hr0 : 0 < r) (hr1 : r < 1) (hunb : ¬ ∃ B, ∀ n, |lam n| * r ^ n ≤ B)
    (C) (ht : ∀ n, |t n| ≤ C * n ^ 2) : Lemma0 (fun n => lam n - t n)
theorem fejer_goodWeights (hη : 0 ≤ η) : GoodWeights fejer η      -- fejer N n = (1 - n/(N+1))₊
```

`Main.lean` bundles the inputs as `Thm1Inputs θ lam t Z` (H_Li, t(0) = 0, H_N, H_w+H_τ) and the finer `Thm1Lemmas θ lam t Z` (H_Li, t(0) = 0, H_N, H7, and H_w+H_tiny+H4 per η; `Thm1Lemmas.toInputs`), and states `li_thm1`, `li_cor1`, `li_cor2`, `li_cor3`, `li_omega`, `li_cor4`, `li_rh_iff`, and, for the actual trend `t = liTrend` with Lemma 0 derived from the radius input, `li_omega_radius`, `li_rh_iff_radius`.

## Hypotheses and sources

| Name | Lean form | Source |
|---|---|---|
| H_Li | `∀ n, HasSum (fun k => 2*(1 - cos (n*θ k))) (lam n)` | Li 1997, Bombieri–Lagarias 1999: λ_n = Σ_ρ[1 − (1 − 1/ρ)ⁿ]; under RH 1 − 1/ρ = e^{iθ_γ}, θ_γ = 2 arctan(1/(2γ)), ρ and ρ̄ paired, zeros with multiplicity (proofs.md lines 3–6) |
| H_N | `RvM Z` | Riemann–von Mangoldt, N(T) = (T/2π) log(T/2πe) + O(log T) (line 57) |
| H_w | `GoodWeights w η`, one family per η > 0 | lines 35–37: φ = ψ∗ψ, 0 ≤ φ ≤ φ(0) = 1, Φ = NΣ_kφ̂(N(x + 2πk)) ≥ 0 (Poisson, φ̂ = ψ̂² ≥ 0), Φ(0) = N∫φ + O(N^{−K}), ∫φ → 1. Satisfiable: `fejer_goodWeights` |
| H_tiny | `TinyTerms w θ lam Z Xs η` | lines 48–49, 52: X = N³; rapid decay \|Φ^{(j)}(x)\| ≪ N^{1+j}(1 + N‖x‖)^{−K} (line 37) at ‖x‖ ≥ θ_γ ≥ N^{−1/2} |
| H4 | `TrendPt w θ t Z` | step 4, lines 42–45: τ(θ) = log(2πθ)/θ² + O(log N) on [N^{−1/2}, 1] ∋ θ_γ for γ ≤ T ≤ √N/2 |
| H7 | `TauT θ Z` | line 51 and line 57: τ_T = −Σ_{γ≤T} log(2πθ_γ)/θ_γ² = T³L²/(6π)(1 + O(1/L)) (partial summation with N(T)) |
| H_2 | `Thm2Upper R` | Theorem 2, lines 63–76 (Lemma A, zero spacing, Schur, Selberg's ∫S²); the weaker bound (log N)^{3/2} |
| H_0 | `¬P → ∃ r ∈ (0,1), ¬∃ B, ∀ n, \|λ_n\| rⁿ ≤ B` (or `Lemma0 R` directly) | lines 9–11: Σλ_n zⁿ has the singularities ρ/(ρ−1) of log ξ(1/(1−z)) in \|z\| < 1 when Re ρ < 1/2; Cauchy–Hadamard for any r between the radius and 1 (the series of log ξ(1/(1−z)) has coefficients λ_n/n; same radius) |
| λ̃ | `∃ C, ∀ n, \|lt n - R n\| ≤ C` | lines 86–87, 97: λ̃_n = R(n) − 3/4 + O(1/n) |

`thm1` takes H_w and H_τ (= `TrendTail`); `thm1_of_lemmas` replaces H_τ by H_tiny, H4, H7 (with H_N). Lemma A is not an input of Theorem 1: proofs.md proves Theorem 1 from Li's formula directly; Lemma A enters only Theorems 2–4 (here: H_2).

## Fidelity (Lean ↔ proofs.md)

- `Phi`, `ip`: line 36 Φ(x) = Σ_{n∈ℤ}φ(n/N)e^{inx}, line 38 ⟨F,G⟩ = Σ_{n∈ℤ}φ(n/N)F(n)G(n), written on ℕ for even sequences (φ even, supp φ ⊂ [−1,1]): `w 0 + 2Σ_{n=1}^N` is the two-sided sum.
- `core_bound`: hypotheses `hw0`, `hw1` line 35, `hΦ` line 36; conclusion line 56, Σ_{n≤N}R(n)² ≥ Φ(0)N(T)(1 − τ'/(Φ(0)N(T)))², with `J` = {γ ≤ T}, `J.card` = N(T), τ any upper bound of the exact residual (0 ≤ τ < Φ(0)N(T) is line 55's 2D > τ').
- `lam_zero`, `ht0`: line 38 (λ_0 = 0, t(0) := 0), so ‖R̃‖² = 2Σ_{1≤n≤N}φ(n/N)R(n)² ≤ 2Σ_{n≤N}R(n)² (line 39).
- `li_pair`, `hasSum_B`: lines 40–41. `trend`: line 42 (real since φ, t even).
- `Wfun`, `ip_W_W`: lines 46, 52. `diag_le`: line 54. `ip_sq_le`: line 55 (Cauchy–Schwarz via the discriminant); the monotonicity of (2x − τ')²/x for x ≥ D: line 55, inside `core_bound`.
- `ip_R_W`, `tauPrime`: lines 47–53. `tauPrime` is the exact residual with ⟨R̃,W⟩ = τ' − 2‖W‖²: 2N(T)Σ_{γ≤T}Φ(θ_γ) + Σ_{γ≤T}Σ_{γ'>T}B(θ_γ,θ_{γ'}) − Σ_{γ≤T}τ(θ_γ) (`tail_eq_tsum_compl`). proofs.md's τ' (line 53) is an upper bound for it: it drops −Σ_{T<γ'≤X}Φ(θ_γ − θ_{γ'}) ≤ 0 (line 50; Lean: `tail_le`), and the cross terms Φ(θ_γ + θ_{γ'}), γ, γ' ≤ T, which line 52 calls tiny, cancel exactly in Lean.
- `thm1`, `step7`, `facts_at`: lines 57–59. T = √(6N/(5 log N)), i.e. T²·(½ log N) = (3/5)N, instead of T² log(T/2π) = (3/5)N∫φ; both give u = T²L/N → 3/5 (L = ½ log N + O(log log N), line 59). ∫φ → 1 is the η-family of `GoodWeights`. `const_bound`: (1 − η)³(1 − v)² ≥ 16/25 − (16/5)η, v = (1 + η)/(5(1 − η)²) → u/3; `c₁_eq`: c₁ = (16/25)√(6/5)/(4π) = √u(1 − u/3)²/(2π√2) at u = 3/5 (line 58); `opt`: (M − τ)²/M ≥ M₀(1 − v)² for M ≥ M₀ > 0, 0 ≤ τ ≤ vM₀, v < 1.
- `TrendTail`'s range √N/log N ≤ T ≤ √N/2: line 46 needs T ≤ √N/2; T ≥ √N/log N gives log N ≤ 3 log T, so N(T) log N = o(T³L²) uniformly (the chosen T lies in the range). `htrend_of` proves this absorption from RvM.
- `Thm1Bound`: line 30, N^{3/2}(log N)^{1/2} = N√N√(log N), "N ≥ N₀(ε)" = `∀ᶠ N in atTop`.
- `cor1`: line 32, limsup |R(n)|/(n log n)^{1/4} ≥ (3c₁/2)^{1/2}, stated equivalently: for every c < (3c₁/2)^{1/2}, |R(n)| ≥ c(n log n)^{1/4} infinitely often. `cor1_num`: 0.289 < (3c₁/2)^{1/2} = 0.28928 (`c₁_gt`: c₁ > 0.05569). Uses Σ_{n≤N}√n ≤ (2/3)N√N + √N (`sum_sqrt_le`).
- `omega_uncond`: line 32. `cor2`: line 33. `cor3`, `rh_iff_exists`, `rh_iff_forall`: lines 65–66 (⇐ by dyadic blocks from H_2: `block_bound`, `summable_of_dyadic`; ⇒ from `cor2`; RH ⇐ by Lemma 0, `not_summable_of_lemma0`).
- `cor4`: lines 95–97. `cor4_rpow`: discoveries/li-remainder.md claim 5 (O(n^ε) fails for ε < 1/4), a further consequence, not in proofs.md's Corollary 4.
- `liTrend`: line 7 (κ = 1 + log 2π − γ_E, Mathlib's `eulerMascheroniConstant`); `liTrend_bound`: |t(n)| ≤ 4n² (line 11, "t(n) is polynomial").
- `lemma0_of_unbounded`: lines 9–11 (c = −½ log r).
- `fejer_goodWeights`: φ(u) = (1 − |u|)₊ = ψ∗ψ, ψ = 1_{(−1/2,1/2)} (line 36's limit): Φ = Fejér kernel, (N+1)Φ(x) = |Σ_{k≤N}e^{ikx}|² ≥ 0, Φ(0) = N + 1. It does not satisfy H_tiny (Φ ≍ 1/(Nx²)); proofs.md needs smooth ψ there.

Review 2 (reviewer subagent, Trend/Fejer parts): FAITHFUL WITH CAVEATS. `TrendPt` and `TauT` are the correct one-sided consequences of step 4 and of τ_T ~ T³L²/(6π); `TinyTerms` is one combined signed remainder (a single one-sided bound, weaker than separate absolute bounds of the tiny terms); `Xs N` is only required to contain `Z T`, the cutoff X = N³ is the intended instantiation, not encoded; `fejer_goodWeights` is a fair non-vacuity witness for `GoodWeights` alone (no joint witness with `TinyTerms`).

Review 1 (reviewer subagent, statement fidelity, before Trend/Fejer were added): FAITHFUL WITH CAVEATS. Caveats: `tauPrime` is the exact residual, proofs.md's τ' an upper envelope (docstrings now say so); `cor4_rpow` goes beyond proofs.md's Corollary 4; `trendTail_of` proves only the reduction, not the analytic estimates; "no vacuity or joint contradiction identified"; the results are conditional, with no instantiation by the actual zeros of ζ.

Not formalized (the hypotheses above): Li's formula, Riemann–von Mangoldt, the smooth ψ∗ψ weights (Poisson summation, φ̂ = ψ̂², rapid decay), the trend transform (step 4), the τ_T asymptotics (partial summation), Theorem 2, ¬RH ⇒ radius < 1, the instantiation by ζ's zeros. Theorem 3, Theorem 4, Lemma A and the barrier are not formalized.

## Build

`.lake` was a real directory (jerbylean's Mathlib cache, 7.4 GB): reused, left in place, never deleted or re-fetched. Each module compiled singly (one process at a time; peak 454 MB with the targeted Mathlib imports, 799 MB in a first attempt with `import Mathlib`), oleans to `.scratch/lilean/olean`, not `.lake/build` (`.scratch/lilean/mk.sh`):

```
O=$PWD/.scratch/lilean/olean
for n in Core Cor Cor3 Asym Reduce Trend Lemma0 Fejer Main; do
  lake env sh -c "LEAN_PATH=\"$O:\$LEAN_PATH\" lean -o $O/Riemann/Deep/LiRate/Lean/$n.olean Riemann/Deep/LiRate/Lean/$n.lean"
done
```

Output: every module exit 0, no errors, no warnings; Main's `#print axioms`:

```
'LiRate.core_bound' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.tail_eq_tsum_compl' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.trendTail_of' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.thm1' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.htrend_of' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.thm1_of_lemmas' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.Thm1Lemmas.toInputs' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.li_thm1' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.li_cor1' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.li_cor2' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.li_cor3' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.li_omega' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.li_cor4' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.li_rh_iff' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.lemma0_of_unbounded' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.li_omega_radius' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.li_rh_iff_radius' depends on axioms: [propext, Classical.choice, Quot.sound]
'LiRate.fejer_goodWeights' depends on axioms: [propext, Classical.choice, Quot.sound]
EXIT 0
```

## FORMALIZATIONS.md entry

- LiRate (our discovery [li-remainder](discoveries/li-remainder.md)) in Lean 4 + Mathlib: Theorem 1, Σ_{n≤N}R(n)² ≥ (c₁ − ε)N^{3/2}√log N (LiRate.thm1_of_lemmas; the Fourier-positivity/Cauchy–Schwarz core and u = 3/5 proved), and Corollaries 1–4 (limsup ≥ 0.289, unconditional Ω((n log n)^{1/4}), RH ⟺ Σ m^δ y_m² < ∞ for some/every δ < 1/2, Merlini–Sala–Sala refuted), sorry-free, standard axioms, 1,700 lines; inputs as named hypotheses: Li's formula, Riemann–von Mangoldt, Theorem 2's bound, ¬RH ⇒ radius < 1, and proofs.md's Fourier estimates for smooth ψ∗ψ weights. Reviewer (two passes): FAITHFUL WITH CAVEATS. [LEAN.md](Riemann/Deep/LiRate/LEAN.md)

Build logs (per module, with #print axioms): logs/lean_*.log, logs/lean_mk.sh.
