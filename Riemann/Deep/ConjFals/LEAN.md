FORMALIZED: lim_{a→∞} s_a = 1, the conjecture of Johnston–Leong–Tudzi (arXiv 2408.04143v3, TeX l. 1142–1147), where s_a = limsup_{x→∞} |W_a(x)|/x^{log₂ a} and W_a(x) = Σ_{n≤x} (−a)^{Ω(n)}. It is stated with the literal limsup. Also proved: L(a) ≤ s_a ≤ U(a) for a ≥ 64, the limsup-free form (i)/(ii), and the exact 2-adic identity. No sorry, standard axioms only, no native_decide. Lean v4.35.0-rc2, Mathlib v4.35.0-rc2 (rev 0653561).

[Discovery](../../../discoveries/omega-power-sums.md) · [REPORT §2](REPORT.md)

## Main theorem (`Lean/Main.lean`)

```lean
theorem OmegaLean.tendsto_s : Filter.Tendsto OmegaLean.s Filter.atTop (nhds 1)
```

## Definitions

```lean
noncomputable def W (a x : ℝ) : ℝ := ∑ n ∈ Icc 1 ⌊x⌋₊, (-a) ^ (ArithmeticFunction.cardFactors n)
noncomputable def s (a : ℝ) : ℝ := limsup (fun x : ℝ => |W a x| / x ^ logb 2 a) atTop
noncomputable def U (a : ℝ) : ℝ := (a / (a - 1)) * (1 + 9 * a ^ (-(1/3 : ℝ)))
noncomputable def L (a : ℝ) : ℝ := a / (a + 1) - (a / (a - 1)) * (9 * a ^ (-(1/3 : ℝ)))
```

How each definition matches the paper:
- Ω is Mathlib's `cardFactors` (prime factors counted with multiplicity).
- The sum runs over 1 ≤ n ≤ ⌊x⌋.
- x^α is `Real.rpow`, with α = `Real.logb 2 a` = log a / log 2.
- For x < 1, W = 0, so the ratio is ≥ 0 for every real x (`ratio_nonneg`).
- For a ≥ 64 the ratio is at most U a once x ≥ 1. So s a is a genuine (bounded, cobounded) limsup, not a junk value.
- For a < 64 Lean's s a may be a junk value. The limit never looks at those a.

A reviewer subagent was given the definitions, the statements and the axiom output. Its verdict: FAITHFUL, "no unbounded-limsup junk value drives the conclusion".

## Other theorems

| Theorem | Statement |
|---|---|
| `s_le`, `le_s` (Main) | 64 ≤ a → s a ≤ U a and L a ≤ s a |
| `explicit_bounds` (Main) | ∀ ε > 0, ∃ a₀, ∀ a ≥ a₀: (∀ x ≥ 1, \|W a x\| ≤ (1+ε) x^α) ∧ (∀ᶠ k : ℕ, (1−ε)(2^k)^α ≤ \|W a 2^k\|). This is stronger than the brief's (i) eventually and (ii) frequently. |
| `W_two_pow_ge` (Main) | 64 ≤ a → L a · a^k − 1 ≤ \|W a 2^k\| |
| `abs_WN_le` (Bounds) | 1 < a, 6 ≤ log₂ a → \|W(N)\| ≤ (a/(a−1))(1 + 9a^{−1/3}) N^α, for every N |
| `abs_WN_two_pow_ge` (Bounds) | same hypotheses → \|W(2^k)\| ≥ (a^{k+1}−1)/(a+1) − (a/(a−1))·9a^{−1/3}·a^k |
| `WN_eq` (Identity) | W(N) = Σ_{m≤N odd} (−a)^{Ω(m)} Σ_{j<N, m2^j≤N} (−a)^j |
| `WN_two_adic` (TwoAdic) | a ≠ −1 → W(N) = (1/(1+a)) Σ_{m≤N odd} (−a)^{Ω(m)} + (a/(1+a)) Σ_{m≤N odd} (−a)^{Ω(m)} (−a)^{Nat.log 2 (N/m)}. This is the brief's identity (1): Nat.log 2 (N/m) is the largest j with m2^j ≤ N (`mul_two_pow_le_iff`), i.e. ⌊log₂(x/m)⌋. |
| `three_pow_cardFactors_le`, `two_pow_three_mul_le` (Identity) | odd m → 3^{Ω(m)} ≤ m and 2^{3Ω(m)} ≤ m² |

## Proof as formalized

This route is simpler than the brief's plan. All sums are finite: there are no tsums, Euler products or dominated-convergence steps.
1. W(N) = O(N) − a·W(⌊N/2⌋), with O the sum over odd n (n = 2n′ and Ω(2n′) = Ω(n′) + 1). Iterating this and swapping the sums (`Finset.sum_comm'`) gives `WN_eq`.
2. Inner sum: |Σ_{j∈F} (−a)^j| ≤ Σ_{j∈F} a^j ≤ a^K·a/(a−1) ≤ (N/m)^α·a/(a−1), with K = max F. No cancellation is needed, because a/(a−1) → 1.
3. For odd m, 3^{Ω(m)} ≤ m, so 2^{3Ω(m)} ≤ 9^{Ω(m)} ≤ m². For odd m ≥ 3 and α ≥ 6 this gives a^{Ω(m)} m^{−α} ≤ 9a^{−1/3} m^{−2} (a linear inequality between logarithms). Then Σ_{m≥2} m^{−2} ≤ 1 (Mathlib `sum_Ioo_inv_sq_le`).
4. Upper bound: |W(x)| ≤ U(a)·x^α for every x ≥ 1.
5. Lower bound at x = 2^k: the m = 1 term is Σ_{j≤k} (−a)^j, of size at least (a^{k+1}−1)/(a+1). The odd m ≥ 3 terms total at most (a/(a−1))·9a^{−1/3}·a^k.
6. Conclude with `limsup_le_of_le` and `le_limsup_of_frequently_le`, then squeeze between L → 1 and U → 1.

## Files (`Riemann/Deep/ConjFals/Lean/`, 655 lines)

| File | Lines | Content |
|---|---|---|
| `Identity.lean` | 118 | WN, ON, E, the recursion and its iteration, `WN_eq`, 3^{Ω(m)} ≤ m |
| `Bounds.lean` | 257 | the explicit bounds on WN (imports Identity) |
| `Main.lean` | 195 | W, s, U, L, the limsup bounds, `tendsto_s`, `explicit_bounds` (imports Bounds) |
| `TwoAdic.lean` | 65 | the closed-form identity (imports Identity; not needed by Main) |
| `Axioms.lean` | 20 | `#check`, `#print`, `#print axioms` |

## Build

Run from the repo root. `.lake` must be a real directory: delete the empty placeholder file, then run `lake exe cache get`.

```
lake build Riemann.Deep.ConjFals.Lean.Main
lake build Riemann.Deep.ConjFals.Lean.TwoAdic
lake env lean Riemann/Deep/ConjFals/Lean/Axioms.lean
```

Output of the clean rebuild (my oleans deleted first), 11 Nov 2026, under `taskpolicy -b`:

```
✔ [2401/2403] Built Riemann.Deep.ConjFals.Lean.Identity (50s)
✔ [2402/2403] Built Riemann.Deep.ConjFals.Lean.Bounds (129s)
✔ [2403/2403] Built Riemann.Deep.ConjFals.Lean.Main (271s)
Build completed successfully (2403 jobs).
✔ [1214/1214] Built Riemann.Deep.ConjFals.Lean.TwoAdic (42s)
Build completed successfully (1214 jobs).
'OmegaLean.tendsto_s' depends on axioms: [propext, Classical.choice, Quot.sound]
'OmegaLean.explicit_bounds' depends on axioms: [propext, Classical.choice, Quot.sound]
'OmegaLean.le_s' depends on axioms: [propext, Classical.choice, Quot.sound]
'OmegaLean.s_le' depends on axioms: [propext, Classical.choice, Quot.sound]
'OmegaLean.WN_eq' depends on axioms: [propext, Classical.choice, Quot.sound]
'OmegaLean.WN_two_adic' depends on axioms: [propext, Classical.choice, Quot.sound]
```

- Warnings and errors: none.
- `rg 'sorry|native_decide|admit|^axiom|implemented_by|@\[extern'` over the Lean files: no match.
- Peak memory footprint (from `/usr/bin/time -l`): at most 173 MB per build.
- Imports are targeted (no `import Mathlib`). Each module loads about 0.25–0.55 GB of oleans.

## Prior formalizations

None found, from a GitHub code search on 11 Nov 2026:
- `2408.04143`: only TauCeti roadmap notes on Johnston–Leong–Tudzi's Möbius bounds (Thms A.1–A.2).
- `Tudzi`: nothing relevant.
- `cardFactors` with `limsup`/log₂: unrelated files.
- `"(-a) ^" cardFactors` in Lean: 0 hits.
- `cardFactors` in openai/math, facebookresearch/atlas-lean and anthropics/formal-math (zeta23): other topics only.

## FORMALIZATIONS.md entry

- Johnston–Leong–Tudzi's conjecture lim_{a→∞} s_a = 1 (our discovery [omega-power-sums](discoveries/omega-power-sums.md)) machine-checked in Lean 4 + Mathlib: OmegaLean.tendsto_s : Tendsto s atTop (𝓝 1), with s a = limsup_{x→∞} |Σ_{n≤x}(−a)^{Ω(n)}|/x^{log₂ a} stated literally. It also proves L(a) ≤ s_a ≤ U(a) for a ≥ 64 (both bounds → 1) and the 2-adic identity. Sorry-free, standard axioms, no native_decide, 655 lines. [LEAN.md](Riemann/Deep/ConjFals/LEAN.md) Reviewer (statement fidelity): FAITHFUL. Online: no prior formalization found.
