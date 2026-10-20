import Mathlib.Analysis.SpecialFunctions.Gamma.Basic
import Mathlib.Analysis.Calculus.LogDeriv
import Mathlib.Analysis.SpecialFunctions.Log.Deriv
import Mathlib.NumberTheory.ArithmeticFunction.VonMangoldt
import Mathlib.NumberTheory.LSeries.RiemannZeta
import Mathlib.MeasureTheory.Integral.Bochner.Basic
import Mathlib.Analysis.SpecialFunctions.Log.Basic
import Mathlib.Analysis.Complex.ExponentialBounds

/-!
# R3.1 deep dive — Weil positivity on a window

Explicit-formula functional (Guinand–Weil), conventions checked numerically against zeta zeros
(`.scratch/weil/check_ef.py`, `check_zeros.py`):

  Σ_γ h(γ) = h(i/2) + h(−i/2) + (1/2π)∫ h(r) Re ψ(1/4 + ir/2) dr − g(0) log π
             − Σ_n Λ(n) n^{-1/2} (g(log n) + g(−log n)),

with `g = f * f̃` (autocorrelation) and `h(z) = ∫ g(x) e^{izx} dx`.
`WeilPos L`: the functional is ≥ 0 for every smooth `f` supported in `[−L/2, L/2]`
(so `g` is supported in `[−L, L]` and only prime powers `n ≤ e^L` enter).

Proved here: monotonicity in the window; only primes `n ≤ e^{2a}` enter; no prime enters below `log 2`.
Cited theorems are `def … : Prop` (never axioms).
-/

namespace Riemann.Deep.Weil

open MeasureTheory Complex ArithmeticFunction

/-- Autocorrelation `g = f * f̃`, `g x = ∫ f y · conj (f (y − x)) dy`. -/
noncomputable def autocorr (f : ℝ → ℂ) (x : ℝ) : ℂ :=
  ∫ y : ℝ, f y * starRingEnd ℂ (f (y - x))

/-- Fourier–Laplace transform `h z = ∫ g x · e^{izx} dx`. -/
noncomputable def ft (g : ℝ → ℂ) (z : ℂ) : ℂ :=
  ∫ x : ℝ, g x * cexp (I * z * x)

/-- Pole terms `h(i/2) + h(−i/2)` (from `s = 1` and `s = 0`). -/
noncomputable def poleTerm (g : ℝ → ℂ) : ℂ :=
  ft g (I / 2) + ft g (-(I / 2))

/-- Archimedean term `(1/2π) ∫ h(r) Re ψ(1/4 + ir/2) dr − g(0) log π`, `ψ = Γ'/Γ`. -/
noncomputable def archTerm (g : ℝ → ℂ) : ℂ :=
  (1 / (2 * Real.pi) : ℂ) *
      ∫ r : ℝ, ft g r * ((logDeriv Complex.Gamma (1 / 4 + I * r / 2)).re : ℂ) -
    g 0 * (Real.log Real.pi : ℂ)

/-- Prime term `Σ_n Λ(n) n^{-1/2} (g(log n) + g(−log n))`. -/
noncomputable def primeTerm (g : ℝ → ℂ) : ℂ :=
  ∑' n : ℕ, ((Λ n / Real.sqrt n : ℝ) : ℂ) * (g (Real.log n) + g (-Real.log n))

/-- The Weil functional `Q(f)` (equals `Σ_ρ h(γ_ρ)` by the explicit formula). -/
noncomputable def weilQ (f : ℝ → ℂ) : ℂ :=
  poleTerm (autocorr f) + archTerm (autocorr f) - primeTerm (autocorr f)

/-- Test functions: smooth (`C^∞`, not analytic), supported in `[−a, a]`. -/
def TestFn (a : ℝ) (f : ℝ → ℂ) : Prop :=
  ContDiff ℝ ((⊤ : ℕ∞) : WithTop ℕ∞) f ∧ Function.support f ⊆ Set.Icc (-a) a

/-- Window positivity: autocorrelation support `[−L, L]`. -/
def WeilPos (L : ℝ) : Prop :=
  ∀ f : ℝ → ℂ, TestFn (L / 2) f → 0 ≤ (weilQ f).re

theorem TestFn.mono {a b : ℝ} {f : ℝ → ℂ} (hab : a ≤ b) (hf : TestFn a f) : TestFn b f :=
  ⟨hf.1, hf.2.trans (Set.Icc_subset_Icc (by linarith) hab)⟩

/-- Monotonicity: positivity on a window gives positivity on every smaller window. -/
theorem WeilPos.mono {L L' : ℝ} (h : L' ≤ L) (hL : WeilPos L) : WeilPos L' :=
  fun f hf => hL f (hf.mono (by linarith))

/-- The autocorrelation of `f` supported in `[−a, a]` vanishes outside `[−2a, 2a]`. -/
theorem autocorr_eq_zero {a x : ℝ} {f : ℝ → ℂ} (hf : Function.support f ⊆ Set.Icc (-a) a)
    (hx : 2 * a < |x|) : autocorr f x = 0 := by
  unfold autocorr
  have : (fun y : ℝ => f y * starRingEnd ℂ (f (y - x))) = fun _ => 0 := by
    funext y
    by_cases h1 : f y = 0
    · simp [h1]
    by_cases h2 : f (y - x) = 0
    · simp [h2]
    exfalso
    have hy := hf (Function.mem_support.mpr h1)
    have hyx := hf (Function.mem_support.mpr h2)
    simp only [Set.mem_Icc] at hy hyx
    rcases abs_cases x with ⟨hx', _⟩ | ⟨hx', _⟩ <;> rw [hx'] at hx <;> linarith [hy.1, hy.2, hyx.1, hyx.2]
  rw [this, integral_zero]

/-- Only prime powers `n ≤ e^{2a}` enter the prime term (window `L = 2a`). -/
theorem primeTerm_eq_finite {a : ℝ} {f : ℝ → ℂ} (hf : Function.support f ⊆ Set.Icc (-a) a) :
    primeTerm (autocorr f) =
      ∑ n ∈ Finset.range (⌊Real.exp (2 * a)⌋₊ + 1),
        ((Λ n / Real.sqrt n : ℝ) : ℂ) * (autocorr f (Real.log n) + autocorr f (-Real.log n)) := by
  unfold primeTerm
  apply tsum_eq_sum
  intro n hn
  simp only [Finset.mem_range, not_lt] at hn
  have hn' : Real.exp (2 * a) < n := by
    have := Nat.lt_floor_add_one (Real.exp (2 * a))
    exact this.trans_le (by exact_mod_cast hn)
  have hpos : (0 : ℝ) < n := (Real.exp_pos _).trans hn'
  have hlog : 2 * a < Real.log n := (Real.lt_log_iff_exp_lt hpos).mpr hn'
  have hlog0 : 0 ≤ Real.log n := by
    have : (1 : ℝ) ≤ n := by
      rcases Nat.eq_zero_or_pos n with h0 | h0
      · subst h0; simp at hpos
      · exact_mod_cast h0
    exact Real.log_nonneg this
  rw [autocorr_eq_zero hf (by rwa [abs_of_nonneg hlog0]),
    autocorr_eq_zero hf (by rwa [abs_neg, abs_of_nonneg hlog0])]
  simp

/-- Below `log 2` no prime enters: the functional is pole + archimedean only. -/
theorem primeTerm_eq_zero {a : ℝ} {f : ℝ → ℂ} (hf : Function.support f ⊆ Set.Icc (-a) a)
    (ha : 2 * a < Real.log 2) : primeTerm (autocorr f) = 0 := by
  rw [primeTerm_eq_finite hf]
  have hfl : ⌊Real.exp (2 * a)⌋₊ ≤ 1 := by
    have h2 : Real.exp (2 * a) < 2 := by
      have := Real.exp_lt_exp.mpr ha
      rwa [Real.exp_log (by norm_num)] at this
    have : ⌊Real.exp (2 * a)⌋₊ < 2 := (Nat.floor_lt (Real.exp_pos _).le).mpr (by exact_mod_cast h2)
    omega
  apply Finset.sum_eq_zero
  intro n hn
  simp only [Finset.mem_range] at hn
  have : n = 0 ∨ n = 1 := by omega
  rcases this with rfl | rfl <;> simp

/-! ## Cited results (statements only) -/

/-- Weil's criterion (Weil 1952; Bombieri, Rend. Lincei 2000): RH ⟺ positivity on all windows. -/
def WeilCriterion : Prop := RiemannHypothesis ↔ ∀ L : ℝ, WeilPos L

/-- Yoshida (1992); Connes–Consani, Selecta Math. 2021 (arXiv 2006.13771): window `log 2`. -/
def YoshidaConnesConsani : Prop := WeilPos (Real.log 2)

/-- Zhu, arXiv 2608.24827 (2026), certified computation: `supp f ⊂ [−0.8, 0.8]`, window 1.6,
`Q(f) ≥ 8.9e-18 ‖f‖²`. -/
def ZhuWindow : Prop := WeilPos (16 / 10)

/-- Liu, alphaXiv 2609 "Certified Weil positivity beyond the unit window" (2026 preprint):
half-width 17/16 (window 17/8, prime power 8 included), `Q(f) ≥ 2^{-49162} ‖f‖²`. -/
def LiuWindow : Prop := WeilPos (17 / 8)

/-- Mori, "Weil positivity on the window of width log 15" (2026 preprint; Arb certificate data
Zenodo 10.5281/zenodo.23281175, C² test functions): `Q(f) ≥ 2^{-2904} ‖f‖²`. -/
def MoriWindow : Prop := WeilPos (Real.log 15)

/-- This repo, `Riemann/Deep/Weil/cert/cert.py` (Arb certificate, not Lean-checked):
`Q(f) ≥ 9.99e-22 ‖f‖²` for every complex `f` supported in `[−7/8, 7/8]`. -/
def CertWindow : Prop := WeilPos (7 / 4)

theorem CertWindow.zhu (h : CertWindow) : ZhuWindow := WeilPos.mono (by norm_num) h

theorem LiuWindow.cert (h : LiuWindow) : CertWindow := WeilPos.mono (by norm_num) h

theorem MoriWindow.liu (h : MoriWindow) : LiuWindow := by
  refine WeilPos.mono ?_ h
  rw [Real.le_log_iff_exp_le (by norm_num)]
  have h1 : Real.exp (17 / 8) ≤ Real.exp (5 / 2) := Real.exp_le_exp.mpr (by norm_num)
  have e5 : Real.exp (5 / 2) ^ 2 = Real.exp 1 ^ 5 := by
    rw [← Real.exp_nat_mul, ← Real.exp_nat_mul]; norm_num
  have h4 : Real.exp 1 ^ 5 < (2.7182818286 : ℝ) ^ 5 := by
    have := Real.exp_one_lt_d9
    have := (Real.exp_pos 1).le
    gcongr
  have h5 : (2.7182818286 : ℝ) ^ 5 < 225 := by norm_num
  nlinarith [Real.exp_pos (5 / 2)]

/-- Zhu's window contains the classical one. -/
theorem ZhuWindow.yoshida (h : ZhuWindow) : YoshidaConnesConsani :=
  WeilPos.mono (le_of_lt (Real.log_two_lt_d9.trans (by norm_num))) h

/-- Under Weil's criterion, RH is exactly "every window is positive"; no finite window suffices
alone, since each `WeilPos L` is implied by RH. -/
theorem rh_of_all_windows (hW : WeilCriterion) (h : ∀ L, WeilPos L) : RiemannHypothesis :=
  hW.mpr h

end Riemann.Deep.Weil
