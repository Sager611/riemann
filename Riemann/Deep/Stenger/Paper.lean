import Riemann.Deep.FirstZero.FirstZero
import Riemann.Deep.DH.ZeroFinal

/-!
# Stenger, arXiv 1708.01209v19 (July 2024): chain laid down, failing step pinned and refuted

Paper: `G(z) = ∫₀^∞ y^{z-1}/(1+e^y) dy = (1 - 2^{1-z}) Γ(z) ζ(z)` (eq. 1.2, 2.1).
Chain:
* Lemma 2.6(i): `G`, `ζ` have the same zeros in `D` (correct; `(1-2^{1-z})Γ(z)` is zero-free there).
* Lemma 2.9/2.10: with `κ(σ,x) = e^{σx}/(1+exp eˣ) > 0`,
  `∂_σ^{2m} G(σ) = ∫ x^{2m} κ > 0` on (0,1) and `∂_σ^{2m+1} G(σ) = -∫₀^∞ x^{2m+1} Δ < 0` on (0,1/2].
  **Correct** (mpmath: `.scratch/stenger/num.py`). Only (i) is used in §4: `EvenPos`.
* §4, Lemmas 4.2/4.3 (the step, `StepHalf`): if `G(σ ± iτ)` vanish with multiplicity `ν+1`, "then `G`
  vanishes on `ℓ[ξ, ξ̄]` with total multiplicity `2ν+2`. Equivalently, `G^{(2ν+2)}(σ)` vanishes with
  multiplicity one, which is not allowed by Lemma 2.10." This is a complex Rolle theorem; it is false
  (`e^z - 1` vanishes at `0, 2πi`, its derivative never does), and Lemma 2.10 only speaks of the real point.
* Lemma 4.2 with `ν = -1` is the same sentence for simple zeros; the step kills EVERY non-real zero,
  including the on-line zeros the paper itself cites (Hardy–Littlewood [16]): `chain`.

Refutations (all sorry-free):
* `not_stepHalf_cosh`: `cosh(z - 3/10)` is real on ℝ, satisfies `EvenPos`, and vanishes at `3/10 + iπ/2`.
* `stenger_refuted`: `EvenPos G → ¬ StepHalf G`, via the certified first zeta zero (14.13 < t < 14.14).
* `dh_test`: Γ·DH (Davenport–Heilbronn) has the same positive Mellin kernel (`theta_pos`) and the same
  sign table numerically, yet a certified off-line zero; the step would "prove" RH for it.
-/

namespace Riemann.Deep.Stenger

open Complex

/-- Lemma 2.10(i): `Re ∂_σ^{2m} H(σ) > 0` for `σ ∈ (0,1)`. True for `G` (`∫ x^{2m} κ(σ,x) dx`, `κ > 0`). -/
def EvenPos (H : ℂ → ℂ) : Prop :=
  ∀ σ : ℝ, 0 < σ → σ < 1 → ∀ m : ℕ, 0 < (iteratedDeriv (2 * m) H σ).re

/-- The §4 step on a range `(0, a]`: a zero at `σ + iτ` (its mirror `σ - iτ` comes by reflection)
forces some even derivative `H^{(2m)}(σ)`, `m ≥ 1`, to vanish. -/
def StepOn (a : ℝ) (H : ℂ → ℂ) : Prop :=
  ∀ σ τ : ℝ, 0 < σ → σ ≤ a → τ ≠ 0 → H (σ + τ * I) = 0 →
    ∃ m : ℕ, 1 ≤ m ∧ iteratedDeriv (2 * m) H σ = 0

/-- As used in the paper: Lemma 4.2 (`σ = 1/2`) and Lemma 4.3 (`σ < 1/2`). -/
abbrev StepHalf := StepOn (1 / 2)

/-- The paper's chain: Lemma 2.10 + the step ⇒ no zero of `H` in `D` with `Re ≤ a`, `Im ≠ 0`.
Proves too much: it excludes the on-line zeros too. -/
theorem chain {a : ℝ} {H : ℂ → ℂ} (ha : a < 1) (hpos : EvenPos H) (hstep : StepOn a H) :
    ∀ z : ℂ, 0 < z.re → z.re ≤ a → z.im ≠ 0 → H z ≠ 0 := by
  intro z h0 ha' him hz
  have hz' : H ((z.re : ℂ) + (z.im : ℂ) * I) = 0 := by rw [re_add_im]; exact hz
  obtain ⟨m, _, hm⟩ := hstep z.re z.im h0 ha' him hz'
  have := hpos z.re h0 (lt_of_le_of_lt ha' ha) m
  rw [hm, zero_re] at this
  exact lt_irrefl _ this

/-! ## Explicit counterexample to the step: `cosh (z - 3/10)` -/

/-- Model function: real on ℝ, all even derivatives `cosh (σ - 3/10) ≥ 1`. -/
noncomputable def Hc (z : ℂ) : ℂ := Complex.cosh (z - 3 / 10)

lemma deriv_Hc : deriv Hc = fun z => Complex.sinh (z - 3 / 10) := by
  funext z
  exact ((Complex.hasDerivAt_cosh (z - 3 / 10)).comp_sub_const z (3 / 10)).deriv

lemma deriv_sinh_shift : deriv (fun z : ℂ => Complex.sinh (z - 3 / 10)) = Hc := by
  funext z
  exact ((Complex.hasDerivAt_sinh (z - 3 / 10)).comp_sub_const z (3 / 10)).deriv

lemma iteratedDeriv_Hc (m : ℕ) : iteratedDeriv (2 * m) Hc = Hc := by
  induction m with
  | zero => simp
  | succ n ih =>
    rw [show 2 * (n + 1) = 2 * n + 1 + 1 by ring, iteratedDeriv_succ', iteratedDeriv_succ',
      deriv_Hc, deriv_sinh_shift, ih]

theorem evenPos_Hc : EvenPos Hc := by
  intro σ _ _ m
  rw [iteratedDeriv_Hc, Hc, show (σ : ℂ) - 3 / 10 = ((σ - 3 / 10 : ℝ) : ℂ) by push_cast; ring,
    ← ofReal_cosh, ofReal_re]
  exact Real.cosh_pos _

theorem Hc_zero : Hc (((3 / 10 : ℝ) : ℂ) + ((Real.pi / 2 : ℝ) : ℂ) * I) = 0 := by
  rw [Hc, show ((3 / 10 : ℝ) : ℂ) + ((Real.pi / 2 : ℝ) : ℂ) * I - 3 / 10 = ((Real.pi / 2 : ℝ) : ℂ) * I by
    push_cast; ring, cosh_mul_I, ← ofReal_cos, Real.cos_pi_div_two, ofReal_zero]

/-- The step is false for a function satisfying Lemma 2.10(i), with a simple zero off the line. -/
theorem not_stepHalf_cosh : ¬ StepHalf Hc := by
  intro h
  exact chain (by norm_num) evenPos_Hc h (((3 / 10 : ℝ) : ℂ) + ((Real.pi / 2 : ℝ) : ℂ) * I)
    (by norm_num) (by norm_num) (by norm_num [Real.pi_ne_zero]) Hc_zero

/-! ## The paper's own `G`: the step contradicts the certified first zero -/

/-- `G(z) = (1 - 2^{1-z}) Γ(z) ζ(z)` (eq. 2.1). -/
noncomputable def G (z : ℂ) : ℂ := (1 - (2 : ℂ) ^ (1 - z)) * Complex.Gamma z * riemannZeta z

/-- `G` vanishes at the first zeta zero `1/2 + it`, `14.13 < t < 14.14`. -/
theorem G_first_zero : ∃ t : ℝ, 14.13 < t ∧ G (((1 / 2 : ℝ) : ℂ) + (t : ℂ) * I) = 0 := by
  obtain ⟨t, ht1, _, hz⟩ := Fayed.FirstZero.first_zero
  refine ⟨t, ht1, ?_⟩
  rw [G, show (((1 / 2 : ℝ) : ℂ)) = (1 / 2 : ℂ) by push_cast; ring, hz, mul_zero]

/-- **Main refutation.** Lemma 2.10(i) (true) and the §4 step cannot both hold for `G`. -/
theorem stenger_refuted (hpos : EvenPos G) : ¬ StepHalf G := by
  intro h
  obtain ⟨t, ht, hz⟩ := G_first_zero
  have hne := chain (by norm_num) hpos h (((1 / 2 : ℝ) : ℂ) + (t : ℂ) * I)
    (by simp) (by simp) (by simp; linarith)
  exact hne hz

/-! ## Davenport–Heilbronn test -/

/-- Mellin kernel of `Γ·DH` in `u = e^{-y}`: `θ = (u + κu² - κu³ - u⁴)/(1 - u⁵)` (coefficients
`1, κ, -κ, -1, 0`, κ ≈ 0.284). Positive, like `1/(1+e^y)` for `G`, so Lemmas 2.9/2.10 carry over. -/
theorem theta_pos {u k : ℝ} (hu0 : 0 < u) (hu1 : u < 1) (hk : |k| < 1) :
    0 < (u + k * u ^ 2 - k * u ^ 3 - u ^ 4) / (1 - u ^ 5) := by
  have hk' := abs_lt.mp hk
  have h5 : u ^ 5 < 1 := pow_lt_one₀ hu0.le hu1 (by norm_num)
  apply div_pos _ (by linarith)
  have hfac : u + k * u ^ 2 - k * u ^ 3 - u ^ 4 = u * (1 - u) * (1 + (1 + k) * u + u ^ 2) := by ring
  rw [hfac]
  have : 0 < 1 + (1 + k) * u + u ^ 2 := by nlinarith
  have : 0 < u * (1 - u) := mul_pos hu0 (by linarith)
  positivity

/-- `Γ(z)·DH(z)`, the DH analogue of `G`. -/
noncomputable def GD (z : ℂ) : ℂ := Complex.Gamma z * Fayed.DH.DH z

/-- DH test: the step (on (0, 9/10], where Lemma 2.10(i) for `GD` also holds) is refuted by the
certified DH zero `0.808517 + 85.699348 i`. -/
theorem dh_test (hpos : EvenPos GD) : ¬ StepOn (9 / 10) GD := by
  intro h
  obtain ⟨s, hs, hre, him⟩ := Fayed.DH.Zero.offLineZero_holds
  rw [abs_lt] at hre him
  have hne := chain (by norm_num) hpos h s (by norm_num at hre ⊢; linarith) (by norm_num at hre ⊢; linarith)
    (by norm_num at him ⊢; linarith)
  exact hne (by rw [GD, hs, mul_zero])

end Riemann.Deep.Stenger

#print axioms Riemann.Deep.Stenger.not_stepHalf_cosh
#print axioms Riemann.Deep.Stenger.stenger_refuted
#print axioms Riemann.Deep.Stenger.dh_test
#print axioms Riemann.Deep.Stenger.theta_pos
