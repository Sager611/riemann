module

public import Riemann.Deep.LeanLi.LogCos
public import Mathlib.Analysis.SpecialFunctions.Integrals.Basic

/-!
# LeanLi: the model trend `I(n) = ∫_0^π 2(1 - cos nv) a₀(v) dv = Σ_{k odd < n} (n-k)/k - (n/2) log 4π`

`a₀(v) = log(cot(v/2)/(4π))/(8π sin²(v/2))`: in `t = ½ cot(v/2)`, `I(n) = π⁻¹ ∫_0^∞ (1 - cos nθ(t)) log(t/2π) dt`.
-/

@[expose] public section

open Real Filter Topology Set MeasureTheory
open scoped Interval

namespace LeanLi

/-- `ℓ(v) = log(cot(v/2)/(4π))` (on `(0, π)`). -/
noncomputable def ell (v : ℝ) : ℝ := log (cos (v / 2)) - log (sin (v / 2)) - log (4 * π)

/-- `a₀(v) = ℓ(v)/(8π sin²(v/2))`. -/
noncomputable def a0 (v : ℝ) : ℝ := ell v / (8 * π * sin (v / 2) ^ 2)

/-- The model trend `I(n)`. -/
noncomputable def Imod (n : ℕ) : ℝ := ∫ v in (0 : ℝ)..π, 2 * (1 - cos (n * v)) * a0 v

/-- Its closed form. -/
noncomputable def Iclosed (n : ℕ) : ℝ :=
  ∑ k ∈ (Finset.range n).filter Odd, ((n : ℝ) - k) / k - (n / 2) * log (4 * π)

/-- Fejér: `sin²(v/2)(n + 2Σ_{j=1}^{n-1}(n-j) cos jv) = sin²(nv/2)`. -/
lemma fejer_identity (v : ℝ) (n : ℕ) :
    sin (v / 2) ^ 2 * (n + 2 * ∑ j ∈ Finset.Ico 1 n, ((n : ℝ) - j) * cos (j * v)) =
      sin (n * v / 2) ^ 2 := by
  induction n with
  | zero => simp
  | succ n ih =>
    have hsplit : ∑ j ∈ Finset.Ico 1 (n + 1), (((n + 1 : ℕ) : ℝ) - j) * cos (j * v) =
        ∑ j ∈ Finset.Ico 1 n, ((n : ℝ) - j) * cos (j * v) + ∑ j ∈ Finset.Ico 1 (n + 1), cos (j * v) := by
      rcases Nat.eq_zero_or_pos n with h0 | hpos
      · subst h0; simp
      · rw [Finset.sum_Ico_succ_top (by omega : 1 ≤ n), Finset.sum_Ico_succ_top (by omega : 1 ≤ n),
          ← add_assoc]
        rw [← Finset.sum_add_distrib]
        push_cast
        congr 1
        · refine Finset.sum_congr rfl fun j _ => by ring
        · ring
    rw [hsplit]
    have hD : sin (v / 2) * (1 + 2 * ∑ j ∈ Finset.Ico 1 (n + 1), cos (j * v)) =
        sin ((2 * ((n + 1 : ℕ) : ℝ) - 1) * v / 2) := dirichlet_sin v (n + 1) (by omega)
    have hsq : sin ((n + 1 : ℕ) * v / 2) ^ 2 - sin (n * v / 2) ^ 2 =
        sin (v / 2) * sin ((2 * ((n + 1 : ℕ) : ℝ) - 1) * v / 2) := by
      push_cast
      have e1 : ((n : ℝ) + 1) * v / 2 = (2 * n + 1) * v / 4 + v / 4 := by ring
      have e2 : (n : ℝ) * v / 2 = (2 * n + 1) * v / 4 - v / 4 := by ring
      rw [e1, e2, sin_add, sin_sub, show (2 * ((n : ℝ) + 1) - 1) * v / 2 = 2 * ((2 * n + 1) * v / 4) by ring,
        show v / 2 = 2 * (v / 4) by ring, sin_two_mul, sin_two_mul]
      ring
    have : sin (v / 2) ^ 2 * (((n + 1 : ℕ) : ℝ) + 2 * (∑ j ∈ Finset.Ico 1 n, ((n : ℝ) - j) * cos (j * v) +
        ∑ j ∈ Finset.Ico 1 (n + 1), cos (j * v))) =
        sin (v / 2) ^ 2 * (n + 2 * ∑ j ∈ Finset.Ico 1 n, ((n : ℝ) - j) * cos (j * v)) +
        sin (v / 2) * (sin (v / 2) * (1 + 2 * ∑ j ∈ Finset.Ico 1 (n + 1), cos (j * v))) := by
      push_cast; ring
    rw [this, ih, hD]
    linarith [hsq]

lemma intervalIntegrable_log_sin_half (a b : ℝ) :
    IntervalIntegrable (fun v => log (sin (v / 2))) volume a b := by
  have := (intervalIntegrable_log_sin (a := a / 2) (b := b / 2)).comp_mul_left (c := 1 / 2)
  simpa [div_eq_mul_inv, mul_comm] using this

lemma intervalIntegrable_log_cos_half (a b : ℝ) :
    IntervalIntegrable (fun v => log (cos (v / 2))) volume a b := by
  have := (intervalIntegrable_log_cos (a := a / 2) (b := b / 2)).comp_mul_left (c := 1 / 2)
  simpa [div_eq_mul_inv, mul_comm] using this

lemma intervalIntegrable_ell (a b : ℝ) : IntervalIntegrable ell volume a b :=
  ((intervalIntegrable_log_cos_half a b).sub (intervalIntegrable_log_sin_half a b)).sub
    intervalIntegrable_const

/-- Reflection `v ↦ π - v`. -/
lemma integral_cos_mul_log_cos_half (j : ℕ) :
    ∫ v in (0 : ℝ)..π, cos (j * v) * log (cos (v / 2)) =
      (-1) ^ j * ∫ v in (0 : ℝ)..π, cos (j * v) * log (sin (v / 2)) := by
  have h := intervalIntegral.integral_comp_sub_left
    (fun u => cos (j * u) * log (sin (u / 2))) (a := 0) (b := π) π
  simp only [sub_self, sub_zero] at h
  rw [← h, ← intervalIntegral.integral_const_mul]
  refine intervalIntegral.integral_congr fun v _ => ?_
  have h1 : cos (j * (π - v)) = (-1) ^ j * cos (j * v) := by
    rw [show (j : ℝ) * (π - v) = j * π - j * v by ring, cos_sub,
      show (j : ℝ) * π = (j : ℕ) * π by rfl, cos_nat_mul_pi, sin_nat_mul_pi]
    ring
  have h2 : sin ((π - v) / 2) = cos (v / 2) := by
    rw [show (π - v) / 2 = π / 2 - v / 2 by ring, sin_pi_div_two_sub]
  rw [h1, h2]
  have : ((-1 : ℝ) ^ j) ^ 2 = 1 := by rw [← pow_mul, mul_comm, pow_mul]; simp
  calc cos (j * v) * log (cos (v / 2)) = ((-1) ^ j) ^ 2 * (cos (j * v) * log (cos (v / 2))) := by
        rw [this, one_mul]
    _ = (-1) ^ j * ((-1) ^ j * cos (j * v) * log (cos (v / 2))) := by ring

lemma integral_cos_mul_ell (j : ℕ) (hj : 1 ≤ j) :
    ∫ v in (0 : ℝ)..π, cos (j * v) * ell v = (π / (2 * j)) * (1 - (-1) ^ j) := by
  have hc : ∀ a b : ℝ, IntervalIntegrable (fun v => cos (j * v)) volume a b := fun a b =>
    (Continuous.continuousOn (by fun_prop)).intervalIntegrable
  have h1 : IntervalIntegrable (fun v => cos (j * v) * log (cos (v / 2))) volume 0 π :=
    (intervalIntegrable_log_cos_half 0 π).continuousOn_mul (Continuous.continuousOn (by fun_prop))
  have h2 : IntervalIntegrable (fun v => cos (j * v) * log (sin (v / 2))) volume 0 π :=
    (intervalIntegrable_log_sin_half 0 π).continuousOn_mul (Continuous.continuousOn (by fun_prop))
  have e : (fun v => cos (j * v) * ell v) = fun v => (cos (j * v) * log (cos (v / 2)) -
      cos (j * v) * log (sin (v / 2))) - log (4 * π) * cos (j * v) := by
    funext v; simp only [ell]; ring
  rw [e, intervalIntegral.integral_sub (h1.sub h2) ((hc 0 π).const_mul _),
    intervalIntegral.integral_sub h1 h2, intervalIntegral.integral_const_mul,
    integral_cos_mul_log_cos_half, integral_cos_mul_log_sin_half j hj]
  have hj0 : (j : ℝ) ≠ 0 := by exact_mod_cast (by omega : j ≠ 0)
  have hcos : ∫ v in (0 : ℝ)..π, cos (j * v) = 0 := by
    rw [intervalIntegral.integral_comp_mul_left (fun v => cos v) hj0, integral_cos]
    simp only [mul_zero, sin_zero, sub_zero]
    rw [show (j : ℝ) * π = (j : ℕ) * π by rfl, sin_nat_mul_pi]; simp
  rw [hcos]
  field_simp
  ring

lemma integral_ell : ∫ v in (0 : ℝ)..π, ell v = -π * log (4 * π) := by
  have h := intervalIntegral.integral_comp_sub_left (fun u => log (sin (u / 2))) (a := 0) (b := π) π
  simp only [sub_self, sub_zero] at h
  have hcs : ∫ v in (0 : ℝ)..π, log (cos (v / 2)) = ∫ v in (0 : ℝ)..π, log (sin (v / 2)) := by
    rw [← h]
    refine intervalIntegral.integral_congr fun v _ => ?_
    rw [show (π - v) / 2 = π / 2 - v / 2 by ring, sin_pi_div_two_sub]
  unfold ell
  rw [intervalIntegral.integral_sub ((intervalIntegrable_log_cos_half 0 π).sub
    (intervalIntegrable_log_sin_half 0 π)) intervalIntegrable_const,
    intervalIntegral.integral_sub (intervalIntegrable_log_cos_half 0 π) (intervalIntegrable_log_sin_half 0 π),
    hcs, intervalIntegral.integral_const]
  simp

/-- On `sin(v/2) ≠ 0`: `2(1 - cos nv) a₀(v) = (1/2π)(n + 2Σ(n-j) cos jv) ℓ(v)`. -/
lemma integrand_eq (n : ℕ) {v : ℝ} (hv : sin (v / 2) ≠ 0) :
    2 * (1 - cos (n * v)) * a0 v =
      (1 / (2 * π)) * ((n + 2 * ∑ j ∈ Finset.Ico 1 n, ((n : ℝ) - j) * cos (j * v)) * ell v) := by
  have hf := fejer_identity v n
  have hcos : 1 - cos (n * v) = 2 * sin (n * v / 2) ^ 2 := by
    have h2 := sin_sq_add_cos_sq (n * v / 2)
    have e1 : cos (n * v) = cos (2 * (n * v / 2)) := by congr 1; ring
    rw [e1, cos_two_mul]
    nlinarith [h2]
  unfold a0
  rw [hcos, ← hf]
  have hπ : π ≠ 0 := pi_ne_zero
  field_simp
  ring

theorem Imod_eq_closed (n : ℕ) : Imod n = Iclosed n := by
  unfold Imod
  have hae : ∀ᵐ v ∂volume, v ∈ Ι (0 : ℝ) π → 2 * (1 - cos (n * v)) * a0 v =
      (1 / (2 * π)) * ((n : ℝ) * ell v +
        2 * ∑ j ∈ Finset.Ico 1 n, ((n : ℝ) - j) * (cos (j * v) * ell v)) := by
    refine Eventually.of_forall fun v hv => ?_
    rw [uIoc_of_le pi_pos.le] at hv
    rw [integrand_eq n (sin_pos_of_pos_of_lt_pi (by linarith [hv.1]) (by linarith [hv.2, pi_pos])).ne']
    congr 1
    rw [add_mul, mul_assoc 2, Finset.sum_mul]
    congr 2
    exact Finset.sum_congr rfl fun j _ => by ring
  rw [intervalIntegral.integral_congr_ae hae, intervalIntegral.integral_const_mul]
  have hint2 : ∀ j ∈ Finset.Ico 1 n,
      IntervalIntegrable (fun v => ((n : ℝ) - j) * (cos (j * v) * ell v)) volume 0 π :=
    fun j _ => ((intervalIntegrable_ell 0 π).continuousOn_mul
      (Continuous.continuousOn (by fun_prop))).const_mul _
  have hsumint : IntervalIntegrable
      (fun v => ∑ j ∈ Finset.Ico 1 n, ((n : ℝ) - j) * (cos (j * v) * ell v)) volume 0 π := by
    have := IntervalIntegrable.sum (Finset.Ico 1 n) hint2
    convert this using 1
    funext v
    simp [Finset.sum_apply]
  rw [intervalIntegral.integral_add ((intervalIntegrable_ell 0 π).const_mul _) (hsumint.const_mul 2),
    intervalIntegral.integral_const_mul, intervalIntegral.integral_const_mul,
    intervalIntegral.integral_finsetSum hint2, integral_ell]
  have hsum : ∑ j ∈ Finset.Ico 1 n, ∫ v in (0 : ℝ)..π, ((n : ℝ) - j) * (cos (j * v) * ell v) =
      ∑ j ∈ Finset.Ico 1 n, ((n : ℝ) - j) * ((π / (2 * j)) * (1 - (-1) ^ j)) := by
    refine Finset.sum_congr rfl fun j hj => ?_
    rw [intervalIntegral.integral_const_mul, integral_cos_mul_ell j (Finset.mem_Ico.1 hj).1]
  rw [hsum]
  have hodd : ∑ k ∈ (Finset.range n).filter Odd, ((n : ℝ) - k) / k =
      ∑ j ∈ Finset.Ico 1 n, ((n : ℝ) - j) * ((1 - (-1) ^ j) / (2 * j)) := by
    rw [Finset.sum_filter, Finset.range_eq_Ico]
    rcases Nat.eq_zero_or_pos n with h0 | hpos
    · subst h0; simp
    rw [Finset.sum_eq_sum_Ico_succ_bot hpos]
    simp only [Nat.not_odd_zero, if_false, zero_add]
    refine Finset.sum_congr rfl fun j hj => ?_
    have hj0 : (j : ℝ) ≠ 0 := by
      have := (Finset.mem_Ico.1 hj).1; exact_mod_cast (by omega : j ≠ 0)
    rcases Nat.even_or_odd j with he | ho
    · rw [if_neg (Nat.not_odd_iff_even.2 he), he.neg_one_pow]; ring
    · rw [if_pos ho, ho.neg_one_pow]; field_simp; ring
  unfold Iclosed
  rw [hodd]
  have h2 : 2 * ∑ j ∈ Finset.Ico 1 n, ((n : ℝ) - j) * ((π / (2 * j)) * (1 - (-1) ^ j)) =
      2 * π * ∑ j ∈ Finset.Ico 1 n, ((n : ℝ) - j) * ((1 - (-1) ^ j) / (2 * j)) := by
    rw [Finset.mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun j hj => ?_
    have hj0 : (j : ℝ) ≠ 0 := by
      have := (Finset.mem_Ico.1 hj).1; exact_mod_cast (by omega : j ≠ 0)
    field_simp
  rw [h2]
  have hπ : π ≠ 0 := pi_ne_zero
  field_simp
  ring

end LeanLi
