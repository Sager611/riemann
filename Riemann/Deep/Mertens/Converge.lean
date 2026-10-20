import Riemann.Deep.Mertens.ZeroFree

/-!
# Companion: a Mertens bound makes Σ μ(n) n^{-s} converge (conditionally) for Re s > θ

For 0 ≤ θ, MertensBound θ and θ < Re s, the partial sums Σ_{n ≤ N} μ(n) n^{-s} tend to
g(s) = s ∫_1^∞ M(t) t^{-s-1} dt (Abel summation; adapted from Mathlib's
LSeries_eq_mul_integral, with the absolute-summability input removed).
-/

open Complex Asymptotics Filter Topology MeasureTheory Set Finset

namespace Fayed.Mertens

theorem tendsto_partialSum {θ : ℝ} (h : Fayed.R2.MertensBound θ) (hr : 0 ≤ θ) {s : ℂ}
    (hs : θ < s.re) :
    Tendsto (fun n : ℕ => ∑ k ∈ Icc 0 n, LSeries.term mu s k) atTop (𝓝 (g s)) := by
  have hf : mu 0 = 0 := by simp [mu]
  have hO : (fun n ↦ ∑ k ∈ Icc 1 n, mu k) =O[atTop] fun n ↦ (n : ℝ) ^ θ := isBigO_Msum h
  have h₁ : (-s - 1).re + θ < -1 := by
    rwa [sub_re, one_re, neg_re, neg_sub_left, neg_add_lt_iff_lt_add, add_neg_cancel_comm]
  have h₂ : s ≠ 0 := ne_zero_of_re_pos (hr.trans_lt hs)
  have h₃ (t : ℝ) (ht : t ∈ Set.Ici 1) : DifferentiableAt ℝ (fun x : ℝ ↦ (x : ℂ) ^ (-s)) t :=
    differentiableAt_id.ofReal_cpow_const (zero_lt_one.trans_le ht).ne' (neg_ne_zero.mpr h₂)
  have h₄ : ∀ n, ∑ k ∈ Icc 0 n, mu k = ∑ k ∈ Icc 1 n, mu k := fun n ↦ by
    rw [← Finset.insert_Icc_add_one_left_eq_Icc n.zero_le, sum_insert (by aesop), hf, zero_add,
      zero_add]
  simp_rw [← h₄] at hO
  rw [g, ← integral_eq_mellin, ← integral_const_mul]
  simp_rw [LSeries.term_def₀ hf, mul_comm (mu _)]
  convert!
    tendsto_sum_mul_atTop_nhds_one_sub_integral₀ (f := fun x ↦ (x : ℂ) ^ (-s)) (l := 0) ?_ hf h₃ ?_
      ?_ ?_ (integrableAtFilter_rpow_atTop_iff.mpr h₁)
  · rw [zero_sub, ← integral_neg]
    refine setIntegral_congr_fun measurableSet_Ioi fun t ht ↦ ?_
    rw [deriv_ofReal_cpow_const (zero_lt_one.trans ht).ne', h₄]
    · ring_nf
    · exact neg_ne_zero.mpr <| ne_zero_of_re_pos (hr.trans_lt hs)
  · refine (Iff.mpr integrableOn_Ici_iff_integrableOn_Ioi <|
      integrableOn_Ioi_deriv_ofReal_cpow zero_lt_one
        (by simpa using! hr.trans_lt hs)).locallyIntegrableOn
  · have hlim : Tendsto (fun n : ℕ ↦ (n : ℝ) ^ (-(s.re - θ))) atTop (𝓝 0) :=
      (tendsto_rpow_neg_atTop (by rwa [sub_pos])).comp tendsto_natCast_atTop_atTop
    refine (IsBigO.mul_atTop_rpow_natCast_of_isBigO_rpow (-s.re) _ _ ?_ hO ?_).trans_tendsto hlim
    · exact isBigO_norm_left.mp <| (norm_ofReal_cpow_eventually_eq_atTop _).isBigO.natCast_atTop
    · linarith
  · refine .mul_atTop_rpow_of_isBigO_rpow (-(s + 1).re) θ _ ?_ ?_ (by rw [← neg_re, neg_add'])
    · simpa [-neg_add_rev, neg_add'] using! isBigO_deriv_ofReal_cpow_const_atTop _
    · exact (hO.comp_tendsto tendsto_nat_floor_atTop).trans <|
        isEquivalent_nat_floor.isBigO.rpow hr (eventually_ge_atTop 0)

/-- Companion statement: under a Mertens bound with 0 ≤ θ, the Dirichlet series of μ converges
(as N → ∞ partial sums) on Re s > θ, to the holomorphic function g. -/
theorem mertens_LSeries_converges {θ : ℝ} (h : Fayed.R2.MertensBound θ) (hr : 0 ≤ θ) {s : ℂ}
    (hs : θ < s.re) :
    ∃ L : ℂ, Tendsto (fun n : ℕ => ∑ k ∈ Icc 1 n, (ArithmeticFunction.moebius k : ℂ) /
      (k : ℂ) ^ s) atTop (𝓝 L) := by
  refine ⟨g s, (tendsto_partialSum h hr hs).congr fun n => ?_⟩
  rw [← Finset.insert_Icc_add_one_left_eq_Icc n.zero_le, sum_insert (by aesop), LSeries.term_zero,
    zero_add]
  refine sum_congr rfl fun k hk => ?_
  rw [LSeries.term_of_ne_zero (by have := (mem_Icc.mp hk).1; omega), mu]

end Fayed.Mertens

#print axioms Fayed.Mertens.mertens_LSeries_converges
