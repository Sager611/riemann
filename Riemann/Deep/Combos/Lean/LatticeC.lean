import Riemann.Deep.Combos.Lean.LatticeB

/-! # Lemma 1, the final arithmetic and statement (NOTE l.67–75) -/

open Finset

namespace RingClass

lemma PQset_finite (a b c : ℤ) (ha : 0 < a) (hb : |b| ≤ a) (hac : a ≤ c) (N : ℝ) (hN : 0 ≤ N) :
    (PQset a b c N).Finite := by
  have haR : (0 : ℝ) < a := by exact_mod_cast ha
  have ha1 : (1 : ℝ) ≤ a := by exact_mod_cast ha
  have hΔ3 := disc_facts (a : ℝ) b c haR (by exact_mod_cast hb) (by exact_mod_cast hac)
  have haN : 0 ≤ (a : ℝ) * N := mul_nonneg haR.le hN
  refine (PQset_subset_card a b c ha hb hac N (a * N + 1) 1 (by positivity) one_pos ?_ ?_).1
  · have e : ((a : ℝ) * N + 1) ^ 2 = ((a : ℝ) * N) ^ 2 + 2 * ((a : ℝ) * N) + 1 := by ring
    have := sq_nonneg ((a : ℝ) * N)
    linarith
  · have h1 : (1 : ℝ) ≤ (a : ℝ) ^ 2 := one_le_pow₀ ha1
    have : ((1 : ℝ)) ^ 2 = 1 := by norm_num
    linarith

/-- 3a² ≤ s⁴ gives α ≤ s when α² = a. -/
lemma root_le (a α s : ℝ) (_hα : 0 < α) (hs : 0 < s) (hα2 : α ^ 2 = a) (has : a ^ 2 ≤ s ^ 4) :
    α ≤ s := by
  by_contra h0
  have h : s < α := lt_of_not_ge h0
  have h1 : s * s < α * α := mul_self_lt_mul_self hs.le h
  have h2 : (s * s) * (s * s) < (α * α) * (α * α) := mul_self_lt_mul_self (by positivity) h1
  have h3 : (α * α) * (α * α) = a ^ 2 := by rw [← hα2]; ring
  have h4 : (s * s) * (s * s) = s ^ 4 := by ring
  linarith

/-- arithmetic of the last line of NOTE l.74: 2V·L = 8N/s², 2V ≤ 4r/s. -/
lemma lemma1_arith (a N α r s : ℝ) (_ha : 0 < a) (hα : 0 < α) (hr : 0 ≤ r) (hs : 0 < s)
    (hα2 : α ^ 2 = a) (hr2 : r ^ 2 = N) (has : a ^ 2 ≤ s ^ 4) :
    2 + 2 * (2 * (α * r) / s ^ 2) * (2 * (α * r) / a + 1) ≤ 2 + 8 * N / s ^ 2 + 4 * r / s := by
  have hs2 : 0 < s ^ 2 := by positivity
  have hmain : 2 * (2 * (α * r) / s ^ 2) * (2 * (α * r) / a) = 8 * N / s ^ 2 := by
    rw [← hr2, ← hα2]; field_simp; ring
  have hαs := root_le a α s hα hs hα2 has
  have hside : 2 * (2 * (α * r) / s ^ 2) ≤ 4 * r / s := by
    have e : 2 * (2 * (α * r) / s ^ 2) = 4 * r * α / s ^ 2 := by ring
    rw [e, div_le_div_iff₀ hs2 hs]
    have := mul_le_mul_of_nonneg_left hαs (by positivity : (0 : ℝ) ≤ 4 * r * s)
    linarith
  linarith

end RingClass
