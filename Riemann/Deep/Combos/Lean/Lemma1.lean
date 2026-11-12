import Riemann.Deep.Combos.Lean.LatticeC

/-! # Lemma 1 (primitive lattice points), NOTE l.67–75: the statement -/

open Finset

namespace RingClass

/-- **Lemma 1** (NOTE l.67–75). For 0 < a, |b| ≤ a ≤ c (so D = b² − 4ac < 0) and N ≥ 1:
P_Q(N) ≤ 2 + 8N|D|^{−1/2} + 4N^{1/2}|D|^{−1/4} (with square roots). -/
theorem lemma1 (a b c : ℤ) (ha : 0 < a) (hb : |b| ≤ a) (hac : a ≤ c) (N : ℝ) (hN : 1 ≤ N) :
    (PQ a b c N : ℝ) ≤ 2 + 8 * N / Real.sqrt |((b : ℝ) ^ 2 - 4 * a * c)|
      + 4 * Real.sqrt N / Real.sqrt (Real.sqrt |((b : ℝ) ^ 2 - 4 * a * c)|) := by
  have haR : (0 : ℝ) < a := by exact_mod_cast ha
  have hΔ3 := disc_facts (a : ℝ) b c haR (by exact_mod_cast hb) (by exact_mod_cast hac)
  have h3a : (0 : ℝ) < 3 * (a : ℝ) ^ 2 := by positivity
  have habs : |((b : ℝ) ^ 2 - 4 * a * c)| = 4 * (a : ℝ) * c - (b : ℝ) ^ 2 := by
    rw [abs_of_neg (by linarith)]; ring
  rw [habs]
  set Δ : ℝ := 4 * (a : ℝ) * c - (b : ℝ) ^ 2 with hΔ
  have hΔpos : 0 < Δ := by linarith
  have hN0 : (0 : ℝ) ≤ N := by linarith
  set s : ℝ := Real.sqrt (Real.sqrt Δ) with hsdef
  have hsqD : 0 < Real.sqrt Δ := Real.sqrt_pos.mpr hΔpos
  have hs0 : 0 < s := Real.sqrt_pos.mpr hsqD
  have hs2 : s ^ 2 = Real.sqrt Δ := Real.sq_sqrt hsqD.le
  have hs4 : s ^ 4 = Δ := by
    rw [show s ^ 4 = (s ^ 2) ^ 2 by ring, hs2, Real.sq_sqrt hΔpos.le]
  set α : ℝ := Real.sqrt a with hαdef
  have hα0 : 0 < α := Real.sqrt_pos.mpr haR
  have hα2 : α ^ 2 = a := Real.sq_sqrt haR.le
  set r : ℝ := Real.sqrt N with hrdef
  have hr0 : 0 ≤ r := Real.sqrt_nonneg N
  have hr2 : r ^ 2 = N := Real.sq_sqrt hN0
  have hp : (a : ℝ) * N ≤ (α * r) ^ 2 := by rw [mul_pow, hα2, hr2]
  have hq : (s ^ 2) ^ 2 ≤ 4 * (a : ℝ) * c - b ^ 2 := by
    rw [show (s ^ 2) ^ 2 = s ^ 4 by ring, hs4]
  have h := (PQset_subset_card a b c ha hb hac N (α * r) (s ^ 2) (mul_nonneg hα0.le hr0)
    (by positivity) hp hq).2
  have has : (a : ℝ) ^ 2 ≤ s ^ 4 := by
    rw [hs4]; have := sq_nonneg (a : ℝ); linarith
  have h2 := lemma1_arith (a : ℝ) N α r s haR hα0 hr0 hs0 hα2 hr2 has
  unfold PQ
  rw [← hs2]
  linarith

end RingClass

