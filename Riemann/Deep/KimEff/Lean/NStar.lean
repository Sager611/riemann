import Riemann.Deep.KimEff.Lean.Arith

/-!
# Explicit Kim–Lee wedge, part 8: n_* = ceil Z(12) (certify.py, u_* = 12)

998616377770 < Z(12) ≤ 998616377771, from π to 20 digits (Mathlib) and e to 18 digits
(Taylor bound Real.exp_bound with 20 terms); hence for u > 0 with Z(u) = n:
u ≥ 12 iff n ≥ n_*. So the hypothesis range n ≥ n_* of (i) is exactly u(n) ≥ u_* = 12.
-/

namespace KimEff

noncomputable section

open Real

lemma exp_one_bounds : (2.718281828459045234 : ℝ) ≤ exp 1 ∧ exp 1 ≤ 2.718281828459045236 := by
  have h := Real.exp_bound (x := 1) (by norm_num) (n := 20) (by norm_num)
  simp only [Finset.sum_range_succ, Finset.sum_range_zero, Nat.factorial] at h
  norm_num at h
  obtain ⟨h1, h2⟩ := abs_le.mp h
  constructor <;> linarith

/-- 998616377770 < Z(12) ≤ 998616377771, i.e. n_* = ceil Z(12). -/
theorem Zfun_12_bounds : (998616377770 : ℝ) < Zfun 12 ∧ Zfun 12 ≤ 998616377771 := by
  obtain ⟨hlo, hhi⟩ := exp_one_bounds
  have he : exp (2 * 12) = exp 1 ^ 24 := by rw [← exp_nat_mul]; norm_num
  have hlo24 : (2.718281828459045234 : ℝ) ^ 24 ≤ exp 1 ^ 24 :=
    pow_le_pow_left₀ (by norm_num) hlo 24
  have hhi24 : exp 1 ^ 24 ≤ (2.718281828459045236 : ℝ) ^ 24 :=
    pow_le_pow_left₀ (exp_pos 1).le hhi 24
  have hp1 := pi_gt_d20
  have hp2 := pi_lt_d20
  unfold Zfun
  rw [he]
  constructor
  · have hm : (3.14159265358979323846 : ℝ) * 2.718281828459045234 ^ 24 ≤ π * exp 1 ^ 24 :=
      mul_le_mul hp1.le hlo24 (by norm_num) pi_pos.le
    have hnum : (998616377770 : ℝ)
        < 12 * (3.14159265358979323846 * 2.718281828459045234 ^ 24 - 9 / 4) := by norm_num
    linarith
  · have hm : π * exp 1 ^ 24 ≤ (3.14159265358979323847 : ℝ) * 2.718281828459045236 ^ 24 :=
      mul_le_mul hp2.le hhi24 (pow_nonneg (exp_pos 1).le 24) (by norm_num)
    have hnum : 12 * ((3.14159265358979323847 : ℝ) * 2.718281828459045236 ^ 24 - 9 / 4)
        ≤ 998616377771 := by norm_num
    linarith

/-- For u > 0 with Z(u) = n: u ≥ u_* = 12 iff n ≥ n_*. -/
theorem twelve_le_u_iff {n : ℕ} {u : ℝ} (hu : 0 < u) (hZ : Zfun u = n) :
    12 ≤ u ↔ nStar ≤ n := by
  obtain ⟨hlo, hhi⟩ := Zfun_12_bounds
  constructor
  · intro h12
    have hm : Zfun 12 ≤ Zfun u :=
      Zfun_strictMonoOn.monotoneOn (by norm_num : (12 : ℝ) ∈ Set.Ici 0)
        (Set.mem_Ici.mpr hu.le) h12
    have h1 : (998616377770 : ℝ) < n := by rw [← hZ]; linarith
    have h2 : 998616377770 < n := by exact_mod_cast h1
    show 998616377771 ≤ n
    omega
  · intro hn
    by_contra h
    have h' := not_le.mp h
    have hm := Zfun_strictMonoOn (Set.mem_Ici.mpr hu.le)
      (by norm_num : (12 : ℝ) ∈ Set.Ici 0) h'
    have hn1 : (998616377771 : ℕ) ≤ n := hn
    have hn2 : (998616377771 : ℝ) ≤ n := by exact_mod_cast hn1
    rw [hZ] at hm
    linarith

end

end KimEff

