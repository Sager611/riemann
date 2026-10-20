import Mathlib

/-!
# Redheffer's theorem: `det R_n = M(n)`

The Redheffer matrix `R_n` has `R i j = 1` iff `i ∣ j` or `j = 1` (1-indexed).
Its determinant is the Mertens function `M(n) = ∑_{k ≤ n} μ(k)` (for `n ≥ 1`).
RH ⟺ `|det R_n| = O(n^{1/2+ε})`.

Proof: `R_n` is the divisibility matrix `D` with column 0 replaced by the all-ones vector.
`det D = 1` (upper unitriangular), and `D⁻¹ = E` with `E i j = [i ∣ j] μ(j/i)`
(Möbius inversion: `∑_{d ∣ c} μ d = [c = 1]`). By Cramer's rule,
`det R = (cramer D 1) 0 = (E *ᵥ 1) 0 = ∑_j μ(j)`.

Note: for `n = 0` the determinant of the empty matrix is `1` while `M(0) = 0`,
so the identity needs `0 < n` (`det_redheffer`); `det_redheffer_zero` records `n = 0`.
-/

open Matrix Finset ArithmeticFunction
open scoped ArithmeticFunction.Moebius

namespace Fayed.Redheffer

/-- The Redheffer matrix (0-indexed). -/
def redheffer (n : ℕ) : Matrix (Fin n) (Fin n) ℤ :=
  fun i j => if (i.1 + 1) ∣ (j.1 + 1) ∨ j.1 = 0 then 1 else 0

/-- The divisibility matrix `D i j = [i+1 ∣ j+1]`. -/
def divMat (n : ℕ) : Matrix (Fin n) (Fin n) ℤ :=
  fun i j => if (i.1 + 1) ∣ (j.1 + 1) then 1 else 0

/-- The Möbius matrix `E i j = [i+1 ∣ j+1] μ((j+1)/(i+1))`. -/
def moebMat (n : ℕ) : Matrix (Fin n) (Fin n) ℤ :=
  fun i j => if (i.1 + 1) ∣ (j.1 + 1) then μ ((j.1 + 1) / (i.1 + 1)) else 0

theorem det_divMat (n : ℕ) : (divMat n).det = 1 := by
  rw [det_of_isUpperTriangular]
  · simp [divMat]
  · intro i j hij
    simp only [divMat]
    rw [ite_eq_right]
    intro h
    have := Nat.le_of_dvd (Nat.succ_pos _) h
    have hij' : j.1 < i.1 := hij
    omega

/-- `∑_{d ∣ c} μ d = [c = 1]`. -/
theorem sum_divisors_moebius (c : ℕ) : ∑ d ∈ c.divisors, μ d = if c = 1 then 1 else 0 := by
  have h := coe_mul_zeta_apply (f := (μ : ArithmeticFunction ℤ)) (x := c)
  rw [moebius_mul_coe_zeta, ArithmeticFunction.one_apply] at h
  exact h.symm

/-- Key arithmetic identity: `∑_{m ∣ b, a ∣ m} μ(m/a) = [a = b]` for `a, b ≥ 1`. -/
theorem sum_divisors_dvd_moebius {a b : ℕ} (ha : 0 < a) (hb : 0 < b) :
    ∑ m ∈ b.divisors, (if a ∣ m then μ (m / a) else 0) = if a = b then 1 else 0 := by
  by_cases hab : a ∣ b
  · obtain ⟨c, rfl⟩ := hab
    have hc : 0 < c := Nat.pos_of_mul_pos_left hb
    rw [← Finset.sum_filter]
    have hset : (a * c).divisors.filter (fun m => a ∣ m) = c.divisors.image (fun d => a * d) := by
      ext m
      simp only [Finset.mem_filter, Nat.mem_divisors, Finset.mem_image]
      constructor
      · rintro ⟨⟨hm, -⟩, ⟨d, rfl⟩⟩
        exact ⟨d, ⟨Nat.dvd_of_mul_dvd_mul_left ha hm, hc.ne'⟩, rfl⟩
      · rintro ⟨d, ⟨hd, -⟩, rfl⟩
        exact ⟨⟨mul_dvd_mul_left a hd, (Nat.mul_pos ha hc).ne'⟩, dvd_mul_right a d⟩
    rw [hset, Finset.sum_image (fun x _ y _ h => Nat.eq_of_mul_eq_mul_left ha h)]
    simp_rw [Nat.mul_div_cancel_left _ ha]
    rw [sum_divisors_moebius]
    congr 1
    apply propext
    constructor
    · rintro rfl; simp
    · intro h
      exact (Nat.mul_eq_left ha.ne').mp h.symm
  · rw [ite_eq_right (show ¬ a = b from fun h => hab (h ▸ dvd_rfl))]
    apply Finset.sum_eq_zero
    intro m hm
    rw [ite_eq_right]
    intro ham
    exact hab (ham.trans (Nat.dvd_of_mem_divisors hm))

theorem moebMat_mul_divMat (n : ℕ) : moebMat n * divMat n = 1 := by
  ext i j
  rw [Matrix.mul_apply, Matrix.one_apply]
  set a := i.1 + 1 with ha_def
  set b := j.1 + 1 with hb_def
  have ha : 0 < a := Nat.succ_pos _
  have hb : 0 < b := Nat.succ_pos _
  have hbn : b ≤ n := j.2
  let g : ℕ → ℤ := fun m => (if a ∣ m then μ (m / a) else 0) * (if m ∣ b then 1 else 0)
  have hsum : ∑ k : Fin n, moebMat n i k * divMat n k j = ∑ k ∈ Finset.range n, g (k + 1) := by
    rw [← Fin.sum_univ_eq_sum_range (fun k => g (k + 1)) n]
    rfl
  rw [hsum]
  have hg0 : g 0 = 0 := by
    simp [g]
  have h1 : ∑ k ∈ Finset.range n, g (k + 1) = ∑ m ∈ Finset.range (n + 1), g m := by
    rw [Finset.sum_range_succ' g n, hg0, add_zero]
  have h2 : ∑ m ∈ Finset.range (n + 1), g m = ∑ m ∈ b.divisors, g m := by
    symm
    apply Finset.sum_subset
    · intro m hm
      have := Nat.divisor_le hm
      rw [Finset.mem_range]; omega
    · intro m _ hm
      have hmb : ¬ m ∣ b := fun hmb => hm (Nat.mem_divisors.mpr ⟨hmb, hb.ne'⟩)
      simp [g, hmb]
  have h3 : ∑ m ∈ b.divisors, g m = ∑ m ∈ b.divisors, (if a ∣ m then μ (m / a) else 0) := by
    apply Finset.sum_congr rfl
    intro m hm
    simp only [g]
    rw [ite_eq_left (Nat.dvd_of_mem_divisors hm), mul_one]
  rw [h1, h2, h3, sum_divisors_dvd_moebius ha hb]
  congr 1
  apply propext
  rw [Fin.ext_iff]
  omega

/-- Mertens sum reindexing. -/
theorem sum_range_succ_eq_sum_Icc (f : ℕ → ℤ) (n : ℕ) :
    ∑ k ∈ Finset.range n, f (k + 1) = ∑ k ∈ Finset.Icc 1 n, f k := by
  induction n with
  | zero => simp
  | succ n ih =>
    rw [Finset.sum_range_succ, ih, Finset.sum_Icc_succ_top (by omega)]

theorem redheffer_eq_updateCol (n : ℕ) (hn : 0 < n) :
    redheffer n = (divMat n).updateCol ⟨0, hn⟩ 1 := by
  ext i j
  rw [updateCol_apply]
  by_cases hj : j = ⟨0, hn⟩
  · subst hj
    simp [redheffer]
  · rw [ite_eq_right hj]
    have hj' : j.1 ≠ 0 := fun h => hj (Fin.ext h)
    simp [redheffer, divMat, hj']

/-- **Redheffer's theorem** (1977): `det R_n = M(n) = ∑_{k=1}^n μ(k)` for `n ≥ 1`. -/
theorem det_redheffer (n : ℕ) (hn : 0 < n) :
    (redheffer n).det = ∑ k ∈ Finset.Icc 1 n, (ArithmeticFunction.moebius k : ℤ) := by
  rw [redheffer_eq_updateCol n hn, ← cramer_apply]
  have hc : cramer (divMat n) 1 = moebMat n *ᵥ 1 := by
    calc cramer (divMat n) 1 = (moebMat n * divMat n) *ᵥ cramer (divMat n) 1 := by
          rw [moebMat_mul_divMat, one_mulVec]
      _ = moebMat n *ᵥ (divMat n *ᵥ cramer (divMat n) 1) := by rw [mulVec_mulVec]
      _ = moebMat n *ᵥ 1 := by rw [mulVec_cramer, det_divMat, one_smul]
  rw [hc, ← sum_range_succ_eq_sum_Icc (fun k => (μ k : ℤ)) n,
    ← Fin.sum_univ_eq_sum_range (fun k => (μ (k + 1) : ℤ)) n]
  simp [mulVec, dotProduct, moebMat]

/-- Degenerate case: the empty matrix has determinant `1`, while `M(0) = 0`. -/
theorem det_redheffer_zero : (redheffer 0).det = 1 := det_isEmpty

end Fayed.Redheffer

#print axioms Fayed.Redheffer.det_redheffer
