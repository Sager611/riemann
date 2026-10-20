import Riemann.Deep.Beta.Reality
import Riemann.Deep.Beta.GammaRatio
import Riemann.Deep.FirstZero.NextLib

/-!
# L(1/2 + it, χ₋₄) through order-8 Euler–Maclaurin for ζ(s, 1/4), ζ(s, 3/4) with J = 3

4^{-s}(M(1/4) - M(3/4)) = 1 - 3^{-s} + 5^{-s} - 7^{-s} + 9^{-s} - 11^{-s} + 13^{-s} K(13/4) - 15^{-s} K(15/4),
since 4^{-s}(j + a)^{-s} = (4j + 4a)^{-s}; only p^{-s} for p = 3, 5, 7, 11, 13 need enclosures.
-/

open Complex Set HurwitzZeta

namespace Fayed.Beta

open Fayed.FirstZero

/-- The coefficient of c^{-s} in the order-8 Euler–Maclaurin main term at c = J + a. -/
noncomputable def KB (c : ℝ) (s : ℂ) : ℂ := -(c : ℂ) / (1 - s) + 1 / 2 +
  s * (1 / 12 * gc s 0 / (c : ℂ) - 1 / 720 * gc s 2 / ((c ^ 3 : ℝ) : ℂ)
    + 1 / 30240 * gc s 4 / ((c ^ 5 : ℝ) : ℂ) - 1 / 1209600 * gc s 6 / ((c ^ 7 : ℝ) : ℂ))

lemma cpow_quarter {x : ℝ} (hx : 0 ≤ x) {n : ℕ} (hn : 4 * x = n) (s : ℂ) :
    ((x : ℝ) : ℂ) ^ (-s) = nps s n / (4 : ℂ) ^ (-s) := by
  have h4 : (4 : ℂ) ^ (-s) ≠ 0 := by rw [Ne, cpow_eq_zero_iff]; norm_num
  rw [eq_div_iff h4, nps, mul_comm, show (4 : ℂ) = ((4 : ℝ) : ℂ) by norm_num,
    ← mul_cpow_ofReal_nonneg (by norm_num) hx, ← ofReal_mul, hn, ofReal_natCast]

lemma cpow_shift_real {x : ℝ} (hx : 0 < x) (s : ℂ) (j : ℕ) :
    ((x : ℝ) : ℂ) ^ (-s - 1 - (j : ℂ)) = (x : ℂ) ^ (-s) * (((x : ℂ)) ^ (j + 1))⁻¹ := by
  have hx' : (x : ℂ) ≠ 0 := by exact_mod_cast hx.ne'
  rw [show -s - 1 - (j : ℂ) = -s + (-((j + 1 : ℕ) : ℂ)) by push_cast; ring,
    cpow_add _ _ hx', cpow_neg (x : ℂ) ((j + 1 : ℕ) : ℂ), cpow_natCast]

/-- The order-8 Euler–Maclaurin main term for ζ(s, a) with J = 3. -/
noncomputable def EMa (a : ℝ) (s : ℂ) : ℂ :=
  ∑ j ∈ Finset.range 3, (((j : ℝ) + a : ℝ) : ℂ) ^ (-s)
    + (-((((3 : ℕ) : ℝ) + a : ℝ) : ℂ) ^ (1 - s) / (1 - s) + ((((3 : ℕ) : ℝ) + a : ℝ) : ℂ) ^ (-s) / 2
      + s * EMPhi s (((3 : ℕ) : ℝ) + a))

lemma EMa_eq {a : ℝ} (ha : 0 < a) (s : ℂ) :
    EMa a s = ((a : ℝ) : ℂ) ^ (-s) + ((1 + a : ℝ) : ℂ) ^ (-s) + ((2 + a : ℝ) : ℂ) ^ (-s)
      + ((3 + a : ℝ) : ℂ) ^ (-s) * KB (3 + a) s := by
  have e3 : (((3 : ℕ) : ℝ) + a : ℝ) = 3 + a := by norm_num
  have hc : (0 : ℝ) < 3 + a := by linarith
  simp only [EMa, KB, e3]
  generalize (3 : ℝ) + a = c at hc ⊢
  have hg : ∀ j : ℕ, gf s j c = gc s j * (c : ℂ) ^ (-s) * ((c : ℂ) ^ (j + 1))⁻¹ := fun j => by
    rw [gf, cpow_shift_real hc]; ring
  have h1 : ((c : ℝ) : ℂ) ^ (1 - s) = (c : ℂ) * (c : ℂ) ^ (-s) := by
    rw [ofReal_cpow_eq_mul hc, show (1 : ℂ) - s - 1 = -s by ring]
  simp only [EMPhi, hg, h1, Finset.sum_range_succ, Finset.sum_range_zero]
  push_cast
  ring_nf

/-- The combined main term: 4^{-s}(M(1/4) - M(3/4)) in terms of p^{-s}. -/
lemma beta_main_eq (s : ℂ) :
    (4 : ℂ) ^ (-s) * (EMa (1 / 4) s - EMa (3 / 4) s) =
      1 - nps s 3 + nps s 5 - nps s 7 + nps s 9 - nps s 11 + nps s 13 * KB (13 / 4) s
        - nps s 15 * KB (15 / 4) s := by
  have h4 : (4 : ℂ) ^ (-s) ≠ 0 := by rw [Ne, cpow_eq_zero_iff]; norm_num
  rw [EMa_eq (by norm_num), EMa_eq (by norm_num)]
  rw [cpow_quarter (x := 1 / 4) (n := 1) (by norm_num) (by norm_num),
    cpow_quarter (x := 1 + 1 / 4) (n := 5) (by norm_num) (by norm_num),
    cpow_quarter (x := 2 + 1 / 4) (n := 9) (by norm_num) (by norm_num),
    cpow_quarter (x := 3 + 1 / 4) (n := 13) (by norm_num) (by norm_num),
    cpow_quarter (x := 3 / 4) (n := 3) (by norm_num) (by norm_num),
    cpow_quarter (x := 1 + 3 / 4) (n := 7) (by norm_num) (by norm_num),
    cpow_quarter (x := 2 + 3 / 4) (n := 11) (by norm_num) (by norm_num),
    cpow_quarter (x := 3 + 3 / 4) (n := 15) (by norm_num) (by norm_num)]
  have hn1 : nps s 1 = 1 := by simp [nps]
  rw [hn1, show (3 : ℝ) + 1 / 4 = 13 / 4 by norm_num, show (3 : ℝ) + 3 / 4 = 15 / 4 by norm_num]
  field_simp
  ring

lemma c_rpow_le {c cN : ℝ} (hc : 0 < c) (hcN : 0 < cN) (h : cN ^ 2 ≤ c) :
    c ^ (-((1 / 2 : ℝ) + 7)) ≤ 1 / (c ^ 7 * cN) := by
  rw [Real.rpow_neg hc.le, Real.rpow_add hc, show (7 : ℝ) = ((7 : ℕ) : ℝ) by norm_num,
    Real.rpow_natCast, ← Real.sqrt_eq_rpow, one_div]
  have hs : cN ≤ Real.sqrt c := Real.le_sqrt_of_sq_le h
  apply inv_anti₀ (by positivity)
  rw [mul_comm]
  exact mul_le_mul_of_nonneg_right hs (by positivity)

lemma norm_four_cpow (t : ℝ) : ‖(4 : ℂ) ^ (-((1 / 2 : ℂ) + (t : ℂ) * I))‖ = 1 / 2 := by
  rw [show (4 : ℂ) = ((4 : ℝ) : ℂ) by norm_num, norm_cpow_eq_rpow_re_of_pos (by norm_num)]
  have hre : (-((1 / 2 : ℂ) + (t : ℂ) * I)).re = -(1 / 2) := by simp
  rw [hre, Real.rpow_neg (by norm_num), ← Real.sqrt_eq_rpow, show (4 : ℝ) = 2 ^ 2 by norm_num,
    Real.sqrt_sq (by norm_num)]
  norm_num

lemma ball_sub {z w : ℂ} {a b r c d q e f R : ℝ} (hz : ‖z - (a + b * I)‖ ≤ r)
    (hw : ‖w - (c + d * I)‖ ≤ q) (h : |a - c - e| + |b - d - f| + r + q ≤ R) :
    ‖z - w - (e + f * I)‖ ≤ R := by
  have hw' : ‖-w - ((-c : ℝ) + (-d : ℝ) * I)‖ ≤ q := by
    rw [show -w - ((-c : ℝ) + (-d : ℝ) * I) = -(w - (c + d * I)) by push_cast; ring, norm_neg]
    exact hw
  have := ball_add hz hw' (e := e) (f := f) (R := R) (by
    rw [show a + -c - e = a - c - e by ring, show b + -d - f = b - d - f by ring]; exact h)
  rwa [← sub_eq_add_neg] at this

/-- Turn a ball for the main term into a ball for L(1/2 + it, χ₋₄). -/
lemma beta_ball (t : ℝ) (ht : t ≠ 0) {A B R1 R2 R3 : ℝ}
    (hmain : ‖1 - nps ((1 / 2 : ℂ) + (t : ℂ) * I) 3 + nps ((1 / 2 : ℂ) + (t : ℂ) * I) 5
      - nps ((1 / 2 : ℂ) + (t : ℂ) * I) 7 + nps ((1 / 2 : ℂ) + (t : ℂ) * I) 9
      - nps ((1 / 2 : ℂ) + (t : ℂ) * I) 11
      + nps ((1 / 2 : ℂ) + (t : ℂ) * I) 13 * KB (13 / 4) ((1 / 2 : ℂ) + (t : ℂ) * I)
      - nps ((1 / 2 : ℂ) + (t : ℂ) * I) 15 * KB (15 / 4) ((1 / 2 : ℂ) + (t : ℂ) * I)
      - (A + B * I)‖ ≤ R1)
    (hem1 : ‖(1 / 2 : ℂ) + (t : ℂ) * I‖ * (‖gc ((1 / 2 : ℂ) + (t : ℂ) * I) 7‖ / 1000000 *
      ((1 / ((13 / 4 : ℝ) ^ 7 * (9 / 5))) / (1 / 2 + 7))) ≤ R2)
    (hem2 : ‖(1 / 2 : ℂ) + (t : ℂ) * I‖ * (‖gc ((1 / 2 : ℂ) + (t : ℂ) * I) 7‖ / 1000000 *
      ((1 / ((15 / 4 : ℝ) ^ 7 * (19 / 10))) / (1 / 2 + 7))) ≤ R3) :
    ‖DirichletCharacter.LFunction chiM4 ((1 / 2 : ℂ) + (t : ℂ) * I) - (A + B * I)‖ ≤
      R1 + (R2 + R3) / 2 := by
  set s : ℂ := (1 / 2 : ℂ) + (t : ℂ) * I with hsdef
  have hre : s.re = 1 / 2 := by simp [hsdef]
  have hs0 : 0 < s.re := by rw [hre]; norm_num
  have hs1 : s ≠ 1 := by
    intro h; have := congrArg Complex.im h; simp [hsdef] at this; exact ht this
  have E1 := hurwitz_EM8 (a := 1 / 4) ⟨by norm_num, by norm_num⟩ 3 hs0 hs1
  have E2 := hurwitz_EM8 (a := 3 / 4) ⟨by norm_num, by norm_num⟩ 3 hs0 hs1
  rw [← EMa] at E1 E2
  have hb1 : ‖s‖ * (‖gc s 7‖ / 1000000 * ((((3 : ℕ) : ℝ) + 1 / 4) ^ (-(s.re + 7)) / (s.re + 7)))
      ≤ R2 := by
    refine le_trans ?_ hem1
    rw [hre, show ((3 : ℕ) : ℝ) + 1 / 4 = 13 / 4 by norm_num]
    gcongr
    exact c_rpow_le (by norm_num) (by norm_num) (by norm_num)
  have hb2 : ‖s‖ * (‖gc s 7‖ / 1000000 * ((((3 : ℕ) : ℝ) + 3 / 4) ^ (-(s.re + 7)) / (s.re + 7)))
      ≤ R3 := by
    refine le_trans ?_ hem2
    rw [hre, show ((3 : ℕ) : ℝ) + 3 / 4 = 15 / 4 by norm_num]
    gcongr
    exact c_rpow_le (by norm_num) (by norm_num) (by norm_num)
  rw [LFunction_chiM4]
  have hm := beta_main_eq s
  have e : (4 : ℂ) ^ (-s) * (hurwitzZeta ((1 / 4 : ℝ) : UnitAddCircle) s
      - hurwitzZeta ((3 / 4 : ℝ) : UnitAddCircle) s) - (A + B * I) =
      (4 : ℂ) ^ (-s) * ((hurwitzZeta ((1 / 4 : ℝ) : UnitAddCircle) s - EMa (1 / 4) s)
        - (hurwitzZeta ((3 / 4 : ℝ) : UnitAddCircle) s - EMa (3 / 4) s))
      + ((4 : ℂ) ^ (-s) * (EMa (1 / 4) s - EMa (3 / 4) s) - (A + B * I)) := by ring
  rw [e, hm]
  have h4 : ‖(4 : ℂ) ^ (-s)‖ = 1 / 2 := norm_four_cpow t
  calc _ ≤ ‖(4 : ℂ) ^ (-s) * ((hurwitzZeta ((1 / 4 : ℝ) : UnitAddCircle) s - EMa (1 / 4) s)
        - (hurwitzZeta ((3 / 4 : ℝ) : UnitAddCircle) s - EMa (3 / 4) s))‖ + R1 :=
        (norm_add_le _ _).trans (by gcongr)
    _ ≤ 1 / 2 * (R2 + R3) + R1 := by
        rw [norm_mul, h4]
        gcongr
        exact (norm_sub_le _ _).trans (add_le_add (E1.trans hb1) (E2.trans hb2))
    _ = R1 + (R2 + R3) / 2 := by ring

end Fayed.Beta
