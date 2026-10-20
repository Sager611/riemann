import Riemann.Deep.FirstZero.SimpleNextLib

/-!
# Zero count up to height 33.4, part 1: an empty off-axis box

Same argument principle as SimpleNextLib.box_count, with the certified arg change around
[σL, 2] × [a, b] below 2π instead of 4π: the zeros in the box then have total order < 1, so there
are none. Zeros left of σL are excluded by the reflection ρ ↦ 1 − conj ρ (σL ≤ 1/2).
-/

open Complex Set

set_option linter.unusedSimpArgs false

namespace Fayed.FirstZero.Count
open Fayed.FirstZero Fayed.FirstZero.Min Fayed.FirstZero.SimpleNext

/-- **Box count**, general bound: a certified arg change < 2πk around [σL, 2] × [a, b] gives
total order < k for the zeros of ζ in the box. -/
theorem box_count_lt {σL a b Ubot Utop Uleft Lra Urb k : ℝ} (hσ0 : 0 < σL) (hσ1 : σL < 1)
    (ha : 0 < a) (hab : a < b)
    (hbot : BotOK a σL 2 Ubot) (htop : TopOK b σL 2 Utop) (hleft : LeftOK σL a b Uleft)
    (hra : Lra ≤ arg (riemannZeta ((2 : ℝ) + a * I)))
    (hrb : arg (riemannZeta ((2 : ℝ) + b * I)) ≤ Urb)
    (hsum : Ubot + Utop + Uleft + (Urb - Lra) < 2 * Real.pi * k) :
    ∃ Z : Finset ℂ,
      (∀ s ∈ Rectangle ((σL : ℂ) + a * I) (2 + b * I), riemannZeta s = 0 → s ∈ Z) ∧
      (∑ ρ ∈ Z, (analyticOrderNatAt riemannZeta ρ : ℝ)) < k := by
  classical
  set z : ℂ := (σL : ℂ) + a * I with hz
  set w : ℂ := 2 + b * I with hw
  have hzre : z.re = σL := by simp [z]
  have hzim : z.im = a := by simp [z]
  have hwre : w.re = 2 := by simp [w]
  have hwim : w.im = b := by simp [w]
  have hσ2 : σL ≤ 2 := by linarith
  have hcomp : IsCompact (Rectangle z w) := isCompact_uIcc.reProdIm isCompact_uIcc
  have hfin := Zeta23.RvM.riemannZeta_zeros_finite_of_isCompact hcomp
  set Z := hfin.toFinset with hZdef
  have hZmem : ∀ s, s ∈ Z ↔ s ∈ Rectangle z w ∧ s ≠ 1 ∧ riemannZeta s = 0 := by
    intro s; rw [hZdef, Set.Finite.mem_toFinset]; rfl
  have hne1 : ∀ s ∈ Rectangle z w, s ≠ 1 := by
    intro s hs h
    have h3 := ((mem_box hσ2 hab.le).mp hs).2.2.1
    rw [h, one_im] at h3; linarith
  have hZ : ∀ s ∈ Rectangle z w, riemannZeta s = 0 ↔ s ∈ Z := by
    intro s hs
    rw [hZmem]
    exact ⟨fun h => ⟨hs, hne1 s hs, h⟩, fun h => h.2.2⟩
  refine ⟨Z, fun s hs h0 => (hZ s hs).mp h0, ?_⟩
  have hstd : ∀ s : ℂ, (s.re : ℂ) + s.im * I = s := fun s => re_add_im s
  have hborder : ∀ s ∈ RectangleBorder z w, riemannZeta s ≠ 0 := by
    intro s hs
    simp only [RectangleBorder, mem_union, mem_reProdIm, mem_singleton_iff, hzre, hzim, hwre,
      hwim, uIcc_of_le hσ2, uIcc_of_le hab.le, mem_Icc] at hs
    rcases hs with ((⟨hr, hi⟩ | ⟨hr, hi⟩) | ⟨hr, hi⟩) | ⟨hr, _⟩
    · have h1 := hbot.1 s.re hr
      rw [← hi, hstd] at h1; exact h1
    · have h1 := hleft.1 s.im hi
      rw [← hr, hstd] at h1; exact h1
    · have h1 := htop.1 s.re hr
      rw [← hi, hstd] at h1; exact h1
    · exact riemannZeta_ne_zero_of_one_le_re (by rw [hr]; norm_num)
  have hf : AnalyticOnNhd ℂ riemannZeta (Rectangle z w) := fun s hs =>
    Zeta23.riemannZeta_analyticOnNhd_compl_one s (hne1 s hs)
  have hg : AnalyticOnNhd ℂ (fun _ : ℂ => (1 : ℂ)) (Rectangle z w) := fun _ _ => analyticAt_const
  have hZsub : (Z : Set ℂ) ⊆ Rectangle z w := fun s hs => ((hZmem s).mp hs).1
  have key := Zeta23.Analytic.rectangleIntegral'_mul_logDeriv (f := riemannZeta)
    (g := fun _ => 1) (by rw [hzre, hwre]; exact hσ2) (by rw [hzim, hwim]; exact hab.le)
    hf hg hborder Z hZ hZsub
  simp only [one_mul, mul_one] at key
  set S : ℝ := ∑ ρ ∈ Z, (analyticOrderNatAt riemannZeta ρ : ℝ) with hS
  have hSc : ((S : ℝ) : ℂ) = ∑ ρ ∈ Z, (analyticOrderNatAt riemannZeta ρ : ℂ) := by
    rw [hS]; push_cast; rfl
  have hpi : (2 * (Real.pi : ℂ) * I) ≠ 0 := by
    simp [Real.pi_ne_zero, I_ne_zero]
  have hRI : RectangleIntegral (fun s => logDeriv riemannZeta s) z w =
      2 * (Real.pi : ℂ) * I * ((S : ℝ) : ℂ) := by
    rw [hSc, ← key, RectangleIntegral', smul_eq_mul]
    field_simp
  have him := congrArg Complex.im hRI
  simp only [RectangleIntegral, HIntegral, VIntegral, hzre, hzim, hwre, hwim, smul_eq_mul,
    sub_im, add_im, mul_im, I_re, I_im] at him
  simp only [mul_im, mul_re, ofReal_re, ofReal_im, I_re, I_im, re_ofNat, im_ofNat] at him
  have hR := right_edge_ab hab.le
  have hb := hbot.2
  have ht := htop.2
  have hl := hleft.2
  have hS2 : 2 * Real.pi * S < 2 * Real.pi * k := by
    nlinarith [hR, hb, ht, hl, him]
  have := Real.pi_pos
  nlinarith

/-- **An empty box**: with the certified arg change around [σL, 2] × [a, b] below 2π, ζ has no zero
in the critical strip with a ≤ Im s ≤ b. -/
theorem no_zero_of_box {σL a b Ubot Utop Uleft Lra Urb : ℝ} (hσ0 : 0 < σL) (hσh : σL ≤ 1 / 2)
    (ha : 0 < a) (hab : a < b)
    (hbot : BotOK a σL 2 Ubot) (htop : TopOK b σL 2 Utop) (hleft : LeftOK σL a b Uleft)
    (hra : Lra ≤ arg (riemannZeta ((2 : ℝ) + a * I)))
    (hrb : arg (riemannZeta ((2 : ℝ) + b * I)) ≤ Urb)
    (hsum : Ubot + Utop + Uleft + (Urb - Lra) < 2 * Real.pi) :
    ∀ s : ℂ, riemannZeta s = 0 → 0 < s.re → s.re < 1 → a ≤ s.im → s.im ≤ b → False := by
  classical
  obtain ⟨Z, hZ, hlt⟩ := box_count_lt (k := 1) hσ0 (by linarith) ha hab hbot htop hleft hra hrb
    (by rw [mul_one]; exact hsum)
  have hσ2 : σL ≤ 2 := by linarith
  have key : ∀ x : ℂ, riemannZeta x = 0 → σL ≤ x.re → x.re < 1 → a ≤ x.im → x.im ≤ b → False := by
    intro x hx h1 h2 h3 h4
    have hxZ : x ∈ Z := hZ x ((mem_box hσ2 hab.le).mpr ⟨h1, by linarith, h3, h4⟩) hx
    have hx1 : 1 ≤ analyticOrderNatAt riemannZeta x := one_le_ord_of_zero hx (by linarith) h2
    have hsub := Finset.single_le_sum (f := fun ρ => (analyticOrderNatAt riemannZeta ρ : ℝ))
      (fun _ _ => Nat.cast_nonneg _) hxZ
    have e1 : (1 : ℝ) ≤ analyticOrderNatAt riemannZeta x := by exact_mod_cast hx1
    linarith
  intro s hs h0 h1 h2 h3
  rcases le_or_gt σL s.re with hσ | hσ
  · exact key s hs hσ h1 h2 h3
  · have hr := Zeta23.zeta_reflect_zero s ⟨hs, h0, h1⟩
    have hre : (Zeta23.reflect s).re = 1 - s.re := by simp [Zeta23.reflect]
    have hrim : (Zeta23.reflect s).im = s.im := by simp [Zeta23.reflect]
    exact key _ hr.1 (by rw [hre]; linarith) (by rw [hre]; linarith) (by rw [hrim]; exact h2)
      (by rw [hrim]; exact h3)

end Fayed.FirstZero.Count
