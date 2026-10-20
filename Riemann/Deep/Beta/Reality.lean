import Riemann.Deep.Beta.HurwitzEM
import Mathlib.NumberTheory.LSeries.DirichletContinuation
import Mathlib.NumberTheory.LegendreSymbol.ZModChar
import Mathlib.Analysis.Fourier.ZMod

/-!
# χ₋₄: the L-function as a Hurwitz difference, conjugation symmetry, and reality on the line

* chiM4 (notation χ₋₄) = Mathlib's ZMod.χ₄ composed with ℤ → ℂ, an odd Dirichlet character mod 4.
* L(s, χ₋₄) = 4^{-s}(ζ(s, 1/4) - ζ(s, 3/4)) (unfolding Mathlib's definition).
* L(s̄) = conj L(s) on Re s > 0 (from the first-order Euler–Maclaurin formula of HurwitzEM.lean),
  hence Λ(s̄) = conj Λ(s) for the completed L-function Λ(s) = L(s) Γℝ(s + 1).
* With the odd functional equation Λ(1 - s) = 4^{s-1}·i·Λ(𝓕χ, s) and 𝓕χ = κ·χ (χ is odd on ZMod 4),
  F(t) = e^{i t log 4 / 2} Λ(1/2 + it) satisfies conj F(t) = ν F(t) for a constant ν.
* zero_of_conj_sign: such an F with Re(F(t₂) conj F(t₁)) < 0 vanishes in (t₁, t₂).
-/

open Complex Set Filter Topology MeasureTheory HurwitzZeta
open scoped ZMod

namespace Fayed.Beta

/-- The primitive Dirichlet character mod 4: Mathlib's ZMod.χ₄, valued in ℂ. -/
noncomputable def chiM4 : DirichletCharacter ℂ 4 := ZMod.χ₄.ringHomComp (Int.castRingHom ℂ)

@[inherit_doc] scoped notation "χ₋₄" => chiM4

lemma chiM4_apply (j : ZMod 4) : chiM4 j = ((ZMod.χ₄ j : ℤ) : ℂ) := rfl

lemma chiM4_vals : chiM4 0 = 0 ∧ chiM4 1 = 1 ∧ chiM4 2 = 0 ∧ chiM4 3 = -1 := by
  have h0 : ZMod.χ₄ 0 = 0 := by decide
  have h1 : ZMod.χ₄ 1 = 1 := by decide
  have h2 : ZMod.χ₄ 2 = 0 := by decide
  have h3 : ZMod.χ₄ 3 = -1 := by decide
  simp [chiM4_apply, h0, h1, h2, h3]

lemma chiM4_odd : chiM4.Odd := by
  show chiM4 (-1) = -1
  have : (-1 : ZMod 4) = 3 := by decide
  rw [this]; exact chiM4_vals.2.2.2

lemma sum_zmod4 (g : ZMod 4 → ℂ) : ∑ j : ZMod 4, g j = g 0 + g 1 + g 2 + g 3 :=
  Fin.sum_univ_four g

/-- **L(s, χ₋₄) = 4^{-s} (ζ(s, 1/4) - ζ(s, 3/4)).** -/
lemma LFunction_chiM4 (s : ℂ) : DirichletCharacter.LFunction chiM4 s =
    (4 : ℂ) ^ (-s) * (hurwitzZeta ((1 / 4 : ℝ) : UnitAddCircle) s
      - hurwitzZeta ((3 / 4 : ℝ) : UnitAddCircle) s) := by
  rw [DirichletCharacter.LFunction, ZMod.LFunction, sum_zmod4]
  obtain ⟨h0, h1, h2, h3⟩ := chiM4_vals
  have v1 : ((1 : ZMod 4).val : ℝ) / ((4 : ℕ) : ℝ) = 1 / 4 := by
    have : (1 : ZMod 4).val = 1 := rfl
    rw [this]; norm_num
  have v3 : ((3 : ZMod 4).val : ℝ) / ((4 : ℕ) : ℝ) = 3 / 4 := by
    have : (3 : ZMod 4).val = 3 := rfl
    rw [this]; norm_num
  rw [h0, h1, h2, h3, ZMod.toAddCircle_apply (N := 4) 1, ZMod.toAddCircle_apply (N := 4) 3, v1, v3]
  push_cast
  ring

lemma ofReal_cpow_conj {x : ℝ} (hx : 0 ≤ x) (w : ℂ) :
    (x : ℂ) ^ (starRingEnd ℂ w) = starRingEnd ℂ ((x : ℂ) ^ w) := by
  rw [cpow_conj _ _ (by rw [arg_ofReal_of_nonneg hx]; exact Real.pi_pos.ne), conj_ofReal]

lemma hurwitzFO_conj {a : ℝ} (ha : 0 < a) (J : ℕ) (s : ℂ) :
    hurwitzFO a J (starRingEnd ℂ s) = starRingEnd ℂ (hurwitzFO a J s) := by
  have hc : 0 ≤ (J : ℝ) + a := by positivity
  have hI : ∫ y in Ioi ((J : ℝ) + a), (saw a y : ℂ) * (y : ℂ) ^ (-(starRingEnd ℂ s + 1)) =
      starRingEnd ℂ (∫ y in Ioi ((J : ℝ) + a), (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1))) := by
    rw [← integral_conj]
    refine setIntegral_congr_fun measurableSet_Ioi fun y hy => ?_
    have hy0 : 0 ≤ y := le_trans hc (le_of_lt hy)
    rw [map_mul, conj_ofReal, ← ofReal_cpow_conj hy0]
    simp [map_neg, map_add, map_one]
  have hS : ∑ j ∈ Finset.range J, (((j : ℝ) + a : ℝ) : ℂ) ^ (-starRingEnd ℂ s) =
      starRingEnd ℂ (∑ j ∈ Finset.range J, (((j : ℝ) + a : ℝ) : ℂ) ^ (-s)) := by
    rw [map_sum]
    refine Finset.sum_congr rfl fun j _ => ?_
    rw [← ofReal_cpow_conj (by positivity), map_neg]
  have h1 : (((J : ℝ) + a : ℝ) : ℂ) ^ (1 - starRingEnd ℂ s) =
      starRingEnd ℂ ((((J : ℝ) + a : ℝ) : ℂ) ^ (1 - s)) := by
    rw [← ofReal_cpow_conj hc, map_sub, map_one]
  have h2 : (((J : ℝ) + a : ℝ) : ℂ) ^ (-starRingEnd ℂ s) =
      starRingEnd ℂ ((((J : ℝ) + a : ℝ) : ℂ) ^ (-s)) := by
    rw [← ofReal_cpow_conj hc, map_neg]
  simp only [hurwitzFO, hI, hS, h1, h2, map_add, map_mul, map_div₀, map_neg, map_sub, map_one,
    map_ofNat]

lemma hurwitz_conj {a : ℝ} (ha : a ∈ Ioc (0 : ℝ) 1) {s : ℂ} (hs : 0 < s.re) (hs1 : s ≠ 1) :
    hurwitzZeta a (starRingEnd ℂ s) = starRingEnd ℂ (hurwitzZeta a s) := by
  have hs' : 0 < (starRingEnd ℂ s).re := by rwa [conj_re]
  have hs1' : starRingEnd ℂ s ≠ 1 := by
    intro h; apply hs1; rw [← conj_conj s, h, map_one]
  rw [hurwitz_first_order ha 0 hs' hs1', hurwitz_first_order ha 0 hs hs1, hurwitzFO_conj ha.1]

lemma LFunction_conj {s : ℂ} (hs : 0 < s.re) (hs1 : s ≠ 1) :
    DirichletCharacter.LFunction chiM4 (starRingEnd ℂ s) =
      starRingEnd ℂ (DirichletCharacter.LFunction chiM4 s) := by
  have h4 : (4 : ℂ) = ((4 : ℝ) : ℂ) := by norm_num
  rw [LFunction_chiM4, LFunction_chiM4,
    hurwitz_conj (a := 1 / 4) ⟨by norm_num, by norm_num⟩ hs hs1,
    hurwitz_conj (a := 3 / 4) ⟨by norm_num, by norm_num⟩ hs hs1, h4,
    show -(starRingEnd ℂ s) = starRingEnd ℂ (-s) by rw [map_neg],
    ofReal_cpow_conj (by norm_num), map_mul, map_sub]

lemma Gammaℝ_conj (s : ℂ) : Gammaℝ (starRingEnd ℂ s) = starRingEnd ℂ (Gammaℝ s) := by
  simp only [Gammaℝ]
  rw [map_mul, ← Gamma_conj, ← ofReal_cpow_conj Real.pi_pos.le]
  congr 2 <;> simp [map_div₀, map_neg, map_ofNat]

lemma completed_eq_mul {s : ℂ} (hs : 0 < (s + 1).re) :
    DirichletCharacter.completedLFunction chiM4 s =
      DirichletCharacter.LFunction chiM4 s * Gammaℝ (s + 1) := by
  rw [DirichletCharacter.LFunction_eq_completed_div_gammaFactor chiM4 s (Or.inr (by norm_num)),
    chiM4_odd.gammaFactor_def, div_mul_cancel₀ _ (Gammaℝ_ne_zero_of_re_pos hs)]

lemma completed_conj {s : ℂ} (hs : 0 < s.re) (hs1 : s ≠ 1) :
    DirichletCharacter.completedLFunction chiM4 (starRingEnd ℂ s) =
      starRingEnd ℂ (DirichletCharacter.completedLFunction chiM4 s) := by
  have h1 : 0 < (s + 1).re := by simp; linarith
  have h2 : 0 < (starRingEnd ℂ s + 1).re := by simp; linarith
  rw [completed_eq_mul h1, completed_eq_mul h2, LFunction_conj hs hs1,
    show starRingEnd ℂ s + 1 = starRingEnd ℂ (s + 1) by rw [map_add, map_one], Gammaℝ_conj, map_mul]

/-- The discrete Fourier transform of χ₋₄ is a multiple of χ₋₄. -/
lemma dft_chiM4 : 𝓕 (fun j : ZMod 4 => chiM4 j) = (𝓕 (fun j : ZMod 4 => chiM4 j) 1) •
    (fun j : ZMod 4 => chiM4 j) := by
  have hΦ : Function.Odd (fun j : ZMod 4 => chiM4 j) := chiM4_odd.to_fun
  have hF : Function.Odd (𝓕 (fun j : ZMod 4 => chiM4 j)) := ZMod.dft_odd_iff.mpr hΦ
  obtain ⟨h0, h1, h2, h3⟩ := chiM4_vals
  funext j
  rw [Pi.smul_apply, smul_eq_mul]
  have hj : j = 0 ∨ j = 1 ∨ j = 2 ∨ j = 3 := by revert j; decide
  rcases hj with rfl | rfl | rfl | rfl
  · rw [hF.map_zero, h0, mul_zero]
  · rw [h1, mul_one]
  · have e := hF 2
    have : (-2 : ZMod 4) = 2 := by decide
    rw [this] at e
    have hz : 𝓕 (fun j : ZMod 4 => chiM4 j) 2 = 0 := by
      have : (2 : ℂ) * 𝓕 (fun j : ZMod 4 => chiM4 j) 2 = 0 := by linear_combination e
      simpa using this
    rw [hz, h2, mul_zero]
  · have e := hF 1
    have : (-1 : ZMod 4) = 3 := by decide
    rw [this] at e
    rw [h3, e]; ring

lemma completed_smul (κ : ℂ) (Φ : ZMod 4 → ℂ) (s : ℂ) :
    ZMod.completedLFunction (κ • Φ) s = κ * ZMod.completedLFunction Φ s := by
  simp only [ZMod.completedLFunction, Pi.smul_apply, smul_eq_mul, mul_assoc, ← Finset.mul_sum]
  ring

/-- **Functional equation for χ₋₄** with κ = 𝓕χ₋₄(1). -/
lemma completed_one_sub (s : ℂ) :
    DirichletCharacter.completedLFunction chiM4 (1 - s) =
      (4 : ℂ) ^ (s - 1) * I * 𝓕 (fun j : ZMod 4 => chiM4 j) 1 *
        DirichletCharacter.completedLFunction chiM4 s := by
  have hΦ : Function.Odd (fun j : ZMod 4 => chiM4 j) := chiM4_odd.to_fun
  have hd : ZMod.completedLFunction (𝓕 (fun j : ZMod 4 => chiM4 j)) s =
      𝓕 (fun j : ZMod 4 => chiM4 j) 1 * ZMod.completedLFunction (fun j : ZMod 4 => chiM4 j) s := by
    conv_lhs => rw [dft_chiM4]
    exact completed_smul _ _ _
  rw [DirichletCharacter.completedLFunction, DirichletCharacter.completedLFunction,
    ZMod.completedLFunction_one_sub_odd hΦ, hd]
  push_cast
  ring

/-- The rotated completed L-function on the critical line. -/
noncomputable def Fb (t : ℝ) : ℂ :=
  cexp (I * ((t * Real.log 4 / 2 : ℝ) : ℂ)) *
    DirichletCharacter.completedLFunction chiM4 ((1 / 2 : ℂ) + (t : ℂ) * I)

/-- conj F(t) = ν F(t) with ν independent of t. -/
lemma Fb_conj (t : ℝ) : starRingEnd ℂ (Fb t) =
    ((Real.exp (-(Real.log 4 / 2)) : ℝ) * I * 𝓕 (fun j : ZMod 4 => chiM4 j) 1) * Fb t := by
  set s : ℂ := (1 / 2 : ℂ) + (t : ℂ) * I with hs
  have hre : 0 < s.re := by simp [hs]
  have hs1 : s ≠ 1 := by
    intro h; have := congrArg Complex.re h; norm_num [hs] at this
  have hcs : starRingEnd ℂ s = 1 - s := by
    apply Complex.ext <;> simp [hs]; norm_num
  have hL : starRingEnd ℂ (DirichletCharacter.completedLFunction chiM4 s) =
      (4 : ℂ) ^ (s - 1) * I * 𝓕 (fun j : ZMod 4 => chiM4 j) 1 *
        DirichletCharacter.completedLFunction chiM4 s := by
    rw [← completed_conj hre hs1, hcs, completed_one_sub]
  have h4 : (4 : ℂ) ^ (s - 1) = cexp (((Real.log 4 : ℝ) : ℂ) * (s - 1)) := by
    rw [cpow_def_of_ne_zero (by norm_num), ofReal_log (by norm_num)]
    norm_num
  simp only [Fb]
  rw [← hs, map_mul, hL, h4, ← exp_conj]
  have e1 : cexp (starRingEnd ℂ (I * ((t * Real.log 4 / 2 : ℝ) : ℂ))) *
      (cexp (((Real.log 4 : ℝ) : ℂ) * (s - 1))) =
      ((Real.exp (-(Real.log 4 / 2)) : ℝ) : ℂ) * cexp (I * ((t * Real.log 4 / 2 : ℝ) : ℂ)) := by
    rw [ofReal_exp, ← Complex.exp_add, ← Complex.exp_add]
    congr 1
    rw [map_mul, conj_I, conj_ofReal, hs]
    push_cast
    ring
  linear_combination (I * 𝓕 (fun j : ZMod 4 => chiM4 j) 1 *
    DirichletCharacter.completedLFunction chiM4 s) * e1

/-- A continuous F with conj F = ν F and Re(F(t₂) conj F(t₁)) < 0 vanishes in (t₁, t₂). -/
lemma zero_of_conj_sign {F : ℝ → ℂ} {t1 t2 : ℝ} (hlt : t1 < t2)
    (hcont : ContinuousOn F (Icc t1 t2)) {ν : ℂ} (hconj : ∀ t, starRingEnd ℂ (F t) = ν * F t)
    (hsign : (F t2 * starRingEnd ℂ (F t1)).re < 0) : ∃ t ∈ Ioo t1 t2, F t = 0 := by
  have hF1 : F t1 ≠ 0 := by intro h; rw [h] at hsign; simp at hsign
  have hν : starRingEnd ℂ ν * ν = 1 := by
    have h := congrArg (starRingEnd ℂ) (hconj t1)
    rw [conj_conj, map_mul, hconj t1] at h
    have h' : (starRingEnd ℂ ν * ν - 1) * F t1 = 0 := by linear_combination -h
    rcases mul_eq_zero.mp h' with h'' | h''
    · linear_combination h''
    · exact absurd h'' hF1
  obtain ⟨ω, hω0, hω⟩ : ∃ ω : ℂ, ω ≠ 0 ∧ starRingEnd ℂ ω * ν = ω := by
    by_cases hν1 : ν = -1
    · exact ⟨I, I_ne_zero, by rw [hν1, conj_I]; ring⟩
    · refine ⟨1 + ν, fun h => hν1 (by linear_combination h), ?_⟩
      rw [map_add, map_one]; linear_combination hν
  have hreal : ∀ t, (ω * F t).im = 0 := by
    intro t
    apply conj_eq_iff_im.mp
    rw [map_mul, hconj t, ← mul_assoc, hω]
  set G : ℝ → ℝ := fun t => (ω * F t).re with hG
  have hGcont : ContinuousOn G (Icc t1 t2) :=
    continuous_re.comp_continuousOn (continuousOn_const.mul hcont)
  have hsign' : G t1 * G t2 < 0 := by
    have e : (ω * F t2 * starRingEnd ℂ (ω * F t1)) =
        ((normSq ω : ℝ) : ℂ) * (F t2 * starRingEnd ℂ (F t1)) := by
      rw [map_mul, ← mul_conj]; ring
    have hre : (ω * F t2 * starRingEnd ℂ (ω * F t1)).re = G t2 * G t1 := by
      simp only [mul_re (ω * F t2), conj_re, conj_im, hreal, hG]; ring
    rw [e, re_ofReal_mul] at hre
    have hn : 0 < normSq ω := normSq_pos.mpr hω0
    nlinarith [mul_neg_of_pos_of_neg hn hsign]
  obtain ⟨t, ht, hGt⟩ : ∃ t ∈ Ioo t1 t2, G t = 0 := by
    rcases lt_or_gt_of_ne (show G t1 ≠ 0 by
        intro h; rw [h, zero_mul] at hsign'; exact lt_irrefl _ hsign') with hneg | hpos
    · have h2 : 0 < G t2 := by by_contra h; push Not at h; nlinarith
      exact intermediate_value_Ioo hlt.le hGcont ⟨hneg, h2⟩
    · have h2 : G t2 < 0 := by by_contra h; push Not at h; nlinarith
      exact intermediate_value_Ioo' hlt.le hGcont ⟨h2, hpos⟩
  refine ⟨t, ht, ?_⟩
  have h0 : ω * F t = 0 := Complex.ext hGt (hreal t)
  rcases mul_eq_zero.mp h0 with h | h
  · exact absurd h hω0
  · exact h

end Fayed.Beta
