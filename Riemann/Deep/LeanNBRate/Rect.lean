/-
Parts of this file are ported from PrimeNumberTheoremAnd (https://github.com/AlexKontorovich/PrimeNumberTheoremAnd),
file ResidueCalcOnRectangles.lean (via the anthropics/zeta-23-lean port, FromPNTPlus); Copyright the
PrimeNumberTheoremAnd contributors; Apache License 2.0 (http://www.apache.org/licenses/LICENSE-2.0).
Ported: the explicit rectangle integral of c/s (ResidueTheoremAtOrigin') and its helper integrals.
-/
import Mathlib.Analysis.Complex.CauchyIntegral
import Mathlib.Analysis.Complex.RemovableSingularity
import Mathlib.Analysis.SpecialFunctions.Integrals.Basic
import Mathlib.Analysis.SpecialFunctions.Trigonometric.Arctan

/-!
# Rectangle integrals, the residue of c/(s - p), Cauchy-Goursat with finitely many simple poles
-/

open Complex Filter Topology Set MeasureTheory intervalIntegral Real
open scoped Interval

noncomputable section

namespace Riemann.Deep.LeanNBRate

/-- The boundary integral of f over the rectangle with corners z (lower left) and w (upper right),
as in Mathlib's Cauchy-Goursat theorem. -/
def RectI (f : ℂ → ℂ) (z w : ℂ) : ℂ :=
  (∫ x : ℝ in z.re..w.re, f (x + z.im * I)) - (∫ x : ℝ in z.re..w.re, f (x + w.im * I)) +
    I * (∫ y : ℝ in z.im..w.im, f (w.re + y * I)) - I * ∫ y : ℝ in z.im..w.im, f (z.re + y * I)

theorem RectI_translate (f : ℂ → ℂ) (z w p : ℂ) :
    RectI (fun s => f (s - p)) z w = RectI f (z - p) (w - p) := by
  simp_rw [RectI, sub_re, sub_im, ← intervalIntegral.integral_comp_sub_right]
  congr <;> ext <;> congr 1 <;> simp [Complex.ext_iff]

lemma inv_re_add_im' {x y : ℝ} : ((x : ℂ) + y * I)⁻¹ = (x - I * y) / (x ^ 2 + y ^ 2) := by
  rw [Complex.inv_def, div_eq_mul_inv]
  congr <;> simp [conj_ofReal, normSq] <;> ring

lemma sq_add_sq_ne_zero' {x y : ℝ} (hy : y ≠ 0) : x ^ 2 + y ^ 2 ≠ 0 := by
  linarith [sq_nonneg x, sq_pos_iff.mpr hy]

lemma continuous_self_div_sq_add_sq' {y : ℝ} (hy : y ≠ 0) :
    Continuous fun x : ℝ => x / (x ^ 2 + y ^ 2) :=
  continuous_id.div (continuous_id.pow 2 |>.add continuous_const) (fun _ => sq_add_sq_ne_zero' hy)

lemma integral_self_div_sq_add_sq' {x₁ x₂ y : ℝ} (hy : y ≠ 0) :
    ∫ x in x₁..x₂, x / (x ^ 2 + y ^ 2) =
    Real.log (x₂ ^ 2 + y ^ 2) / 2 - Real.log (x₁ ^ 2 + y ^ 2) / 2 := by
  let f (x : ℝ) : ℝ := Real.log (x ^ 2 + y ^ 2) / 2
  have e1 {x} := HasDerivAt.add_const (y ^ 2) (by simpa using hasDerivAt_pow 2 x)
  have e2 {x} : HasDerivAt f (x / (x ^ 2 + y ^ 2)) x := by
    convert (e1.log (sq_add_sq_ne_zero' hy)).div_const 2 using 1
    field_simp
  have e3 : deriv f = fun x => x / (x ^ 2 + y ^ 2) := funext (fun _ => e2.deriv)
  have e4 : Continuous (deriv f) := by simpa only [e3] using continuous_self_div_sq_add_sq' hy
  simp_rw [← e2.deriv]
  exact integral_deriv_eq_sub (fun _ _ => e2.differentiableAt) (e4.intervalIntegrable _ _)

lemma integral_const_div_sq_add_sq' {x₁ x₂ y : ℝ} (hy : y ≠ 0) :
    ∫ x in x₁..x₂, y / (x ^ 2 + y ^ 2) = arctan (x₂ / y) - arctan (x₁ / y) := by
  nth_rewrite 1 [← div_mul_cancel₀ x₁ hy, ← div_mul_cancel₀ x₂ hy]
  simp_rw [← mul_integral_comp_mul_right, ← intervalIntegral.integral_const_mul,
    ← integral_one_div_one_add_sq]
  exact integral_congr fun x _ => by
    field_simp
    ring

lemma integral_const_div_self_add_im' {A : ℂ} {x₁ x₂ y : ℝ} (hy : y ≠ 0) :
    ∫ x : ℝ in x₁..x₂, A / (x + y * I) =
    A * (Real.log (x₂ ^ 2 + y ^ 2) / 2 - Real.log (x₁ ^ 2 + y ^ 2) / 2) -
    A * I * (arctan (x₂ / y) - arctan (x₁ / y)) := by
  have e1 {x : ℝ} : A / (x + y * I) = A * x / (x ^ 2 + y ^ 2) - A * I * y / (x ^ 2 + y ^ 2) := by
    ring_nf
    simp_rw [inv_re_add_im']
    ring
  have e2 : IntervalIntegrable (fun x : ℝ ↦ A * x / (x ^ 2 + y ^ 2)) volume x₁ x₂ := by
    apply Continuous.intervalIntegrable
    simp_rw [mul_div_assoc]
    norm_cast
    exact continuous_const.mul (continuous_ofReal.comp (continuous_self_div_sq_add_sq' hy))
  have e3 : IntervalIntegrable (fun x : ℝ ↦ A * I * y / (x ^ 2 + y ^ 2)) volume x₁ x₂ := by
    apply Continuous.intervalIntegrable
    refine continuous_const.div (by fun_prop) (fun x => ?_)
    norm_cast
    exact sq_add_sq_ne_zero' hy
  simp_rw [integral_congr (fun _ _ => e1), integral_sub e2 e3, mul_div_assoc]
  norm_cast
  simp_rw [intervalIntegral.integral_const_mul, intervalIntegral.integral_ofReal,
    integral_self_div_sq_add_sq' hy, integral_const_div_sq_add_sq' hy]

lemma integral_const_div_re_add_self' {A : ℂ} {x y₁ y₂ : ℝ} (hx : x ≠ 0) :
    ∫ y : ℝ in y₁..y₂, A / (x + y * I) =
    A / I * (Real.log (y₂ ^ 2 + (-x) ^ 2) / 2 - Real.log (y₁ ^ 2 + (-x) ^ 2) / 2) -
    A / I * I * (arctan (y₂ / -x) - arctan (y₁ / -x)) := by
  have l1 {y : ℝ} : A / (x + y * I) = A / I / (y + ↑(-x) * I) := by
    have e1 : (x : ℂ) + y * I ≠ 0 := by
      contrapose! hx
      simpa using congr_arg re hx
    have e2 : (y : ℂ) + I * ↑(-x) ≠ 0 := by
      contrapose! hx
      simpa using congr_arg im hx
    field_simp [*]
    push_cast
    ring_nf
    simp
  have l2 : -x ≠ 0 := by rwa [neg_ne_zero]
  simp_rw [l1, integral_const_div_self_add_im' l2]

/-- The rectangle integral of c/s around the origin is 2 pi i c. -/
lemma rectI_const_div_origin {z w c : ℂ}
    (h1 : z.re < 0) (h2 : z.im < 0) (h3 : 0 < w.re) (h4 : 0 < w.im) :
    RectI (fun s => c / s) z w = 2 * I * π * c := by
  simp only [RectI]
  rw [integral_const_div_re_add_self' h1.ne, integral_const_div_re_add_self' h3.ne.symm]
  rw [integral_const_div_self_add_im' h2.ne, integral_const_div_self_add_im' h4.ne.symm]
  have l1 : z.im * w.re⁻¹ = (w.re * z.im⁻¹)⁻¹ := by group
  have l3 := arctan_inv_of_neg <| mul_neg_of_pos_of_neg h3 <| inv_lt_zero.mpr h2
  have l4 : w.im * z.re⁻¹ = (z.re * w.im⁻¹)⁻¹ := by group
  have l6 := arctan_inv_of_neg <| mul_neg_of_neg_of_pos h1 <| inv_pos.mpr h4
  have r1 : z.im * z.re⁻¹ = (z.re * z.im⁻¹)⁻¹ := by group
  have r3 := arctan_inv_of_pos <| mul_pos_of_neg_of_neg h1 <| inv_lt_zero.mpr h2
  have r4 : w.im * w.re⁻¹ = (w.re * w.im⁻¹)⁻¹ := by group
  have r6 := arctan_inv_of_pos <| mul_pos h3 <| inv_pos.mpr h4
  ring_nf
  simp only [one_div, inv_I, mul_neg, neg_mul, I_sq, neg_neg, arctan_neg, ofReal_neg,
    sub_neg_eq_add]
  rw [l1, l3, l4, l6, r1, r3, r4, r6]
  ring_nf
  simp only [I_sq, ofReal_sub, ofReal_mul, ofReal_ofNat, ofReal_div, ofReal_neg, ofReal_one]
  ring_nf

/-- The rectangle integral of c/(s - p) is 2 pi i c for p in the open rectangle. -/
theorem rectI_const_div_sub {z w p c : ℂ} (h1 : z.re < p.re) (h2 : p.re < w.re) (h3 : z.im < p.im)
    (h4 : p.im < w.im) : RectI (fun s => c / (s - p)) z w = 2 * π * I * c := by
  rw [RectI_translate (fun s => c / s), rectI_const_div_origin (by simp; linarith) (by simp; linarith)
    (by simp; linarith) (by simp; linarith)]
  ring


theorem continuousAt_finset_sum' {ι : Type*} (P : Finset ι) {g : ι → ℂ → ℂ} {x : ℂ}
    (h : ∀ i ∈ P, ContinuousAt (g i) x) : ContinuousAt (fun s => ∑ i ∈ P, g i s) x :=
  tendsto_finsetSum P h

theorem RectI_eq_zero_of_off_countable {f : ℂ → ℂ} {z w : ℂ} (S : Set ℂ) (hS : S.Countable)
    (hc : ContinuousOn f ([[z.re, w.re]] ×ℂ [[z.im, w.im]]))
    (hd : ∀ x ∈ Ioo (min z.re w.re) (max z.re w.re) ×ℂ Ioo (min z.im w.im) (max z.im w.im) \ S,
      DifferentiableAt ℂ f x) : RectI f z w = 0 := by
  have h := Complex.integral_boundary_rect_eq_zero_of_differentiable_on_off_countable f z w S hS hc hd
  simpa [RectI, smul_eq_mul] using h

/-- The four side-integrability conditions for RectI. -/
structure RectIntegrable (f : ℂ → ℂ) (z w : ℂ) : Prop where
  bot : IntervalIntegrable (fun x : ℝ => f (x + z.im * I)) volume z.re w.re
  top : IntervalIntegrable (fun x : ℝ => f (x + w.im * I)) volume z.re w.re
  rgt : IntervalIntegrable (fun y : ℝ => f (w.re + y * I)) volume z.im w.im
  lft : IntervalIntegrable (fun y : ℝ => f (z.re + y * I)) volume z.im w.im

theorem RectI_sub {f g : ℂ → ℂ} {z w : ℂ} (hf : RectIntegrable f z w) (hg : RectIntegrable g z w) :
    RectI (fun s => f s - g s) z w = RectI f z w - RectI g z w := by
  simp only [RectI]
  rw [intervalIntegral.integral_sub hf.bot hg.bot, intervalIntegral.integral_sub hf.top hg.top,
    intervalIntegral.integral_sub hf.rgt hg.rgt, intervalIntegral.integral_sub hf.lft hg.lft]
  ring

theorem RectIntegrable.sub {f g : ℂ → ℂ} {z w : ℂ} (hf : RectIntegrable f z w)
    (hg : RectIntegrable g z w) : RectIntegrable (fun s => f s - g s) z w :=
  ⟨hf.bot.sub hg.bot, hf.top.sub hg.top, hf.rgt.sub hg.rgt, hf.lft.sub hg.lft⟩

theorem RectIntegrable.of_continuousOn {f : ℂ → ℂ} {z w : ℂ} (hzw1 : z.re ≤ w.re) (hzw2 : z.im ≤ w.im)
    (hf : ∀ x ∈ ([[z.re, w.re]] ×ℂ [[z.im, w.im]]), (x.re = z.re ∨ x.re = w.re ∨ x.im = z.im ∨ x.im = w.im) →
      ContinuousAt f x) : RectIntegrable f z w := by
  have hb : ∀ y ∈ ({z.im, w.im} : Set ℝ), IntervalIntegrable (fun x : ℝ => f (x + y * I)) volume z.re w.re := by
    intro y hy
    refine ContinuousOn.intervalIntegrable fun x hx => ?_
    have hm : (x : ℂ) + y * I ∈ ([[z.re, w.re]] ×ℂ [[z.im, w.im]]) := by
      rw [Complex.mem_reProdIm]
      have e1 : ((x : ℂ) + y * I).re = x := by simp
      have e2 : ((x : ℂ) + y * I).im = y := by simp
      rw [e1, e2]
      refine ⟨hx, ?_⟩
      rcases hy with rfl | rfl
      · exact Set.left_mem_uIcc
      · exact Set.right_mem_uIcc
    have := hf _ hm (by rcases hy with rfl | rfl <;> simp)
    have hcont : Continuous fun t : ℝ => (t : ℂ) + y * I := by fun_prop
    exact (ContinuousAt.comp (f := fun t : ℝ => (t : ℂ) + y * I) (x := x) this
      hcont.continuousAt).continuousWithinAt
  have hv : ∀ x ∈ ({z.re, w.re} : Set ℝ), IntervalIntegrable (fun y : ℝ => f (x + y * I)) volume z.im w.im := by
    intro x hx
    refine ContinuousOn.intervalIntegrable fun y hy => ?_
    have hm : (x : ℂ) + y * I ∈ ([[z.re, w.re]] ×ℂ [[z.im, w.im]]) := by
      rw [Complex.mem_reProdIm]
      have e1 : ((x : ℂ) + y * I).re = x := by simp
      have e2 : ((x : ℂ) + y * I).im = y := by simp
      rw [e1, e2]
      refine ⟨?_, hy⟩
      rcases hx with rfl | rfl
      · exact Set.left_mem_uIcc
      · exact Set.right_mem_uIcc
    have := hf _ hm (by rcases hx with rfl | rfl <;> simp)
    have hcont : Continuous fun t : ℝ => (x : ℂ) + t * I := by fun_prop
    exact (ContinuousAt.comp (f := fun t : ℝ => (x : ℂ) + t * I) (x := y) this
      hcont.continuousAt).continuousWithinAt
  exact ⟨hb _ (by simp), hb _ (by simp), hv _ (by simp), hv _ (by simp)⟩

/-- **Cauchy-Goursat with finitely many simple poles** in the open rectangle: if near each pole p,
f = A p/(s - p) + g with g continuous at p, then RectI f = 2 pi i sum A p. -/
theorem RectI_eq_sum_residues {f : ℂ → ℂ} {z w : ℂ} (hzw1 : z.re < w.re) (hzw2 : z.im < w.im)
    (P : Finset ℂ) (A : ℂ → ℂ)
    (hP : ∀ p ∈ P, z.re < p.re ∧ p.re < w.re ∧ z.im < p.im ∧ p.im < w.im)
    (hc : ∀ x ∈ ([[z.re, w.re]] ×ℂ [[z.im, w.im]]), x ∉ P → ContinuousAt f x)
    (hd : ∀ x ∈ Ioo z.re w.re ×ℂ Ioo z.im w.im, x ∉ P → DifferentiableAt ℂ f x)
    (hpole : ∀ p ∈ P, ∃ g : ℂ → ℂ, ContinuousAt g p ∧
      ∀ᶠ s in 𝓝[≠] p, f s = A p / (s - p) + g s) :
    RectI f z w = 2 * π * I * ∑ p ∈ P, A p := by
  classical
  set h : ℂ → ℂ := fun s => ∑ p ∈ P, A p / (s - p) with hh
  set G₀ : ℂ → ℂ := fun s => f s - h s with hG₀
  set G : ℂ → ℂ := fun s => if s ∈ P then limUnder (𝓝[≠] s) G₀ else G₀ s with hG
  have hPc : IsClosed (P : Set ℂ) := P.finite_toSet.isClosed
  -- local form of G near each point
  have hloc_out : ∀ x, x ∉ P → G =ᶠ[𝓝 x] G₀ := by
    intro x hx
    filter_upwards [hPc.isOpen_compl.mem_nhds hx] with y hy
    simp only [hG, if_neg (show y ∉ P from hy)]
  have hcont_h : ∀ x, x ∉ P → ContinuousAt h x := by
    intro x hx
    refine continuousAt_finset_sum' P fun p hp => ?_
    exact continuousAt_const.div (continuousAt_id.sub continuousAt_const)
      (sub_ne_zero.mpr fun e => hx (e ▸ hp))
  have hpole' : ∀ p ∈ P, ∃ L : ℂ, Tendsto G₀ (𝓝[≠] p) (𝓝 L) ∧ G p = L ∧
      ∃ K : ℂ → ℂ, ContinuousAt K p ∧ K p = L ∧ G₀ =ᶠ[𝓝[≠] p] K := by
    intro p hp
    obtain ⟨g, hgc, hgeq⟩ := hpole p hp
    set K : ℂ → ℂ := fun s => g s - ∑ q ∈ P.erase p, A q / (s - q) with hK
    have hKc : ContinuousAt K p := by
      refine hgc.sub (continuousAt_finset_sum' (P.erase p) fun q hq => ?_)
      exact continuousAt_const.div (continuousAt_id.sub continuousAt_const)
        (sub_ne_zero.mpr fun e => (Finset.mem_erase.mp hq).1 e.symm)
    have hev : G₀ =ᶠ[𝓝[≠] p] K := by
      filter_upwards [hgeq] with s hs
      simp only [hG₀, hh, hK, hs, ← Finset.add_sum_erase P _ hp]
      ring
    have ht : Tendsto G₀ (𝓝[≠] p) (𝓝 (K p)) :=
      (hKc.tendsto.mono_left nhdsWithin_le_nhds).congr' hev.symm
    refine ⟨K p, ht, ?_, K, hKc, rfl, hev⟩
    simp only [hG, if_pos hp]
    exact ht.limUnder_eq
  have hGcont : ContinuousOn G ([[z.re, w.re]] ×ℂ [[z.im, w.im]]) := by
    intro x hx
    refine ContinuousAt.continuousWithinAt ?_
    by_cases hxP : x ∈ P
    · obtain ⟨L, hL, hGL, K, hKc, hKL, hev⟩ := hpole' x hxP
      -- G = K near x
      have hGK : G =ᶠ[𝓝 x] K := by
        have h1 : ∀ᶠ y in 𝓝 x, y ≠ x → y ∉ P := by
          have : ((P.erase x : Finset ℂ) : Set ℂ)ᶜ ∈ 𝓝 x :=
            (P.erase x).finite_toSet.isClosed.isOpen_compl.mem_nhds (by simp)
          filter_upwards [this] with y hy hyx
          intro hyP; exact hy (Finset.mem_erase.mpr ⟨hyx, hyP⟩)
        have h2 : ∀ᶠ y in 𝓝 x, y ≠ x → G₀ y = K y := eventually_nhdsWithin_iff.mp hev
        filter_upwards [h1, h2] with y hy1 hy2
        by_cases hyx : y = x
        · rw [hyx, hGL, hKL]
        · simp only [hG, if_neg (hy1 hyx)]; exact hy2 hyx
      exact hKc.congr hGK.symm
    · have hfc := hc x hx hxP
      exact (hfc.sub (hcont_h x hxP)).congr (hloc_out x hxP).symm
  have hGd : ∀ x ∈ Ioo (min z.re w.re) (max z.re w.re) ×ℂ Ioo (min z.im w.im) (max z.im w.im) \ (P : Set ℂ),
      DifferentiableAt ℂ G x := by
    intro x ⟨hx, hxP⟩
    rw [min_eq_left hzw1.le, max_eq_right hzw1.le, min_eq_left hzw2.le, max_eq_right hzw2.le] at hx
    have hfd := hd x hx hxP
    have hhd : DifferentiableAt ℂ h x := by
      refine DifferentiableAt.fun_sum fun p hp => ?_
      exact (differentiableAt_const _).div (differentiableAt_id.sub (differentiableAt_const _))
        (sub_ne_zero.mpr fun e => hxP (e ▸ hp))
    exact (hfd.sub hhd).congr_of_eventuallyEq (hloc_out x hxP)
  have hG0 := RectI_eq_zero_of_off_countable (P : Set ℂ) P.finite_toSet.countable hGcont hGd
  -- on the boundary G = f - h
  have hbd : ∀ x ∈ ([[z.re, w.re]] ×ℂ [[z.im, w.im]]),
      (x.re = z.re ∨ x.re = w.re ∨ x.im = z.im ∨ x.im = w.im) → x ∉ P := by
    intro x _ hx hxP
    obtain ⟨a1, a2, a3, a4⟩ := hP x hxP
    rcases hx with h | h | h | h <;> linarith
  have hfi : RectIntegrable f z w := RectIntegrable.of_continuousOn hzw1.le hzw2.le
    fun x hx hb => hc x hx (hbd x hx hb)
  have hhi : ∀ p ∈ P, RectIntegrable (fun s => A p / (s - p)) z w := by
    intro p hp
    refine RectIntegrable.of_continuousOn hzw1.le hzw2.le fun x hx hb => ?_
    exact continuousAt_const.div (continuousAt_id.sub continuousAt_const)
      (sub_ne_zero.mpr fun e => hbd x hx hb (e ▸ hp))
  have hGeq : RectI G z w = RectI (fun s => f s - h s) z w := by
    simp only [RectI]
    have e : ∀ x ∈ ([[z.re, w.re]] ×ℂ [[z.im, w.im]]),
        (x.re = z.re ∨ x.re = w.re ∨ x.im = z.im ∨ x.im = w.im) → G x = f x - h x := by
      intro x hx hb
      simp only [hG, if_neg (hbd x hx hb), hG₀]
    congr 1
    · congr 1
      · congr 1
        · refine intervalIntegral.integral_congr fun x hx => e _ ?_ (by simp)
          exact ⟨by simpa using hx, by simp⟩
        · refine intervalIntegral.integral_congr fun x hx => e _ ?_ (by simp)
          exact ⟨by simpa using hx, by simp⟩
      · congr 1
        refine intervalIntegral.integral_congr fun y hy => e _ ?_ (by simp)
        exact ⟨by simp, by simpa using hy⟩
    · congr 1
      refine intervalIntegral.integral_congr fun y hy => e _ ?_ (by simp)
      exact ⟨by simp, by simpa using hy⟩
  have hsum : ∀ (Q : Finset ℂ), Q ⊆ P → RectIntegrable (fun s => ∑ p ∈ Q, A p / (s - p)) z w ∧
      RectI (fun s => ∑ p ∈ Q, A p / (s - p)) z w = ∑ p ∈ Q, 2 * π * I * A p := by
    intro Q
    induction Q using Finset.induction_on with
    | empty =>
      intro _
      refine ⟨⟨by simp, by simp, by simp, by simp⟩, by simp [RectI]⟩
    | insert q Q hq ih =>
      intro hQ
      have hqP : q ∈ P := hQ (Finset.mem_insert_self q Q)
      obtain ⟨ih1, ih2⟩ := ih (fun x hx => hQ (Finset.mem_insert_of_mem hx))
      have hq1 := hhi q hqP
      have e1 : (fun s => ∑ p ∈ insert q Q, A p / (s - p)) =
          fun s => A q / (s - q) - (-(∑ p ∈ Q, A p / (s - p))) := by
        funext s; rw [Finset.sum_insert hq]; ring
      have hneg : RectIntegrable (fun s => -(∑ p ∈ Q, A p / (s - p))) z w :=
        ⟨ih1.bot.neg, ih1.top.neg, ih1.rgt.neg, ih1.lft.neg⟩
      have hnegI : RectI (fun s => -(∑ p ∈ Q, A p / (s - p))) z w =
          -RectI (fun s => ∑ p ∈ Q, A p / (s - p)) z w := by
        simp only [RectI, intervalIntegral.integral_neg]; ring
      refine ⟨by rw [e1]; exact hq1.sub hneg, ?_⟩
      rw [e1, RectI_sub hq1 hneg, hnegI, ih2, Finset.sum_insert hq,
        rectI_const_div_sub (hP q hqP).1 (hP q hqP).2.1 (hP q hqP).2.2.1 (hP q hqP).2.2.2]
      ring
  obtain ⟨hsi, hsI⟩ := hsum P subset_rfl
  rw [hGeq, RectI_sub hfi hsi, hsI, ← Finset.mul_sum] at hG0
  exact sub_eq_zero.mp hG0


end Riemann.Deep.LeanNBRate
