import Mathlib

/-!
# Davenport–Heilbronn zero, part 3: a Newton–Kantorovich existence lemma

newton_zero: if F is entire, ‖F z0‖ is small, a difference quotient of F at z0 is close to d ≠ 0
and F'' is bounded near z0, then z ↦ z − F z / d contracts the closed ball of radius r around z0,
so F has a zero there (Banach fixed point). Cauchy estimates turn sup bounds into derivative
bounds (deriv2_le, diff_quot_le).
-/

open Complex Set Metric Filter Topology

namespace Fayed.DH.Zero

/-- Cauchy estimate for F'' from a bound on F. -/
lemma deriv2_le {F : ℂ → ℂ} (hF : Differentiable ℂ F) {z0 : ℂ} {R ρ M : ℝ} (hρ : 0 < ρ)
    (hM : ∀ z, ‖z - z0‖ ≤ R + 2 * ρ → ‖F z‖ ≤ M) :
    ∀ u, ‖u - z0‖ ≤ R → ‖deriv (deriv F) u‖ ≤ M / ρ ^ 2 := by
  have hF' : Differentiable ℂ (deriv F) := fun z =>
    ((hF.differentiableOn.deriv isOpen_univ) z (mem_univ z)).differentiableAt
      (isOpen_univ.mem_nhds (mem_univ z))
  have h1 : ∀ v, ‖v - z0‖ ≤ R + ρ → ‖deriv F v‖ ≤ M / ρ := by
    intro v hv
    refine Complex.norm_deriv_le_of_forall_mem_sphere_norm_le hρ hF.diffContOnCl fun z hz => ?_
    refine hM z ?_
    rw [mem_sphere_iff_norm] at hz
    calc ‖z - z0‖ = ‖(z - v) + (v - z0)‖ := by ring_nf
      _ ≤ ‖z - v‖ + ‖v - z0‖ := norm_add_le _ _
      _ ≤ R + 2 * ρ := by linarith
  intro u hu
  have := Complex.norm_deriv_le_of_forall_mem_sphere_norm_le hρ hF'.diffContOnCl fun z hz => by
    refine h1 z ?_
    rw [mem_sphere_iff_norm] at hz
    calc ‖z - z0‖ = ‖(z - u) + (u - z0)‖ := by ring_nf
      _ ≤ ‖z - u‖ + ‖u - z0‖ := norm_add_le _ _
      _ ≤ R + ρ := by linarith
  rwa [div_div, ← sq] at this

/-- Difference quotient of a function bounded on a neighbourhood (Cauchy + mean value). -/
lemma diff_quot_le {E : ℂ → ℂ} {z0 k : ℂ} {ρ eE : ℝ} (hρ : 0 < ρ)
    (hE : ∀ z, ‖z - z0‖ ≤ ‖k‖ + ρ → DifferentiableAt ℂ E z)
    (hEb : ∀ z, ‖z - z0‖ ≤ ‖k‖ + ρ → ‖E z‖ ≤ eE) :
    ‖E (z0 + k) - E z0‖ ≤ eE / ρ * ‖k‖ := by
  set h := ‖k‖ with hhdef
  have hh : 0 ≤ h := norm_nonneg k
  have hD : ∀ u ∈ closedBall z0 h, ‖deriv E u‖ ≤ eE / ρ := by
    intro u hu
    rw [mem_closedBall, dist_eq_norm] at hu
    have hsub : closedBall u ρ ⊆ closedBall z0 (h + ρ) := by
      intro z hz
      rw [mem_closedBall, dist_eq_norm] at hz ⊢
      calc ‖z - z0‖ = ‖(z - u) + (u - z0)‖ := by ring_nf
        _ ≤ ‖z - u‖ + ‖u - z0‖ := norm_add_le _ _
        _ ≤ h + ρ := by linarith
    have hdo : DifferentiableOn ℂ E (closedBall z0 (h + ρ)) := fun z hz =>
      (hE z (by rwa [mem_closedBall, dist_eq_norm] at hz)).differentiableWithinAt
    refine Complex.norm_deriv_le_of_forall_mem_sphere_norm_le hρ (hdo.diffContOnCl_ball hsub)
      fun z hz => hEb z ?_
    rw [mem_sphere_iff_norm] at hz
    calc ‖z - z0‖ = ‖(z - u) + (u - z0)‖ := by ring_nf
      _ ≤ ‖z - u‖ + ‖u - z0‖ := norm_add_le _ _
      _ ≤ h + ρ := by linarith
  have hdiff : ∀ u ∈ closedBall z0 h, DifferentiableAt ℂ E u := fun u hu => hE u (by
    rw [mem_closedBall, dist_eq_norm] at hu; linarith)
  have := (convex_closedBall z0 h).norm_image_sub_le_of_norm_deriv_le hdiff hD
    (mem_closedBall_self hh) (by rw [mem_closedBall, dist_eq_norm, add_sub_cancel_left])
  rwa [add_sub_cancel_left] at this

/-- **Newton–Kantorovich zero existence.** -/
theorem newton_zero {F : ℂ → ℂ} (hF : Differentiable ℂ F) {z0 d k : ℂ} {r M2 e0 ed q : ℝ}
    (hr : 0 ≤ r) (hk : k ≠ 0)
    (hM2 : ∀ u, ‖u - z0‖ ≤ max r ‖k‖ → ‖deriv (deriv F) u‖ ≤ M2)
    (hF0 : ‖F z0‖ ≤ e0)
    (hD : ‖(F (z0 + k) - F z0) / k - d‖ ≤ ed)
    (hd : 0 < ‖d‖)
    (hq : ((r + ‖k‖) * M2 + ed) / ‖d‖ ≤ q) (hq1 : q < 1)
    (hc : e0 / ‖d‖ ≤ (1 - q) * r) :
    ∃ z, ‖z - z0‖ ≤ r ∧ F z = 0 := by
  set h := ‖k‖ with hhdef
  have hh : 0 < h := norm_pos_iff.mpr hk
  have hd0 : d ≠ 0 := norm_pos_iff.mp hd
  have hF' : Differentiable ℂ (deriv F) := fun z =>
    ((hF.differentiableOn.deriv isOpen_univ) z (mem_univ z)).differentiableAt
      (isOpen_univ.mem_nhds (mem_univ z))
  set B := closedBall z0 (max r h) with hB
  have hM2' : ∀ u ∈ B, ‖deriv (deriv F) u‖ ≤ M2 := fun u hu => hM2 u (by
    rwa [hB, mem_closedBall, dist_eq_norm] at hu)
  have hM2nn : 0 ≤ M2 := (norm_nonneg _).trans (hM2' z0 (mem_closedBall_self (by positivity)))
  -- (1) Lipschitz bound for F'
  have hL1 : ∀ w ∈ B, ‖deriv F w - deriv F z0‖ ≤ M2 * ‖w - z0‖ := fun w hw =>
    (convex_closedBall z0 (max r h)).norm_image_sub_le_of_norm_deriv_le
      (fun u _ => hF' u) hM2' (mem_closedBall_self (by positivity)) hw
  -- (2) the difference quotient is close to F'(z0)
  have hT : ‖(F (z0 + k) - F z0) / k - deriv F z0‖ ≤ M2 * h := by
    set φ : ℂ → ℂ := fun z => F z - z * deriv F z0 with hφ
    have hφd : ∀ u, HasDerivAt φ (deriv F u - deriv F z0) u := fun u => by
      have := ((hF u).hasDerivAt).sub ((hasDerivAt_id u).mul_const (deriv F z0))
      rw [one_mul] at this
      exact this
    have hbd : ∀ u ∈ closedBall z0 h, ‖deriv φ u‖ ≤ M2 * h := by
      intro u hu
      rw [(hφd u).deriv]
      have hu' : u ∈ B := closedBall_subset_closedBall (le_max_right r h) hu
      refine (hL1 u hu').trans ?_
      rw [mem_closedBall, dist_eq_norm] at hu
      exact mul_le_mul_of_nonneg_left hu hM2nn
    have hmv := (convex_closedBall z0 h).norm_image_sub_le_of_norm_deriv_le
      (fun u _ => (hφd u).differentiableAt) hbd (mem_closedBall_self hh.le)
      (by rw [mem_closedBall, dist_eq_norm, add_sub_cancel_left])
    have e : (F (z0 + k) - F z0) / k - deriv F z0 = (φ (z0 + k) - φ z0) / k := by
      simp only [hφ]; field_simp; ring
    rw [e, norm_div, ← hhdef, div_le_iff₀ hh]
    rwa [add_sub_cancel_left, ← hhdef] at hmv
  -- (3) derivative of the Newton map
  have hgd : ∀ w ∈ closedBall z0 r, ‖1 - deriv F w / d‖ ≤ q := by
    intro w hw
    have hw' : w ∈ B := closedBall_subset_closedBall (le_max_left r h) hw
    rw [mem_closedBall, dist_eq_norm] at hw
    have e : 1 - deriv F w / d = ((d - (F (z0 + k) - F z0) / k) +
        ((F (z0 + k) - F z0) / k - deriv F z0) + (deriv F z0 - deriv F w)) / d := by
      field_simp; ring
    rw [e, norm_div, div_le_iff₀ hd]
    have h1 := hL1 w hw'
    rw [norm_sub_rev] at h1
    have h0 : ‖d - (F (z0 + k) - F z0) / k‖ ≤ ed := by rw [norm_sub_rev]; exact hD
    have hq' := (div_le_iff₀ hd).mp hq
    calc ‖(d - (F (z0 + k) - F z0) / k) + ((F (z0 + k) - F z0) / k - deriv F z0) +
          (deriv F z0 - deriv F w)‖
        ≤ ed + M2 * h + M2 * ‖w - z0‖ := by
          refine (norm_add_le _ _).trans ?_
          refine add_le_add ((norm_add_le _ _).trans (add_le_add h0 hT)) h1
      _ ≤ ((r + h) * M2 + ed) := by nlinarith
      _ ≤ q * ‖d‖ := hq'
  set g : ℂ → ℂ := fun z => z - F z / d with hg
  have hgD : ∀ w, HasDerivAt g (1 - deriv F w / d) w := fun w => by
    exact (hasDerivAt_id w).sub ((hF w).hasDerivAt.div_const d)
  have hLip : ∀ x ∈ closedBall z0 r, ∀ y ∈ closedBall z0 r, ‖g y - g x‖ ≤ q * ‖y - x‖ :=
    fun x hx y hy => (convex_closedBall z0 r).norm_image_sub_le_of_norm_deriv_le
      (fun u _ => (hgD u).differentiableAt) (fun u hu => by rw [(hgD u).deriv]; exact hgd u hu) hx hy
  have hq0 : 0 ≤ q := by
    have := hgd z0 (mem_closedBall_self hr)
    linarith [norm_nonneg (1 - deriv F z0 / d)]
  have hmaps : MapsTo g (closedBall z0 r) (closedBall z0 r) := by
    intro z hz
    rw [mem_closedBall, dist_eq_norm]
    have h1 := hLip z0 (mem_closedBall_self hr) z hz
    rw [mem_closedBall, dist_eq_norm] at hz
    have h2 : ‖g z0 - z0‖ ≤ e0 / ‖d‖ := by
      simp only [hg, sub_sub_cancel_left, norm_neg, norm_div]
      exact div_le_div_of_nonneg_right hF0 hd.le
    calc ‖g z - z0‖ = ‖(g z - g z0) + (g z0 - z0)‖ := by ring_nf
      _ ≤ ‖g z - g z0‖ + ‖g z0 - z0‖ := norm_add_le _ _
      _ ≤ q * r + e0 / ‖d‖ := add_le_add (h1.trans (mul_le_mul_of_nonneg_left hz hq0)) h2
      _ ≤ r := by linarith
  set K : NNReal := ⟨q, hq0⟩ with hK
  have hcontr : ContractingWith K (hmaps.restrict g (closedBall z0 r) (closedBall z0 r)) := by
    refine ⟨by rw [hK]; exact_mod_cast hq1, LipschitzWith.of_dist_le_mul fun x y => ?_⟩
    simp only [Subtype.dist_eq, MapsTo.val_restrict_apply, dist_eq_norm, hK]
    exact hLip y y.2 x x.2
  obtain ⟨y, hy, hfix, -⟩ := hcontr.exists_fixedPoint' isClosed_closedBall.isComplete hmaps
    (mem_closedBall_self hr) (edist_ne_top _ _)
  refine ⟨y, by rwa [mem_closedBall, dist_eq_norm] at hy, ?_⟩
  have : F y / d = 0 := by
    have := hfix.eq
    simp only [hg] at this
    linear_combination -this
  rwa [div_eq_zero_iff, or_iff_left hd0] at this

end Fayed.DH.Zero
