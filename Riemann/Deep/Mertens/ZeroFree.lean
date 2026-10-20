import Riemann.Tree.R002

/-!
# Mertens bound ⇒ zero-free half-plane (Littlewood's direction)

If |M(N)| ≤ C N^θ then ζ(s) ≠ 0 for Re s > θ.

Proof: g(s) := s · mellin M(⌊·⌋) (-s) is holomorphic on Re s > θ
(mellin_differentiableAt_of_isBigO_rpow), equals L(μ, s) = 1/ζ(s) on Re s > max θ 1
(LSeries_eq_mul_integral), so ζ·g = 1 on the connected set {Re s > θ} \ {1} by the identity
theorem; s = 1 is riemannZeta_one_ne_zero.
-/

open Complex Asymptotics Filter Topology MeasureTheory Set

namespace Fayed.Mertens

/-- Möbius coefficients as complex numbers (defeq to the LSeries coercion of μ). -/
noncomputable def mu : ℕ → ℂ := fun n => (ArithmeticFunction.moebius n : ℂ)

/-- Mertens function M(n) = Σ_{k ≤ n} μ(k), complex-valued. -/
noncomputable def Msum (n : ℕ) : ℂ := ∑ k ∈ Finset.Icc 1 n, mu k

/-- Step function t ↦ M(⌊t⌋). -/
noncomputable def Mstep (t : ℝ) : ℂ := Msum ⌊t⌋₊

/-- g(s) = s ∫_1^∞ M(t) t^{-s-1} dt, written as a Mellin transform. -/
noncomputable def g (s : ℂ) : ℂ := s * mellin Mstep (-s)

lemma Msum_eq (n : ℕ) :
    Msum n = ((∑ k ∈ Finset.Icc 1 n, (ArithmeticFunction.moebius k : ℝ) : ℝ) : ℂ) := by
  simp [Msum, mu]

lemma isBigO_Msum {θ : ℝ} (h : Fayed.R2.MertensBound θ) :
    (fun n : ℕ => Msum n) =O[atTop] (fun n : ℕ => (n : ℝ) ^ θ) := by
  obtain ⟨C, hC⟩ := h
  refine IsBigO.of_bound C (Eventually.of_forall fun n => ?_)
  rw [Msum_eq, Complex.norm_real, Real.norm_eq_abs,
    Real.norm_of_nonneg (Real.rpow_nonneg (Nat.cast_nonneg _) _)]
  exact hC n

lemma norm_mu_le (k : ℕ) : ‖mu k‖ ≤ 1 := by
  have h := ArithmeticFunction.abs_moebius_le_one (n := k)
  rw [mu, Complex.norm_intCast]
  exact_mod_cast h

lemma norm_Msum_le (n : ℕ) : ‖Msum n‖ ≤ n := by
  calc ‖Msum n‖ ≤ ∑ k ∈ Finset.Icc 1 n, ‖mu k‖ := norm_sum_le _ _
    _ ≤ ∑ k ∈ Finset.Icc 1 n, (1 : ℝ) := Finset.sum_le_sum fun k _ => norm_mu_le k
    _ = n := by simp

lemma Mstep_eq_zero {t : ℝ} (ht : t < 1) : Mstep t = 0 := by
  have : ⌊t⌋₊ = 0 := Nat.floor_eq_zero.mpr ht
  simp [Mstep, Msum, this]

lemma isBigO_Mstep_atTop {θ : ℝ} (h : Fayed.R2.MertensBound θ) :
    Mstep =O[atTop] (fun t : ℝ => t ^ (-(-θ))) := by
  have h1 : (fun t : ℝ => (⌊t⌋₊ : ℝ) ^ θ) =O[atTop] (fun t : ℝ => t ^ θ) :=
    ((isEquivalent_nat_floor (R := ℝ)).symm.rpow (r := θ)
      (fun t => Nat.cast_nonneg ⌊t⌋₊)).symm.isBigO
  rw [neg_neg]
  exact ((isBigO_Msum h).comp_tendsto tendsto_nat_floor_atTop).trans h1

lemma isBigO_Mstep_zero (b : ℝ) : Mstep =O[𝓝[>] 0] (fun t : ℝ => t ^ (-b)) := by
  have : Mstep =ᶠ[𝓝[>] (0 : ℝ)] (fun _ => (0 : ℂ)) := by
    filter_upwards [Ioo_mem_nhdsGT (zero_lt_one' ℝ)] with t ht using Mstep_eq_zero ht.2
  exact (isBigO_zero _ _).congr' this.symm EventuallyEq.rfl

lemma measurable_Mstep : Measurable Mstep :=
  show Measurable (Msum ∘ Nat.floor) from (measurable_from_nat (f := Msum)).comp Nat.measurable_floor

lemma locallyIntegrableOn_Mstep : LocallyIntegrableOn Mstep (Ioi 0) := by
  refine ((continuous_id (X := ℝ)).locallyIntegrable.locallyIntegrableOn (Ioi 0)).mono
    measurable_Mstep.aestronglyMeasurable (Eventually.of_forall fun t => ?_)
  refine (norm_Msum_le _).trans ?_
  rw [id, Real.norm_eq_abs]
  rcases le_or_gt 0 t with ht | ht
  · exact (Nat.floor_le ht).trans (le_abs_self t)
  · rw [Nat.floor_of_nonpos ht.le, Nat.cast_zero]; exact abs_nonneg t

lemma integral_eq_mellin (s : ℂ) :
    ∫ t in Ioi (1 : ℝ), (∑ k ∈ Finset.Icc 1 ⌊t⌋₊, mu k) * (t : ℂ) ^ (-(s + 1)) =
      mellin Mstep (-s) := by
  rw [mellin, ← integral_Ici_eq_integral_Ioi,
    setIntegral_eq_of_subset_of_ae_sdiff_eq_zero measurableSet_Ioi.nullMeasurableSet
      (fun t (ht : 1 ≤ t) => (zero_lt_one.trans_le ht : (0 : ℝ) < t))
      (Eventually.of_forall fun t ht => ?_)]
  · refine setIntegral_congr_fun measurableSet_Ici fun t _ => ?_
    rw [smul_eq_mul, mul_comm, Mstep, Msum, show -s - 1 = -(s + 1) by ring]
  · have : t < 1 := not_le.mp ht.2
    rw [Mstep_eq_zero this, smul_zero]

lemma LSeries_eq_g {θ : ℝ} (h : Fayed.R2.MertensBound θ) {s : ℂ} (hs : max θ 1 < s.re) :
    LSeries mu s = g s := by
  have hs1 : 1 < s.re := lt_of_le_of_lt (le_max_right _ _) hs
  have hr : max θ 0 < s.re := max_lt (by linarith [le_max_left θ 1]) (by linarith)
  have hO : (fun n : ℕ => ∑ k ∈ Finset.Icc 1 n, mu k) =O[atTop]
      (fun n : ℕ => (n : ℝ) ^ (max θ 0)) := by
    refine (isBigO_Msum h).trans (IsBigO.of_bound 1 ?_)
    filter_upwards [eventually_ge_atTop 1] with n hn
    rw [one_mul, Real.norm_of_nonneg (by positivity), Real.norm_of_nonneg (by positivity)]
    exact Real.rpow_le_rpow_of_exponent_le (by exact_mod_cast hn) (le_max_left _ _)
  have hS : LSeriesSummable mu s := ArithmeticFunction.LSeriesSummable_moebius_iff.mpr hs1
  rw [LSeries_eq_mul_integral mu (le_max_right θ 0) hr hS hO, g, integral_eq_mellin]

lemma differentiableAt_g {θ : ℝ} (h : Fayed.R2.MertensBound θ) {s : ℂ} (hs : θ < s.re) :
    DifferentiableAt ℂ g s := by
  have hm : DifferentiableAt ℂ (mellin Mstep) (-s) :=
    mellin_differentiableAt_of_isBigO_rpow (a := -θ) (b := -s.re - 1) locallyIntegrableOn_Mstep
      (isBigO_Mstep_atTop h) (by rw [neg_re]; linarith) (isBigO_Mstep_zero _)
      (by rw [neg_re]; linarith)
  exact differentiableAt_id.mul (hm.comp s differentiableAt_id.neg)

lemma isPreconnected_U (θ : ℝ) : IsPreconnected ({s : ℂ | θ < s.re} \ {1}) := by
  let c : ℂ := ((max θ 1 + 1 : ℝ) : ℂ)
  have hθ1 := le_max_left θ 1
  have h11 := le_max_right θ 1
  have cP1 : IsPreconnected ({s : ℂ | θ < s.re} ∩ {s | 0 < s.im}) :=
    ((convex_halfSpace_re_gt θ).inter (convex_halfSpace_im_gt 0)).isPreconnected
  have cP2 : IsPreconnected ({s : ℂ | θ < s.re} ∩ {s | s.im < 0}) :=
    ((convex_halfSpace_re_gt θ).inter (convex_halfSpace_im_lt 0)).isPreconnected
  have cP3 : IsPreconnected {s : ℂ | max θ 1 < s.re} :=
    (convex_halfSpace_re_gt _).isPreconnected
  have cP4 : IsPreconnected ({s : ℂ | θ < s.re} ∩ {s | s.re < 1}) :=
    ((convex_halfSpace_re_gt θ).inter (convex_halfSpace_re_lt 1)).isPreconnected
  have sP1 : {s : ℂ | θ < s.re} ∩ {s | 0 < s.im} ⊆ {s : ℂ | θ < s.re} \ {1} := by
    rintro s ⟨h1, h2⟩
    refine ⟨h1, fun hs => ?_⟩
    rw [mem_singleton_iff] at hs; subst hs; simp at h2
  have sP2 : {s : ℂ | θ < s.re} ∩ {s | s.im < 0} ⊆ {s : ℂ | θ < s.re} \ {1} := by
    rintro s ⟨h1, h2⟩
    refine ⟨h1, fun hs => ?_⟩
    rw [mem_singleton_iff] at hs; subst hs; simp at h2
  have sP3 : {s : ℂ | max θ 1 < s.re} ⊆ {s : ℂ | θ < s.re} \ {1} := by
    intro s (h1 : max θ 1 < s.re)
    refine ⟨show θ < s.re by linarith, fun hs => ?_⟩
    rw [mem_singleton_iff] at hs; subst hs; simp at h1
  have sP4 : {s : ℂ | θ < s.re} ∩ {s | s.re < 1} ⊆ {s : ℂ | θ < s.re} \ {1} := by
    rintro s ⟨h1, h2⟩
    refine ⟨h1, fun hs => ?_⟩
    rw [mem_singleton_iff] at hs; subst hs; simp at h2
  have hc3 : c ∈ {s : ℂ | max θ 1 < s.re} := by
    show max θ 1 < c.re
    simp only [c, ofReal_re]; linarith
  have hcI3 : c + I ∈ {s : ℂ | max θ 1 < s.re} := by
    show max θ 1 < (c + I).re
    simp only [c, add_re, ofReal_re, I_re]; linarith
  have hcI1 : c + I ∈ {s : ℂ | θ < s.re} ∩ {s | 0 < s.im} := by
    refine ⟨show θ < (c + I).re by simp only [c, add_re, ofReal_re, I_re]; linarith,
      show 0 < (c + I).im by simp [c]⟩
  have hcJ3 : c - I ∈ {s : ℂ | max θ 1 < s.re} := by
    show max θ 1 < (c - I).re
    simp only [c, sub_re, ofReal_re, I_re]; linarith
  have hcJ2 : c - I ∈ {s : ℂ | θ < s.re} ∩ {s | s.im < 0} := by
    refine ⟨show θ < (c - I).re by simp only [c, sub_re, ofReal_re, I_re]; linarith,
      show (c - I).im < 0 by simp [c]⟩
  refine isPreconnected_of_forall c fun y hy => ?_
  obtain ⟨hyθ, hy1⟩ := hy
  have hyθ : θ < y.re := hyθ
  by_cases h3 : max θ 1 < y.re
  · exact ⟨_, sP3, hc3, h3, cP3⟩
  by_cases hpos : 0 < y.im
  · exact ⟨_, union_subset sP1 sP3, Or.inr hc3, Or.inl ⟨hyθ, hpos⟩,
      cP1.union (c + I) hcI1 hcI3 cP3⟩
  by_cases hneg : y.im < 0
  · exact ⟨_, union_subset sP2 sP3, Or.inr hc3, Or.inl ⟨hyθ, hneg⟩,
      cP2.union (c - I) hcJ2 hcJ3 cP3⟩
  have him : y.im = 0 := by linarith [not_lt.mp hpos, not_lt.mp hneg]
  have hre : y.re < 1 := by
    by_contra hge
    push Not at hge h3
    have : y.re = 1 := by
      rcases max_cases θ 1 with ⟨hm, _⟩ | ⟨hm, _⟩ <;> rw [hm] at h3 <;> linarith
    exact hy1 (Complex.ext (by simpa using this) (by simpa using him))
  let w : ℂ := ⟨y.re, 1⟩
  have hw4 : w ∈ {s : ℂ | θ < s.re} ∩ {s | s.re < 1} := ⟨hyθ, hre⟩
  have hw1 : w ∈ {s : ℂ | θ < s.re} ∩ {s | 0 < s.im} := ⟨hyθ, show (0 : ℝ) < 1 from one_pos⟩
  exact ⟨_, union_subset (union_subset sP4 sP1) sP3, Or.inr hc3, Or.inl (Or.inl ⟨hyθ, hre⟩),
    (cP4.union w hw4 hw1 cP1).union (c + I) (Or.inr hcI1) hcI3 cP3⟩

theorem mertensZeroFree : Fayed.R2.MertensZeroFree := by
  intro θ h s hs
  by_cases hs1 : s = 1
  · subst hs1; exact riemannZeta_one_ne_zero
  have hUo : IsOpen ({s : ℂ | θ < s.re} \ {1}) :=
    (isOpen_lt continuous_const Complex.continuous_re).sdiff isClosed_singleton
  have hd : DifferentiableOn ℂ (fun z => riemannZeta z * g z) ({s : ℂ | θ < s.re} \ {1}) :=
    fun z hz => ((differentiableAt_riemannZeta hz.2).mul
      (differentiableAt_g h hz.1)).differentiableWithinAt
  let c : ℂ := ((max θ 1 + 1 : ℝ) : ℂ)
  have hc : c ∈ {s : ℂ | θ < s.re} \ {1} := by
    refine ⟨show θ < c.re by simp only [c, ofReal_re]; linarith [le_max_left θ 1], fun hc1 => ?_⟩
    have := congrArg Complex.re (mem_singleton_iff.mp hc1)
    simp only [c, ofReal_re, one_re] at this
    linarith [le_max_right θ 1]
  have hev : (fun z => riemannZeta z * g z) =ᶠ[𝓝 c] (fun _ => 1) := by
    have hmem : {z : ℂ | max θ 1 < z.re} ∈ 𝓝 c :=
      (isOpen_lt continuous_const Complex.continuous_re).mem_nhds
        (show max θ 1 < c.re by simp only [c, ofReal_re]; linarith)
    filter_upwards [hmem] with z hz
    have hz1 : 1 < z.re := lt_of_le_of_lt (le_max_right _ _) hz
    rw [← LSeries_eq_g h hz, ← ArithmeticFunction.LSeries_zeta_eq_riemannZeta hz1]
    exact ArithmeticFunction.LSeries_zeta_mul_Lseries_moebius hz1
  have key := (hd.analyticOnNhd hUo).eqOn_of_preconnected_of_eventuallyEq analyticOnNhd_const
    (isPreconnected_U θ) hc hev ⟨hs, hs1⟩
  intro h0
  simp [h0] at key

end Fayed.Mertens

/-- The R2.3 formalization target, proved. -/
theorem Fayed.R2.mertensZeroFree : Fayed.R2.MertensZeroFree := Fayed.Mertens.mertensZeroFree

#print axioms Fayed.R2.mertensZeroFree
