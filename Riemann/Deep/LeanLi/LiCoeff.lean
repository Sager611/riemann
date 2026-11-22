module

public import Riemann.Deep.LeanLi.ZeroSums
public import Riemann.Deep.LeanLi.Statement

/-!
# LeanLi: convergence and reality of `λ_n` (Part 0); Li's formula under RH
-/

@[expose] public section

open Complex Filter Topology Set ComplexConjugate MeasureTheory

namespace LeanLi

/-- `1 - wⁿ = n(1 - w) - (1 - w)² Q_n(w)`, `Q_n = Σ_{k<n} Σ_{j<k} wʲ`. -/
lemma one_sub_pow_eq (w : ℂ) (n : ℕ) :
    1 - w ^ n = n * (1 - w) - (1 - w) ^ 2 * ∑ k ∈ Finset.range n, ∑ j ∈ Finset.range k, w ^ j := by
  induction n with
  | zero => simp
  | succ n ih =>
    have hg := geom_sum_mul w n
    rw [Finset.sum_range_succ, pow_succ]
    have : 1 - w ^ n * w = (1 - w ^ n) + w ^ n * (1 - w) := by ring
    rw [this, ih]
    push_cast
    linear_combination (-(1 - w)) * hg

lemma norm_geom_le (w : ℂ) (k : ℕ) :
    ‖∑ j ∈ Finset.range k, w ^ j‖ ≤ k * (max 1 ‖w‖) ^ k := by
  calc ‖∑ j ∈ Finset.range k, w ^ j‖ ≤ ∑ j ∈ Finset.range k, ‖w‖ ^ j := by
        refine (norm_sum_le _ _).trans (le_of_eq ?_)
        simp [norm_pow]
    _ ≤ ∑ j ∈ Finset.range k, (max 1 ‖w‖) ^ k := by
        refine Finset.sum_le_sum fun j hj => ?_
        calc ‖w‖ ^ j ≤ (max 1 ‖w‖) ^ j := pow_le_pow_left₀ (norm_nonneg _) (le_max_right _ _) _
          _ ≤ (max 1 ‖w‖) ^ k := pow_le_pow_right₀ (le_max_left _ _) (Finset.mem_range.1 hj).le
    _ = k * (max 1 ‖w‖) ^ k := by simp

lemma norm_one_sub_pow_sub_le (w : ℂ) (n : ℕ) :
    ‖(1 - w ^ n) - n * (1 - w)‖ ≤ n ^ 2 * ‖1 - w‖ ^ 2 * (max 1 ‖w‖) ^ n := by
  rw [one_sub_pow_eq, show (n : ℂ) * (1 - w) - (1 - w) ^ 2 * _ - n * (1 - w) =
    -((1 - w) ^ 2 * ∑ k ∈ Finset.range n, ∑ j ∈ Finset.range k, w ^ j) by ring, norm_neg, norm_mul,
    norm_pow]
  have hM : 1 ≤ max 1 ‖w‖ := le_max_left _ _
  have hQ : ‖∑ k ∈ Finset.range n, ∑ j ∈ Finset.range k, w ^ j‖ ≤ n ^ 2 * (max 1 ‖w‖) ^ n := by
    calc ‖∑ k ∈ Finset.range n, ∑ j ∈ Finset.range k, w ^ j‖
        ≤ ∑ k ∈ Finset.range n, ‖∑ j ∈ Finset.range k, w ^ j‖ := norm_sum_le _ _
      _ ≤ ∑ k ∈ Finset.range n, (n : ℝ) * (max 1 ‖w‖) ^ n := by
          refine Finset.sum_le_sum fun k hk => (norm_geom_le w k).trans ?_
          have hk' : k ≤ n := (Finset.mem_range.1 hk).le
          exact mul_le_mul (by exact_mod_cast hk') (pow_le_pow_right₀ hM hk') (by positivity)
            (Nat.cast_nonneg _)
      _ = n ^ 2 * (max 1 ‖w‖) ^ n := by simp; ring
  calc ‖1 - w‖ ^ 2 * ‖∑ k ∈ Finset.range n, ∑ j ∈ Finset.range k, w ^ j‖
      ≤ ‖1 - w‖ ^ 2 * (n ^ 2 * (max 1 ‖w‖) ^ n) := mul_le_mul_of_nonneg_left hQ (by positivity)
    _ = n ^ 2 * ‖1 - w‖ ^ 2 * (max 1 ‖w‖) ^ n := by ring

/-- The paired term of an upper zero `ρ`: `m_ρ([1 - (1-1/ρ)ⁿ] + [1 - (1-1/ρ̄)ⁿ]) = 2 m_ρ Re(1 - (1-1/ρ)ⁿ)`. -/
noncomputable def upperTerm (n : ℕ) (ρ : ℂ) : ℝ := 2 * mult ρ * (1 - (1 - 1 / ρ) ^ n).re

lemma abs_upperTerm_le {n : ℕ} {c : ℝ} (hc : 0 < c) {ρ : ℂ} (hρ : ρ ∈ NontrivialZeros)
    (hcρ : c ≤ ρ.im) :
    |upperTerm n ρ| ≤ 2 * (n + n ^ 2 * (1 + 1 / c) ^ n) * ((mult ρ : ℝ) / ρ.im ^ 2) := by
  have hre0 := re_pos_of_mem hρ
  have hre1 := re_lt_one_of_mem hρ
  have him : 0 < ρ.im := hc.trans_le hcρ
  have hρ0 : ρ ≠ 0 := fun h => by rw [h] at him; simp at him
  have hnorm : ρ.im ≤ ‖ρ‖ := (le_abs_self _).trans (Complex.abs_im_le_norm ρ)
  have hnpos : 0 < ‖ρ‖ := norm_pos_iff.2 hρ0
  set w : ℂ := 1 - 1 / ρ with hw
  have h1w : 1 - w = 1 / ρ := by rw [hw]; ring
  have hre_inv : (1 / ρ).re = ρ.re / ‖ρ‖ ^ 2 := by
    rw [one_div, Complex.inv_re, Complex.normSq_eq_norm_sq]
  have hnorm_inv : ‖1 / ρ‖ = 1 / ‖ρ‖ := by rw [norm_div, norm_one]
  have hwle : max 1 ‖w‖ ≤ 1 + 1 / c := by
    refine max_le (by have : 0 < 1 / c := by positivity
                      linarith) ?_
    calc ‖w‖ ≤ ‖(1 : ℂ)‖ + ‖1 / ρ‖ := norm_sub_le _ _
      _ = 1 + 1 / ‖ρ‖ := by rw [norm_one, hnorm_inv]
      _ ≤ 1 + 1 / c := by gcongr; exact hcρ.trans hnorm
  have hmain := norm_one_sub_pow_sub_le w n
  rw [h1w, hnorm_inv] at hmain
  have hre : |(1 - w ^ n).re| ≤ n * (1 / ‖ρ‖ ^ 2) + n ^ 2 * (1 / ‖ρ‖) ^ 2 * (1 + 1 / c) ^ n := by
    have h2 : |((1 - w ^ n) - n * (1 / ρ)).re| ≤ n ^ 2 * (1 / ‖ρ‖) ^ 2 * (1 + 1 / c) ^ n := by
      refine (Complex.abs_re_le_norm _).trans (hmain.trans ?_)
      gcongr
    have h3 : |((n : ℂ) * (1 / ρ)).re| ≤ n * (1 / ‖ρ‖ ^ 2) := by
      have : ((n : ℂ) * (1 / ρ)).re = n * (ρ.re / ‖ρ‖ ^ 2) := by
        rw [← hre_inv]; simp [Complex.mul_re]
      rw [this, abs_mul, Nat.abs_cast]
      gcongr
      rw [abs_of_pos (by positivity)]
      exact div_le_div_of_nonneg_right hre1.le (by positivity)
    have := abs_add_le ((1 - w ^ n) - n * (1 / ρ)).re ((n : ℂ) * (1 / ρ)).re
    rw [← Complex.add_re, sub_add_cancel] at this
    linarith
  have hγ : 1 / ‖ρ‖ ^ 2 ≤ 1 / ρ.im ^ 2 := by
    gcongr
  unfold upperTerm
  rw [← hw, abs_mul, abs_of_nonneg (by positivity : (0 : ℝ) ≤ 2 * mult ρ)]
  calc 2 * (mult ρ : ℝ) * |(1 - w ^ n).re|
      ≤ 2 * (mult ρ : ℝ) * (n * (1 / ρ.im ^ 2) + n ^ 2 * (1 / ρ.im ^ 2) * (1 + 1 / c) ^ n) := by
        gcongr
        refine hre.trans ?_
        have : (1 / ‖ρ‖) ^ 2 = 1 / ‖ρ‖ ^ 2 := by ring
        rw [this]
        gcongr
    _ = 2 * (n + n ^ 2 * (1 + 1 / c) ^ n) * ((mult ρ : ℝ) / ρ.im ^ 2) := by ring

lemma summable_upperTerm (L : LitCore) (n : ℕ) :
    Summable (fun ρ : UpperZeros => upperTerm n ρ) := by
  obtain ⟨c, hc0, _, hc⟩ := exists_im_lower
  refine Summable.of_norm_bounded ((summable_mult_div_sq L).mul_left (2 * (n + n ^ 2 * (1 + 1 / c) ^ n)))
    fun ρ => ?_
  rw [Real.norm_eq_abs]
  exact abs_upperTerm_le hc0 ρ.2.1 (hc _ ρ.2.1 ρ.2.2)

/-- The upper zeros with `|ρ| ≤ T`, as a finset of the subtype. -/
noncomputable def upperFin (T : ℝ) : Finset UpperZeros :=
  ((finite_norm_le T).preimage (f := fun ρ : UpperZeros => (ρ : ℂ)) Subtype.val_injective.injOn).toFinset

lemma mem_upperFin {T : ℝ} {ρ : UpperZeros} : ρ ∈ upperFin T ↔ ‖(ρ : ℂ)‖ ≤ T := by
  simp [upperFin, ρ.2.1]

/-- Pairing `ρ ↔ ρ̄`: the partial Li sum is real and equals the sum of the paired terms. -/
lemma liPartial_eq (L : LitCore) (n : ℕ) (T : ℝ) :
    liPartial n T = ((∑ ρ ∈ upperFin T, upperTerm n ρ : ℝ) : ℂ) := by
  set F : ℂ → ℂ := fun ρ => (mult ρ : ℂ) * (1 - (1 - 1 / ρ) ^ n) with hF
  set S := (finite_norm_le T).toFinset with hS
  have h0 : liPartial n T = ∑ ρ ∈ S, F ρ := by
    unfold liPartial
    rw [finsum_mem_eq_finite_toFinset_sum _ (finite_norm_le T)]
  have hsplit := Finset.sum_filter_add_sum_filter_not S (fun ρ : ℂ => 0 < ρ.im) F
  set Sp := S.filter (fun ρ : ℂ => 0 < ρ.im) with hSp
  have hneg : S.filter (fun ρ : ℂ => ¬ 0 < ρ.im) = Sp.image conj := by
    ext z
    simp only [hSp, hS, Finset.mem_filter, Finset.mem_image, Set.Finite.mem_toFinset, mem_setOf_eq]
    constructor
    · rintro ⟨⟨hz, hzT⟩, hzim⟩
      refine ⟨conj z, ⟨⟨conj_mem hz, by rwa [Complex.norm_conj]⟩, ?_⟩, conj_conj z⟩
      have := im_ne_zero L.zeta_unit_interval hz
      simp only [Complex.conj_im]
      push Not at hzim
      exact neg_pos.2 (lt_of_le_of_ne hzim this)
    · rintro ⟨y, ⟨⟨hy, hyT⟩, hyim⟩, rfl⟩
      refine ⟨⟨conj_mem hy, by rwa [Complex.norm_conj]⟩, ?_⟩
      simp only [Complex.conj_im, not_lt]
      linarith
  have hinj : Set.InjOn conj (Sp : Set ℂ) := fun x _ y _ h => by simpa using congrArg conj h
  rw [hneg, Finset.sum_image hinj] at hsplit
  have hconj : ∀ ρ ∈ Sp, F (conj ρ) = conj (F ρ) := by
    intro ρ hρ
    simp only [hSp, hS, Finset.mem_filter, Set.Finite.mem_toFinset, mem_setOf_eq] at hρ
    simp only [hF, map_mul, map_sub, map_one, map_pow, map_div₀, map_natCast, mult_conj hρ.1.1.2.2]
  rw [h0, ← hsplit, ← Finset.sum_add_distrib, Finset.sum_congr rfl (fun ρ hρ => by rw [hconj ρ hρ])]
  have hre : ∀ ρ, F ρ + conj (F ρ) = ((2 * (F ρ).re : ℝ) : ℂ) := by
    intro ρ; rw [Complex.add_conj]
  simp only [hre]
  push_cast
  -- reindex Sp as upperFin T
  have hmap : Sp = (upperFin T).map (Function.Embedding.subtype _) := by
    ext z
    simp only [hSp, hS, Finset.mem_filter, Set.Finite.mem_toFinset, mem_setOf_eq, Finset.mem_map,
      Function.Embedding.coe_subtype]
    constructor
    · rintro ⟨⟨hz, hzT⟩, hzim⟩
      exact ⟨⟨z, hz, hzim⟩, mem_upperFin.2 hzT, rfl⟩
    · rintro ⟨ρ, hρ, rfl⟩
      exact ⟨⟨ρ.2.1, mem_upperFin.1 hρ⟩, ρ.2.2⟩
  rw [hmap, Finset.sum_map]
  refine Finset.sum_congr rfl fun ρ _ => ?_
  simp only [hF, Function.Embedding.coe_subtype, upperTerm]
  rw [show ((mult (ρ : ℂ) : ℕ) : ℂ) = ((mult (ρ : ℂ) : ℝ) : ℂ) from (Complex.ofReal_natCast _).symm,
    Complex.re_ofReal_mul]
  push_cast; ring

lemma tendsto_upperFin : Tendsto upperFin atTop atTop := by
  refine tendsto_atTop_finset_of_monotone (fun S T h ρ hρ => ?_) fun ρ => ⟨‖(ρ : ℂ)‖, ?_⟩
  · rw [mem_upperFin] at hρ ⊢; linarith
  · exact mem_upperFin.2 le_rfl

/-- **Part 0**: the partial Li sums converge, to a real number: `λ_n = Σ_{γ>0} 2 m Re(1-(1-1/ρ)ⁿ)`. -/
theorem tendsto_liPartial (L : LitCore) (n : ℕ) :
    Tendsto (liPartial n) atTop (𝓝 ((∑' ρ : UpperZeros, upperTerm n ρ : ℝ) : ℂ)) := by
  have h := ((summable_upperTerm L n).hasSum.comp tendsto_upperFin)
  have h2 := (Complex.continuous_ofReal.tendsto _).comp h
  refine h2.congr' (Eventually.of_forall fun T => ?_)
  simp only [Function.comp]
  rw [liPartial_eq L n T]

theorem liCoeff_eq (L : LitCore) (n : ℕ) :
    liCoeff n = ((∑' ρ : UpperZeros, upperTerm n ρ : ℝ) : ℂ) :=
  (tendsto_liPartial L n).limUnder_eq

theorem lam_eq (L : LitCore) (n : ℕ) : lam n = ∑' ρ : UpperZeros, upperTerm n ρ := by
  unfold lam; rw [liCoeff_eq L n, Complex.ofReal_re]

theorem part0 (L : LitCore) : Part0 := by
  intro n
  refine ⟨?_, ?_⟩
  · rw [liCoeff_eq L n]; exact tendsto_liPartial L n
  · rw [liCoeff_eq L n, Complex.ofReal_im]

end LeanLi
