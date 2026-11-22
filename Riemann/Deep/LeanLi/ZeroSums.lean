module

public import Riemann.Deep.LeanLi.CountBounds

/-!
# LeanLi: sums over the zeros (Abel summation against `N(t)`), `Σ_{γ>0} m/γ² < ∞`
-/

@[expose] public section

open Complex Filter Topology Set ComplexConjugate MeasureTheory

namespace LeanLi

lemma zeroCount_intervalIntegrable (a b : ℝ) :
    IntervalIntegrable (fun t => (zeroCount t : ℝ)) volume a b :=
  zeroCount_mono.intervalIntegrable

/-- **Abel summation over the zeros**: for `0 ≤ a ≤ b`,
`Σ_{a<γ≤b} m f(γ) = f(b)(N(b) - N(a)) - ∫_a^b f'(t)(N(t) - N(a)) dt`. -/
theorem zsum_abel {a b : ℝ} (hab : a ≤ b) {f f' : ℝ → ℝ}
    (hf : ∀ t ∈ Icc a b, HasDerivAt f (f' t) t) (hf' : IntervalIntegrable f' volume a b) :
    ∑ ρ ∈ (zset b).filter (fun ρ => a < ρ.im), (mult ρ : ℝ) * f ρ.im =
      f b * ((zeroCount b : ℝ) - zeroCount a) -
        ∫ t in a..b, f' t * ((zeroCount t : ℝ) - zeroCount a) := by
  have h := abel_sum ((zset b).filter (fun ρ => a < ρ.im)) (fun ρ => (mult ρ : ℝ)) (fun ρ => ρ.im)
    (a := a) (b := b) (fun ρ hρ => by
      rw [Finset.mem_filter, mem_zset] at hρ
      exact ⟨hρ.2, hρ.1.2.2⟩) hf hf'
  have hsplit : ∀ t, a ≤ t → t ≤ b →
      cnt ((zset b).filter (fun ρ => a < ρ.im)) (fun ρ => (mult ρ : ℝ)) (fun ρ => ρ.im) t =
        (zeroCount t : ℝ) - zeroCount a := by
    intro t hat htb
    have e1 := cnt_zset (b := b) (t := t) htb
    have e2 := cnt_zset (b := b) (t := a) hab
    unfold cnt at e1 e2 ⊢
    rw [Finset.sum_filter, ← e1, ← e2, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun ρ _ => ?_
    by_cases h1 : a < ρ.im
    · have : ¬ ρ.im ≤ a := not_le.2 h1
      by_cases h2 : ρ.im ≤ t <;> simp [h1, h2, this]
    · have h1' : ρ.im ≤ a := not_lt.1 h1
      have h2 : ρ.im ≤ t := h1'.trans hat
      simp [h1, h1', h2]
  have hsum : ∑ ρ ∈ (zset b).filter (fun ρ => a < ρ.im), (mult ρ : ℝ) = (zeroCount b : ℝ) - zeroCount a := by
    have := hsplit b hab le_rfl
    unfold cnt at this
    rw [← this]
    refine Finset.sum_congr rfl fun ρ hρ => ?_
    rw [Finset.mem_filter, mem_zset] at hρ
    simp [hρ.1.2.2]
  rw [h, hsum]
  congr 1
  refine intervalIntegral.integral_congr fun t ht => ?_
  rw [uIcc_of_le hab] at ht
  rw [hsplit t ht.1 ht.2]

/-- `N(t) ≤ K + 2 t √t` for all `t ≥ 0`. -/
lemma count_le_poly (L : LitCore) : ∃ K : ℝ, 0 ≤ K ∧ ∀ t, 0 ≤ t → (zeroCount t : ℝ) ≤ K + 2 * t * √t := by
  obtain ⟨T₀, hT₀, h⟩ := count_le L
  refine ⟨zeroCount T₀, Nat.cast_nonneg _, fun t ht => ?_⟩
  by_cases htT : t ≤ T₀
  · have := zeroCount_mono htT
    have : 0 ≤ 2 * t * √t := by positivity
    linarith
  · push Not at htT
    have h1 := h t htT.le
    have hlog : Real.log t ≤ 2 * √t := by
      have ht1 : 0 < √t := Real.sqrt_pos.2 (by linarith)
      have := Real.log_le_sub_one_of_pos ht1
      rw [Real.log_sqrt (by linarith)] at this
      linarith
    have : t * Real.log t ≤ t * (2 * √t) := mul_le_mul_of_nonneg_left hlog (by linarith)
    have : (0 : ℝ) ≤ zeroCount T₀ := Nat.cast_nonneg _
    nlinarith

lemma hasDerivAt_inv_sq {t : ℝ} (ht : t ≠ 0) : HasDerivAt (fun t : ℝ => 1 / t ^ 2) (-2 / t ^ 3) t := by
  have := (hasDerivAt_pow 2 t).inv (pow_ne_zero 2 ht)
  convert this using 1
  · funext x; simp [one_div]
  · field_simp; ring

/-- The positive ordinates: `Σ_{0<γ≤b} m/γ² ≤ B` uniformly in `b`. -/
lemma sum_inv_sq_le (L : LitCore) : ∃ B : ℝ, ∀ b : ℝ,
    ∑ ρ ∈ zset b, (mult ρ : ℝ) / ρ.im ^ 2 ≤ B := by
  obtain ⟨c, hc0, hc1, hc⟩ := exists_im_lower
  obtain ⟨K, hK0, hK⟩ := count_le_poly L
  set a := c / 2 with ha
  have ha0 : 0 < a := by positivity
  refine ⟨2 * K / a ^ 2 + 10 / √a, fun b => ?_⟩
  by_cases hb : b < c
  · have : zset b = ∅ := by
      ext ρ; simp only [mem_zset, Finset.notMem_empty, iff_false]
      rintro ⟨hρ, hpos, hle⟩
      have := hc ρ hρ hpos; linarith
    rw [this, Finset.sum_empty]; positivity
  push Not at hb
  have hab : a ≤ b := by linarith
  have hNa : (zeroCount a : ℝ) = 0 := by
    rw [zeroCount_eq]
    have : zset a = ∅ := by
      ext ρ; simp only [mem_zset, Finset.notMem_empty, iff_false]
      rintro ⟨hρ, hpos, hle⟩
      have := hc ρ hρ hpos; linarith
    rw [this, Finset.sum_empty]
  have hfilt : (zset b).filter (fun ρ => a < ρ.im) = zset b := by
    ext ρ; simp only [Finset.mem_filter, mem_zset]
    constructor
    · exact fun h => h.1
    · intro h; exact ⟨h, by have := hc ρ h.1 h.2.1; linarith⟩
  have hd : ∀ t ∈ Icc a b, HasDerivAt (fun t : ℝ => 1 / t ^ 2) (-2 / t ^ 3) t :=
    fun t ht => hasDerivAt_inv_sq (by linarith [ht.1])
  have hcont : ContinuousOn (fun t : ℝ => -2 / t ^ 3) (uIcc a b) := by
    rw [uIcc_of_le hab]
    exact ContinuousOn.div continuousOn_const (continuousOn_pow 3) fun t ht => by
      have : 0 < t := by linarith [ht.1]
      positivity
  have hint : IntervalIntegrable (fun t : ℝ => -2 / t ^ 3) volume a b := hcont.intervalIntegrable
  have h := zsum_abel hab hd hint
  rw [hfilt, hNa] at h
  simp only [sub_zero] at h
  have hsum : ∑ ρ ∈ zset b, (mult ρ : ℝ) / ρ.im ^ 2 = ∑ ρ ∈ zset b, (mult ρ : ℝ) * (1 / ρ.im ^ 2) := by
    refine Finset.sum_congr rfl fun ρ _ => by ring
  rw [hsum, h]
  -- bound the integral
  have hNint : IntervalIntegrable (fun t => (2 / t ^ 3) * (zeroCount t : ℝ)) volume a b := by
    refine (zeroCount_intervalIntegrable a b).continuousOn_mul ?_
    rw [uIcc_of_le hab]
    exact ContinuousOn.div continuousOn_const (continuousOn_pow 3) fun t ht => by
      have : 0 < t := by linarith [ht.1]
      positivity
  have hcmp : ∀ t ∈ Icc a b, (2 / t ^ 3) * (zeroCount t : ℝ) ≤ 2 * K / t ^ 3 + 4 / (t * √t) := by
    intro t ht
    have ht0 : 0 < t := by linarith [ht.1]
    have h1 := hK t ht0.le
    have hs : 0 < √t := Real.sqrt_pos.2 ht0
    have hss : √t * √t = t := Real.mul_self_sqrt ht0.le
    rw [div_mul_eq_mul_div, div_le_iff₀ (by positivity)]
    have : 4 / (t * √t) * t ^ 3 = 4 * t * √t := by
      field_simp
      nlinarith [hss]
    have e : (2 * K / t ^ 3 + 4 / (t * √t)) * t ^ 3 = 2 * K + 4 * t * √t := by
      rw [add_mul, this]; field_simp
    rw [e]; nlinarith
  have hFTC : ∫ t in a..b, (2 * K / t ^ 3 + 4 / (t * √t)) = (K / a ^ 2 - K / b ^ 2) + (8 / √a - 8 / √b) := by
    have hd2 : ∀ t ∈ uIcc a b, HasDerivAt (fun t : ℝ => -K / t ^ 2 - 8 / √t) (2 * K / t ^ 3 + 4 / (t * √t)) t := by
      intro t ht
      rw [uIcc_of_le hab] at ht
      have ht0 : 0 < t := by linarith [ht.1]
      have hs : 0 < √t := Real.sqrt_pos.2 ht0
      have h1 := (hasDerivAt_inv_sq ht0.ne').const_mul (-K)
      have h2 := ((Real.hasDerivAt_sqrt ht0.ne').inv hs.ne').const_mul (-8)
      have := h1.add h2
      convert this using 1
      · funext x; simp [one_div]; ring
      · have hss : √t * √t = t := Real.mul_self_sqrt ht0.le
        field_simp
        rw [show t ^ 3 = t * t * t by ring]
        nlinarith [hss]
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hd2]
    · ring
    · refine ContinuousOn.intervalIntegrable ?_
      rw [uIcc_of_le hab]
      refine ContinuousOn.add (ContinuousOn.div continuousOn_const (continuousOn_pow 3) fun t ht => ?_)
        (ContinuousOn.div continuousOn_const (continuousOn_id.mul (Real.continuous_sqrt.continuousOn)) fun t ht => ?_)
      · have : 0 < t := by linarith [ht.1]
        positivity
      · have : 0 < t := by linarith [ht.1]
        have : 0 < √t := Real.sqrt_pos.2 this
        positivity
  have hI : -∫ t in a..b, (-2 / t ^ 3) * (zeroCount t : ℝ) ≤ (K / a ^ 2 - K / b ^ 2) + (8 / √a - 8 / √b) := by
    rw [← hFTC, ← intervalIntegral.integral_neg]
    refine intervalIntegral.integral_mono_on hab ?_ ?_ ?_
    · have e : (fun t => -(-2 / t ^ 3 * (zeroCount t : ℝ))) = fun t => 2 / t ^ 3 * (zeroCount t : ℝ) := by
        funext t; ring
      rw [e]; exact hNint
    · refine ContinuousOn.intervalIntegrable ?_
      rw [uIcc_of_le hab]
      refine ContinuousOn.add (ContinuousOn.div continuousOn_const (continuousOn_pow 3) fun t ht => ?_)
        (ContinuousOn.div continuousOn_const (continuousOn_id.mul (Real.continuous_sqrt.continuousOn)) fun t ht => ?_)
      · have : 0 < t := by linarith [ht.1]
        positivity
      · have : 0 < t := by linarith [ht.1]
        have : 0 < √t := Real.sqrt_pos.2 this
        positivity
    · intro t ht
      have := hcmp t ht
      have e : -(-2 / t ^ 3 * (zeroCount t : ℝ)) = 2 / t ^ 3 * (zeroCount t : ℝ) := by ring
      rw [e]; exact this
  have hb0 : 0 < b := by linarith
  have hNb : 1 / b ^ 2 * (zeroCount b : ℝ) ≤ K / a ^ 2 + 2 / √a := by
    have h1 := hK b hb0.le
    have hsb : 0 < √b := Real.sqrt_pos.2 hb0
    have hsa : 0 < √a := Real.sqrt_pos.2 ha0
    have hsab : √a ≤ √b := Real.sqrt_le_sqrt hab
    have hbb : √b * √b = b := Real.mul_self_sqrt hb0.le
    calc 1 / b ^ 2 * (zeroCount b : ℝ) ≤ 1 / b ^ 2 * (K + 2 * b * √b) :=
          mul_le_mul_of_nonneg_left h1 (by positivity)
      _ = K / b ^ 2 + 2 / √b := by field_simp; nlinarith [hbb]
      _ ≤ K / a ^ 2 + 2 / √a := by
          gcongr
  have hKb : 0 ≤ K / b ^ 2 := by positivity
  have h8b : 0 ≤ 8 / √b := by positivity
  set I := ∫ t in a..b, -2 / t ^ 3 * (zeroCount t : ℝ) with hIdef
  have h10 : 2 * K / a ^ 2 + 10 / √a = (K / a ^ 2 + 2 / √a) + (K / a ^ 2 + 8 / √a) := by ring
  rw [h10]
  have hI' : -I ≤ K / a ^ 2 + 8 / √a := by linarith
  linarith

/-- The zeros in the upper half plane. -/
def UpperZeros : Set ℂ := {ρ | ρ ∈ NontrivialZeros ∧ 0 < ρ.im}

lemma summable_mult_div_sq (L : LitCore) :
    Summable (fun ρ : UpperZeros => (mult (ρ : ℂ) : ℝ) / (ρ : ℂ).im ^ 2) := by
  obtain ⟨B, hB⟩ := sum_inv_sq_le L
  refine summable_of_sum_le (c := B) (fun ρ => by positivity) fun u => ?_
  set v : Finset ℂ := u.map (Function.Embedding.subtype _) with hv
  set b : ℝ := ∑ z ∈ v, |z.im| with hb
  have hsub : v ⊆ zset b := by
    intro z hz
    rw [mem_zset]
    rw [hv, Finset.mem_map] at hz
    obtain ⟨ρ, _, rfl⟩ := hz
    refine ⟨ρ.2.1, ρ.2.2, ?_⟩
    have : |(ρ : ℂ).im| ≤ b := by
      rw [hb]
      exact Finset.single_le_sum (f := fun z : ℂ => |z.im|) (fun _ _ => abs_nonneg _)
        (Finset.mem_map_of_mem _ (by assumption))
    exact (le_abs_self _).trans this
  calc ∑ ρ ∈ u, (mult (ρ : ℂ) : ℝ) / (ρ : ℂ).im ^ 2 = ∑ z ∈ v, (mult z : ℝ) / z.im ^ 2 := by
        rw [hv, Finset.sum_map]; rfl
    _ ≤ ∑ z ∈ zset b, (mult z : ℝ) / z.im ^ 2 :=
        Finset.sum_le_sum_of_subset_of_nonneg hsub (fun _ _ _ => by positivity)
    _ ≤ B := hB b

end LeanLi
