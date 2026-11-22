module

public import Riemann.Deep.LeanLi.TauT

/-!
# LeanLi: tails `Σ_{γ > X} 1/γ² ≤ K/X² + 10/√X` and the tail of Li's sum
-/

@[expose] public section

open Real Finset Filter Set MeasureTheory

namespace LeanLi

section
variable (L : LitCore)
include L

/-- `Σ_{a < γ ≤ b} m/γ² ≤ K/a² + 10/√a` uniformly in `b`. -/
lemma sum_inv_sq_tail : ∃ K : ℝ, 0 ≤ K ∧ ∀ a b : ℝ, 0 < a → a ≤ b →
    ∑ ρ ∈ (zset b).filter (fun ρ => a < ρ.im), (mult ρ : ℝ) / ρ.im ^ 2 ≤ K / a ^ 2 + 10 / √a := by
  obtain ⟨K, hK0, hK⟩ := count_le_poly L
  refine ⟨K, hK0, fun a b ha0 hab => ?_⟩
  have hd : ∀ t ∈ Icc a b, HasDerivAt (fun t : ℝ => 1 / t ^ 2) (-2 / t ^ 3) t :=
    fun t ht => hasDerivAt_inv_sq (by linarith [ht.1])
  have hcont : ContinuousOn (fun t : ℝ => -2 / t ^ 3) (uIcc a b) := by
    rw [uIcc_of_le hab]
    exact ContinuousOn.div continuousOn_const (continuousOn_pow 3) fun t ht => by
      have : 0 < t := by linarith [ht.1]
      positivity
  have hint : IntervalIntegrable (fun t : ℝ => -2 / t ^ 3) volume a b := hcont.intervalIntegrable
  have h := zsum_abel hab hd hint
  have hsum : ∑ ρ ∈ (zset b).filter (fun ρ => a < ρ.im), (mult ρ : ℝ) / ρ.im ^ 2 =
      ∑ ρ ∈ (zset b).filter (fun ρ => a < ρ.im), (mult ρ : ℝ) * (1 / ρ.im ^ 2) := by
    refine Finset.sum_congr rfl fun ρ _ => by ring
  rw [hsum, h]
  have hb0 : 0 < b := by linarith
  have hNa0 : (0 : ℝ) ≤ zeroCount a := Nat.cast_nonneg _
  have hcont3 : ContinuousOn (fun t : ℝ => 2 / t ^ 3) (uIcc a b) := by
    rw [uIcc_of_le hab]
    exact ContinuousOn.div continuousOn_const (continuousOn_pow 3) fun t ht => by
      have : 0 < t := by linarith [ht.1]
      positivity
  have hNint : IntervalIntegrable (fun t => (2 / t ^ 3) * (zeroCount t : ℝ)) volume a b :=
    (zeroCount_intervalIntegrable a b).continuousOn_mul hcont3
  have hNint2 : IntervalIntegrable (fun t => -2 / t ^ 3 * ((zeroCount t : ℝ) - zeroCount a)) volume a b :=
    ((zeroCount_intervalIntegrable a b).sub intervalIntegrable_const).continuousOn_mul hcont
  -- -∫ (-2/t³)(N(t) - N(a)) ≤ ∫ (2/t³) N(t)
  have hI1 : -∫ t in a..b, -2 / t ^ 3 * ((zeroCount t : ℝ) - zeroCount a) ≤
      ∫ t in a..b, (2 / t ^ 3) * (zeroCount t : ℝ) := by
    rw [← intervalIntegral.integral_neg]
    apply intervalIntegral.integral_mono_on hab hNint2.neg hNint
    intro t ht
    have ht0 : 0 < t := by linarith [ht.1]
    have hmono : (zeroCount a : ℝ) ≤ zeroCount t := zeroCount_mono ht.1
    show -(-2 / t ^ 3 * ((zeroCount t : ℝ) - zeroCount a)) ≤ 2 / t ^ 3 * (zeroCount t : ℝ)
    have hp : 0 < 2 / t ^ 3 := by positivity
    have e : -(-2 / t ^ 3 * ((zeroCount t : ℝ) - zeroCount a)) =
        2 / t ^ 3 * ((zeroCount t : ℝ) - zeroCount a) := by ring
    rw [e]
    exact mul_le_mul_of_nonneg_left (by linarith) hp.le
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
  have hI2 : ∫ t in a..b, (2 / t ^ 3) * (zeroCount t : ℝ) ≤ (K / a ^ 2 - K / b ^ 2) + (8 / √a - 8 / √b) := by
    rw [← hFTC]
    refine intervalIntegral.integral_mono_on hab hNint ?_ hcmp
    refine ContinuousOn.intervalIntegrable ?_
    rw [uIcc_of_le hab]
    refine ContinuousOn.add (ContinuousOn.div continuousOn_const (continuousOn_pow 3) fun t ht => ?_)
      (ContinuousOn.div continuousOn_const (continuousOn_id.mul (Real.continuous_sqrt.continuousOn)) fun t ht => ?_)
    · have : 0 < t := by linarith [ht.1]
      positivity
    · have : 0 < t := by linarith [ht.1]
      have : 0 < √t := Real.sqrt_pos.2 this
      positivity
  have hNb : 1 / b ^ 2 * ((zeroCount b : ℝ) - zeroCount a) ≤ K / b ^ 2 + 2 / √b := by
    have h1 := hK b hb0.le
    have hsb : 0 < √b := Real.sqrt_pos.2 hb0
    have hbb : √b * √b = b := Real.mul_self_sqrt hb0.le
    calc 1 / b ^ 2 * ((zeroCount b : ℝ) - zeroCount a) ≤ 1 / b ^ 2 * (K + 2 * b * √b) :=
          mul_le_mul_of_nonneg_left (by linarith) (by positivity)
      _ = K / b ^ 2 + 2 / √b := by field_simp; nlinarith [hbb]
  have hsab : 2 / √b ≤ 2 / √a := by
    have hsa : 0 < √a := Real.sqrt_pos.2 ha0
    exact div_le_div_of_nonneg_left (by norm_num) hsa (Real.sqrt_le_sqrt hab)
  have h8b : 0 ≤ 8 / √b := by positivity
  have e10 : 10 / √a = 8 / √a + 2 / √a := by ring
  linarith

lemma summable_inv_ordK_sq : Summable (fun k : ZeroIdx => 1 / ordK k ^ 2) := by
  have hnn : ∀ k : ZeroIdx, 0 ≤ 1 / ordK k ^ 2 := fun k => by positivity
  rw [summable_sigma_of_nonneg hnn]
  refine ⟨fun ρ => (hasSum_fintype _).summable, ?_⟩
  refine (summable_mult_div_sq L).congr fun ρ => ?_
  rw [tsum_fintype]
  simp [ordK, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul, div_eq_mul_inv]

/-- **Tail over the zero index**: ‵Σ_{γ_k > X} 1/γ_k² ≤ K/X² + 10/√X‵. -/
lemma tsum_tail_le : ∃ K : ℝ, 0 ≤ K ∧ ∀ X : ℝ, 0 < X →
    ∑' k : ↑((ZK X : Set ZeroIdx)ᶜ), 1 / ordK k ^ 2 ≤ K / X ^ 2 + 10 / √X := by
  obtain ⟨K, hK0, hK⟩ := sum_inv_sq_tail L
  refine ⟨K, hK0, fun X hX => ?_⟩
  have hsum := (summable_inv_ordK_sq L).subtype ((ZK X : Set ZeroIdx)ᶜ)
  refine hsum.tsum_le_of_sum_le fun s => ?_
  set t : Finset ZeroIdx := s.map (Function.Embedding.subtype _) with ht
  set b : ℝ := X + ∑ k ∈ t, |ordK k| with hb
  have hXb : X ≤ b := by
    have : 0 ≤ ∑ k ∈ t, |ordK k| := Finset.sum_nonneg fun _ _ => abs_nonneg _
    linarith
  have hsub : t ⊆ (ZK b).filter (fun k => X < ordK k) := by
    intro k hk
    rw [Finset.mem_filter, mem_ZK]
    rw [ht, Finset.mem_map] at hk
    obtain ⟨k', _, rfl⟩ := hk
    have hk1 : (k' : ZeroIdx) ∉ ZK X := k'.2
    rw [mem_ZK, not_le] at hk1
    refine ⟨?_, hk1⟩
    have h1 : |ordK (k' : ZeroIdx)| ≤ ∑ k ∈ t, |ordK k| :=
      Finset.single_le_sum (f := fun k => |ordK k|) (fun _ _ => abs_nonneg _)
        (by rw [ht]; exact Finset.mem_map_of_mem _ (by assumption))
    have h2 : 0 ≤ X := hX.le
    show ordK (k' : ZeroIdx) ≤ b
    linarith [le_abs_self (ordK (k' : ZeroIdx))]
  have e1 : ∑ i ∈ s, ((fun k => 1 / ordK k ^ 2) ∘ Subtype.val) i = ∑ k ∈ t, 1 / ordK k ^ 2 := by
    rw [ht, Finset.sum_map]; rfl
  rw [e1]
  have e2 : ∑ k ∈ (ZK b).filter (fun k => X < ordK k), 1 / ordK k ^ 2 =
      ∑ ρ ∈ (zset b).filter (fun ρ => X < ρ.im), (mult ρ : ℝ) / ρ.im ^ 2 := by
    rw [Finset.sum_filter, Finset.sum_filter]
    have := sum_ZK b (fun γ => if X < γ then 1 / γ ^ 2 else 0)
    rw [this]
    refine Finset.sum_congr rfl fun ρ _ => ?_
    split_ifs <;> simp [div_eq_mul_inv]
  calc ∑ k ∈ t, 1 / ordK k ^ 2 ≤ ∑ k ∈ (ZK b).filter (fun k => X < ordK k), 1 / ordK k ^ 2 :=
        Finset.sum_le_sum_of_subset_of_nonneg hsub fun k _ _ => by positivity
    _ = _ := e2
    _ ≤ K / X ^ 2 + 10 / √X := hK X b hX hXb

end

end LeanLi
