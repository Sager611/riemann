module

public import Mathlib.MeasureTheory.Integral.IntervalIntegral.FundThmCalculus

/-!
# LeanLi: Abel (Stieltjes) summation for finite families of points

`Σ_i m_i f(x_i) = f(b) Σ_i m_i - ∫_a^b f'(t) M(t) dt`, `M(t) = Σ_{x_i ≤ t} m_i`, all `x_i ∈ (a, b]`.
-/

@[expose] public section

open MeasureTheory Set Filter Topology

namespace LeanLi

/-- The counting function of a weighted finite family: `Σ_{i ∈ s, x_i ≤ t} m_i`. -/
noncomputable def cnt {ι : Type*} (s : Finset ι) (m : ι → ℝ) (x : ι → ℝ) (t : ℝ) : ℝ :=
  ∑ i ∈ s, if x i ≤ t then m i else 0

lemma integral_mul_ite_lt {g : ℝ → ℝ} {a b c : ℝ} (hac : a ≤ c) (hcb : c ≤ b) :
    ∫ t in a..b, g t * (if c < t then (1 : ℝ) else 0) = ∫ t in c..b, g t := by
  rw [intervalIntegral.integral_of_le (hac.trans hcb), intervalIntegral.integral_of_le hcb]
  have h1 : (fun t => g t * (if c < t then (1 : ℝ) else 0)) = (Ioi c).indicator g := by
    funext t
    by_cases h : c < t
    · rw [if_pos h, mul_one, Set.indicator_of_mem (show t ∈ Ioi c from h)]
    · rw [if_neg h, mul_zero, Set.indicator_of_notMem (show t ∉ Ioi c from h)]
  rw [h1, integral_indicator measurableSet_Ioi, Measure.restrict_restrict measurableSet_Ioi]
  have h2 : Ioi c ∩ Ioc a b = Ioc c b := by
    ext t; simp only [mem_inter_iff, mem_Ioi, mem_Ioc]
    constructor
    · rintro ⟨h1, _, h3⟩; exact ⟨h1, h3⟩
    · rintro ⟨h1, h2⟩; exact ⟨h1, by linarith, h2⟩
  rw [h2]

lemma integral_mul_ite_le {g : ℝ → ℝ} {a b c : ℝ} (hac : a ≤ c) (hcb : c ≤ b) :
    ∫ t in a..b, g t * (if c ≤ t then (1 : ℝ) else 0) = ∫ t in c..b, g t := by
  rw [← integral_mul_ite_lt (g := g) hac hcb]
  refine intervalIntegral.integral_congr_ae ?_
  filter_upwards [(Set.countable_singleton c).ae_notMem volume] with t ht _
  have htc : t ≠ c := fun h => ht (by simp [h])
  by_cases h : c < t
  · rw [if_pos h, if_pos h.le]
  · rw [if_neg h, if_neg (fun h' => h (lt_of_le_of_ne h' (Ne.symm htc)))]

lemma intervalIntegrable_mul_ite {g : ℝ → ℝ} {a b c : ℝ} (hab : a ≤ b)
    (hg : IntervalIntegrable g volume a b) :
    IntervalIntegrable (fun t => g t * (if c ≤ t then (1 : ℝ) else 0)) volume a b := by
  have h1 : (fun t => g t * (if c ≤ t then (1 : ℝ) else 0)) = (Ici c).indicator g := by
    funext t
    by_cases h : c ≤ t
    · rw [if_pos h, mul_one, Set.indicator_of_mem (show t ∈ Ici c from h)]
    · rw [if_neg h, mul_zero, Set.indicator_of_notMem (show t ∉ Ici c from h)]
  rw [h1]
  rw [intervalIntegrable_iff_integrableOn_Ioc_of_le hab] at hg ⊢
  exact hg.indicator measurableSet_Ici

/-- **Abel summation** over a finite weighted family of points in `(a, b]`. -/
theorem abel_sum {ι : Type*} (s : Finset ι) (m x : ι → ℝ) {a b : ℝ}
    (hx : ∀ i ∈ s, a < x i ∧ x i ≤ b) {f f' : ℝ → ℝ}
    (hf : ∀ t ∈ Icc a b, HasDerivAt f (f' t) t) (hf' : IntervalIntegrable f' volume a b) :
    ∑ i ∈ s, m i * f (x i) = f b * ∑ i ∈ s, m i - ∫ t in a..b, f' t * cnt s m x t := by
  by_cases hs : s = ∅
  · subst hs; simp [cnt]
  obtain ⟨i₀, hi₀⟩ := Finset.nonempty_iff_ne_empty.2 hs
  have hab : a ≤ b := ((hx i₀ hi₀).1.trans_le (hx i₀ hi₀).2).le
  have hpt : ∀ i ∈ s, f (x i) = f b - ∫ t in a..b, f' t * (if x i ≤ t then (1 : ℝ) else 0) := by
    intro i hi
    obtain ⟨h1, h2⟩ := hx i hi
    rw [integral_mul_ite_le h1.le h2]
    have hd : ∀ t ∈ uIcc (x i) b, HasDerivAt f (f' t) t := by
      intro t ht
      rw [uIcc_of_le h2] at ht
      exact hf t ⟨by linarith [ht.1], ht.2⟩
    have hint : IntervalIntegrable f' volume (x i) b :=
      hf'.mono_set (by rw [uIcc_of_le h2, uIcc_of_le hab]; exact Icc_subset_Icc h1.le le_rfl)
    rw [intervalIntegral.integral_eq_sub_of_hasDerivAt hd hint]
    ring
  have hcnt : (fun t => f' t * cnt s m x t) =
      fun t => ∑ i ∈ s, m i * (f' t * (if x i ≤ t then (1 : ℝ) else 0)) := by
    funext t
    unfold cnt
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun i _ => ?_
    by_cases h : x i ≤ t <;> simp [h]; ring
  rw [hcnt, intervalIntegral.integral_finsetSum]
  · rw [Finset.sum_congr rfl (fun i hi => by rw [hpt i hi]), Finset.mul_sum, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [intervalIntegral.integral_const_mul]
    ring
  · intro i _
    exact (intervalIntegrable_mul_ite hab hf').const_mul _

end LeanLi

