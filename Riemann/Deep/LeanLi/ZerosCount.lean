module

public import Riemann.Deep.LeanLi.ZerosBasic
public import Riemann.Deep.LeanLi.Abel

/-!
# LeanLi: counting the zeros (finite sets, `N(T)`, Titchmarsh's averaged `N(T)`)
-/

@[expose] public section

open Complex Filter Topology Set ComplexConjugate

namespace LeanLi

/-- The finite set of zeros with `0 < γ ≤ T`. -/
noncomputable def zset (T : ℝ) : Finset ℂ := (finite_zerosUpTo T).toFinset

/-- The finite set of zeros with `0 < γ < T`. -/
noncomputable def zsetBelow (T : ℝ) : Finset ℂ := (finite_zerosBelow T).toFinset

lemma mem_zset {T : ℝ} {ρ : ℂ} : ρ ∈ zset T ↔ ρ ∈ zerosUpTo T := by
  simp [zset]

lemma mem_zsetBelow {T : ℝ} {ρ : ℂ} : ρ ∈ zsetBelow T ↔ ρ ∈ zerosBelow T := by
  simp [zsetBelow]

lemma zeroCount_eq (T : ℝ) : (zeroCount T : ℝ) = ∑ ρ ∈ zset T, (mult ρ : ℝ) := by
  unfold zeroCount zset
  rw [finsum_mem_eq_finite_toFinset_sum _ (finite_zerosUpTo T)]
  push_cast; rfl

lemma zeroCountBelow_eq (T : ℝ) :
    ((∑ᶠ ρ ∈ zerosBelow T, mult ρ : ℕ) : ℝ) = ∑ ρ ∈ zsetBelow T, (mult ρ : ℝ) := by
  unfold zsetBelow
  rw [finsum_mem_eq_finite_toFinset_sum _ (finite_zerosBelow T)]
  push_cast; rfl

lemma zset_mono {S T : ℝ} (h : S ≤ T) : zset S ⊆ zset T := by
  intro ρ hρ
  rw [mem_zset] at hρ ⊢
  exact ⟨hρ.1, hρ.2.1, hρ.2.2.trans h⟩

lemma zeroCount_mono : Monotone (fun T => (zeroCount T : ℝ)) := by
  intro S T h
  simp only [zeroCount_eq]
  exact Finset.sum_le_sum_of_subset_of_nonneg (zset_mono h) (fun _ _ _ => Nat.cast_nonneg _)

lemma zsetBelow_subset (T : ℝ) : zsetBelow T ⊆ zset T := by
  intro ρ hρ
  rw [mem_zsetBelow] at hρ; rw [mem_zset]
  exact ⟨hρ.1, hρ.2.1, hρ.2.2.le⟩

lemma zeroCount_le_two_avg (T : ℝ) : (zeroCount T : ℝ) ≤ 2 * zeroCountAvg T := by
  unfold zeroCountAvg
  rw [zeroCountBelow_eq]
  have : 0 ≤ ∑ ρ ∈ zsetBelow T, (mult ρ : ℝ) := Finset.sum_nonneg fun _ _ => Nat.cast_nonneg _
  linarith

lemma avg_le_zeroCount (T : ℝ) : zeroCountAvg T ≤ zeroCount T := by
  unfold zeroCountAvg
  rw [zeroCountBelow_eq, zeroCount_eq]
  have := Finset.sum_le_sum_of_subset_of_nonneg (zsetBelow_subset T)
    (f := fun ρ => (mult ρ : ℝ)) (fun _ _ _ => Nat.cast_nonneg _)
  linarith

/-- The counting function of `zset b` is `N(t)` for `t ≤ b`. -/
lemma cnt_zset {b t : ℝ} (ht : t ≤ b) :
    cnt (zset b) (fun ρ => (mult ρ : ℝ)) (fun ρ => ρ.im) t = zeroCount t := by
  unfold cnt
  rw [zeroCount_eq, ← Finset.sum_filter]
  congr 1
  ext ρ
  simp only [Finset.mem_filter, mem_zset, zerosUpTo, mem_setOf_eq]
  constructor
  · rintro ⟨⟨h1, h2, _⟩, h4⟩; exact ⟨h1, h2, h4⟩
  · rintro ⟨h1, h2, h3⟩; exact ⟨⟨h1, h2, h3.trans ht⟩, h3⟩

/-- A positive lower bound for the positive ordinates. -/
lemma exists_im_lower : ∃ c : ℝ, 0 < c ∧ c ≤ 1 ∧
    ∀ ρ : ℂ, ρ ∈ NontrivialZeros → 0 < ρ.im → c ≤ ρ.im := by
  by_cases h : (zset 1).Nonempty
  · obtain ⟨ρ₀, hρ₀, hmin⟩ := (zset 1).exists_min_image (fun ρ => ρ.im) h
    have hρ₀' := (mem_zset.1 hρ₀)
    refine ⟨ρ₀.im, hρ₀'.2.1, hρ₀'.2.2, fun ρ hρ hpos => ?_⟩
    by_cases h1 : ρ.im ≤ 1
    · exact hmin ρ (mem_zset.2 ⟨hρ, hpos, h1⟩)
    · linarith [hρ₀'.2.2]
  · refine ⟨1, one_pos, le_rfl, fun ρ hρ hpos => ?_⟩
    by_contra h1
    push Not at h1
    exact h ⟨ρ, mem_zset.2 ⟨hρ, hpos, h1.le⟩⟩

/-- In every interval of positive reals there is a point that is not an ordinate. -/
lemma exists_not_ordinate {a b : ℝ} (ha : 0 ≤ a) (hab : a < b) :
    ∃ t ∈ Ioo a b, ∀ ρ ∈ NontrivialZeros, ρ.im ≠ t := by
  have hfin : ((fun ρ : ℂ => ρ.im) '' (zerosUpTo b)).Finite := (finite_zerosUpTo b).image _
  obtain ⟨t, ht, htn⟩ := ((Set.Ioo_infinite hab).diff hfin).nonempty
  refine ⟨t, ht, fun ρ hρ heq => htn ?_⟩
  exact ⟨ρ, ⟨hρ, by rw [heq]; linarith [ht.1], heq ▸ ht.2.le⟩, heq⟩

/-- At a non-ordinate, Titchmarsh's averaged count is the count. -/
lemma avg_eq_of_not_ordinate {t : ℝ} (ht : ∀ ρ ∈ NontrivialZeros, ρ.im ≠ t) :
    zeroCountAvg t = zeroCount t := by
  unfold zeroCountAvg
  have : zsetBelow t = zset t := by
    ext ρ; rw [mem_zsetBelow, mem_zset]
    constructor
    · rintro ⟨h1, h2, h3⟩; exact ⟨h1, h2, h3.le⟩
    · rintro ⟨h1, h2, h3⟩; exact ⟨h1, h2, lt_of_le_of_ne h3 (ht ρ h1)⟩
  rw [zeroCountBelow_eq, this, ← zeroCount_eq]; ring

end LeanLi
