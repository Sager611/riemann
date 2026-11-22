module

public import Riemann.Deep.LeanLi.ThetaRS
public import Riemann.Deep.LeanLi.ZeroSums

/-!
# LeanLi: Abel summation over zeros against `dN = dA + dS`

For `0 ≤ a ≤ X` and `f ∈ C¹[a, X]`:
`Σ_{a<γ≤X} m f(γ) = f(X)(N(X) - A(X)) - f(a)(N(a) - A(a)) + ∫_a^X f A' - ∫_a^X f' S`.
-/

@[expose] public section

open Real Filter Set MeasureTheory
open scoped Interval

namespace LeanLi

lemma zeroCountAvg_mono : Monotone zeroCountAvg := by
  intro s t h
  unfold zeroCountAvg
  rw [zeroCountBelow_eq, zeroCountBelow_eq]
  have h1 : ∑ ρ ∈ zsetBelow s, (mult ρ : ℝ) ≤ ∑ ρ ∈ zsetBelow t, (mult ρ : ℝ) := by
    apply Finset.sum_le_sum_of_subset_of_nonneg _ fun _ _ _ => Nat.cast_nonneg _
    intro ρ hρ
    rw [mem_zsetBelow] at hρ ⊢
    exact ⟨hρ.1, hρ.2.1, lt_of_lt_of_le hρ.2.2 h⟩
  have h2 : (zeroCount s : ℝ) ≤ zeroCount t := zeroCount_mono h
  linarith

lemma zetaS_eq (t : ℝ) : zetaS t = zeroCountAvg t - Afun t := by
  unfold zetaS Afun; ring

lemma intervalIntegrable_zetaS (a b : ℝ) : IntervalIntegrable zetaS volume a b := by
  have h1 : IntervalIntegrable zeroCountAvg volume a b := zeroCountAvg_mono.intervalIntegrable
  have h2 : IntervalIntegrable Afun volume a b :=
    (continuous_iff_continuousAt.2 fun t => (hasDerivAt_Afun t).continuousAt).intervalIntegrable _ _
  exact (h1.sub h2).congr_uIoo fun t _ => (zetaS_eq t).symm

/-- `N = A + S` almost everywhere on `[a, X]` (off the finitely many ordinates). -/
lemma zeroCount_ae_eq {a X : ℝ} :
    ∀ᵐ t ∂(volume.restrict (Ι a X)), (zeroCount t : ℝ) = Afun t + zetaS t := by
  have hfin : ((fun ρ : ℂ => ρ.im) '' (zerosUpTo (max |a| |X|))).Finite :=
    (finite_zerosUpTo _).image _
  have hnull : volume ((fun ρ : ℂ => ρ.im) '' (zerosUpTo (max |a| |X|))) = 0 := hfin.measure_zero _
  rw [ae_restrict_iff' measurableSet_uIoc]
  have : ∀ᵐ t ∂volume, t ∉ (fun ρ : ℂ => ρ.im) '' (zerosUpTo (max |a| |X|)) :=
    measure_eq_zero_iff_ae_notMem.1 hnull
  filter_upwards [this] with t ht htI
  have hbd : |t| ≤ max |a| |X| := by
    rcases htI with ⟨h1, h2⟩
    rw [abs_le]
    constructor
    · rcases le_total a X with h | h
      · rw [min_eq_left h] at h1; have := neg_abs_le a; have := le_max_left |a| |X|; linarith
      · rw [min_eq_right h] at h1; have := neg_abs_le X; have := le_max_right |a| |X|; linarith
    · rcases le_total a X with h | h
      · rw [max_eq_right h] at h2; have := le_abs_self X; have := le_max_right |a| |X|; linarith
      · rw [max_eq_left h] at h2; have := le_abs_self a; have := le_max_left |a| |X|; linarith
  rw [zetaS_eq]
  rcases le_or_gt t 0 with ht0 | ht0
  · -- no zeros with 0 < γ ≤ t
    have h1 : zset t = ∅ := by
      ext ρ; simp only [mem_zset, zerosUpTo, Set.mem_setOf_eq, Finset.notMem_empty, iff_false]
      rintro ⟨_, h2, h3⟩; linarith
    have h2 : zsetBelow t = ∅ := by
      ext ρ; simp only [mem_zsetBelow, zerosBelow, Set.mem_setOf_eq, Finset.notMem_empty, iff_false]
      rintro ⟨_, h2, h3⟩; linarith
    have hc : (zeroCount t : ℝ) = 0 := by rw [zeroCount_eq, h1, Finset.sum_empty]
    have ha : zeroCountAvg t = 0 := by
      unfold zeroCountAvg; rw [zeroCountBelow_eq, h2, Finset.sum_empty, hc]; ring
    rw [hc, ha]; ring
  · have hno : ∀ ρ ∈ NontrivialZeros, ρ.im ≠ t := by
      intro ρ hρ heq
      exact ht ⟨ρ, ⟨hρ, by rw [heq]; exact ht0, by rw [heq]; exact (le_abs_self t).trans hbd⟩, heq⟩
    rw [avg_eq_of_not_ordinate hno]; ring

/-- **Abel summation against `dA + dS`.** -/
theorem zsum_abel_AS {a X : ℝ} (haX : a ≤ X) {f f' : ℝ → ℝ}
    (hf : ∀ t ∈ Icc a X, HasDerivAt f (f' t) t) (hf' : ContinuousOn f' (Icc a X)) :
    ∑ ρ ∈ (zset X).filter (fun ρ => a < ρ.im), (mult ρ : ℝ) * f ρ.im =
      f X * ((zeroCount X : ℝ) - Afun X) - f a * ((zeroCount a : ℝ) - Afun a) +
        (∫ t in a..X, f t * Ad t) - ∫ t in a..X, f' t * zetaS t := by
  have hui : uIcc a X = Icc a X := uIcc_of_le haX
  have hf'i : IntervalIntegrable f' volume a X := (hui ▸ hf').intervalIntegrable
  rw [zsum_abel haX hf hf'i]
  have hae : ∀ᵐ t ∂volume, t ∈ Ι a X →
      f' t * ((zeroCount t : ℝ) - zeroCount a) = f' t * (Afun t + zetaS t) - (zeroCount a : ℝ) * f' t := by
    have h := (ae_restrict_iff' measurableSet_uIoc).1 (zeroCount_ae_eq (a := a) (X := X))
    filter_upwards [h] with t ht hmem
    rw [ht hmem]; ring
  rw [intervalIntegral.integral_congr_ae hae]
  have hA : ContinuousOn Afun (uIcc a X) := fun t _ => (hasDerivAt_Afun t).continuousAt.continuousWithinAt
  have i1 : IntervalIntegrable (fun t => f' t * Afun t) volume a X := (hf'i.mul_continuousOn hA)
  have i2 : IntervalIntegrable (fun t => f' t * zetaS t) volume a X :=
    (intervalIntegrable_zetaS a X).continuousOn_mul (hui ▸ hf')
  have i3 : IntervalIntegrable (fun t => (zeroCount a : ℝ) * f' t) volume a X := hf'i.const_mul _
  have e1 : ∫ t in a..X, (f' t * (Afun t + zetaS t) - (zeroCount a : ℝ) * f' t) =
      (∫ t in a..X, f' t * Afun t) + (∫ t in a..X, f' t * zetaS t) - (zeroCount a : ℝ) * ∫ t in a..X, f' t := by
    rw [← intervalIntegral.integral_const_mul, ← intervalIntegral.integral_add i1 i2,
      ← intervalIntegral.integral_sub (i1.add i2) i3]
    refine intervalIntegral.integral_congr fun t _ => ?_
    ring
  rw [e1]
  have hparts := intervalIntegral.integral_mul_deriv_eq_deriv_mul (a := a) (b := X)
    (fun t _ => hasDerivAt_Afun t) (fun t ht => hf t (hui ▸ ht))
    (continuous_Ad.intervalIntegrable _ _) hf'i
  have e2 : ∫ t in a..X, f' t * Afun t = Afun X * f X - Afun a * f a - ∫ t in a..X, f t * Ad t := by
    have : ∫ t in a..X, f' t * Afun t = ∫ t in a..X, Afun t * f' t :=
      intervalIntegral.integral_congr fun t _ => by ring
    rw [this, hparts]
    congr 1
    exact intervalIntegral.integral_congr fun t _ => by ring
  have e3 : ∫ t in a..X, f' t = f X - f a :=
    intervalIntegral.integral_eq_sub_of_hasDerivAt (fun t ht => hf t (hui ▸ ht)) hf'i
  rw [e2, e3]
  ring

end LeanLi
