import Riemann.Deep.LeanGale.Gale
import Riemann.Deep.LeanGale.Splx
import Riemann.Deep.LeanGale.BallAlg
import Mathlib.LinearAlgebra.FiniteDimensional.Lemmas

/-!
# LeanGale T6 (2/3): the regular-simplex partition of the ball

* `IsGram U`: `U : Fin (n + 1) → ℝⁿ` has the Gram matrix of a regular simplex with unit circumradius centred at 0,
  ⟪U j, U k⟫ = 1 (j = k), −1/n (j ≠ k); `isGram_splxE`: the simplex `splx n` of Splx.lean.
* `frame`: then Σⱼ U j = 0, the first n vectors are linearly independent, hence span ℝⁿ, and
  x = n/(n+1) Σⱼ ⟪U j, x⟫ U j; so ⟪x, y⟫ = n/(n+1) Σⱼ ⟪U j, x⟫⟪U j, y⟫ (`inner_eq_frame`).
* `cell U x`: an index maximising ⟪U j, x⟫ (the colouring; its classes are the cones over the facets of −U).
* `cone_ineq`: same cell i ⇒ ⟪x, y⟫ + κₙ/2 (‖x‖² + ‖y‖²) ≥ 0 (BallAlg `psi_nonneg` with aⱼ = ⟪U i, x⟫ − ⟪U j, x⟫,
  bⱼ = ⟪U i, y⟫ − ⟪U j, y⟫, common zero at i); then ‖x − y‖² ≤ (1 + κₙ)(‖x‖² + ‖y‖²) ≤ (1 + κₙ)/2 = σₙ² for
  ‖x‖, ‖y‖ ≤ 1/2 (`dist_le_sigma`).
* `beta_closedBall_le_sigma`: β_{n+1}(closedBall 0 (1/2)) ≤ σₙ for every n ≥ 1.
-/

namespace LeanGale

open Gale7 Finset
open scoped InnerProductSpace

/-- σₙ = √((1 + κₙ)/2) -/
theorem sigma_eq_sqrt_kappa (n : ℕ) : sigma n = Real.sqrt ((1 + kappa n) / 2) := by
  have h2 : (n : ℝ) + 2 ≠ 0 := by positivity
  unfold sigma kappa
  split_ifs with h
  · congr 1
    field_simp
    ring
  · congr 1
    ring

theorem sum_mul_gram {ι : Type*} [Fintype ι] [DecidableEq ι] (m : ℝ) (j : ι) (c : ι → ℝ) :
    ∑ k, c k * (if j = k then (1 : ℝ) else -1 / m) = c j - (∑ k, c k - c j) / m := by
  have h : ∀ k, c k * (if j = k then (1 : ℝ) else -1 / m) =
      (if j = k then c k + c k / m else 0) - c k / m := by
    intro k
    split_ifs <;> ring
  rw [Finset.sum_congr rfl fun k _ => h k, Finset.sum_sub_distrib, ← Finset.sum_div]
  simp only [Finset.sum_ite_eq, Finset.mem_univ, ite_true]
  ring

section Frame

variable {n : ℕ}

/-- `U` has the Gram matrix of a regular simplex with unit circumradius: ⟪U j, U k⟫ = 1 (j = k), −1/n (j ≠ k) -/
def IsGram (U : Fin (n + 1) → E n) : Prop :=
  ∀ j k, ⟪U j, U k⟫_ℝ = if j = k then (1 : ℝ) else -1 / n

theorem sum_eq_zero_of_isGram {U : Fin (n + 1) → E n} (hn : n ≠ 0) (hU : IsGram U) : ∑ j, U j = 0 := by
  have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast hn
  rw [← inner_self_eq_zero (𝕜 := ℝ), sum_inner]
  have h : ∀ j : Fin (n + 1), ⟪U j, ∑ k, U k⟫_ℝ = 0 := by
    intro j
    rw [inner_sum]
    have := sum_mul_gram (n : ℝ) j (fun _ => (1 : ℝ))
    simp only [one_mul, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul,
      mul_one] at this
    simp only [hU j]
    rw [this]
    push_cast
    rw [add_sub_cancel_right, div_self hn', sub_self]
  simp [h]

theorem inner_sum_eq_zero_of_isGram {U : Fin (n + 1) → E n} (hn : n ≠ 0) (hU : IsGram U) (x : E n) :
    ∑ j, ⟪U j, x⟫_ℝ = 0 := by
  rw [← sum_inner, sum_eq_zero_of_isGram hn hU, inner_zero_left]

theorem linearIndependent_of_isGram {U : Fin (n + 1) → E n} (hn : n ≠ 0) (hU : IsGram U) :
    LinearIndependent ℝ (fun i : Fin n => U (Fin.castSucc i)) := by
  have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast hn
  rw [Fintype.linearIndependent_iff]
  intro g hg
  have hk : ∀ k : Fin n, ((n : ℝ) + 1) * g k = ∑ i, g i := by
    intro k
    have h := congrArg (fun v => ⟪U (Fin.castSucc k), v⟫_ℝ) hg
    simp only [inner_sum, inner_smul_right, inner_zero_right, hU (Fin.castSucc k),
      Fin.castSucc_inj] at h
    rw [sum_mul_gram (n : ℝ) k g] at h
    field_simp at h
    linarith
  have hsum : ∑ i, g i = 0 := by
    have h := Finset.sum_congr (s₁ := (univ : Finset (Fin n))) rfl fun k _ => hk k
    rw [← Finset.mul_sum, Finset.sum_const, Finset.card_univ, Fintype.card_fin, nsmul_eq_mul] at h
    linarith
  intro i
  have h := hk i
  rw [hsum] at h
  have hpos : (0 : ℝ) < (n : ℝ) + 1 := by positivity
  exact (mul_eq_zero.mp h).resolve_left hpos.ne'

theorem inner_residual {U : Fin (n + 1) → E n} (hn : n ≠ 0) (hU : IsGram U) (x : E n) (k : Fin (n + 1)) :
    ⟪U k, x - ((n : ℝ) / (n + 1)) • ∑ j, ⟪U j, x⟫_ℝ • U j⟫_ℝ = 0 := by
  have hn' : (n : ℝ) ≠ 0 := by exact_mod_cast hn
  have hn1 : (n : ℝ) + 1 ≠ 0 := by positivity
  rw [inner_sub_right, inner_smul_right, inner_sum]
  simp only [inner_smul_right, hU k]
  rw [sum_mul_gram (n : ℝ) k (fun j => ⟪U j, x⟫_ℝ), inner_sum_eq_zero_of_isGram hn hU x]
  field_simp
  ring

/-- **the frame identity**: x = n/(n+1) Σⱼ ⟪U j, x⟫ U j -/
theorem frame {U : Fin (n + 1) → E n} (hn : n ≠ 0) (hU : IsGram U) (x : E n) :
    x = ((n : ℝ) / (n + 1)) • ∑ j, ⟪U j, x⟫_ℝ • U j := by
  have horth : ∀ k, ⟪U k, x - ((n : ℝ) / (n + 1)) • ∑ j, ⟪U j, x⟫_ℝ • U j⟫_ℝ = 0 :=
    inner_residual hn hU x
  have hspan := (linearIndependent_of_isGram hn hU).span_eq_top_of_card_eq_finrank' (by simp)
  have hmem : x - ((n : ℝ) / (n + 1)) • ∑ j, ⟪U j, x⟫_ℝ • U j ∈
      Submodule.span ℝ (Set.range fun i : Fin n => U (Fin.castSucc i)) := by
    rw [hspan]; exact Submodule.mem_top
  obtain ⟨c, hc⟩ := (Submodule.mem_span_range_iff_exists_fun ℝ).mp hmem
  have hrr : ⟪x - ((n : ℝ) / (n + 1)) • ∑ j, ⟪U j, x⟫_ℝ • U j,
      x - ((n : ℝ) / (n + 1)) • ∑ j, ⟪U j, x⟫_ℝ • U j⟫_ℝ = 0 := by
    rw [← hc, sum_inner]
    exact Finset.sum_eq_zero fun i _ => by rw [real_inner_smul_left, hc, horth, mul_zero]
  exact sub_eq_zero.mp (inner_self_eq_zero.mp hrr)

/-- ⟪x, y⟫ = n/(n+1) Σⱼ ⟪U j, x⟫⟪U j, y⟫ -/
theorem inner_eq_frame {U : Fin (n + 1) → E n} (hn : n ≠ 0) (hU : IsGram U) (x y : E n) :
    ⟪x, y⟫_ℝ = (n : ℝ) / (n + 1) * ∑ j, ⟪U j, x⟫_ℝ * ⟪U j, y⟫_ℝ := by
  conv_lhs => rw [frame hn hU x]
  rw [real_inner_smul_left, sum_inner]
  simp only [real_inner_smul_left]

/-- the colouring: an index maximising ⟪U j, x⟫ -/
noncomputable def cell (U : Fin (n + 1) → E n) (x : E n) : Fin (n + 1) :=
  Classical.choose (Finset.exists_max_image Finset.univ (fun j => ⟪U j, x⟫_ℝ) Finset.univ_nonempty)

theorem le_cell (U : Fin (n + 1) → E n) (x : E n) (j : Fin (n + 1)) : ⟪U j, x⟫_ℝ ≤ ⟪U (cell U x), x⟫_ℝ :=
  (Classical.choose_spec (Finset.exists_max_image Finset.univ (fun j => ⟪U j, x⟫_ℝ)
    Finset.univ_nonempty)).2 j (Finset.mem_univ j)

/-- **the cone inequality**: x, y in the cell i ⇒ ⟪x, y⟫ + κₙ/2 (‖x‖² + ‖y‖²) ≥ 0 -/
theorem cone_ineq {U : Fin (n + 1) → E n} (hn : n ≠ 0) (hU : IsGram U) (i : Fin (n + 1)) (x y : E n)
    (hx : ∀ j, ⟪U j, x⟫_ℝ ≤ ⟪U i, x⟫_ℝ) (hy : ∀ j, ⟪U j, y⟫_ℝ ≤ ⟪U i, y⟫_ℝ) :
    0 ≤ ⟪x, y⟫_ℝ + kappa n / 2 * (‖x‖ ^ 2 + ‖y‖ ^ 2) := by
  have hint : ∀ a b : ℕ, a + b + 1 ≤ Fintype.card (Fin (n + 1)) →
      (a : ℝ) * b ≤ kappa n ^ 2 * ((Fintype.card (Fin (n + 1)) : ℝ) - a) *
        ((Fintype.card (Fin (n + 1)) : ℝ) - b) := by
    intro a b hab
    rw [Fintype.card_fin] at hab ⊢
    push_cast
    exact kappa_int n a b (by omega)
  have hpsi := psi_nonneg (kappa n) (kappa_nonneg n) (kappa_le_one n) hint i
    (fun j => ⟪U i, x⟫_ℝ - ⟪U j, x⟫_ℝ) (fun j => ⟪U i, y⟫_ℝ - ⟪U j, y⟫_ℝ)
    (fun j => sub_nonneg.mpr (hx j)) (fun j => sub_nonneg.mpr (hy j)) (sub_self _) (sub_self _)
  rw [psi_shift (kappa n) (fun j => ⟪U j, x⟫_ℝ) (fun j => ⟪U j, y⟫_ℝ) i
    (inner_sum_eq_zero_of_isGram hn hU x) (inner_sum_eq_zero_of_isGram hn hU y), Fintype.card_fin] at hpsi
  have hin := (mul_nonneg_iff_of_pos_left (by positivity)).mp hpsi
  rw [← real_inner_self_eq_norm_sq, ← real_inner_self_eq_norm_sq, inner_eq_frame hn hU x y,
    inner_eq_frame hn hU x x, inner_eq_frame hn hU y y]
  simp only [← sq]
  have hq : 0 ≤ (n : ℝ) / (n + 1) := by positivity
  nlinarith [mul_nonneg hq hin]

/-- same cell, ‖x‖, ‖y‖ ≤ 1/2 ⇒ dist x y ≤ σₙ -/
theorem dist_le_sigma {U : Fin (n + 1) → E n} (hn : n ≠ 0) (hU : IsGram U) (i : Fin (n + 1)) {x y : E n}
    (hxB : ‖x‖ ≤ 1 / 2) (hyB : ‖y‖ ≤ 1 / 2) (hx : ∀ j, ⟪U j, x⟫_ℝ ≤ ⟪U i, x⟫_ℝ)
    (hy : ∀ j, ⟪U j, y⟫_ℝ ≤ ⟪U i, y⟫_ℝ) : dist x y ≤ sigma n := by
  have hc := cone_ineq hn hU i x y hx hy
  have hκ0 := kappa_nonneg n
  have hx2 : ‖x‖ ^ 2 ≤ 1 / 4 := by nlinarith [norm_nonneg x]
  have hy2 : ‖y‖ ^ 2 ≤ 1 / 4 := by nlinarith [norm_nonneg y]
  have h : ‖x - y‖ ^ 2 ≤ (1 + kappa n) / 2 := by
    rw [norm_sub_sq_real]
    have := mul_le_mul_of_nonneg_left (add_le_add hx2 hy2) (by linarith : (0 : ℝ) ≤ 1 + kappa n)
    nlinarith
  rw [dist_eq_norm, sigma_eq_sqrt_kappa]
  calc ‖x - y‖ = Real.sqrt (‖x - y‖ ^ 2) := (Real.sqrt_sq (norm_nonneg _)).symm
    _ ≤ Real.sqrt ((1 + kappa n) / 2) := Real.sqrt_le_sqrt h

/-- the cells of `U` cut the ball of diameter 1 into n + 1 parts of diameter ≤ σₙ -/
theorem beta_ball_le_sigma {U : Fin (n + 1) → E n} (hn : n ≠ 0) (hU : IsGram U) :
    beta (n + 1) (Metric.closedBall (0 : E n) (1 / 2)) ≤ sigma n := by
  refine (ciInf_le (bddBelow_parts _) (cell U)).trans (ciSup_le fun i => ?_)
  refine Metric.diam_le_of_forall_dist_le (sigma_nonneg n) ?_
  rintro x ⟨hxB, hxi⟩ y ⟨hyB, hyi⟩
  rw [Set.mem_preimage, Set.mem_singleton_iff] at hxi hyi
  rw [Metric.mem_closedBall, dist_zero_right] at hxB hyB
  refine dist_le_sigma hn hU i hxB hyB (fun j => ?_) (fun j => ?_)
  · rw [← hxi]; exact le_cell U x j
  · rw [← hyi]; exact le_cell U y j

end Frame

/-- the regular simplex `splx n` (Splx.lean) as points of ℝⁿ -/
noncomputable def splxE (n : ℕ) : Fin (n + 1) → E n := fun j => toE (splx n j)

theorem isGram_splxE {n : ℕ} (hn : 1 ≤ n) : IsGram (splxE n) := by
  intro j k
  rw [← splx_simplex hn j k]
  simp only [splxE, toE, EuclideanSpace.inner_toLp_toLp, dotProduct, Pi.star_apply, star_trivial]
  exact Finset.sum_congr rfl fun t _ => mul_comm _ _

/-- **β_{n+1}(B) ≤ σₙ** for the closed ball B of diameter 1 in ℝⁿ, every n ≥ 1 -/
theorem beta_closedBall_le_sigma {n : ℕ} (hn : 1 ≤ n) :
    beta (n + 1) (Metric.closedBall (0 : E n) (1 / 2)) ≤ sigma n :=
  beta_ball_le_sigma (by omega) (isGram_splxE hn)

end LeanGale
