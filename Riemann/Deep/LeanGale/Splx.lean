import Riemann.Deep.Gale7.Lean.Fam
import Mathlib.Analysis.Real.Sqrt

/-!
# LeanGale: a regular simplex of unit vectors in every dimension

`splx m`: m + 1 unit vectors in ℝ^m with pairwise inner products −1/m (the vertices of a regular simplex with
circumradius 1, centred at 0), built recursively: in ℝ^(m+1) the last vertex is (0, …, 0, −1) and the others
are (α u_j, 1/(m+1)) with u_j the vertices in ℝ^m and α = √(m(m+2))/(m+1).
`splx_simplex`: `IsUnitSimplex (splx m)` (Gale7's definition) for every m ≥ 1.
-/

namespace LeanGale

open Gale7

/-- unit vertex vectors of a regular simplex in ℝ^m (junk for m = 0) -/
noncomputable def splx : (m : ℕ) → Fin (m + 1) → Fin m → ℝ
  | 0 => fun _ => Fin.elim0
  | m + 1 => fun i => Fin.lastCases (motive := fun _ => Fin (m + 1) → ℝ) (Fin.snoc 0 (-1))
      (fun j => Fin.snoc ((Real.sqrt ((m : ℝ) * (m + 2)) / (m + 1)) • splx m j) (1 / ((m : ℝ) + 1))) i

theorem splx_last (m : ℕ) : splx (m + 1) (Fin.last (m + 1)) = Fin.snoc 0 (-1) := by
  simp only [splx, Fin.lastCases_last]

theorem splx_castSucc (m : ℕ) (j : Fin (m + 1)) :
    splx (m + 1) (Fin.castSucc j) =
      Fin.snoc ((Real.sqrt ((m : ℝ) * (m + 2)) / (m + 1)) • splx m j) (1 / ((m : ℝ) + 1)) := by
  simp only [splx, Fin.lastCases_castSucc]

theorem dot_snoc {n : ℕ} (x y : Fin n → ℝ) (a b : ℝ) :
    ∑ t, (Fin.snoc x a : Fin (n + 1) → ℝ) t * (Fin.snoc y b : Fin (n + 1) → ℝ) t = ∑ t, x t * y t + a * b := by
  rw [Fin.sum_univ_castSucc]
  simp only [Fin.snoc_castSucc, Fin.snoc_last]

theorem dot_smul {n : ℕ} (α : ℝ) (x y : Fin n → ℝ) :
    ∑ t, (α • x) t * (α • y) t = α ^ 2 * ∑ t, x t * y t := by
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun t _ => ?_
  simp only [Pi.smul_apply, smul_eq_mul]
  ring

/-- the Gram matrix of `splx (m + 1)` -/
theorem splx_simplex_succ : ∀ m : ℕ, IsUnitSimplex (splx (m + 1)) := by
  intro m
  induction m with
  | zero =>
      intro i k
      rcases Fin.eq_castSucc_or_eq_last i with ⟨j, rfl⟩ | rfl <;>
        rcases Fin.eq_castSucc_or_eq_last k with ⟨l, rfl⟩ | rfl
      · rw [splx_castSucc, splx_castSucc, dot_snoc, Fin.fin_one_eq_zero j, Fin.fin_one_eq_zero l]
        simp
      · rw [splx_castSucc, splx_last, dot_snoc, if_neg (Fin.castSucc_ne_last j)]
        simp
      · rw [splx_castSucc, splx_last, dot_snoc, if_neg (Fin.castSucc_ne_last l).symm]
        simp
      · rw [splx_last, dot_snoc, if_pos rfl]
        simp
  | succ m ih =>
      have hα : Real.sqrt (((m + 1 : ℕ) : ℝ) * ((m + 1 : ℕ) + 2)) ^ 2 = ((m + 1 : ℕ) : ℝ) * ((m + 1 : ℕ) + 2) :=
        Real.sq_sqrt (by positivity)
      have hm : ((m + 1 : ℕ) : ℝ) ≠ 0 := by positivity
      have hm1 : ((m + 1 : ℕ) : ℝ) + 1 ≠ 0 := by positivity
      intro i k
      rcases Fin.eq_castSucc_or_eq_last i with ⟨j, rfl⟩ | rfl <;>
        rcases Fin.eq_castSucc_or_eq_last k with ⟨l, rfl⟩ | rfl
      · rw [splx_castSucc, splx_castSucc, dot_snoc, dot_smul, ih j l, div_pow, hα]
        by_cases hjl : j = l
        · subst hjl
          rw [if_pos rfl, if_pos rfl]
          field_simp
          push_cast
          ring
        · rw [if_neg hjl, if_neg (fun h => hjl (Fin.castSucc_injective _ h))]
          field_simp
          push_cast
          ring
      · rw [splx_castSucc, splx_last, dot_snoc, if_neg (Fin.castSucc_ne_last j)]
        simp only [Pi.zero_apply, mul_zero, Finset.sum_const_zero, zero_add]
        push_cast
        field_simp
      · rw [splx_castSucc, splx_last, dot_snoc, if_neg (Fin.castSucc_ne_last l).symm]
        simp only [Pi.zero_apply, zero_mul, Finset.sum_const_zero, zero_add]
        push_cast
        field_simp
      · rw [splx_last, dot_snoc, if_pos rfl]
        simp

/-- **Regular simplices exist in every dimension m ≥ 1.** -/
theorem splx_simplex {m : ℕ} (hm : 1 ≤ m) : IsUnitSimplex (splx m) := by
  obtain ⟨k, rfl⟩ : ∃ k, m = k + 1 := ⟨m - 1, by omega⟩
  exact splx_simplex_succ k

end LeanGale
