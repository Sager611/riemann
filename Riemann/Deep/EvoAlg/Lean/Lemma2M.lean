import Riemann.Deep.EvoAlg.Lean.Lemma2
import Mathlib.Data.Matrix.Mul

/-! # Lemma 2 (NOTE §4) with Mathlib matrices -/

namespace EvoAlg

open Matrix

variable {𝕜 : Type*} [Field 𝕜]

/-- **Lemma 2** (NOTE §4), matrix form: `P ∈ M_k(𝕜)`, `PᵀP = αI`, `P𝟙 = γ𝟙` ⟹ `k(γ² − α) = 0`. -/
theorem lemma2_matrix {ι : Type*} [Fintype ι] [DecidableEq ι] (P : Matrix ι ι 𝕜) (α γ : 𝕜)
    (h1 : Pᵀ * P = α • (1 : Matrix ι ι 𝕜)) (h2 : P *ᵥ (fun _ => (1 : 𝕜)) = fun _ => γ) :
    (Fintype.card ι : 𝕜) * (γ ^ 2 - α) = 0 := by
  refine lemma2 (fun w a => P w a) α γ (fun a b => ?_) (fun w => ?_)
  · have := congrFun (congrFun h1 a) b
    simp only [Matrix.mul_apply, Matrix.transpose_apply, Matrix.smul_apply, Matrix.one_apply,
      smul_eq_mul, mul_ite, mul_one, mul_zero] at this
    exact this
  · have := congrFun h2 w
    simp only [Matrix.mulVec, dotProduct, mul_one] at this
    exact this

end EvoAlg

