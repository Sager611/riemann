module

public import Mathlib.Basic.Real.Basic
public import Mathlib.Tactic.Linarith
public import Mathlib.Tactic.Ring

@[expose] public section

/-! # Theorem 1′(a): power-free arithmetic of the final combination (RCApps2 l.33, 57–61)

Kept out of the modules that import `Real.rpow` (their import alone is close to the memory cap). -/

namespace RingClass

lemma eps_facts (ε : ℝ) (hε : 0 < ε) (hε1 : ε ≤ 1) : 0 < ε / 10 ∧ ε / 10 ≤ 1 / 10 := by
  constructor <;> linarith

lemma eps5 (e ε : ℝ) : e + 5 * (ε / 10) = e + ε / 2 := by ring

lemma T_ge (t : ℝ) : 1 ≤ |t| + 3 := by
  have := abs_nonneg t
  linarith

lemma one_le_of_three_le (x : ℝ) (h : 3 ≤ x) : 1 ≤ x := by linarith

lemma W_pos (D T : ℝ) (hD : 1 ≤ D) (hT : 1 ≤ T) : 0 < 9 * D * T := by
  have h9 : (0 : ℝ) < 9 * D := by linarith
  exact mul_pos h9 (by linarith)

lemma A_nonneg (c ε : ℝ) (hc : 0 ≤ c) (hε : 0 < ε) : 0 ≤ c * (2 / ε) :=
  mul_nonneg hc (div_nonneg zero_le_two hε.le)

lemma KI_nonneg (CI ε : ℝ) (h : 0 ≤ CI) : 0 ≤ CI * (8 / ε ^ 2) :=
  mul_nonneg h (div_nonneg (by norm_num) (sq_nonneg _))

/-- the case dispatch: from either case bound and P ≤ 81 R, N ≤ 81 (A (K_I + K_II + 1 + |K_S|)) R -/
lemma case_combine (N A KI KII KS P R : ℝ) (hA : 0 ≤ A) (hP : 0 ≤ P) (hKI : 0 ≤ KI)
    (hKII : 0 ≤ KII) (hc : N ≤ A * (KI + KII + 1) * P ∨ N ≤ A * |KS| * P) (hPR : P ≤ 81 * R) :
    N ≤ 81 * (A * (KI + KII + 1 + |KS|)) * R := by
  have hab := abs_nonneg KS
  have hKt : 0 ≤ KI + KII + 1 + |KS| := by linarith
  have hK : N ≤ A * (KI + KII + 1 + |KS|) * P := by
    rcases hc with hc | hc
    · have g : A * (KI + KII + 1) * P ≤ A * (KI + KII + 1 + |KS|) * P :=
        mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left (by linarith) hA) hP
      linarith
    · have g : A * |KS| * P ≤ A * (KI + KII + 1 + |KS|) * P :=
        mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left (by linarith) hA) hP
      linarith
  have g : A * (KI + KII + 1 + |KS|) * P ≤ A * (KI + KII + 1 + |KS|) * (81 * R) :=
    mul_le_mul_of_nonneg_left hPR (mul_nonneg hA hKt)
  have e : A * (KI + KII + 1 + |KS|) * (81 * R) = 81 * (A * (KI + KII + 1 + |KS|)) * R := by ring
  linarith

end RingClass
