module

public import Mathlib.Basic.Real.Basic
public import Mathlib.Tactic.Linarith
public import Mathlib.Tactic.Ring

@[expose] public section

/-! # Theorem 1′(a): final arithmetic of the large-𝔔 case (RCApps2 l.57), power-free:
E = W^{ε₀}, U₁ = W^{e+4ε₀}, U₂ = W^{e+ε₀}, P = W^{e+5ε₀} enter only through EU₁ = P, EU₂ ≤ P, E ≤ P. -/

namespace RingClass

/-- N ≤ B(x+1), x ≤ B_I + B_II, B ≤ c(2/ε₀)E, B_I ≤ K_I U₁, B_II ≤ K_II U₂ ⟹ N ≤ c(2/ε₀)(K_I+K_II+1)P -/
lemma large_arith (N B x BI BII c ε₀ KI KII E U1 U2 P : ℝ) (hε₀ : 0 < ε₀)
    (hc : 0 ≤ c) (hKI : 0 ≤ KI) (hKII : 0 ≤ KII) (hU1 : 0 ≤ U1) (hU2 : 0 ≤ U2)
    (f1 : E * U1 = P) (f2 : E * U2 ≤ P) (f3 : E ≤ P)
    (hN : N ≤ B * (x + 1)) (hx : x ≤ BI + BII) (hB0 : 0 ≤ B) (hB : B ≤ c * (2 / ε₀) * E)
    (hBI : BI ≤ KI * U1) (hBII : BII ≤ KII * U2) :
    N ≤ c * (2 / ε₀) * (KI + KII + 1) * P := by
  have hC0 : 0 ≤ c * (2 / ε₀) := mul_nonneg hc (div_nonneg zero_le_two hε₀.le)
  have hR0 : 0 ≤ KI * U1 + KII * U2 + 1 := by
    have h1 := mul_nonneg hKI hU1
    have h2 := mul_nonneg hKII hU2
    linarith only [h1, h2]
  have g1 : B * (x + 1) ≤ B * (KI * U1 + KII * U2 + 1) :=
    mul_le_mul_of_nonneg_left (by linarith only [hx, hBI, hBII]) hB0
  have g2 : B * (KI * U1 + KII * U2 + 1) ≤ c * (2 / ε₀) * E * (KI * U1 + KII * U2 + 1) :=
    mul_le_mul_of_nonneg_right hB hR0
  have g3 : E * (KI * U1 + KII * U2 + 1) ≤ (KI + KII + 1) * P := by
    have e1 : E * (KI * U1 + KII * U2 + 1) = KI * (E * U1) + KII * (E * U2) + E := by ring
    rw [e1, f1]
    have k2 := mul_le_mul_of_nonneg_left f2 hKII
    linarith only [k2, f3]
  have g4 := mul_le_mul_of_nonneg_left g3 hC0
  have e2 : c * (2 / ε₀) * E * (KI * U1 + KII * U2 + 1)
      = c * (2 / ε₀) * (E * (KI * U1 + KII * U2 + 1)) := by ring
  have e3 : c * (2 / ε₀) * ((KI + KII + 1) * P) = c * (2 / ε₀) * (KI + KII + 1) * P := by ring
  linarith only [hN, g1, g2, g4, e2, e3]

/-- the small-𝔔 case (RCApps2 l.33): N ≤ B(x+1), x ≤ h, h + 1 ≤ K_S, B ≤ c(2/ε₀)E, E ≤ P -/
lemma small_arith (N B x h c ε₀ KS E P : ℝ) (hε₀ : 0 < ε₀) (hc : 0 ≤ c) (_hE : 0 ≤ E)
    (hN : N ≤ B * (x + 1)) (hx : x ≤ h) (hKS : h + 1 ≤ KS) (hB0 : 0 ≤ B)
    (hB : B ≤ c * (2 / ε₀) * E) (f3 : E ≤ P) :
    N ≤ c * (2 / ε₀) * |KS| * P := by
  have hC0 : 0 ≤ c * (2 / ε₀) := mul_nonneg hc (div_nonneg zero_le_two hε₀.le)
  have hKS' : KS ≤ |KS| := le_abs_self KS
  have g1 : B * (x + 1) ≤ B * |KS| := mul_le_mul_of_nonneg_left (by linarith only [hx, hKS, hKS']) hB0
  have g2 : B * |KS| ≤ c * (2 / ε₀) * E * |KS| := mul_le_mul_of_nonneg_right hB (abs_nonneg _)
  have g3 : c * (2 / ε₀) * E * |KS| ≤ c * (2 / ε₀) * P * |KS| :=
    mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left f3 hC0) (abs_nonneg _)
  have e : c * (2 / ε₀) * P * |KS| = c * (2 / ε₀) * |KS| * P := by ring
  linarith only [hN, g1, g2, g3, e]

end RingClass

