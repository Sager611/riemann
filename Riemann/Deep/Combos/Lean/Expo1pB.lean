module

public import Riemann.Deep.Combos.Lean.Expo1pA

@[expose] public section

/-! # Theorem 1′(a): class II exponent (RCApps2 l.47–48, 59) -/

open Real

namespace RingClass

/-- class II (RCApps2 l.47–48, 59): Y^{1/2−σ} 𝔔^{1/2+ε₀} T^{1/2+ε₀} = (𝔔T)^{e(σ)+ε₀} -/
lemma classII_exp (Q T σ ε₀ : ℝ) (hQ : 1 ≤ Q) (hT : 1 ≤ T) (h1 : σ ≤ 1) :
    ((Q * T) ^ (1 / (3 - 2 * σ))) ^ (1 / 2 - σ) * Q ^ (1 / 2 + ε₀) * T ^ (1 / 2 + ε₀)
      = (Q * T) ^ (ee σ + ε₀) := by
  have hW : 0 < Q * T := by positivity
  rw [mul_assoc, ← Real.mul_rpow (by linarith) (by linarith), ← Real.rpow_mul hW.le,
    ← Real.rpow_add hW, ← balance σ h1]
  congr 1; ring

end RingClass
