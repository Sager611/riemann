import Riemann.Deep.SunPrim.Lean.F9Defs

/-! Field axioms of `F9`, part Mul (`decide`). -/

namespace SunLean

namespace F9

theorem mul_assoc' : ∀ a b c : F9, a * b * c = a * (b * c) := by decide +kernel
theorem zero_mul' : ∀ a : F9, 0 * a = 0 := by decide +kernel
theorem mul_zero' : ∀ a : F9, a * 0 = 0 := by decide +kernel
theorem one_mul' : ∀ a : F9, 1 * a = a := by decide +kernel
theorem mul_one' : ∀ a : F9, a * 1 = a := by decide +kernel
theorem mul_comm' : ∀ a b : F9, a * b = b * a := by decide +kernel
theorem mul_inv_cancel' : ∀ a : F9, a ≠ 0 → a * a⁻¹ = 1 := by decide +kernel
theorem inv_zero' : (0 : F9)⁻¹ = 0 := by decide +kernel

end F9

end SunLean
