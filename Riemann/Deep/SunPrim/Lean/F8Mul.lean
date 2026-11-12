import Riemann.Deep.SunPrim.Lean.F8Defs

/-! Field axioms of `F8`, part Mul (`decide`). -/

namespace SunLean

namespace F8

theorem mul_assoc' : ∀ a b c : F8, a * b * c = a * (b * c) := by decide +kernel
theorem zero_mul' : ∀ a : F8, 0 * a = 0 := by decide +kernel
theorem mul_zero' : ∀ a : F8, a * 0 = 0 := by decide +kernel
theorem one_mul' : ∀ a : F8, 1 * a = a := by decide +kernel
theorem mul_one' : ∀ a : F8, a * 1 = a := by decide +kernel
theorem mul_comm' : ∀ a b : F8, a * b = b * a := by decide +kernel
theorem mul_inv_cancel' : ∀ a : F8, a ≠ 0 → a * a⁻¹ = 1 := by decide +kernel
theorem inv_zero' : (0 : F8)⁻¹ = 0 := by decide +kernel

end F8

end SunLean
