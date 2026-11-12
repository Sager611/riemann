import Riemann.Deep.SunPrim.Lean.F9Defs

/-! Field axioms of `F9`, part Add (`decide`). -/

namespace SunLean

namespace F9

theorem add_assoc' : ∀ a b c : F9, a + b + c = a + (b + c) := by decide +kernel
theorem zero_add' : ∀ a : F9, 0 + a = a := by decide +kernel
theorem add_zero' : ∀ a : F9, a + 0 = a := by decide +kernel
theorem neg_add_cancel' : ∀ a : F9, -a + a = 0 := by decide +kernel
theorem add_comm' : ∀ a b : F9, a + b = b + a := by decide +kernel
theorem zero_ne_one' : (0 : F9) ≠ 1 := by decide +kernel

end F9

end SunLean
