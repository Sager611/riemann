import Riemann.Deep.SunPrim.Lean.F8Defs

/-! Field axioms of `F8`, part Add (`decide`). -/

namespace SunLean

namespace F8

theorem add_assoc' : ∀ a b c : F8, a + b + c = a + (b + c) := by decide +kernel
theorem zero_add' : ∀ a : F8, 0 + a = a := by decide +kernel
theorem add_zero' : ∀ a : F8, a + 0 = a := by decide +kernel
theorem neg_add_cancel' : ∀ a : F8, -a + a = 0 := by decide +kernel
theorem add_comm' : ∀ a b : F8, a + b = b + a := by decide +kernel
theorem zero_ne_one' : (0 : F8) ≠ 1 := by decide +kernel

end F8

end SunLean
