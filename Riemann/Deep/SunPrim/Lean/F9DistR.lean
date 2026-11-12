import Riemann.Deep.SunPrim.Lean.F9Defs

/-! Field axioms of `F9`, part DistR (`decide`). -/

namespace SunLean

namespace F9

theorem right_distrib' : ∀ a b c : F9, (a + b) * c = a * c + b * c := by decide +kernel

end F9

end SunLean
