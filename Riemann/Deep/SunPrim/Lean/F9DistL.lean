import Riemann.Deep.SunPrim.Lean.F9Defs

/-! Field axioms of `F9`, part DistL (`decide`). -/

namespace SunLean

namespace F9

theorem left_distrib' : ∀ a b c : F9, a * (b + c) = a * b + a * c := by decide +kernel

end F9

end SunLean
