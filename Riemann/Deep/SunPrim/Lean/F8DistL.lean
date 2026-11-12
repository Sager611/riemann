import Riemann.Deep.SunPrim.Lean.F8Defs

/-! Field axioms of `F8`, part DistL (`decide`). -/

namespace SunLean

namespace F8

theorem left_distrib' : ∀ a b c : F8, a * (b + c) = a * b + a * c := by decide +kernel

end F8

end SunLean
