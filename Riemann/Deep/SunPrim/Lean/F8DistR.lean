import Riemann.Deep.SunPrim.Lean.F8Defs

/-! Field axioms of `F8`, part DistR (`decide`). -/

namespace SunLean

namespace F8

theorem right_distrib' : ∀ a b c : F8, (a + b) * c = a * c + b * c := by decide +kernel

end F8

end SunLean
