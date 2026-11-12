import Riemann.Deep.SunPrim.Lean.F8Search

/-! Search leaves, part 3: `K k 0 x1 … xk` by `decide +kernel` (full subtrees at level 2; vacuous ones where the last guard fails). -/

namespace SunLean

namespace F8

theorem k2_6_4 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨4, by decide⟩) := by decide +kernel
theorem k2_6_5 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨5, by decide⟩) := by decide +kernel
theorem k2_6_6 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨6, by decide⟩) := by decide +kernel
theorem k2_6_7 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨7, by decide⟩) := by decide +kernel
theorem k2_7_0 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_7_1 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨1, by decide⟩) := by decide +kernel
theorem k2_7_2 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨2, by decide⟩) := by decide +kernel
theorem k2_7_3 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨3, by decide⟩) := by decide +kernel
theorem k2_7_4 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨4, by decide⟩) := by decide +kernel
theorem k2_7_5 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨5, by decide⟩) := by decide +kernel
theorem k2_7_6 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨6, by decide⟩) := by decide +kernel
theorem k2_7_7 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨7, by decide⟩) := by decide +kernel

end F8

end SunLean
