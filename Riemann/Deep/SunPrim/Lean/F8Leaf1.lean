import Riemann.Deep.SunPrim.Lean.F8Search

/-! Search leaves, part 1: `K k 0 x1 … xk` by `decide +kernel` (full subtrees at level 2; vacuous ones where the last guard fails). -/

namespace SunLean

namespace F8

theorem k2_2_6 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨6, by decide⟩) := by decide +kernel
theorem k2_2_7 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨7, by decide⟩) := by decide +kernel
theorem k2_3_0 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_3_1 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨1, by decide⟩) := by decide +kernel
theorem k2_3_2 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨2, by decide⟩) := by decide +kernel
theorem k2_3_3 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨3, by decide⟩) := by decide +kernel
theorem k2_3_4 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨4, by decide⟩) := by decide +kernel
theorem k2_3_5 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨5, by decide⟩) := by decide +kernel
theorem k2_3_6 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨6, by decide⟩) := by decide +kernel
theorem k2_3_7 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨7, by decide⟩) := by decide +kernel
theorem k2_4_0 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_4_1 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨1, by decide⟩) := by decide +kernel
theorem k2_4_2 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨2, by decide⟩) := by decide +kernel
theorem k2_4_3 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨3, by decide⟩) := by decide +kernel
theorem k2_4_4 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨4, by decide⟩) := by decide +kernel

end F8

end SunLean
