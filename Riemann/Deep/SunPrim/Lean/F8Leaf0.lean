import Riemann.Deep.SunPrim.Lean.F8Search

/-! Search leaves, part 0: `K k 0 x1 … xk` by `decide +kernel` (full subtrees at level 2; vacuous ones where the last guard fails). -/

namespace SunLean

namespace F8

theorem k1_0 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_1_0 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_1_1 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨1, by decide⟩) := by decide +kernel
theorem k2_1_2 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨2, by decide⟩) := by decide +kernel
theorem k2_1_3 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨3, by decide⟩) := by decide +kernel
theorem k2_1_4 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨4, by decide⟩) := by decide +kernel
theorem k2_1_5 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨5, by decide⟩) := by decide +kernel
theorem k2_1_6 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨6, by decide⟩) := by decide +kernel
theorem k2_1_7 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨7, by decide⟩) := by decide +kernel
theorem k2_2_0 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_2_1 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨1, by decide⟩) := by decide +kernel
theorem k2_2_2 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨2, by decide⟩) := by decide +kernel
theorem k2_2_3 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨3, by decide⟩) := by decide +kernel
theorem k2_2_4 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨4, by decide⟩) := by decide +kernel
theorem k2_2_5 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨5, by decide⟩) := by decide +kernel

end F8

end SunLean
