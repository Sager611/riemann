import Riemann.Deep.SunPrim.Lean.F9Search

/-! Search leaves, part 0: `K k 0 x1 … xk` by `decide +kernel` (full subtrees at level 3; vacuous ones where the last guard fails). -/

namespace SunLean

namespace F9

theorem k1_0 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_1_0 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_1_1 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨1, by decide⟩) := by decide +kernel
theorem k2_2_0 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_2_2 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨2, by decide⟩) := by decide +kernel
theorem k2_3_0 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_3_3 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨3, by decide⟩) := by decide +kernel
theorem k2_3_6 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨6, by decide⟩) := by decide +kernel
theorem k2_4_0 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_4_4 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨4, by decide⟩) := by decide +kernel
theorem k2_4_6 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨6, by decide⟩) := by decide +kernel
theorem k2_5_0 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨5, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_5_5 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨5, by decide⟩) (ofFin ⟨5, by decide⟩) := by decide +kernel
theorem k2_5_6 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨5, by decide⟩) (ofFin ⟨6, by decide⟩) := by decide +kernel
theorem k2_6_0 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_6_3 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨3, by decide⟩) := by decide +kernel
theorem k2_6_6 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨6, by decide⟩) := by decide +kernel
theorem k2_7_0 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_7_3 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨3, by decide⟩) := by decide +kernel
theorem k2_7_7 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨7, by decide⟩) := by decide +kernel
theorem k2_8_0 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k2_8_3 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨3, by decide⟩) := by decide +kernel
theorem k2_8_8 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨8, by decide⟩) := by decide +kernel
theorem k3_1_2_0 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k3_1_2_1 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨1, by decide⟩) := by decide +kernel
theorem k3_1_2_2 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨2, by decide⟩) := by decide +kernel
theorem k3_1_2_3 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨3, by decide⟩) := by decide +kernel

end F9

end SunLean
