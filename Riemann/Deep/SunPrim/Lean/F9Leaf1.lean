import Riemann.Deep.SunPrim.Lean.F9Search

/-! Search leaves, part 1: `K k 0 x1 … xk` by `decide +kernel` (full subtrees at level 3; vacuous ones where the last guard fails). -/

namespace SunLean

namespace F9

theorem k3_1_2_4 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨4, by decide⟩) := by decide +kernel
theorem k3_1_2_5 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨5, by decide⟩) := by decide +kernel
theorem k3_1_2_6 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨6, by decide⟩) := by decide +kernel
theorem k3_1_2_7 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨7, by decide⟩) := by decide +kernel
theorem k3_1_2_8 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨8, by decide⟩) := by decide +kernel
theorem k3_1_3_0 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k3_1_3_1 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨1, by decide⟩) := by decide +kernel
theorem k3_1_3_2 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨2, by decide⟩) := by decide +kernel
theorem k3_1_3_3 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨3, by decide⟩) := by decide +kernel
theorem k3_1_3_4 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨4, by decide⟩) := by decide +kernel
theorem k3_1_3_5 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨5, by decide⟩) := by decide +kernel
theorem k3_1_3_6 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨6, by decide⟩) := by decide +kernel
theorem k3_1_3_7 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨7, by decide⟩) := by decide +kernel
theorem k3_1_3_8 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨8, by decide⟩) := by decide +kernel
theorem k3_1_4_0 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k3_1_4_1 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨1, by decide⟩) := by decide +kernel
theorem k3_1_4_2 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨2, by decide⟩) := by decide +kernel
theorem k3_1_4_3 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨3, by decide⟩) := by decide +kernel
theorem k3_1_4_4 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨4, by decide⟩) := by decide +kernel
theorem k3_1_4_5 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨5, by decide⟩) := by decide +kernel
theorem k3_1_4_6 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨6, by decide⟩) := by decide +kernel
theorem k3_1_4_7 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨7, by decide⟩) := by decide +kernel
theorem k3_1_4_8 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨8, by decide⟩) := by decide +kernel
theorem k3_1_5_0 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨5, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k3_1_5_1 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨5, by decide⟩) (ofFin ⟨1, by decide⟩) := by decide +kernel
theorem k3_1_5_2 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨5, by decide⟩) (ofFin ⟨2, by decide⟩) := by decide +kernel
theorem k3_1_5_3 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨5, by decide⟩) (ofFin ⟨3, by decide⟩) := by decide +kernel

end F9

end SunLean
