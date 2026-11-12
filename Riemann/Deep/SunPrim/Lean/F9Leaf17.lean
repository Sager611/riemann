import Riemann.Deep.SunPrim.Lean.F9Search

/-! Search leaves, part 17: `K k 0 x1 … xk` by `decide +kernel` (full subtrees at level 3; vacuous ones where the last guard fails). -/

namespace SunLean

namespace F9

theorem k3_8_6_4 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨4, by decide⟩) := by decide +kernel
theorem k3_8_6_5 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨5, by decide⟩) := by decide +kernel
theorem k3_8_6_6 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨6, by decide⟩) := by decide +kernel
theorem k3_8_6_7 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨7, by decide⟩) := by decide +kernel
theorem k3_8_6_8 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨8, by decide⟩) := by decide +kernel
theorem k3_8_7_0 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨0, by decide⟩) := by decide +kernel
theorem k3_8_7_1 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨1, by decide⟩) := by decide +kernel
theorem k3_8_7_2 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨2, by decide⟩) := by decide +kernel
theorem k3_8_7_3 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨3, by decide⟩) := by decide +kernel
theorem k3_8_7_4 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨4, by decide⟩) := by decide +kernel
theorem k3_8_7_5 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨5, by decide⟩) := by decide +kernel
theorem k3_8_7_6 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨6, by decide⟩) := by decide +kernel
theorem k3_8_7_7 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨7, by decide⟩) := by decide +kernel
theorem k3_8_7_8 : K3 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨8, by decide⟩) := by decide +kernel

end F9

end SunLean
