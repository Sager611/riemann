import Riemann.Deep.SunPrim.Lean.F8Leaf0
import Riemann.Deep.SunPrim.Lean.F8Leaf1
import Riemann.Deep.SunPrim.Lean.F8Leaf2
import Riemann.Deep.SunPrim.Lean.F8Leaf3

/-! Assembly of the leaves: no circular order of F₈ starting at 0 has pairwise distinct sums. -/

namespace SunLean

namespace F8

theorem k1_1 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) := fun _ y _ =>
  match (y : Fin 8) with
  | ⟨0, _⟩ => k2_1_0
  | ⟨1, _⟩ => k2_1_1
  | ⟨2, _⟩ => k2_1_2
  | ⟨3, _⟩ => k2_1_3
  | ⟨4, _⟩ => k2_1_4
  | ⟨5, _⟩ => k2_1_5
  | ⟨6, _⟩ => k2_1_6
  | ⟨7, _⟩ => k2_1_7
  | ⟨n + 8, h⟩ => absurd h (by omega)

theorem k1_2 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) := fun _ y _ =>
  match (y : Fin 8) with
  | ⟨0, _⟩ => k2_2_0
  | ⟨1, _⟩ => k2_2_1
  | ⟨2, _⟩ => k2_2_2
  | ⟨3, _⟩ => k2_2_3
  | ⟨4, _⟩ => k2_2_4
  | ⟨5, _⟩ => k2_2_5
  | ⟨6, _⟩ => k2_2_6
  | ⟨7, _⟩ => k2_2_7
  | ⟨n + 8, h⟩ => absurd h (by omega)

theorem k1_3 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) := fun _ y _ =>
  match (y : Fin 8) with
  | ⟨0, _⟩ => k2_3_0
  | ⟨1, _⟩ => k2_3_1
  | ⟨2, _⟩ => k2_3_2
  | ⟨3, _⟩ => k2_3_3
  | ⟨4, _⟩ => k2_3_4
  | ⟨5, _⟩ => k2_3_5
  | ⟨6, _⟩ => k2_3_6
  | ⟨7, _⟩ => k2_3_7
  | ⟨n + 8, h⟩ => absurd h (by omega)

theorem k1_4 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) := fun _ y _ =>
  match (y : Fin 8) with
  | ⟨0, _⟩ => k2_4_0
  | ⟨1, _⟩ => k2_4_1
  | ⟨2, _⟩ => k2_4_2
  | ⟨3, _⟩ => k2_4_3
  | ⟨4, _⟩ => k2_4_4
  | ⟨5, _⟩ => k2_4_5
  | ⟨6, _⟩ => k2_4_6
  | ⟨7, _⟩ => k2_4_7
  | ⟨n + 8, h⟩ => absurd h (by omega)

theorem k1_5 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨5, by decide⟩) := fun _ y _ =>
  match (y : Fin 8) with
  | ⟨0, _⟩ => k2_5_0
  | ⟨1, _⟩ => k2_5_1
  | ⟨2, _⟩ => k2_5_2
  | ⟨3, _⟩ => k2_5_3
  | ⟨4, _⟩ => k2_5_4
  | ⟨5, _⟩ => k2_5_5
  | ⟨6, _⟩ => k2_5_6
  | ⟨7, _⟩ => k2_5_7
  | ⟨n + 8, h⟩ => absurd h (by omega)

theorem k1_6 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) := fun _ y _ =>
  match (y : Fin 8) with
  | ⟨0, _⟩ => k2_6_0
  | ⟨1, _⟩ => k2_6_1
  | ⟨2, _⟩ => k2_6_2
  | ⟨3, _⟩ => k2_6_3
  | ⟨4, _⟩ => k2_6_4
  | ⟨5, _⟩ => k2_6_5
  | ⟨6, _⟩ => k2_6_6
  | ⟨7, _⟩ => k2_6_7
  | ⟨n + 8, h⟩ => absurd h (by omega)

theorem k1_7 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) := fun _ y _ =>
  match (y : Fin 8) with
  | ⟨0, _⟩ => k2_7_0
  | ⟨1, _⟩ => k2_7_1
  | ⟨2, _⟩ => k2_7_2
  | ⟨3, _⟩ => k2_7_3
  | ⟨4, _⟩ => k2_7_4
  | ⟨5, _⟩ => k2_7_5
  | ⟨6, _⟩ => k2_7_6
  | ⟨7, _⟩ => k2_7_7
  | ⟨n + 8, h⟩ => absurd h (by omega)

theorem no_order_zero : ∀ x1 ∈ LL, K1 (ofFin ⟨0, by decide⟩) x1 := fun y _ =>
  match (y : Fin 8) with
  | ⟨0, _⟩ => k1_0
  | ⟨1, _⟩ => k1_1
  | ⟨2, _⟩ => k1_2
  | ⟨3, _⟩ => k1_3
  | ⟨4, _⟩ => k1_4
  | ⟨5, _⟩ => k1_5
  | ⟨6, _⟩ => k1_6
  | ⟨7, _⟩ => k1_7
  | ⟨n + 8, h⟩ => absurd h (by omega)

end F8

end SunLean

#print axioms SunLean.F8.no_order_zero
