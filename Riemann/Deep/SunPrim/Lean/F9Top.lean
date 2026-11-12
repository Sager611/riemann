import Riemann.Deep.SunPrim.Lean.F9Leaf0
import Riemann.Deep.SunPrim.Lean.F9Leaf1
import Riemann.Deep.SunPrim.Lean.F9Leaf2
import Riemann.Deep.SunPrim.Lean.F9Leaf3
import Riemann.Deep.SunPrim.Lean.F9Leaf4
import Riemann.Deep.SunPrim.Lean.F9Leaf5
import Riemann.Deep.SunPrim.Lean.F9Leaf6
import Riemann.Deep.SunPrim.Lean.F9Leaf7
import Riemann.Deep.SunPrim.Lean.F9Leaf8
import Riemann.Deep.SunPrim.Lean.F9Leaf9
import Riemann.Deep.SunPrim.Lean.F9Leaf10
import Riemann.Deep.SunPrim.Lean.F9Leaf11
import Riemann.Deep.SunPrim.Lean.F9Leaf12
import Riemann.Deep.SunPrim.Lean.F9Leaf13
import Riemann.Deep.SunPrim.Lean.F9Leaf14
import Riemann.Deep.SunPrim.Lean.F9Leaf15
import Riemann.Deep.SunPrim.Lean.F9Leaf16
import Riemann.Deep.SunPrim.Lean.F9Leaf17

/-! Assembly of the leaves: no circular order of F₉ starting at 0 has pairwise distinct sums. -/

namespace SunLean

namespace F9

theorem k2_1_2 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨2, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_1_2_0
  | ⟨1, _⟩ => k3_1_2_1
  | ⟨2, _⟩ => k3_1_2_2
  | ⟨3, _⟩ => k3_1_2_3
  | ⟨4, _⟩ => k3_1_2_4
  | ⟨5, _⟩ => k3_1_2_5
  | ⟨6, _⟩ => k3_1_2_6
  | ⟨7, _⟩ => k3_1_2_7
  | ⟨8, _⟩ => k3_1_2_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_1_3 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨3, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_1_3_0
  | ⟨1, _⟩ => k3_1_3_1
  | ⟨2, _⟩ => k3_1_3_2
  | ⟨3, _⟩ => k3_1_3_3
  | ⟨4, _⟩ => k3_1_3_4
  | ⟨5, _⟩ => k3_1_3_5
  | ⟨6, _⟩ => k3_1_3_6
  | ⟨7, _⟩ => k3_1_3_7
  | ⟨8, _⟩ => k3_1_3_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_1_4 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨4, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_1_4_0
  | ⟨1, _⟩ => k3_1_4_1
  | ⟨2, _⟩ => k3_1_4_2
  | ⟨3, _⟩ => k3_1_4_3
  | ⟨4, _⟩ => k3_1_4_4
  | ⟨5, _⟩ => k3_1_4_5
  | ⟨6, _⟩ => k3_1_4_6
  | ⟨7, _⟩ => k3_1_4_7
  | ⟨8, _⟩ => k3_1_4_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_1_5 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨5, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_1_5_0
  | ⟨1, _⟩ => k3_1_5_1
  | ⟨2, _⟩ => k3_1_5_2
  | ⟨3, _⟩ => k3_1_5_3
  | ⟨4, _⟩ => k3_1_5_4
  | ⟨5, _⟩ => k3_1_5_5
  | ⟨6, _⟩ => k3_1_5_6
  | ⟨7, _⟩ => k3_1_5_7
  | ⟨8, _⟩ => k3_1_5_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_1_6 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨6, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_1_6_0
  | ⟨1, _⟩ => k3_1_6_1
  | ⟨2, _⟩ => k3_1_6_2
  | ⟨3, _⟩ => k3_1_6_3
  | ⟨4, _⟩ => k3_1_6_4
  | ⟨5, _⟩ => k3_1_6_5
  | ⟨6, _⟩ => k3_1_6_6
  | ⟨7, _⟩ => k3_1_6_7
  | ⟨8, _⟩ => k3_1_6_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_1_7 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨7, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_1_7_0
  | ⟨1, _⟩ => k3_1_7_1
  | ⟨2, _⟩ => k3_1_7_2
  | ⟨3, _⟩ => k3_1_7_3
  | ⟨4, _⟩ => k3_1_7_4
  | ⟨5, _⟩ => k3_1_7_5
  | ⟨6, _⟩ => k3_1_7_6
  | ⟨7, _⟩ => k3_1_7_7
  | ⟨8, _⟩ => k3_1_7_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_1_8 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) (ofFin ⟨8, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_1_8_0
  | ⟨1, _⟩ => k3_1_8_1
  | ⟨2, _⟩ => k3_1_8_2
  | ⟨3, _⟩ => k3_1_8_3
  | ⟨4, _⟩ => k3_1_8_4
  | ⟨5, _⟩ => k3_1_8_5
  | ⟨6, _⟩ => k3_1_8_6
  | ⟨7, _⟩ => k3_1_8_7
  | ⟨8, _⟩ => k3_1_8_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_2_1 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨1, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_2_1_0
  | ⟨1, _⟩ => k3_2_1_1
  | ⟨2, _⟩ => k3_2_1_2
  | ⟨3, _⟩ => k3_2_1_3
  | ⟨4, _⟩ => k3_2_1_4
  | ⟨5, _⟩ => k3_2_1_5
  | ⟨6, _⟩ => k3_2_1_6
  | ⟨7, _⟩ => k3_2_1_7
  | ⟨8, _⟩ => k3_2_1_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_2_3 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨3, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_2_3_0
  | ⟨1, _⟩ => k3_2_3_1
  | ⟨2, _⟩ => k3_2_3_2
  | ⟨3, _⟩ => k3_2_3_3
  | ⟨4, _⟩ => k3_2_3_4
  | ⟨5, _⟩ => k3_2_3_5
  | ⟨6, _⟩ => k3_2_3_6
  | ⟨7, _⟩ => k3_2_3_7
  | ⟨8, _⟩ => k3_2_3_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_2_4 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨4, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_2_4_0
  | ⟨1, _⟩ => k3_2_4_1
  | ⟨2, _⟩ => k3_2_4_2
  | ⟨3, _⟩ => k3_2_4_3
  | ⟨4, _⟩ => k3_2_4_4
  | ⟨5, _⟩ => k3_2_4_5
  | ⟨6, _⟩ => k3_2_4_6
  | ⟨7, _⟩ => k3_2_4_7
  | ⟨8, _⟩ => k3_2_4_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_2_5 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨5, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_2_5_0
  | ⟨1, _⟩ => k3_2_5_1
  | ⟨2, _⟩ => k3_2_5_2
  | ⟨3, _⟩ => k3_2_5_3
  | ⟨4, _⟩ => k3_2_5_4
  | ⟨5, _⟩ => k3_2_5_5
  | ⟨6, _⟩ => k3_2_5_6
  | ⟨7, _⟩ => k3_2_5_7
  | ⟨8, _⟩ => k3_2_5_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_2_6 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨6, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_2_6_0
  | ⟨1, _⟩ => k3_2_6_1
  | ⟨2, _⟩ => k3_2_6_2
  | ⟨3, _⟩ => k3_2_6_3
  | ⟨4, _⟩ => k3_2_6_4
  | ⟨5, _⟩ => k3_2_6_5
  | ⟨6, _⟩ => k3_2_6_6
  | ⟨7, _⟩ => k3_2_6_7
  | ⟨8, _⟩ => k3_2_6_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_2_7 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨7, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_2_7_0
  | ⟨1, _⟩ => k3_2_7_1
  | ⟨2, _⟩ => k3_2_7_2
  | ⟨3, _⟩ => k3_2_7_3
  | ⟨4, _⟩ => k3_2_7_4
  | ⟨5, _⟩ => k3_2_7_5
  | ⟨6, _⟩ => k3_2_7_6
  | ⟨7, _⟩ => k3_2_7_7
  | ⟨8, _⟩ => k3_2_7_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_2_8 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) (ofFin ⟨8, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_2_8_0
  | ⟨1, _⟩ => k3_2_8_1
  | ⟨2, _⟩ => k3_2_8_2
  | ⟨3, _⟩ => k3_2_8_3
  | ⟨4, _⟩ => k3_2_8_4
  | ⟨5, _⟩ => k3_2_8_5
  | ⟨6, _⟩ => k3_2_8_6
  | ⟨7, _⟩ => k3_2_8_7
  | ⟨8, _⟩ => k3_2_8_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_3_1 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨1, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_3_1_0
  | ⟨1, _⟩ => k3_3_1_1
  | ⟨2, _⟩ => k3_3_1_2
  | ⟨3, _⟩ => k3_3_1_3
  | ⟨4, _⟩ => k3_3_1_4
  | ⟨5, _⟩ => k3_3_1_5
  | ⟨6, _⟩ => k3_3_1_6
  | ⟨7, _⟩ => k3_3_1_7
  | ⟨8, _⟩ => k3_3_1_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_3_2 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨2, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_3_2_0
  | ⟨1, _⟩ => k3_3_2_1
  | ⟨2, _⟩ => k3_3_2_2
  | ⟨3, _⟩ => k3_3_2_3
  | ⟨4, _⟩ => k3_3_2_4
  | ⟨5, _⟩ => k3_3_2_5
  | ⟨6, _⟩ => k3_3_2_6
  | ⟨7, _⟩ => k3_3_2_7
  | ⟨8, _⟩ => k3_3_2_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_3_4 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨4, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_3_4_0
  | ⟨1, _⟩ => k3_3_4_1
  | ⟨2, _⟩ => k3_3_4_2
  | ⟨3, _⟩ => k3_3_4_3
  | ⟨4, _⟩ => k3_3_4_4
  | ⟨5, _⟩ => k3_3_4_5
  | ⟨6, _⟩ => k3_3_4_6
  | ⟨7, _⟩ => k3_3_4_7
  | ⟨8, _⟩ => k3_3_4_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_3_5 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨5, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_3_5_0
  | ⟨1, _⟩ => k3_3_5_1
  | ⟨2, _⟩ => k3_3_5_2
  | ⟨3, _⟩ => k3_3_5_3
  | ⟨4, _⟩ => k3_3_5_4
  | ⟨5, _⟩ => k3_3_5_5
  | ⟨6, _⟩ => k3_3_5_6
  | ⟨7, _⟩ => k3_3_5_7
  | ⟨8, _⟩ => k3_3_5_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_3_7 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨7, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_3_7_0
  | ⟨1, _⟩ => k3_3_7_1
  | ⟨2, _⟩ => k3_3_7_2
  | ⟨3, _⟩ => k3_3_7_3
  | ⟨4, _⟩ => k3_3_7_4
  | ⟨5, _⟩ => k3_3_7_5
  | ⟨6, _⟩ => k3_3_7_6
  | ⟨7, _⟩ => k3_3_7_7
  | ⟨8, _⟩ => k3_3_7_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_3_8 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) (ofFin ⟨8, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_3_8_0
  | ⟨1, _⟩ => k3_3_8_1
  | ⟨2, _⟩ => k3_3_8_2
  | ⟨3, _⟩ => k3_3_8_3
  | ⟨4, _⟩ => k3_3_8_4
  | ⟨5, _⟩ => k3_3_8_5
  | ⟨6, _⟩ => k3_3_8_6
  | ⟨7, _⟩ => k3_3_8_7
  | ⟨8, _⟩ => k3_3_8_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_4_1 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨1, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_4_1_0
  | ⟨1, _⟩ => k3_4_1_1
  | ⟨2, _⟩ => k3_4_1_2
  | ⟨3, _⟩ => k3_4_1_3
  | ⟨4, _⟩ => k3_4_1_4
  | ⟨5, _⟩ => k3_4_1_5
  | ⟨6, _⟩ => k3_4_1_6
  | ⟨7, _⟩ => k3_4_1_7
  | ⟨8, _⟩ => k3_4_1_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_4_2 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨2, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_4_2_0
  | ⟨1, _⟩ => k3_4_2_1
  | ⟨2, _⟩ => k3_4_2_2
  | ⟨3, _⟩ => k3_4_2_3
  | ⟨4, _⟩ => k3_4_2_4
  | ⟨5, _⟩ => k3_4_2_5
  | ⟨6, _⟩ => k3_4_2_6
  | ⟨7, _⟩ => k3_4_2_7
  | ⟨8, _⟩ => k3_4_2_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_4_3 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨3, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_4_3_0
  | ⟨1, _⟩ => k3_4_3_1
  | ⟨2, _⟩ => k3_4_3_2
  | ⟨3, _⟩ => k3_4_3_3
  | ⟨4, _⟩ => k3_4_3_4
  | ⟨5, _⟩ => k3_4_3_5
  | ⟨6, _⟩ => k3_4_3_6
  | ⟨7, _⟩ => k3_4_3_7
  | ⟨8, _⟩ => k3_4_3_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_4_5 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨5, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_4_5_0
  | ⟨1, _⟩ => k3_4_5_1
  | ⟨2, _⟩ => k3_4_5_2
  | ⟨3, _⟩ => k3_4_5_3
  | ⟨4, _⟩ => k3_4_5_4
  | ⟨5, _⟩ => k3_4_5_5
  | ⟨6, _⟩ => k3_4_5_6
  | ⟨7, _⟩ => k3_4_5_7
  | ⟨8, _⟩ => k3_4_5_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_4_7 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨7, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_4_7_0
  | ⟨1, _⟩ => k3_4_7_1
  | ⟨2, _⟩ => k3_4_7_2
  | ⟨3, _⟩ => k3_4_7_3
  | ⟨4, _⟩ => k3_4_7_4
  | ⟨5, _⟩ => k3_4_7_5
  | ⟨6, _⟩ => k3_4_7_6
  | ⟨7, _⟩ => k3_4_7_7
  | ⟨8, _⟩ => k3_4_7_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_4_8 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) (ofFin ⟨8, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_4_8_0
  | ⟨1, _⟩ => k3_4_8_1
  | ⟨2, _⟩ => k3_4_8_2
  | ⟨3, _⟩ => k3_4_8_3
  | ⟨4, _⟩ => k3_4_8_4
  | ⟨5, _⟩ => k3_4_8_5
  | ⟨6, _⟩ => k3_4_8_6
  | ⟨7, _⟩ => k3_4_8_7
  | ⟨8, _⟩ => k3_4_8_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_5_1 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨5, by decide⟩) (ofFin ⟨1, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_5_1_0
  | ⟨1, _⟩ => k3_5_1_1
  | ⟨2, _⟩ => k3_5_1_2
  | ⟨3, _⟩ => k3_5_1_3
  | ⟨4, _⟩ => k3_5_1_4
  | ⟨5, _⟩ => k3_5_1_5
  | ⟨6, _⟩ => k3_5_1_6
  | ⟨7, _⟩ => k3_5_1_7
  | ⟨8, _⟩ => k3_5_1_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_5_2 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨5, by decide⟩) (ofFin ⟨2, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_5_2_0
  | ⟨1, _⟩ => k3_5_2_1
  | ⟨2, _⟩ => k3_5_2_2
  | ⟨3, _⟩ => k3_5_2_3
  | ⟨4, _⟩ => k3_5_2_4
  | ⟨5, _⟩ => k3_5_2_5
  | ⟨6, _⟩ => k3_5_2_6
  | ⟨7, _⟩ => k3_5_2_7
  | ⟨8, _⟩ => k3_5_2_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_5_3 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨5, by decide⟩) (ofFin ⟨3, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_5_3_0
  | ⟨1, _⟩ => k3_5_3_1
  | ⟨2, _⟩ => k3_5_3_2
  | ⟨3, _⟩ => k3_5_3_3
  | ⟨4, _⟩ => k3_5_3_4
  | ⟨5, _⟩ => k3_5_3_5
  | ⟨6, _⟩ => k3_5_3_6
  | ⟨7, _⟩ => k3_5_3_7
  | ⟨8, _⟩ => k3_5_3_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_5_4 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨5, by decide⟩) (ofFin ⟨4, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_5_4_0
  | ⟨1, _⟩ => k3_5_4_1
  | ⟨2, _⟩ => k3_5_4_2
  | ⟨3, _⟩ => k3_5_4_3
  | ⟨4, _⟩ => k3_5_4_4
  | ⟨5, _⟩ => k3_5_4_5
  | ⟨6, _⟩ => k3_5_4_6
  | ⟨7, _⟩ => k3_5_4_7
  | ⟨8, _⟩ => k3_5_4_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_5_7 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨5, by decide⟩) (ofFin ⟨7, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_5_7_0
  | ⟨1, _⟩ => k3_5_7_1
  | ⟨2, _⟩ => k3_5_7_2
  | ⟨3, _⟩ => k3_5_7_3
  | ⟨4, _⟩ => k3_5_7_4
  | ⟨5, _⟩ => k3_5_7_5
  | ⟨6, _⟩ => k3_5_7_6
  | ⟨7, _⟩ => k3_5_7_7
  | ⟨8, _⟩ => k3_5_7_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_5_8 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨5, by decide⟩) (ofFin ⟨8, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_5_8_0
  | ⟨1, _⟩ => k3_5_8_1
  | ⟨2, _⟩ => k3_5_8_2
  | ⟨3, _⟩ => k3_5_8_3
  | ⟨4, _⟩ => k3_5_8_4
  | ⟨5, _⟩ => k3_5_8_5
  | ⟨6, _⟩ => k3_5_8_6
  | ⟨7, _⟩ => k3_5_8_7
  | ⟨8, _⟩ => k3_5_8_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_6_1 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨1, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_6_1_0
  | ⟨1, _⟩ => k3_6_1_1
  | ⟨2, _⟩ => k3_6_1_2
  | ⟨3, _⟩ => k3_6_1_3
  | ⟨4, _⟩ => k3_6_1_4
  | ⟨5, _⟩ => k3_6_1_5
  | ⟨6, _⟩ => k3_6_1_6
  | ⟨7, _⟩ => k3_6_1_7
  | ⟨8, _⟩ => k3_6_1_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_6_2 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨2, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_6_2_0
  | ⟨1, _⟩ => k3_6_2_1
  | ⟨2, _⟩ => k3_6_2_2
  | ⟨3, _⟩ => k3_6_2_3
  | ⟨4, _⟩ => k3_6_2_4
  | ⟨5, _⟩ => k3_6_2_5
  | ⟨6, _⟩ => k3_6_2_6
  | ⟨7, _⟩ => k3_6_2_7
  | ⟨8, _⟩ => k3_6_2_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_6_4 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨4, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_6_4_0
  | ⟨1, _⟩ => k3_6_4_1
  | ⟨2, _⟩ => k3_6_4_2
  | ⟨3, _⟩ => k3_6_4_3
  | ⟨4, _⟩ => k3_6_4_4
  | ⟨5, _⟩ => k3_6_4_5
  | ⟨6, _⟩ => k3_6_4_6
  | ⟨7, _⟩ => k3_6_4_7
  | ⟨8, _⟩ => k3_6_4_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_6_5 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨5, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_6_5_0
  | ⟨1, _⟩ => k3_6_5_1
  | ⟨2, _⟩ => k3_6_5_2
  | ⟨3, _⟩ => k3_6_5_3
  | ⟨4, _⟩ => k3_6_5_4
  | ⟨5, _⟩ => k3_6_5_5
  | ⟨6, _⟩ => k3_6_5_6
  | ⟨7, _⟩ => k3_6_5_7
  | ⟨8, _⟩ => k3_6_5_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_6_7 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨7, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_6_7_0
  | ⟨1, _⟩ => k3_6_7_1
  | ⟨2, _⟩ => k3_6_7_2
  | ⟨3, _⟩ => k3_6_7_3
  | ⟨4, _⟩ => k3_6_7_4
  | ⟨5, _⟩ => k3_6_7_5
  | ⟨6, _⟩ => k3_6_7_6
  | ⟨7, _⟩ => k3_6_7_7
  | ⟨8, _⟩ => k3_6_7_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_6_8 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) (ofFin ⟨8, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_6_8_0
  | ⟨1, _⟩ => k3_6_8_1
  | ⟨2, _⟩ => k3_6_8_2
  | ⟨3, _⟩ => k3_6_8_3
  | ⟨4, _⟩ => k3_6_8_4
  | ⟨5, _⟩ => k3_6_8_5
  | ⟨6, _⟩ => k3_6_8_6
  | ⟨7, _⟩ => k3_6_8_7
  | ⟨8, _⟩ => k3_6_8_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_7_1 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨1, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_7_1_0
  | ⟨1, _⟩ => k3_7_1_1
  | ⟨2, _⟩ => k3_7_1_2
  | ⟨3, _⟩ => k3_7_1_3
  | ⟨4, _⟩ => k3_7_1_4
  | ⟨5, _⟩ => k3_7_1_5
  | ⟨6, _⟩ => k3_7_1_6
  | ⟨7, _⟩ => k3_7_1_7
  | ⟨8, _⟩ => k3_7_1_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_7_2 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨2, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_7_2_0
  | ⟨1, _⟩ => k3_7_2_1
  | ⟨2, _⟩ => k3_7_2_2
  | ⟨3, _⟩ => k3_7_2_3
  | ⟨4, _⟩ => k3_7_2_4
  | ⟨5, _⟩ => k3_7_2_5
  | ⟨6, _⟩ => k3_7_2_6
  | ⟨7, _⟩ => k3_7_2_7
  | ⟨8, _⟩ => k3_7_2_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_7_4 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨4, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_7_4_0
  | ⟨1, _⟩ => k3_7_4_1
  | ⟨2, _⟩ => k3_7_4_2
  | ⟨3, _⟩ => k3_7_4_3
  | ⟨4, _⟩ => k3_7_4_4
  | ⟨5, _⟩ => k3_7_4_5
  | ⟨6, _⟩ => k3_7_4_6
  | ⟨7, _⟩ => k3_7_4_7
  | ⟨8, _⟩ => k3_7_4_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_7_5 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨5, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_7_5_0
  | ⟨1, _⟩ => k3_7_5_1
  | ⟨2, _⟩ => k3_7_5_2
  | ⟨3, _⟩ => k3_7_5_3
  | ⟨4, _⟩ => k3_7_5_4
  | ⟨5, _⟩ => k3_7_5_5
  | ⟨6, _⟩ => k3_7_5_6
  | ⟨7, _⟩ => k3_7_5_7
  | ⟨8, _⟩ => k3_7_5_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_7_6 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨6, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_7_6_0
  | ⟨1, _⟩ => k3_7_6_1
  | ⟨2, _⟩ => k3_7_6_2
  | ⟨3, _⟩ => k3_7_6_3
  | ⟨4, _⟩ => k3_7_6_4
  | ⟨5, _⟩ => k3_7_6_5
  | ⟨6, _⟩ => k3_7_6_6
  | ⟨7, _⟩ => k3_7_6_7
  | ⟨8, _⟩ => k3_7_6_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_7_8 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) (ofFin ⟨8, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_7_8_0
  | ⟨1, _⟩ => k3_7_8_1
  | ⟨2, _⟩ => k3_7_8_2
  | ⟨3, _⟩ => k3_7_8_3
  | ⟨4, _⟩ => k3_7_8_4
  | ⟨5, _⟩ => k3_7_8_5
  | ⟨6, _⟩ => k3_7_8_6
  | ⟨7, _⟩ => k3_7_8_7
  | ⟨8, _⟩ => k3_7_8_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_8_1 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨1, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_8_1_0
  | ⟨1, _⟩ => k3_8_1_1
  | ⟨2, _⟩ => k3_8_1_2
  | ⟨3, _⟩ => k3_8_1_3
  | ⟨4, _⟩ => k3_8_1_4
  | ⟨5, _⟩ => k3_8_1_5
  | ⟨6, _⟩ => k3_8_1_6
  | ⟨7, _⟩ => k3_8_1_7
  | ⟨8, _⟩ => k3_8_1_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_8_2 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨2, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_8_2_0
  | ⟨1, _⟩ => k3_8_2_1
  | ⟨2, _⟩ => k3_8_2_2
  | ⟨3, _⟩ => k3_8_2_3
  | ⟨4, _⟩ => k3_8_2_4
  | ⟨5, _⟩ => k3_8_2_5
  | ⟨6, _⟩ => k3_8_2_6
  | ⟨7, _⟩ => k3_8_2_7
  | ⟨8, _⟩ => k3_8_2_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_8_4 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨4, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_8_4_0
  | ⟨1, _⟩ => k3_8_4_1
  | ⟨2, _⟩ => k3_8_4_2
  | ⟨3, _⟩ => k3_8_4_3
  | ⟨4, _⟩ => k3_8_4_4
  | ⟨5, _⟩ => k3_8_4_5
  | ⟨6, _⟩ => k3_8_4_6
  | ⟨7, _⟩ => k3_8_4_7
  | ⟨8, _⟩ => k3_8_4_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_8_5 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨5, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_8_5_0
  | ⟨1, _⟩ => k3_8_5_1
  | ⟨2, _⟩ => k3_8_5_2
  | ⟨3, _⟩ => k3_8_5_3
  | ⟨4, _⟩ => k3_8_5_4
  | ⟨5, _⟩ => k3_8_5_5
  | ⟨6, _⟩ => k3_8_5_6
  | ⟨7, _⟩ => k3_8_5_7
  | ⟨8, _⟩ => k3_8_5_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_8_6 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨6, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_8_6_0
  | ⟨1, _⟩ => k3_8_6_1
  | ⟨2, _⟩ => k3_8_6_2
  | ⟨3, _⟩ => k3_8_6_3
  | ⟨4, _⟩ => k3_8_6_4
  | ⟨5, _⟩ => k3_8_6_5
  | ⟨6, _⟩ => k3_8_6_6
  | ⟨7, _⟩ => k3_8_6_7
  | ⟨8, _⟩ => k3_8_6_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k2_8_7 : K2 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) (ofFin ⟨7, by decide⟩) := fun _ _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k3_8_7_0
  | ⟨1, _⟩ => k3_8_7_1
  | ⟨2, _⟩ => k3_8_7_2
  | ⟨3, _⟩ => k3_8_7_3
  | ⟨4, _⟩ => k3_8_7_4
  | ⟨5, _⟩ => k3_8_7_5
  | ⟨6, _⟩ => k3_8_7_6
  | ⟨7, _⟩ => k3_8_7_7
  | ⟨8, _⟩ => k3_8_7_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k1_1 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨1, by decide⟩) := fun _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k2_1_0
  | ⟨1, _⟩ => k2_1_1
  | ⟨2, _⟩ => k2_1_2
  | ⟨3, _⟩ => k2_1_3
  | ⟨4, _⟩ => k2_1_4
  | ⟨5, _⟩ => k2_1_5
  | ⟨6, _⟩ => k2_1_6
  | ⟨7, _⟩ => k2_1_7
  | ⟨8, _⟩ => k2_1_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k1_2 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨2, by decide⟩) := fun _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k2_2_0
  | ⟨1, _⟩ => k2_2_1
  | ⟨2, _⟩ => k2_2_2
  | ⟨3, _⟩ => k2_2_3
  | ⟨4, _⟩ => k2_2_4
  | ⟨5, _⟩ => k2_2_5
  | ⟨6, _⟩ => k2_2_6
  | ⟨7, _⟩ => k2_2_7
  | ⟨8, _⟩ => k2_2_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k1_3 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨3, by decide⟩) := fun _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k2_3_0
  | ⟨1, _⟩ => k2_3_1
  | ⟨2, _⟩ => k2_3_2
  | ⟨3, _⟩ => k2_3_3
  | ⟨4, _⟩ => k2_3_4
  | ⟨5, _⟩ => k2_3_5
  | ⟨6, _⟩ => k2_3_6
  | ⟨7, _⟩ => k2_3_7
  | ⟨8, _⟩ => k2_3_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k1_4 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨4, by decide⟩) := fun _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k2_4_0
  | ⟨1, _⟩ => k2_4_1
  | ⟨2, _⟩ => k2_4_2
  | ⟨3, _⟩ => k2_4_3
  | ⟨4, _⟩ => k2_4_4
  | ⟨5, _⟩ => k2_4_5
  | ⟨6, _⟩ => k2_4_6
  | ⟨7, _⟩ => k2_4_7
  | ⟨8, _⟩ => k2_4_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k1_5 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨5, by decide⟩) := fun _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k2_5_0
  | ⟨1, _⟩ => k2_5_1
  | ⟨2, _⟩ => k2_5_2
  | ⟨3, _⟩ => k2_5_3
  | ⟨4, _⟩ => k2_5_4
  | ⟨5, _⟩ => k2_5_5
  | ⟨6, _⟩ => k2_5_6
  | ⟨7, _⟩ => k2_5_7
  | ⟨8, _⟩ => k2_5_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k1_6 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨6, by decide⟩) := fun _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k2_6_0
  | ⟨1, _⟩ => k2_6_1
  | ⟨2, _⟩ => k2_6_2
  | ⟨3, _⟩ => k2_6_3
  | ⟨4, _⟩ => k2_6_4
  | ⟨5, _⟩ => k2_6_5
  | ⟨6, _⟩ => k2_6_6
  | ⟨7, _⟩ => k2_6_7
  | ⟨8, _⟩ => k2_6_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k1_7 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨7, by decide⟩) := fun _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k2_7_0
  | ⟨1, _⟩ => k2_7_1
  | ⟨2, _⟩ => k2_7_2
  | ⟨3, _⟩ => k2_7_3
  | ⟨4, _⟩ => k2_7_4
  | ⟨5, _⟩ => k2_7_5
  | ⟨6, _⟩ => k2_7_6
  | ⟨7, _⟩ => k2_7_7
  | ⟨8, _⟩ => k2_7_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem k1_8 : K1 (ofFin ⟨0, by decide⟩) (ofFin ⟨8, by decide⟩) := fun _ y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k2_8_0
  | ⟨1, _⟩ => k2_8_1
  | ⟨2, _⟩ => k2_8_2
  | ⟨3, _⟩ => k2_8_3
  | ⟨4, _⟩ => k2_8_4
  | ⟨5, _⟩ => k2_8_5
  | ⟨6, _⟩ => k2_8_6
  | ⟨7, _⟩ => k2_8_7
  | ⟨8, _⟩ => k2_8_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

theorem no_order_zero : ∀ x1 ∈ LL, K1 (ofFin ⟨0, by decide⟩) x1 := fun y _ =>
  match (y : Fin 9) with
  | ⟨0, _⟩ => k1_0
  | ⟨1, _⟩ => k1_1
  | ⟨2, _⟩ => k1_2
  | ⟨3, _⟩ => k1_3
  | ⟨4, _⟩ => k1_4
  | ⟨5, _⟩ => k1_5
  | ⟨6, _⟩ => k1_6
  | ⟨7, _⟩ => k1_7
  | ⟨8, _⟩ => k1_8
  | ⟨n + 9, h⟩ => absurd h (by omega)

end F9

end SunLean

#print axioms SunLean.F9.no_order_zero
