import Riemann.Deep.SunPrim.Lean.F9

/-!
# Search tree for circular orders of F₉ with the sums `x ^ 3 + y`

`K k x0 … xk`: after choosing `x0, …, xk` (guards: `xk` is new, its sum `x(k-1) ^ 3 + xk` is new), every way
to continue fails; at the last level the closing sum `x8 ^ 3 + x0` repeats. `K1` with `x0 = 0` says no circular
order starting at 0 has pairwise distinct sums. Each level has its own `Decidable` instance.
-/

namespace SunLean

namespace F9

/-- All elements. -/
def LL : List F9 := [(ofFin ⟨0, by decide⟩), (ofFin ⟨1, by decide⟩), (ofFin ⟨2, by decide⟩), (ofFin ⟨3, by decide⟩), (ofFin ⟨4, by decide⟩), (ofFin ⟨5, by decide⟩), (ofFin ⟨6, by decide⟩), (ofFin ⟨7, by decide⟩), (ofFin ⟨8, by decide⟩)]

theorem mem_LL : ∀ x : F9, x ∈ LL := by decide

/-- The table of `x ^ 3 + y`, packed as in `F9Defs`. -/
def STAB : ℕ := 2156537190471354285207355117605577518562014144151612983279107790287310684937668227712976881463824

/-- `x ^ 3 + y` by table lookup. -/
def sT (a b : F9) : F9 := ofFin ⟨(STAB >>> (4 * (9 * (a : Fin 9).val + (b : Fin 9).val))) % 16 % 9, Nat.mod_lt _ (by decide)⟩

theorem sT_eq : ∀ a b : F9, sT a b = a ^ 3 + b := by decide +kernel

/-- Level 8. -/
def K8 (x0 x1 x2 x3 x4 x5 x6 x7 x8 : F9) : Prop :=
  x8 ∉ [x7, x6, x5, x4, x3, x2, x1, x0] → sT x7 x8 ∉ [sT x6 x7, sT x5 x6, sT x4 x5, sT x3 x4, sT x2 x3, sT x1 x2, sT x0 x1] → sT x8 x0 ∉ [sT x7 x8, sT x6 x7, sT x5 x6, sT x4 x5, sT x3 x4, sT x2 x3, sT x1 x2, sT x0 x1] → False

instance (x0 x1 x2 x3 x4 x5 x6 x7 x8 : F9) : Decidable (K8 x0 x1 x2 x3 x4 x5 x6 x7 x8) := by unfold K8; infer_instance

/-- Level 7. -/
def K7 (x0 x1 x2 x3 x4 x5 x6 x7 : F9) : Prop :=
  x7 ∉ [x6, x5, x4, x3, x2, x1, x0] → sT x6 x7 ∉ [sT x5 x6, sT x4 x5, sT x3 x4, sT x2 x3, sT x1 x2, sT x0 x1] → ∀ x8 ∈ LL, K8 x0 x1 x2 x3 x4 x5 x6 x7 x8

instance (x0 x1 x2 x3 x4 x5 x6 x7 : F9) : Decidable (K7 x0 x1 x2 x3 x4 x5 x6 x7) := by unfold K7; infer_instance

/-- Level 6. -/
def K6 (x0 x1 x2 x3 x4 x5 x6 : F9) : Prop :=
  x6 ∉ [x5, x4, x3, x2, x1, x0] → sT x5 x6 ∉ [sT x4 x5, sT x3 x4, sT x2 x3, sT x1 x2, sT x0 x1] → ∀ x7 ∈ LL, K7 x0 x1 x2 x3 x4 x5 x6 x7

instance (x0 x1 x2 x3 x4 x5 x6 : F9) : Decidable (K6 x0 x1 x2 x3 x4 x5 x6) := by unfold K6; infer_instance

/-- Level 5. -/
def K5 (x0 x1 x2 x3 x4 x5 : F9) : Prop :=
  x5 ∉ [x4, x3, x2, x1, x0] → sT x4 x5 ∉ [sT x3 x4, sT x2 x3, sT x1 x2, sT x0 x1] → ∀ x6 ∈ LL, K6 x0 x1 x2 x3 x4 x5 x6

instance (x0 x1 x2 x3 x4 x5 : F9) : Decidable (K5 x0 x1 x2 x3 x4 x5) := by unfold K5; infer_instance

/-- Level 4. -/
def K4 (x0 x1 x2 x3 x4 : F9) : Prop :=
  x4 ∉ [x3, x2, x1, x0] → sT x3 x4 ∉ [sT x2 x3, sT x1 x2, sT x0 x1] → ∀ x5 ∈ LL, K5 x0 x1 x2 x3 x4 x5

instance (x0 x1 x2 x3 x4 : F9) : Decidable (K4 x0 x1 x2 x3 x4) := by unfold K4; infer_instance

/-- Level 3. -/
def K3 (x0 x1 x2 x3 : F9) : Prop :=
  x3 ∉ [x2, x1, x0] → sT x2 x3 ∉ [sT x1 x2, sT x0 x1] → ∀ x4 ∈ LL, K4 x0 x1 x2 x3 x4

instance (x0 x1 x2 x3 : F9) : Decidable (K3 x0 x1 x2 x3) := by unfold K3; infer_instance

/-- Level 2. -/
def K2 (x0 x1 x2 : F9) : Prop :=
  x2 ∉ [x1, x0] → sT x1 x2 ∉ [sT x0 x1] → ∀ x3 ∈ LL, K3 x0 x1 x2 x3

instance (x0 x1 x2 : F9) : Decidable (K2 x0 x1 x2) := by unfold K2; infer_instance

/-- Level 1. -/
def K1 (x0 x1 : F9) : Prop :=
  x1 ∉ [x0] → ∀ x2 ∈ LL, K2 x0 x1 x2

instance (x0 x1 : F9) : Decidable (K1 x0 x1) := by unfold K1; infer_instance

end F9

end SunLean
