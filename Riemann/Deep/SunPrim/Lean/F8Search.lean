import Riemann.Deep.SunPrim.Lean.F8

/-!
# Search tree for circular orders of F₈ with the sums `x ^ 2 + y`

`K k x0 … xk`: after choosing `x0, …, xk` (guards: `xk` is new, its sum `x(k-1) ^ 2 + xk` is new), every way
to continue fails; at the last level the closing sum `x7 ^ 2 + x0` repeats. `K1` with `x0 = 0` says no circular
order starting at 0 has pairwise distinct sums. Each level has its own `Decidable` instance.
-/

namespace SunLean

namespace F8

/-- All elements. -/
def LL : List F8 := [(ofFin ⟨0, by decide⟩), (ofFin ⟨1, by decide⟩), (ofFin ⟨2, by decide⟩), (ofFin ⟨3, by decide⟩), (ofFin ⟨4, by decide⟩), (ofFin ⟨5, by decide⟩), (ofFin ⟨6, by decide⟩), (ofFin ⟨7, by decide⟩)]

theorem mem_LL : ∀ x : F8, x ∈ LL := by decide

/-- The table of `x ^ 2 + y`, packed as in `F8Defs`. -/
def STAB : ℕ := 31391579654156835434034099347089511431904655154088029213908744818857637851664

/-- `x ^ 2 + y` by table lookup. -/
def sT (a b : F8) : F8 := ofFin ⟨(STAB >>> (4 * (8 * (a : Fin 8).val + (b : Fin 8).val))) % 16 % 8, Nat.mod_lt _ (by decide)⟩

theorem sT_eq : ∀ a b : F8, sT a b = a ^ 2 + b := by decide +kernel

/-- Level 7. -/
def K7 (x0 x1 x2 x3 x4 x5 x6 x7 : F8) : Prop :=
  x7 ∉ [x6, x5, x4, x3, x2, x1, x0] → sT x6 x7 ∉ [sT x5 x6, sT x4 x5, sT x3 x4, sT x2 x3, sT x1 x2, sT x0 x1] → sT x7 x0 ∉ [sT x6 x7, sT x5 x6, sT x4 x5, sT x3 x4, sT x2 x3, sT x1 x2, sT x0 x1] → False

instance (x0 x1 x2 x3 x4 x5 x6 x7 : F8) : Decidable (K7 x0 x1 x2 x3 x4 x5 x6 x7) := by unfold K7; infer_instance

/-- Level 6. -/
def K6 (x0 x1 x2 x3 x4 x5 x6 : F8) : Prop :=
  x6 ∉ [x5, x4, x3, x2, x1, x0] → sT x5 x6 ∉ [sT x4 x5, sT x3 x4, sT x2 x3, sT x1 x2, sT x0 x1] → ∀ x7 ∈ LL, K7 x0 x1 x2 x3 x4 x5 x6 x7

instance (x0 x1 x2 x3 x4 x5 x6 : F8) : Decidable (K6 x0 x1 x2 x3 x4 x5 x6) := by unfold K6; infer_instance

/-- Level 5. -/
def K5 (x0 x1 x2 x3 x4 x5 : F8) : Prop :=
  x5 ∉ [x4, x3, x2, x1, x0] → sT x4 x5 ∉ [sT x3 x4, sT x2 x3, sT x1 x2, sT x0 x1] → ∀ x6 ∈ LL, K6 x0 x1 x2 x3 x4 x5 x6

instance (x0 x1 x2 x3 x4 x5 : F8) : Decidable (K5 x0 x1 x2 x3 x4 x5) := by unfold K5; infer_instance

/-- Level 4. -/
def K4 (x0 x1 x2 x3 x4 : F8) : Prop :=
  x4 ∉ [x3, x2, x1, x0] → sT x3 x4 ∉ [sT x2 x3, sT x1 x2, sT x0 x1] → ∀ x5 ∈ LL, K5 x0 x1 x2 x3 x4 x5

instance (x0 x1 x2 x3 x4 : F8) : Decidable (K4 x0 x1 x2 x3 x4) := by unfold K4; infer_instance

/-- Level 3. -/
def K3 (x0 x1 x2 x3 : F8) : Prop :=
  x3 ∉ [x2, x1, x0] → sT x2 x3 ∉ [sT x1 x2, sT x0 x1] → ∀ x4 ∈ LL, K4 x0 x1 x2 x3 x4

instance (x0 x1 x2 x3 : F8) : Decidable (K3 x0 x1 x2 x3) := by unfold K3; infer_instance

/-- Level 2. -/
def K2 (x0 x1 x2 : F8) : Prop :=
  x2 ∉ [x1, x0] → sT x1 x2 ∉ [sT x0 x1] → ∀ x3 ∈ LL, K3 x0 x1 x2 x3

instance (x0 x1 x2 : F8) : Decidable (K2 x0 x1 x2) := by unfold K2; infer_instance

/-- Level 1. -/
def K1 (x0 x1 : F8) : Prop :=
  x1 ∉ [x0] → ∀ x2 ∈ LL, K2 x0 x1 x2

instance (x0 x1 : F8) : Decidable (K1 x0 x1) := by unfold K1; infer_instance

end F8

end SunLean
