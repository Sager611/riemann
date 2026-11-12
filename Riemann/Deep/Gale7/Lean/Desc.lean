import Riemann.Deep.Gale7.Lean.Lemma
import Mathlib.Algebra.Order.Field.Rat

/-!
# Gale7: descriptors (rational, no reals)

A pair of vertices of `H(m₁+1, m₂+1)` has the descriptor `(orb v, orb w, [ix1 v = ix1 w], [ix2 v = ix2 w])`.
`clsQ`: the class value (squared distance of Vₙ, NOTE §2) of a descriptor; `adjD`: the 17 edge classes of `H`
(NOTE §3); `neD`: descriptors of distinct vertices. `adj_adjD`, `ne_neD` connect them to `Adj` and `≠`.
-/

namespace Gale7

open Vx

variable {m₁ m₂ : ℕ}

/-- first index (0 for `Q j`, whose first-block coefficient is 0) -/
def ix1 : Vx (m₁ + 1) (m₂ + 1) → Fin (m₁ + 1)
  | P i => i
  | Q _ => 0
  | A i _ => i
  | B i _ => i
  | C i _ => i

/-- second index (0 for `P i`) -/
def ix2 : Vx (m₁ + 1) (m₂ + 1) → Fin (m₂ + 1)
  | P _ => 0
  | Q j => j
  | A _ j => j
  | B _ j => j
  | C _ j => j

/-- the five orbits of Vₙ -/
inductive Orb : Type where
  | P | Q | A | B | C
  deriving DecidableEq

instance : Fintype Orb := Fintype.ofList [.P, .Q, .A, .B, .C] (by intro x; cases x <;> simp)

/-- orbit of a vertex -/
def orb {p q : ℕ} : Vx p q → Orb
  | P _ => .P
  | Q _ => .Q
  | A _ _ => .A
  | B _ _ => .B
  | C _ _ => .C

/-- first-block coefficient of an orbit: (r_A, 0, a, −c, e) -/
def k1Q (rA a c e : ℚ) : Orb → ℚ
  | .P => rA
  | .Q => 0
  | .A => a
  | .B => -c
  | .C => e

/-- second-block coefficient of an orbit: (0, r_B, −b, d, f) -/
def k2Q (rB b d f : ℚ) : Orb → ℚ
  | .P => 0
  | .Q => rB
  | .A => -b
  | .B => d
  | .C => f

/-- class value: squared distance between points of orbits `o`, `o'` whose first indices agree iff `e₁`
and second indices agree iff `e₂` (NOTE §2; family.py `class_d2`) -/
def clsQ (m₁ m₂ : ℕ) (rA rB a b c d e f : ℚ) (o o' : Orb) (e₁ e₂ : Bool) : ℚ :=
  (k1Q rA a c e o ^ 2 + k1Q rA a c e o' ^ 2
      - 2 * k1Q rA a c e o * k1Q rA a c e o' * (if e₁ then 1 else -1 / (m₁ : ℚ)))
    + (k2Q rB b d f o ^ 2 + k2Q rB b d f o' ^ 2
      - 2 * k2Q rB b d f o * k2Q rB b d f o' * (if e₂ then 1 else -1 / (m₂ : ℚ)))

/-- the 17 edge classes of `H` (NOTE §3; family.py `EDGES`), descriptor form -/
def adjD : Orb → Orb → Bool → Bool → Bool
  | .P, .P, e, _ => !e
  | .P, .A, e, _ => !e
  | .P, .B, e, _ => e
  | .P, .C, e, _ => !e
  | .Q, .Q, _, e => !e
  | .Q, .A, _, e => e
  | .Q, .B, _, e => !e
  | .Q, .C, _, e => !e
  | .A, .P, e, _ => !e
  | .A, .Q, _, e => e
  | .A, .A, e, _ => !e
  | .A, .C, e, _ => !e
  | .B, .P, e, _ => e
  | .B, .Q, _, e => !e
  | .B, .B, _, e => !e
  | .B, .C, _, e => !e
  | .C, .P, e, _ => !e
  | .C, .Q, _, e => !e
  | .C, .A, e, _ => !e
  | .C, .B, _, e => !e
  | .C, .C, e, e' => !e && !e'
  | _, _, _, _ => false

theorem adj_adjD {v w : Vx (m₁ + 1) (m₂ + 1)} (h : Adj v w) :
    adjD (orb v) (orb w) (decide (ix1 v = ix1 w)) (decide (ix2 v = ix2 w)) = true := by
  cases v <;> cases w <;> simp only [Adj] at h <;> simp only [adjD, orb, ix1, ix2] <;>
    first
    | (simp [h]; done)
    | (simp [Ne.symm h]; done)

/-- descriptors of pairs of distinct vertices -/
def neD : Orb → Orb → Bool → Bool → Bool
  | .P, .P, e, _ => !e
  | .Q, .Q, _, e => !e
  | .A, .A, e, e' => !(e && e')
  | .B, .B, e, e' => !(e && e')
  | .C, .C, e, e' => !(e && e')
  | _, _, _, _ => true

theorem eq_of_neD {m₁ m₂ : ℕ} {v w : Vx (m₁ + 1) (m₂ + 1)}
    (h : neD (orb v) (orb w) (decide (ix1 v = ix1 w)) (decide (ix2 v = ix2 w)) = false) : v = w := by
  cases v <;> cases w <;> simp_all [neD, orb, ix1, ix2] <;>
    first
    | exact of_decide_eq_true h
    | exact ⟨of_decide_eq_true h.1, of_decide_eq_true h.2⟩

theorem ne_neD {m₁ m₂ : ℕ} {v w : Vx (m₁ + 1) (m₂ + 1)} (h : v ≠ w) :
    neD (orb v) (orb w) (decide (ix1 v = ix1 w)) (decide (ix2 v = ix2 w)) = true := by
  cases hb : neD (orb v) (orb w) (decide (ix1 v = ix1 w)) (decide (ix2 v = ix2 w))
  · exact absurd (eq_of_neD hb) h
  · rfl

end Gale7
