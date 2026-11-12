/-!
# Orders of + × − ÷: rules, perfect sequences, a Bool checker

arXiv 2610.10095 (v2 l. 70–86 and 101 = v1 l. 110–129 and 151–153). Start with the number 2.
Move m = 1, 2, 3, ... applies σ_k with k ≡ m (mod 4) to two numbers already made (possibly the
same one twice); a division must give a positive integer; the result joins the collection.
A sequence is perfect for N if it makes all of 1..N in exactly N − 1 moves; then every move makes
a previously unmade number in [1, N].
-/
namespace Orders

/-- The four operations. -/
inductive Op where
  | add | mul | sub | div
  deriving DecidableEq, Repr

/-- `o.Rel a b v`: v is the result of a o b (operands in this order). Numbers are naturals:
a − b = v means v + b = a, a ÷ b = v means b > 0 and v · b = a (exact). Positivity of v is
part of `Valid` (every result lies in [1, N]). -/
def Op.Rel : Op → Nat → Nat → Nat → Prop
  | .add, a, b, v => v = a + b
  | .mul, a, b, v => v = a * b
  | .sub, a, b, v => v + b = a
  | .div, a, b, v => 0 < b ∧ v * b = a

/-- An order σ = σ1 σ2 σ3 σ4 of the operations. -/
structure Order where
  s1 : Op
  s2 : Op
  s3 : Op
  s4 : Op
  deriving DecidableEq, Repr

/-- σ is a permutation of + × − ÷ (the four entries are distinct). -/
def Order.IsPerm (σ : Order) : Prop :=
  σ.s1 ≠ σ.s2 ∧ σ.s1 ≠ σ.s3 ∧ σ.s1 ≠ σ.s4 ∧ σ.s2 ≠ σ.s3 ∧ σ.s2 ≠ σ.s4 ∧ σ.s3 ≠ σ.s4

/-- σ_(k+1) for k = 0, 1, 2, 3. -/
def Order.get (σ : Order) : Nat → Op
  | 0 => σ.s1
  | 1 => σ.s2
  | 2 => σ.s3
  | _ => σ.s4

/-- The operation of move m = i + 1 (0-based index i): σ_k with k ≡ m (mod 4). -/
def Order.op (σ : Order) (i : Nat) : Op := σ.get (i % 4)

/-- A move: operands a, b (already made) and result v. -/
structure Move where
  a : Nat
  b : Nat
  v : Nat
  deriving DecidableEq, Repr

/-- `Valid σ N S i ms`: S lists the numbers made so far and the first move of ms has 0-based
index i; every move uses two numbers already made and the operation of its index, and makes a
new number in [1, N]. -/
def Valid (σ : Order) (N : Nat) : List Nat → Nat → List Move → Prop
  | _, _, [] => True
  | S, i, x :: xs => x.a ∈ S ∧ x.b ∈ S ∧ (σ.op i).Rel x.a x.b x.v ∧ x.v ∉ S ∧
      1 ≤ x.v ∧ x.v ≤ N ∧ Valid σ N (x.v :: S) (i + 1) xs

/-- The numbers made by ms from the start {2}. -/
def made (ms : List Move) : List Nat := 2 :: ms.map Move.v

/-- ms is perfect for N under σ: N − 1 moves (so N ≥ 1) from {2}, each legal and new in
[1, N], and every number 1..N is made. -/
def Perfect (σ : Order) (N : Nat) (ms : List Move) : Prop :=
  ms.length + 1 = N ∧ Valid σ N [2] 0 ms ∧ ∀ x, 1 ≤ x → x ≤ N → x ∈ made ms


/-! ## A computable checker -/

/-- Membership test. -/
def memB (x : Nat) : List Nat → Bool
  | [] => false
  | y :: ys => (x == y) || memB x ys

theorem memB_iff (x : Nat) : ∀ S : List Nat, memB x S = true ↔ x ∈ S
  | [] => by simp [memB]
  | y :: ys => by simp [memB, memB_iff x ys]

theorem memB_false (x : Nat) (S : List Nat) : memB x S = false ↔ x ∉ S := by
  rw [← memB_iff]; simp

def Op.relB : Op → Nat → Nat → Nat → Bool
  | .add, a, b, v => v == a + b
  | .mul, a, b, v => v == a * b
  | .sub, a, b, v => v + b == a
  | .div, a, b, v => decide (0 < b) && v * b == a

theorem Op.relB_iff (o : Op) (a b v : Nat) : o.relB a b v = true ↔ o.Rel a b v := by
  cases o <;> simp [Op.relB, Op.Rel]

def validB (σ : Order) (N : Nat) : List Nat → Nat → List Move → Bool
  | _, _, [] => true
  | S, i, x :: xs => memB x.a S && memB x.b S && (σ.op i).relB x.a x.b x.v && !memB x.v S &&
      decide (1 ≤ x.v) && decide (x.v ≤ N) && validB σ N (x.v :: S) (i + 1) xs

theorem validB_iff (σ : Order) (N : Nat) :
    ∀ (S : List Nat) (i : Nat) (ms : List Move), validB σ N S i ms = true ↔ Valid σ N S i ms
  | _, _, [] => by simp [validB, Valid]
  | S, i, x :: xs => by
      simp only [validB, Valid, Bool.and_eq_true, memB_iff, Op.relB_iff, Bool.not_eq_true',
        memB_false, decide_eq_true_eq, validB_iff σ N (x.v :: S) (i + 1) xs, and_assoc]

/-- f 1 && f 2 && ... && f n. -/
def allB (f : Nat → Bool) : Nat → Bool
  | 0 => true
  | n + 1 => f (n + 1) && allB f n

theorem allB_iff (f : Nat → Bool) : ∀ n, allB f n = true ↔ ∀ x, 1 ≤ x → x ≤ n → f x = true
  | 0 => by
      simp only [allB, true_iff]
      intro x h1 h2; omega
  | n + 1 => by
      simp only [allB, Bool.and_eq_true, allB_iff f n]
      constructor
      · rintro ⟨h1, h2⟩ x hx1 hx2
        by_cases hx : x = n + 1
        · subst hx; exact h1
        · exact h2 x hx1 (by omega)
      · intro h
        exact ⟨h (n + 1) (by omega) (by omega), fun x hx1 hx2 => h x hx1 (by omega)⟩

/-- The checker: exactly N − 1 moves, all legal and new in [1, N], and 1..N all made. -/
def check (σ : Order) (N : Nat) (ms : List Move) : Bool :=
  (ms.length + 1 == N) && validB σ N [2] 0 ms && allB (fun x => memB x (made ms)) N

/-- The checker is exact: `check σ N ms = true` iff ms is a perfect sequence for N. -/
theorem check_iff (σ : Order) (N : Nat) (ms : List Move) :
    check σ N ms = true ↔ Perfect σ N ms := by
  simp only [check, Perfect, Bool.and_eq_true, beq_iff_eq, validB_iff, allB_iff, memB_iff, and_assoc]

theorem Order.op_congr (σ : Order) {i j : Nat} (h : i % 4 = j % 4) : σ.op i = σ.op j := by
  simp [Order.op, h]

end Orders

