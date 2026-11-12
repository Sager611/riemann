import Riemann.Deep.Orders.Lean.Defs

/-! # Exhaustive search over reachable collections, proved sound -/
namespace Orders

/-- The result of a o b if it is a natural number (− positive, ÷ exact). -/
def Op.res : Op → Nat → Nat → Option Nat
  | .add, a, b => some (a + b)
  | .mul, a, b => some (a * b)
  | .sub, a, b => if b < a then some (a - b) else none
  | .div, a, b => if 0 < b ∧ a % b = 0 then some (a / b) else none

theorem Op.res_of_rel {o : Op} {a b v : Nat} (h : o.Rel a b v) (hv : 1 ≤ v) : o.res a b = some v := by
  cases o with
  | add => simp only [Op.Rel] at h; simp [Op.res, h]
  | mul => simp only [Op.Rel] at h; simp [Op.res, h]
  | sub =>
    simp only [Op.Rel] at h
    simp only [Op.res, if_pos (show b < a by omega)]
    congr 1; omega
  | div =>
    simp only [Op.Rel] at h
    obtain ⟨hb, rfl⟩ := h
    have h1 : v * b % b = 0 := Nat.mul_mod_left v b
    simp only [Op.res, if_pos (And.intro hb h1), Nat.mul_div_cancel v hb]

/-- All results of o on pairs from S. -/
def cands (o : Op) (S : List Nat) : List Nat := S.flatMap (fun a => S.filterMap (fun b => o.res a b))

theorem mem_cands {o : Op} {S : List Nat} {a b v : Nat} (ha : a ∈ S) (hb : b ∈ S)
    (h : o.res a b = some v) : v ∈ cands o S := by
  simp only [cands, List.mem_flatMap, List.mem_filterMap]
  exact ⟨a, ha, b, hb, h⟩

/-- Sorted insertion (a canonical form, so that dedup merges equal collections). -/
def ins (v : Nat) : List Nat → List Nat
  | [] => [v]
  | y :: ys => if v ≤ y then v :: y :: ys else y :: ins v ys

theorem mem_ins (v y : Nat) : ∀ S : List Nat, y ∈ ins v S ↔ y = v ∨ y ∈ S
  | [] => by simp [ins]
  | z :: zs => by
    unfold ins
    split
    · simp
    · simp only [List.mem_cons, mem_ins v y zs]
      constructor
      · rintro (h | h | h)
        · exact Or.inr (Or.inl h)
        · exact Or.inl h
        · exact Or.inr (Or.inr h)
      · rintro (h | h | h)
        · exact Or.inr (Or.inl h)
        · exact Or.inl h
        · exact Or.inr (Or.inr h)

/-- Remove repeated collections. -/
def dedup : List (List Nat) → List (List Nat)
  | [] => []
  | x :: xs => if (dedup xs).contains x then dedup xs else x :: dedup xs

theorem mem_dedup {x : List Nat} : ∀ {l : List (List Nat)}, x ∈ l → x ∈ dedup l
  | [], h => by simp at h
  | y :: ys, h => by
    unfold dedup
    rcases List.mem_cons.1 h with rfl | h
    · split
      · next hc => exact List.contains_iff_mem.1 hc
      · exact List.mem_cons_self ..
    · split
      · exact mem_dedup h
      · exact List.mem_cons_of_mem _ (mem_dedup h)

/-- One move: from each collection S, every new result v with 1 ≤ v and ok v. -/
def step (σ : Order) (ok : Nat → Bool) (i : Nat) (L : List (List Nat)) : List (List Nat) :=
  dedup (L.flatMap (fun S => ((cands (σ.op i) S).filter
    (fun v => decide (1 ≤ v) && ok v && !memB v S)).map (fun v => ins v S)))

/-- Collections reachable in n more moves, the first with index i. -/
def reach (σ : Order) (ok : Nat → Bool) : Nat → List (List Nat) → Nat → List (List Nat)
  | _, L, 0 => L
  | i, L, n + 1 => reach σ ok (i + 1) (step σ ok i L) n


theorem Valid.prefix {σ : Order} {N : Nat} {l1 : List Move} :
    ∀ {l2 : List Move} {S : List Nat} {i : Nat}, Valid σ N S i (l1 ++ l2) → Valid σ N S i l1 := by
  induction l1 with
  | nil => intros; trivial
  | cons x xs ih =>
    intro l2 S i h
    obtain ⟨ha, hb, hr, hn, h1, h2, h⟩ := h
    exact ⟨ha, hb, hr, hn, h1, h2, ih h⟩

/-- Soundness of the search: the collection made by any legal prefix is (up to order) reached. -/
theorem reach_sound (σ : Order) (ok : Nat → Bool) (N : Nat)
    (hok : ∀ v, 1 ≤ v → v ≤ N → ok v = true) :
    ∀ (ms : List Move) (i : Nat) (L : List (List Nat)) (S0 S : List Nat),
    S ∈ L → (∀ y, y ∈ S ↔ y ∈ S0) → Valid σ N S0 i ms →
    ∃ S' ∈ reach σ ok i L ms.length, ∀ y, y ∈ S' ↔ y ∈ (ms.map Move.v).reverse ++ S0
  | [], i, L, S0, S, hS, hSS, _ => ⟨S, hS, by simpa using hSS⟩
  | m :: ms, i, L, S0, S, hS, hSS, hv => by
    obtain ⟨ha, hb, hr, hn, h1, h2, hv⟩ := hv
    have hstep : ins m.v S ∈ step σ ok i L := by
      apply mem_dedup
      simp only [List.mem_flatMap, List.mem_map, List.mem_filter]
      refine ⟨S, hS, m.v, ⟨mem_cands ((hSS _).2 ha) ((hSS _).2 hb) (Op.res_of_rel hr h1), ?_⟩, rfl⟩
      simp only [Bool.and_eq_true, decide_eq_true_eq, Bool.not_eq_true', memB_false]
      exact ⟨⟨h1, hok _ h1 h2⟩, fun h => hn ((hSS _).1 h)⟩
    obtain ⟨S', hS', hS'eq⟩ := reach_sound σ ok N hok ms (i + 1) (step σ ok i L) (m.v :: S0)
      (ins m.v S) hstep (fun y => by rw [mem_ins, List.mem_cons, hSS]) hv
    refine ⟨S', hS', fun y => ?_⟩
    rw [hS'eq]
    simp only [List.map_cons, List.reverse_cons, List.mem_append, List.mem_reverse, List.mem_cons,
      List.mem_singleton, List.not_mem_nil, or_false]
    exact or_assoc.symm

/-- No perfect sequence for N: no collection reachable in N − 1 moves with results ≤ N
contains all of 1..N. -/
def noPerfect (σ : Order) (N : Nat) : Bool :=
  !((reach σ (fun v => decide (v ≤ N)) 0 [[2]] (N - 1)).any (fun S => allB (fun x => memB x S) N))

theorem noPerfect_sound {σ : Order} {N : Nat} (h : noPerfect σ N = true) :
    ∀ ms, ¬ Perfect σ N ms := by
  rintro ms ⟨hl, hv, hc⟩
  obtain ⟨S', hS', hS'eq⟩ := reach_sound σ (fun v => decide (v ≤ N)) N
    (fun v _ h2 => decide_eq_true h2) ms 0 [[2]] [2] [2] (List.mem_singleton_self _)
    (fun _ => Iff.rfl) hv
  have hlen : ms.length = N - 1 := by omega
  rw [hlen] at hS'
  simp only [noPerfect, Bool.not_eq_true', List.any_eq_false] at h
  apply h S' hS'
  rw [allB_iff]
  intro x h1 h2
  rw [memB_iff, hS'eq]
  have := hc x h1 h2
  simp only [made, List.mem_cons] at this
  simp only [List.mem_append, List.mem_reverse, List.mem_singleton]
  rcases this with h | h
  · exact Or.inr h
  · exact Or.inl h

/-- If no collection is reachable in k moves (no bound on the numbers), a perfect sequence
for N has fewer than k moves, so N ≤ k. -/
theorem bound_of_reach {σ : Order} {k : Nat} (h : reach σ (fun _ => true) 0 [[2]] k = []) :
    ∀ N ms, Perfect σ N ms → N ≤ k := by
  intro N ms hp
  obtain ⟨hl, hv, _⟩ := hp
  apply Classical.byContradiction
  intro hN
  have hk : k ≤ ms.length := by omega
  rw [← List.take_append_drop k ms] at hv
  obtain ⟨S', hS', _⟩ := reach_sound σ (fun _ => true) N (fun _ _ _ => rfl) (ms.take k) 0 [[2]]
    [2] [2] (List.mem_singleton_self _) (fun _ => Iff.rfl) (Valid.prefix hv)
  rw [List.length_take, Nat.min_eq_left hk, h] at hS'
  simp at hS'

end Orders

