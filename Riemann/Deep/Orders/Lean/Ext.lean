import Riemann.Deep.Orders.Lean.Defs

/-! # Extending a perfect sequence: a list of moves whose results are exactly P+1..N -/
namespace Orders

/-- Operand and operation conditions only (no newness, no range), relative to a predicate S
of numbers made so far. -/
def Steps (σ : Order) : (Nat → Prop) → Nat → List Move → Prop
  | _, _, [] => True
  | S, i, x :: xs => S x.a ∧ S x.b ∧ (σ.op i).Rel x.a x.b x.v ∧
      Steps σ (fun y => y = x.v ∨ S y) (i + 1) xs

theorem Steps.mono {σ : Order} {l : List Move} : ∀ {S S' : Nat → Prop} {i : Nat},
    (∀ y, S y → S' y) → Steps σ S i l → Steps σ S' i l := by
  induction l with
  | nil => intros; trivial
  | cons x xs ih =>
    intro S S' i h hs
    obtain ⟨ha, hb, hr, hs⟩ := hs
    exact ⟨h _ ha, h _ hb, hr, ih (fun y hy => hy.elim Or.inl (fun h' => Or.inr (h y h'))) hs⟩

theorem Steps.index {σ : Order} {l : List Move} : ∀ {S : Nat → Prop} {i j : Nat},
    i % 4 = j % 4 → Steps σ S i l → Steps σ S j l := by
  induction l with
  | nil => intros; trivial
  | cons x xs ih =>
    intro S i j h hs
    obtain ⟨ha, hb, hr, hs⟩ := hs
    refine ⟨ha, hb, ?_, ih (by omega) hs⟩
    rw [← σ.op_congr h]; exact hr

theorem Steps.append {σ : Order} {l1 : List Move} : ∀ {l2 : List Move} {S : Nat → Prop} {i : Nat},
    Steps σ S i l1 → Steps σ (fun y => y ∈ l1.map Move.v ∨ S y) (i + l1.length) l2 →
    Steps σ S i (l1 ++ l2) := by
  induction l1 with
  | nil =>
    intro l2 S i _ h2
    simp only [List.nil_append, List.length_nil, Nat.add_zero] at h2 ⊢
    exact Steps.mono (fun y hy => hy.elim (fun h => by simp at h) id) h2
  | cons x xs ih =>
    intro l2 S i h1 h2
    obtain ⟨ha, hb, hr, hs⟩ := h1
    refine ⟨ha, hb, hr, ih hs ?_⟩
    have e : i + (x :: xs).length = i + 1 + xs.length := by simp only [List.length_cons]; omega
    rw [e] at h2
    refine Steps.mono (fun y hy => ?_) h2
    simp only [List.map_cons, List.mem_cons] at hy
    rcases hy with (hy | hy) | hy
    · exact Or.inr (Or.inl hy)
    · exact Or.inl hy
    · exact Or.inr (Or.inr hy)

theorem Valid.results {σ : Order} {N : Nat} {l : List Move} : ∀ {S : List Nat} {i : Nat},
    Valid σ N S i l → ∀ y ∈ l.map Move.v, 1 ≤ y ∧ y ≤ N := by
  induction l with
  | nil => intro S i _ y hy; simp at hy
  | cons x xs ih =>
    intro S i h y hy
    obtain ⟨_, _, _, _, h1, h2, h⟩ := h
    simp only [List.map_cons, List.mem_cons] at hy
    rcases hy with hy | hy
    · subst hy; exact ⟨h1, h2⟩
    · exact ih h y hy

theorem Valid.mono_N {σ : Order} {P N : Nat} (hPN : P ≤ N) {l : List Move} :
    ∀ {S : List Nat} {i : Nat}, Valid σ P S i l → Valid σ N S i l := by
  induction l with
  | nil => intros; trivial
  | cons x xs ih =>
    intro S i h
    obtain ⟨ha, hb, hr, hn, h1, h2, h⟩ := h
    exact ⟨ha, hb, hr, hn, h1, by omega, ih h⟩


theorem Valid.congr {σ : Order} {N : Nat} {l : List Move} : ∀ {S S' : List Nat} {i : Nat},
    (∀ y, y ∈ S ↔ y ∈ S') → Valid σ N S i l → Valid σ N S' i l := by
  induction l with
  | nil => intros; trivial
  | cons x xs ih =>
    intro S S' i hS h
    obtain ⟨ha, hb, hr, hn, h1, h2, h⟩ := h
    refine ⟨(hS _).1 ha, (hS _).1 hb, hr, fun h' => hn ((hS _).2 h'), h1, h2, ih ?_ h⟩
    intro y; simp only [List.mem_cons, hS]

theorem Valid.append {σ : Order} {N : Nat} {l1 : List Move} : ∀ {l2 : List Move} {S : List Nat} {i : Nat},
    Valid σ N S i l1 → Valid σ N ((l1.map Move.v).reverse ++ S) (i + l1.length) l2 →
    Valid σ N S i (l1 ++ l2) := by
  induction l1 with
  | nil => intro l2 S i _ h2; simpa using h2
  | cons x xs ih =>
    intro l2 S i h1 h2
    obtain ⟨ha, hb, hr, hn, hv1, hv2, hs⟩ := h1
    refine ⟨ha, hb, hr, hn, hv1, hv2, ih hs ?_⟩
    have e : i + (x :: xs).length = i + 1 + xs.length := by simp only [List.length_cons]; omega
    rw [e] at h2
    refine Valid.congr (fun y => ?_) h2
    simp only [List.map_cons, List.reverse_cons, List.mem_append, List.mem_reverse, List.mem_cons,
      List.mem_singleton, List.not_mem_nil, or_false]
    exact or_assoc

theorem valid_of_steps {σ : Order} {N : Nat} {l : List Move} :
    ∀ {L : List Nat} {S : Nat → Prop} {i : Nat}, (∀ y, S y ↔ y ∈ L) → Steps σ S i l →
    (l.map Move.v).Nodup → (∀ y ∈ l.map Move.v, y ∉ L ∧ 1 ≤ y ∧ y ≤ N) → Valid σ N L i l := by
  induction l with
  | nil => intros; trivial
  | cons x xs ih =>
    intro L S i hS hs hnd hnew
    obtain ⟨ha, hb, hr, hs⟩ := hs
    simp only [List.map_cons, List.nodup_cons] at hnd
    have hx := hnew x.v (by simp)
    refine ⟨(hS _).1 ha, (hS _).1 hb, hr, hx.1, hx.2.1, hx.2.2, ?_⟩
    refine ih (S := fun y => y = x.v ∨ S y) (fun y => ?_) hs hnd.2 (fun y hy => ?_)
    · simp only [List.mem_cons, hS]
    · have h1 := hnew y (by simp only [List.map_cons, List.mem_cons]; exact Or.inr hy)
      refine ⟨?_, h1.2⟩
      simp only [List.mem_cons, not_or]
      exact ⟨fun h => hnd.1 (h ▸ hy), h1.1⟩

/-- **Extension.** A perfect sequence for P followed by moves whose operands are in [1, P] or
made earlier in the extension, with the right operations, whose results are N − P distinct
numbers exactly filling (P, N], is a perfect sequence for N. -/
theorem perfect_extend {σ : Order} {P N : Nat} {pre ext : List Move}
    (hpre : Perfect σ P pre) (hPN : P ≤ N)
    (hs : Steps σ (fun y => 1 ≤ y ∧ y ≤ P) (P - 1) ext)
    (hnd : (ext.map Move.v).Nodup)
    (hmem : ∀ y, y ∈ ext.map Move.v ↔ P < y ∧ y ≤ N)
    (hlen : ext.length = N - P) :
    Perfect σ N (pre ++ ext) := by
  obtain ⟨hl, hv, hc⟩ := hpre
  have hP2 : 2 ≤ P := by
    apply Classical.byContradiction
    intro h
    have hP1 : P = 1 := by omega
    subst hP1
    have hnil : pre = [] := List.eq_nil_of_length_eq_zero (by omega)
    subst hnil
    have := hc 1 (Nat.le_refl 1) (Nat.le_refl 1)
    simp [made] at this
  have hL : ∀ y, (1 ≤ y ∧ y ≤ P) ↔ y ∈ (pre.map Move.v).reverse ++ [2] := by
    intro y
    simp only [List.mem_append, List.mem_reverse, List.mem_singleton]
    constructor
    · rintro ⟨h1, h2⟩
      have := hc y h1 h2
      simp only [made, List.mem_cons] at this
      exact this.symm
    · rintro (hy | hy)
      · exact hv.results y hy
      · omega
  refine ⟨?_, ?_, ?_⟩
  · simp only [List.length_append]; omega
  · apply Valid.append (hv.mono_N hPN)
    have hidx : 0 + pre.length = P - 1 := by omega
    rw [hidx]
    refine valid_of_steps hL hs hnd (fun y hy => ?_)
    have h := (hmem y).1 hy
    exact ⟨fun h' => by have := (hL y).2 h'; omega, by omega, h.2⟩
  · intro x h1 h2
    simp only [made, List.map_append, List.mem_cons, List.mem_append]
    by_cases hx : x ≤ P
    · have := hc x h1 hx
      simp only [made, List.mem_cons] at this
      rcases this with h | h
      · exact Or.inl h
      · exact Or.inr (Or.inl h)
    · exact Or.inr (Or.inr ((hmem x).2 ⟨by omega, h2⟩))

end Orders

