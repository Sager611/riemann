import Riemann.Deep.Orders.Lean.Ext

/-! # The bulk rounds (NOTE.md §4, Lemma 2) -/
namespace Orders

theorem Op.forall_iff {p : Op → Prop} : (∀ o, p o) ↔ p .add ∧ p .mul ∧ p .sub ∧ p .div := by
  constructor
  · intro h; exact ⟨h _, h _, h _, h _⟩
  · rintro ⟨h1, h2, h3, h4⟩ o; cases o <;> assumption

theorem cover_aux : ∀ a b c d o : Op, a ≠ b → a ≠ c → a ≠ d → b ≠ c → b ≠ d → c ≠ d →
    o = a ∨ o = b ∨ o = c ∨ o = d := by
  intro a b c d o
  cases a <;> cases b <;> cases c <;> cases d <;> cases o <;> simp

theorem Order.get_inj {σ : Order} (hσ : σ.IsPerm) {a b : Nat} (ha : a < 4) (hb : b < 4)
    (h : σ.get a = σ.get b) : a = b := by
  obtain ⟨h12, h13, h14, h23, h24, h34⟩ := hσ
  rcases (by omega : a = 0 ∨ a = 1 ∨ a = 2 ∨ a = 3) with rfl | rfl | rfl | rfl <;>
  rcases (by omega : b = 0 ∨ b = 1 ∨ b = 2 ∨ b = 3) with rfl | rfl | rfl | rfl <;>
  simp_all [Order.get]

theorem Order.op_ne {σ : Order} (hσ : σ.IsPerm) {i j : Nat} (h : i % 4 ≠ j % 4) :
    σ.op i ≠ σ.op j := fun h' =>
  h (Order.get_inj hσ (Nat.mod_lt _ (by decide)) (Nat.mod_lt _ (by decide)) h')

/-- Four consecutive moves use all four operations. -/
theorem Order.cover {σ : Order} (hσ : σ.IsPerm) (i : Nat) (o : Op) :
    o = σ.op i ∨ o = σ.op (i + 1) ∨ o = σ.op (i + 2) ∨ o = σ.op (i + 3) := by
  have hc := cover_aux σ.s1 σ.s2 σ.s3 σ.s4 o hσ.1 hσ.2.1 hσ.2.2.1 hσ.2.2.2.1 hσ.2.2.2.2.1
    hσ.2.2.2.2.2
  simp only [Order.op]
  rcases (by omega : i % 4 = 0 ∨ i % 4 = 1 ∨ i % 4 = 2 ∨ i % 4 = 3) with h | h | h | h
  · rw [show (i + 1) % 4 = 1 by omega, show (i + 2) % 4 = 2 by omega,
      show (i + 3) % 4 = 3 by omega, h]
    simp only [Order.get]; exact hc
  · rw [show (i + 1) % 4 = 2 by omega, show (i + 2) % 4 = 3 by omega,
      show (i + 3) % 4 = 0 by omega, h]
    simp only [Order.get]; rcases hc with h' | h' | h' | h' <;> simp [h']
  · rw [show (i + 1) % 4 = 3 by omega, show (i + 2) % 4 = 0 by omega,
      show (i + 3) % 4 = 1 by omega, h]
    simp only [Order.get]; rcases hc with h' | h' | h' | h' <;> simp [h']
  · rw [show (i + 1) % 4 = 0 by omega, show (i + 2) % 4 = 1 by omega,
      show (i + 3) % 4 = 2 by omega, h]
    simp only [Order.get]; rcases hc with h' | h' | h' | h' <;> simp [h']

/-- n_j = 3P − (P mod 2) − 2j, even; n_j / 2 = ⌊3P/2⌋ − j. -/
def nj (P j : Nat) : Nat := 3 * P - P % 2 - 2 * j

/-- Round j of the bulk, move with operation o: × 3·(P−j) = C_j, + C_j + 1, − C_j − 1, ÷ n_j/2. -/
def bulkMove (P j : Nat) : Op → Move
  | .mul => ⟨3, P - j, 3 * (P - j)⟩
  | .add => ⟨3 * (P - j), 1, 3 * (P - j) + 1⟩
  | .sub => ⟨3 * (P - j), 1, 3 * (P - j) - 1⟩
  | .div => ⟨nj P j, 2, nj P j / 2⟩

/-- Round j, its four moves having the operations of indices i, i+1, i+2, i+3. -/
def round (σ : Order) (i P j : Nat) : List Move :=
  [bulkMove P j (σ.op i), bulkMove P j (σ.op (i + 1)), bulkMove P j (σ.op (i + 2)),
    bulkMove P j (σ.op (i + 3))]

/-- Rounds j, j+1, ..., j+k−1. -/
def rounds (σ : Order) (i P : Nat) : Nat → Nat → List Move
  | _, 0 => []
  | j, k + 1 => round σ i P j ++ rounds σ i P (j + 1) k

/-- J = ⌊P/2⌋ − 1 − d, the last bulk round. -/
def lastJ (P d : Nat) : Nat := P / 2 - 1 - d

/-- The bulk: rounds j = 3..J. -/
def bulk (σ : Order) (i P d : Nat) : List Move := rounds σ i P 3 (lastJ P d - 2)

theorem rounds_length (σ : Order) (i P : Nat) : ∀ k j, (rounds σ i P j k).length = 4 * k
  | 0, _ => rfl
  | k + 1, j => by
      simp only [rounds, List.length_append, rounds_length σ i P k (j + 1)]
      simp [round]; omega


theorem Steps.cons' {σ : Order} {S : Nat → Prop} {i : Nat} {x : Move} {xs : List Move}
    (ha : S x.a) (hb : S x.b) (hr : (σ.op i).Rel x.a x.b x.v)
    (h : Steps σ (fun y => y = x.v ∨ S y) (i + 1) xs) : Steps σ S i (x :: xs) := ⟨ha, hb, hr, h⟩

/-- A bulk move other than × is legal once C_j, n_j, 1, 2 are made. -/
theorem bulk_legal {P j : Nat} {o o' : Op} (hoo : o' = o) (ho : o ≠ .mul) {S : Nat → Prop}
    (h1 : S 1) (h2 : S 2) (hC : S (3 * (P - j))) (hn : S (nj P j)) (hev : nj P j % 2 = 0)
    (hpos : 1 ≤ 3 * (P - j)) :
    S (bulkMove P j o).a ∧ S (bulkMove P j o).b ∧
      o'.Rel (bulkMove P j o).a (bulkMove P j o).b (bulkMove P j o).v := by
  subst o'
  cases o with
  | add => exact ⟨hC, h1, rfl⟩
  | mul => exact absurd rfl ho
  | sub => exact ⟨hC, h1, by simp only [bulkMove, Op.Rel]; omega⟩
  | div => exact ⟨hn, h2, by simp only [bulkMove, Op.Rel]; omega⟩

theorem round_steps {σ : Order} (hσ : σ.IsPerm) {i P j : Nat} (hmul : σ.op i = .mul)
    {S : Nat → Prop} {i' : Nat} (hi : i' % 4 = i % 4) (h1 : S 1) (h2 : S 2) (h3 : S 3)
    (hPj : S (P - j)) (hn : S (nj P j)) (hj : j + 1 ≤ P) (hev : nj P j % 2 = 0) :
    Steps σ S i' (round σ i P j) := by
  have hne : ∀ t, 1 ≤ t → t < 4 → σ.op (i + t) ≠ .mul := by
    intro t ht1 ht2 h; rw [← hmul] at h; exact Order.op_ne hσ (by omega) h
  have e0 : σ.op i' = .mul := by rw [σ.op_congr hi]; exact hmul
  have e1 : σ.op (i' + 1) = σ.op (i + 1) := σ.op_congr (by omega)
  have e2 : σ.op (i' + 1 + 1) = σ.op (i + 2) := σ.op_congr (by omega)
  have e3 : σ.op (i' + 1 + 1 + 1) = σ.op (i + 3) := σ.op_congr (by omega)
  have hC : 1 ≤ 3 * (P - j) := by omega
  unfold round
  rw [hmul]
  refine Steps.cons' h3 hPj (by rw [e0]; rfl) ?_
  obtain ⟨a1, b1, r1⟩ := bulk_legal (S := fun y => y = 3 * (P - j) ∨ S y) e1
    (hne 1 (by omega) (by omega)) (Or.inr h1) (Or.inr h2) (Or.inl rfl) (Or.inr hn) hev hC
  refine Steps.cons' a1 b1 r1 ?_
  obtain ⟨a2, b2, r2⟩ := bulk_legal
    (S := fun y => y = (bulkMove P j (σ.op (i + 1))).v ∨ y = 3 * (P - j) ∨ S y) e2
    (hne 2 (by omega) (by omega)) (Or.inr (Or.inr h1)) (Or.inr (Or.inr h2)) (Or.inr (Or.inl rfl))
    (Or.inr (Or.inr hn)) hev hC
  refine Steps.cons' a2 b2 r2 ?_
  obtain ⟨a3, b3, r3⟩ := bulk_legal
    (S := fun y => y = (bulkMove P j (σ.op (i + 2))).v ∨ y = (bulkMove P j (σ.op (i + 1))).v ∨
      y = 3 * (P - j) ∨ S y) e3
    (hne 3 (by omega) (by omega)) (Or.inr (Or.inr (Or.inr h1))) (Or.inr (Or.inr (Or.inr h2)))
    (Or.inr (Or.inr (Or.inl rfl))) (Or.inr (Or.inr (Or.inr hn))) hev hC
  exact Steps.cons' a3 b3 r3 trivial

theorem round_mem {σ : Order} (hσ : σ.IsPerm) (i P j : Nat) (o : Op) :
    (bulkMove P j o).v ∈ (round σ i P j).map Move.v := by
  rcases Order.cover hσ i o with h | h | h | h <;> simp [round, h]

theorem round_vals {σ : Order} (hσ : σ.IsPerm) (i P j y : Nat) :
    y ∈ (round σ i P j).map Move.v ↔
      y = 3 * (P - j) ∨ y = 3 * (P - j) + 1 ∨ y = 3 * (P - j) - 1 ∨ y = nj P j / 2 := by
  constructor
  · intro hy
    have key : ∀ o : Op, y = (bulkMove P j o).v →
        y = 3 * (P - j) ∨ y = 3 * (P - j) + 1 ∨ y = 3 * (P - j) - 1 ∨ y = nj P j / 2 := by
      intro o h; cases o <;> simp only [bulkMove] at h <;> simp [h]
    simp only [round, List.map_cons, List.map_nil, List.mem_cons, List.not_mem_nil,
      or_false] at hy
    rcases hy with h | h | h | h <;> exact key _ h
  · rintro (h | h | h | h) <;> subst h
    · exact round_mem hσ i P j .mul
    · exact round_mem hσ i P j .add
    · exact round_mem hσ i P j .sub
    · exact round_mem hσ i P j .div


theorem rounds_steps {σ : Order} (hσ : σ.IsPerm) {i P : Nat} (hmul : σ.op i = .mul) :
    ∀ (k j : Nat) (S : Nat → Prop) (i' : Nat), i' % 4 = i % 4 → 3 ≤ j → 2 * (j + k) ≤ P →
    (∀ y, 1 ≤ y → y ≤ P → S y) → S (nj P 3) →
    (∀ j', 3 ≤ j' → j' < j → S (3 * (P - j') - 1) ∧ S (3 * (P - j')) ∧ S (3 * (P - j') + 1)) →
    Steps σ S i' (rounds σ i P j k)
  | 0, _, _, _, _, _, _, _, _, _ => trivial
  | k + 1, j, S, i', hi, hj, hk, hS, hn3, hprev => by
    simp only [rounds]
    apply Steps.append
    · refine round_steps hσ hmul hi (hS 1 (by omega) (by omega)) (hS 2 (by omega) (by omega))
        (hS 3 (by omega) (by omega)) (hS _ (by omega) (by omega)) ?_ (by omega)
        (by simp only [nj]; omega)
      by_cases hj3 : j = 3
      · subst hj3; exact hn3
      · obtain ⟨ha, hb, hc⟩ := hprev ((2 * j + 1 + P % 2) / 3) (by omega) (by omega)
        rcases (by simp only [nj]; omega : nj P j = 3 * (P - (2 * j + 1 + P % 2) / 3) - 1 ∨
            nj P j = 3 * (P - (2 * j + 1 + P % 2) / 3) ∨
            nj P j = 3 * (P - (2 * j + 1 + P % 2) / 3) + 1) with h | h | h <;> rw [h] <;> assumption
    · have hlen : (round σ i P j).length = 4 := rfl
      rw [hlen]
      refine rounds_steps hσ hmul k (j + 1) _ (i' + 4) (by omega) (by omega) (by omega)
        (fun y h1 h2 => Or.inr (hS y h1 h2)) (Or.inr hn3) ?_
      intro j' h1 h2
      by_cases hjj : j' = j
      · subst hjj
        refine ⟨Or.inl ?_, Or.inl ?_, Or.inl ?_⟩ <;> rw [round_vals hσ] <;> simp
      · obtain ⟨a, b, c⟩ := hprev j' h1 (by omega)
        exact ⟨Or.inr a, Or.inr b, Or.inr c⟩

theorem rounds_vals {σ : Order} (hσ : σ.IsPerm) (i P : Nat) : ∀ (k j y : Nat),
    y ∈ (rounds σ i P j k).map Move.v ↔ ∃ j', j ≤ j' ∧ j' < j + k ∧
      (y = 3 * (P - j') ∨ y = 3 * (P - j') + 1 ∨ y = 3 * (P - j') - 1 ∨ y = nj P j' / 2)
  | 0, j, y => by
    simp only [rounds, List.map_nil, List.not_mem_nil, false_iff, not_exists, not_and]
    intro j' h1 h2; omega
  | k + 1, j, y => by
    simp only [rounds, List.map_append, List.mem_append, round_vals hσ,
      rounds_vals hσ i P k (j + 1) y]
    constructor
    · rintro (h | ⟨j', h1, h2, h3⟩)
      · exact ⟨j, Nat.le_refl _, by omega, h⟩
      · exact ⟨j', by omega, by omega, h3⟩
    · rintro ⟨j', h1, h2, h3⟩
      by_cases hj : j' = j
      · subst hj; exact Or.inl h3
      · exact Or.inr ⟨j', by omega, by omega, h3⟩

theorem bulkV_inj {P j : Nat} (hsep : nj P j / 2 + 2 ≤ 3 * (P - j)) {o o' : Op}
    (h : (bulkMove P j o).v = (bulkMove P j o').v) : o = o' := by
  cases o <;> cases o' <;> simp only [bulkMove] at h <;> first | rfl | omega

theorem round_nodup {σ : Order} (hσ : σ.IsPerm) (i P j : Nat)
    (hsep : nj P j / 2 + 2 ≤ 3 * (P - j)) : ((round σ i P j).map Move.v).Nodup := by
  have h01 : (bulkMove P j (σ.op i)).v ≠ (bulkMove P j (σ.op (i + 1))).v := fun h =>
    Order.op_ne hσ (i := i) (j := i + 1) (by omega) (bulkV_inj hsep h)
  have h02 : (bulkMove P j (σ.op i)).v ≠ (bulkMove P j (σ.op (i + 2))).v := fun h =>
    Order.op_ne hσ (i := i) (j := i + 2) (by omega) (bulkV_inj hsep h)
  have h03 : (bulkMove P j (σ.op i)).v ≠ (bulkMove P j (σ.op (i + 3))).v := fun h =>
    Order.op_ne hσ (i := i) (j := i + 3) (by omega) (bulkV_inj hsep h)
  have h12 : (bulkMove P j (σ.op (i + 1))).v ≠ (bulkMove P j (σ.op (i + 2))).v := fun h =>
    Order.op_ne hσ (i := i + 1) (j := i + 2) (by omega) (bulkV_inj hsep h)
  have h13 : (bulkMove P j (σ.op (i + 1))).v ≠ (bulkMove P j (σ.op (i + 3))).v := fun h =>
    Order.op_ne hσ (i := i + 1) (j := i + 3) (by omega) (bulkV_inj hsep h)
  have h23 : (bulkMove P j (σ.op (i + 2))).v ≠ (bulkMove P j (σ.op (i + 3))).v := fun h =>
    Order.op_ne hσ (i := i + 2) (j := i + 3) (by omega) (bulkV_inj hsep h)
  simp only [round, List.map_cons, List.map_nil, List.nodup_cons, List.mem_cons,
    List.not_mem_nil, or_false, not_or, List.nodup_nil, and_true]
  exact ⟨⟨h01, h02, h03⟩, ⟨h12, h13⟩, h23, not_false⟩


theorem rounds_nodup {σ : Order} (hσ : σ.IsPerm) (i P : Nat) : ∀ k j, 3 ≤ j → 2 * (j + k) ≤ P →
    ((rounds σ i P j k).map Move.v).Nodup
  | 0, _, _, _ => by simp [rounds]
  | k + 1, j, hj, hk => by
    simp only [rounds, List.map_append]
    rw [List.nodup_append]
    refine ⟨round_nodup hσ i P j (by simp only [nj]; omega),
      rounds_nodup hσ i P k (j + 1) (by omega) (by omega), ?_⟩
    intro a ha b hb hab
    rw [round_vals hσ] at ha
    rw [rounds_vals hσ] at hb
    obtain ⟨j', h1, h2, h3⟩ := hb
    subst hab
    simp only [nj] at ha h3
    omega

/-- Lemma 2 (values): the bulk makes exactly X = [P+1+d, H−3] and
Y = [H+e+1+3d, 3P−8], H = ⌊3P/2⌋, e = 1 + 2·(P mod 2). -/
theorem bulk_vals {σ : Order} (hσ : σ.IsPerm) (i P d y : Nat) (hJ : 6 + 2 * d ≤ P) :
    y ∈ (bulk σ i P d).map Move.v ↔
      (P + 1 + d ≤ y ∧ y ≤ 3 * P / 2 - 3) ∨
      (3 * P / 2 + 2 + 2 * (P % 2) + 3 * d ≤ y ∧ y ≤ 3 * P - 8) := by
  unfold bulk lastJ
  rw [rounds_vals hσ]
  constructor
  · rintro ⟨j', h1, h2, h3⟩
    simp only [nj] at h3
    omega
  · rintro (⟨h1, h2⟩ | ⟨h1, h2⟩)
    · exact ⟨3 * P / 2 - y, by omega, by omega,
        Or.inr (Or.inr (Or.inr (by simp only [nj]; omega)))⟩
    · refine ⟨(3 * P - y + 1) / 3, by omega, by omega, ?_⟩
      omega

theorem bulk_nodup {σ : Order} (hσ : σ.IsPerm) (i P d : Nat) (hJ : 6 + 2 * d ≤ P) :
    ((bulk σ i P d).map Move.v).Nodup :=
  rounds_nodup hσ i P _ 3 (Nat.le_refl 3) (by unfold lastJ; omega)

theorem bulk_length (σ : Order) (i P d : Nat) :
    (bulk σ i P d).length = 4 * (lastJ P d - 2) := rounds_length σ i P _ 3

/-- Lemma 2 (legality): if move i is a multiplication and [1, P] and n_3 are made, every bulk
move is legal. -/
theorem bulk_steps {σ : Order} (hσ : σ.IsPerm) {i P d : Nat} (hmul : σ.op i = .mul)
    (hJ : 6 + 2 * d ≤ P) {S : Nat → Prop} {i' : Nat} (hi : i' % 4 = i % 4)
    (hS : ∀ y, 1 ≤ y → y ≤ P → S y) (hn3 : S (nj P 3)) : Steps σ S i' (bulk σ i P d) :=
  rounds_steps hσ hmul _ 3 S i' hi (Nat.le_refl 3) (by unfold lastJ; omega) hS hn3
    (fun j' h1 h2 => absurd h2 (by omega))

end Orders

