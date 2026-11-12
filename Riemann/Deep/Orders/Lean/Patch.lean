import Riemann.Deep.Orders.Lean.Bulk

/-! # Patches: symbolic checker (NOTE.md §4, Lemma 3) and its soundness -/
namespace Orders

/-- The value a·q + b (P = 4q + r). -/
structure Lin where
  a : Int
  b : Int
  deriving DecidableEq, Repr

def Lin.ev (x : Lin) (q : Nat) : Int := x.a * q + x.b

def Lin.val (x : Lin) (q : Nat) : Nat := (x.ev q).toNat

/-- A head or tail move of a patch: result t = u ∘ w. -/
structure PMove where
  op : Op
  t : Lin
  u : Lin
  w : Lin
  deriving DecidableEq, Repr

def PMove.toMove (q : Nat) (m : PMove) : Move := ⟨m.u.val q, m.w.val q, m.t.val q⟩

/-- t = u ∘ w as an identity in q (÷: w a constant k ≥ 2 and u = k·t; ×: u a constant). -/
def PMove.identB : PMove → Bool
  | ⟨.add, t, u, w⟩ => t == ⟨u.a + w.a, u.b + w.b⟩
  | ⟨.sub, t, u, w⟩ => u == ⟨t.a + w.a, t.b + w.b⟩
  | ⟨.mul, t, u, w⟩ => u.a == 0 && t == ⟨u.b * w.a, u.b * w.b⟩
  | ⟨.div, t, u, w⟩ => w.a == 0 && decide (2 ≤ w.b) && u == ⟨w.b * t.a, w.b * t.b⟩

/-- u ≤ v for every q ≥ q0 (slope and value at q0). -/
def leB (q0 : Nat) (u v : Lin) : Bool := decide (u.a ≤ v.a) && decide (u.ev q0 ≤ v.ev q0)

def inB (q0 : Nat) (x lo hi : Lin) : Bool := leB q0 lo x && leB q0 x hi

/-- Moves k, k+1, ... of a segment: operation σ.op (k0 + k), identity, and each operand an
earlier result (done) or witnessed by okB. -/
def segB (σ : Order) (k0 : Nat) (okB : Lin → Bool) : List Lin → Nat → List PMove → Bool
  | _, _, [] => true
  | done, k, x :: xs => (x.op == σ.op (k0 + k)) && x.identB &&
      (done.contains x.u || okB x.u) && (done.contains x.w || okB x.w) &&
      segB σ k0 okB (x.t :: done) (k + 1) xs

theorem leB_sound {q0 : Nat} {u v : Lin} (h : leB q0 u v = true) {q : Nat} (hq : q0 ≤ q) :
    u.ev q ≤ v.ev q := by
  simp only [leB, Bool.and_eq_true, decide_eq_true_eq] at h
  obtain ⟨h1, h2⟩ := h
  simp only [Lin.ev] at *
  have hm : 0 ≤ (v.a - u.a) * ((q : Int) - q0) := Int.mul_nonneg (by omega) (by omega)
  have e : (v.a - u.a) * ((q : Int) - q0) = v.a * q - u.a * q - (v.a * q0 - u.a * q0) := by
    rw [Int.sub_mul, Int.mul_sub, Int.mul_sub]; omega
  omega

theorem inB_sound {q0 : Nat} {x lo hi : Lin} (h : inB q0 x lo hi = true) {q : Nat} (hq : q0 ≤ q) :
    lo.ev q ≤ x.ev q ∧ x.ev q ≤ hi.ev q := by
  simp only [inB, Bool.and_eq_true] at h
  exact ⟨leB_sound h.1 hq, leB_sound h.2 hq⟩

theorem identB_sound {m : PMove} (h : m.identB = true) {q : Nat} (hu : 1 ≤ m.u.ev q)
    (hw : 1 ≤ m.w.ev q) (ht : 1 ≤ m.t.ev q) : m.op.Rel (m.u.val q) (m.w.val q) (m.t.val q) := by
  obtain ⟨op, t, u, w⟩ := m
  simp only at hu hw ht ⊢
  cases op with
  | add =>
    simp only [PMove.identB, beq_iff_eq] at h
    subst h
    simp only [Op.Rel, Lin.val, Lin.ev] at *
    rw [Int.add_mul]; omega
  | sub =>
    simp only [PMove.identB, beq_iff_eq] at h
    subst h
    simp only [Op.Rel, Lin.val, Lin.ev] at *
    rw [Int.add_mul] at hu ⊢; omega
  | mul =>
    simp only [PMove.identB, Bool.and_eq_true, beq_iff_eq] at h
    obtain ⟨h1, h2⟩ := h
    subst h2
    simp only [Op.Rel, Lin.val, Lin.ev, h1, Int.zero_mul, Int.zero_add] at *
    rw [Int.mul_assoc, ← Int.mul_add, Int.toNat_mul (by omega) (by omega)]
  | div =>
    simp only [PMove.identB, Bool.and_eq_true, beq_iff_eq, decide_eq_true_eq] at h
    obtain ⟨⟨h1, h2⟩, h3⟩ := h
    subst h3
    simp only [Op.Rel, Lin.val, Lin.ev, h1, Int.zero_mul, Int.zero_add] at *
    have e : w.b * t.a * (q : Int) + w.b * t.b = (t.a * q + t.b) * w.b := by
      rw [Int.mul_assoc, ← Int.mul_add]; exact Int.mul_comm _ _
    refine ⟨by omega, ?_⟩
    rw [e, Int.toNat_mul (by omega) (by omega)]


theorem segB_sound (σ : Order) (q k0 : Nat) (okB : Lin → Bool) :
    ∀ (l : List PMove) (done : List Lin) (k i : Nat) (S : Nat → Prop),
    segB σ k0 okB done k l = true → i % 4 = (k0 + k) % 4 →
    (∀ x, okB x = true → 1 ≤ x.ev q ∧ S (x.val q)) →
    (∀ x ∈ done, 1 ≤ x.ev q ∧ S (x.val q)) →
    (∀ m ∈ l, 1 ≤ m.t.ev q) →
    Steps σ S i (l.map (PMove.toMove q))
  | [], _, _, _, _, _, _, _, _, _ => trivial
  | x :: xs, done, k, i, S, h, hi, hok, hdone, hpos => by
    simp only [segB, Bool.and_eq_true, beq_iff_eq, Bool.or_eq_true, List.contains_iff_mem] at h
    obtain ⟨⟨⟨⟨hop, hid⟩, hu⟩, hw⟩, hrest⟩ := h
    have wit : ∀ y, y ∈ done ∨ okB y = true → 1 ≤ y.ev q ∧ S (y.val q) := by
      rintro y (hy | hy)
      · exact hdone y hy
      · exact hok y hy
    have hu' := wit _ hu
    have hw' := wit _ hw
    have ht := hpos x (List.mem_cons_self ..)
    refine ⟨hu'.2, hw'.2, ?_, ?_⟩
    · rw [σ.op_congr hi, ← hop]
      exact identB_sound hid hu'.1 hw'.1 ht
    · refine segB_sound σ q k0 okB xs (x.t :: done) (k + 1) (i + 1) _ hrest (by omega)
        (fun y hy => ⟨(hok y hy).1, Or.inr (hok y hy).2⟩) ?_
        (fun m hm => hpos m (List.mem_cons_of_mem _ hm))
      intro y hy
      rcases List.mem_cons.1 hy with rfl | hy
      · exact ⟨ht, Or.inl rfl⟩
      · exact ⟨(hdone y hy).1, Or.inr (hdone y hy).2⟩

/-- ⌊3r/2⌋ and e = 1 + 2·(r mod 2) (H = 6q + ⌊3r/2⌋ = ⌊3P/2⌋). -/
def H1 (r : Nat) : Int := ((3 * r / 2 : Nat) : Int)

def eE (r : Nat) : Nat := 1 + 2 * (r % 2)

/-- The pool L ∪ M ∪ U made by head and tail (NOTE.md §4), as linear forms in q. -/
def poolForms (r c d : Nat) : List Lin :=
  (List.range d).map (fun k => ⟨4, (r : Int) + 1 + k⟩) ++
  (List.range (eE r + 3 * d + 3)).map (fun k => ⟨6, H1 r - 2 + k⟩) ++
  (List.range (c + 8)).map (fun k => ⟨12, 3 * (r : Int) - 7 + k⟩)

theorem pool_vals (r c d q : Nat) (hJ : 6 + 2 * d ≤ 4 * q + r) (y : Nat) :
    y ∈ (poolForms r c d).map (fun x => x.val q) ↔
      (4 * q + r < y ∧ y ≤ 4 * q + r + d) ∨
      (3 * (4 * q + r) / 2 - 2 ≤ y ∧ y ≤ 3 * (4 * q + r) / 2 + eE r + 3 * d) ∨
      (3 * (4 * q + r) - 7 ≤ y ∧ y ≤ 3 * (4 * q + r) + c) := by
  simp only [poolForms, List.map_append, List.mem_append, List.mem_map, List.mem_range,
    Lin.val, Lin.ev, H1, eE]
  constructor
  · rintro ((⟨x, ⟨k, hk, rfl⟩, rfl⟩ | ⟨x, ⟨k, hk, rfl⟩, rfl⟩) | ⟨x, ⟨k, hk, rfl⟩, rfl⟩) <;>
      dsimp only <;> omega
  · rintro (h | h | h)
    · exact Or.inl (Or.inl ⟨_, ⟨y - (4 * q + r + 1), by omega, rfl⟩, by dsimp only; omega⟩)
    · exact Or.inl (Or.inr ⟨_, ⟨y - (3 * (4 * q + r) / 2 - 2), by omega, rfl⟩,
        by dsimp only; omega⟩)
    · exact Or.inr ⟨_, ⟨y - (3 * (4 * q + r) - 7), by omega, rfl⟩, by dsimp only; omega⟩

theorem pool_pos (r c d q : Nat) (hJ : 6 + 2 * d ≤ 4 * q + r) :
    ∀ x ∈ poolForms r c d, 1 ≤ x.ev q := by
  intro x hx
  simp only [poolForms, List.mem_append, List.mem_map, List.mem_range] at hx
  rcases hx with (⟨k, hk, rfl⟩ | ⟨k, hk, rfl⟩) | ⟨k, hk, rfl⟩ <;> simp only [Lin.ev, H1] <;>
    omega


theorem block_nodup (A B : Int) (n q : Nat) (hpos : 0 ≤ A * q + B) :
    (((List.range n).map (fun k : Nat => (⟨A, B + (k : Int)⟩ : Lin))).map (fun x => x.val q)).Nodup := by
  rw [List.map_map]
  apply List.pairwise_map.2
  apply List.pairwise_lt_range.imp
  intro a b hab
  simp only [Function.comp, Lin.val, Lin.ev]
  omega

theorem block_mem (A B : Int) (n q y : Nat) (hpos : 0 ≤ A * q + B) :
    y ∈ ((List.range n).map (fun k : Nat => (⟨A, B + (k : Int)⟩ : Lin))).map (fun x => x.val q) ↔
      (A * q + B ≤ (y : Int) ∧ (y : Int) < A * q + B + n) := by
  simp only [List.map_map, List.mem_map, List.mem_range, Function.comp, Lin.val, Lin.ev]
  constructor
  · rintro ⟨k, hk, rfl⟩; omega
  · intro h; exact ⟨(y - (A * q + B)).toNat, by omega, by omega⟩

theorem pool_nodup (r c d q : Nat) (hJ : 6 + 2 * d ≤ 4 * q + r) :
    ((poolForms r c d).map (fun x => x.val q)).Nodup := by
  unfold poolForms
  simp only [List.map_append]
  rw [List.nodup_append, List.nodup_append]
  have p1 : (0 : Int) ≤ 4 * q + ((r : Int) + 1) := by omega
  have p2 : (0 : Int) ≤ 6 * q + (H1 r - 2) := by simp only [H1]; omega
  have p3 : (0 : Int) ≤ 12 * q + (3 * (r : Int) - 7) := by omega
  refine ⟨⟨block_nodup _ _ _ _ p1, block_nodup _ _ _ _ p2, ?_⟩, block_nodup _ _ _ _ p3, ?_⟩
  · intro a ha b hb hab
    subst hab
    rw [block_mem _ _ _ _ _ p1] at ha
    rw [block_mem _ _ _ _ _ p2] at hb
    simp only [H1, eE] at hb
    omega
  · intro a ha b hb hab
    subst hab
    rw [block_mem _ _ _ _ _ p3] at hb
    rcases List.mem_append.1 ha with ha | ha
    · rw [block_mem _ _ _ _ _ p1] at ha; omega
    · rw [block_mem _ _ _ _ _ p2] at ha; simp only [H1, eE] at ha; omega

/-- n_3 = 3P − (P mod 2) − 6 as a form. -/
def n3Form (r : Nat) : Lin := ⟨12, 3 * (r : Int) - ((r % 2 : Nat) : Int) - 6⟩

/-- Head witness: the value lies in [1, P] for every q ≥ q0. -/
def okHead (r q0 : Nat) (x : Lin) : Bool := inB q0 x ⟨0, 1⟩ ⟨4, r⟩

/-- Tail witness: in [1, P], in X = [P+1+d, H−3] or in Y = [H+e+1+3d, 3P−8]. -/
def okTail (r d q0 : Nat) (x : Lin) : Bool :=
  okHead r q0 x || inB q0 x ⟨4, (r : Int) + 1 + d⟩ ⟨6, H1 r - 3⟩ ||
    inB q0 x ⟨6, H1 r + 2 + 2 * ((r % 2 : Nat) : Int) + 3 * d⟩ ⟨12, 3 * (r : Int) - 8⟩

/-- A patch for offset c: d, threshold q*, 12 head moves (rounds 0–2) and the tail. -/
structure Patch where
  d : Nat
  qstar : Nat
  head : List PMove
  tail : List PMove

/-- The conditions of Lemma 3, as a Bool. -/
def patchOK (σ : Order) (r c : Nat) (p : Patch) : Bool :=
  (σ.op (r + 3) == .mul) && (p.head.length == 12) &&
  (p.tail.length == c + 4 * p.d + 2 * (r % 2)) &&
  decide (6 + 2 * p.d ≤ 4 * p.qstar + r) &&
  segB σ (r + 3) (okHead r p.qstar) [] 0 p.head &&
  segB σ (r + 3) (okTail r p.d p.qstar) (p.head.map PMove.t) 0 p.tail &&
  (p.head.map PMove.t).contains (n3Form r) &&
  ((p.head ++ p.tail).map PMove.t).isPerm (poolForms r c p.d)

/-- The extension of a sequence for P = 4q + r: head, bulk rounds 3..J, tail. -/
def ext (σ : Order) (r : Nat) (p : Patch) (q : Nat) : List Move :=
  p.head.map (PMove.toMove q) ++ bulk σ (4 * q + r + 11) (4 * q + r) p.d ++
    p.tail.map (PMove.toMove q)


theorem okHead_sound {r q0 q : Nat} (hq : q0 ≤ q) {x : Lin} (h : okHead r q0 x = true) :
    1 ≤ x.ev q ∧ 1 ≤ x.val q ∧ x.val q ≤ 4 * q + r := by
  have := inB_sound h hq
  simp only [Lin.ev, Lin.val] at this ⊢
  omega

/-- **Lemma 3 (patch).** A checked patch extends every perfect sequence for P = 4q + r,
q ≥ q*, to a perfect sequence for 3P + c. -/
theorem patch_sound {σ : Order} (hσ : σ.IsPerm) {r c : Nat} {p : Patch}
    (hok : patchOK σ r c p = true) {q : Nat} (hq : p.qstar ≤ q) {pre : List Move}
    (hpre : Perfect σ (4 * q + r) pre) :
    Perfect σ (3 * (4 * q + r) + c) (pre ++ ext σ r p q) := by
  simp only [patchOK, Bool.and_eq_true, beq_iff_eq, decide_eq_true_eq, List.contains_iff_mem,
    List.isPerm_iff] at hok
  obtain ⟨⟨⟨⟨⟨⟨⟨hmul0, hhl⟩, htl⟩, hJ0⟩, hsegH⟩, hsegT⟩, hn3⟩, hperm⟩ := hok
  have hJ : 6 + 2 * p.d ≤ 4 * q + r := by omega
  have e1 : (4 * q + r) % 2 = r % 2 := by omega
  have e2 : 3 * (4 * q + r) / 2 = 6 * q + 3 * r / 2 := by omega
  have e3 : (4 * q + r) / 2 = 2 * q + r / 2 := by omega
  have hres : ∀ m ∈ p.head ++ p.tail, 1 ≤ m.t.ev q := fun m hm =>
    pool_pos r c p.d q hJ _ ((hperm.mem_iff).1 (List.mem_map_of_mem hm))
  have eH : (p.head.map (PMove.toMove q)).map Move.v = (p.head.map PMove.t).map (fun x => x.val q) := by
    simp only [List.map_map]; rfl
  have eT : (p.tail.map (PMove.toMove q)).map Move.v = (p.tail.map PMove.t).map (fun x => x.val q) := by
    simp only [List.map_map]; rfl
  have hHT : ((p.head.map (PMove.toMove q)).map Move.v ++ (p.tail.map (PMove.toMove q)).map Move.v).Perm
      ((poolForms r c p.d).map (fun x => x.val q)) := by
    rw [eH, eT, ← List.map_append, ← List.map_append]; exact hperm.map _
  have hlenH : (p.head.map (PMove.toMove q)).length = 12 := by simp [hhl]
  have hlenB := bulk_length σ (4 * q + r + 11) (4 * q + r) p.d
  have hmulB : σ.op (4 * q + r + 11) = .mul := by rw [σ.op_congr (j := r + 3) (by omega)]; exact hmul0
  -- head
  have hH : Steps σ (fun y => 1 ≤ y ∧ y ≤ 4 * q + r) (4 * q + r - 1) (p.head.map (PMove.toMove q)) :=
    segB_sound σ q (r + 3) _ p.head [] 0 _ _ hsegH (by omega)
      (fun x hx => have h := okHead_sound hq hx; ⟨h.1, h.2.1, h.2.2⟩) (by simp)
      (fun m hm => hres m (List.mem_append_left _ hm))
  -- bulk
  have hn3v : nj (4 * q + r) 3 ∈ (p.head.map (PMove.toMove q)).map Move.v := by
    rw [eH]
    have h := List.mem_map_of_mem (f := fun x => x.val q) hn3
    have e : (n3Form r).val q = nj (4 * q + r) 3 := by
      simp only [n3Form, Lin.val, Lin.ev, nj]; omega
    rw [e] at h; exact h
  have hB : Steps σ (fun y => y ∈ (p.head.map (PMove.toMove q)).map Move.v ∨ (1 ≤ y ∧ y ≤ 4 * q + r))
      (4 * q + r - 1 + (p.head.map (PMove.toMove q)).length) (bulk σ (4 * q + r + 11) (4 * q + r) p.d) :=
    bulk_steps hσ hmulB hJ (by rw [hlenH]; omega) (fun y h1 h2 => Or.inr ⟨h1, h2⟩) (Or.inl hn3v)

  -- tail
  have hposH : ∀ x ∈ p.head.map PMove.t, 1 ≤ x.ev q := fun x hx =>
    pool_pos r c p.d q hJ _ ((hperm.mem_iff).1 (by rw [List.map_append]; exact List.mem_append_left _ hx))
  have hT : Steps σ (fun y => y ∈ (p.head.map (PMove.toMove q) ++
        bulk σ (4 * q + r + 11) (4 * q + r) p.d).map Move.v ∨ (1 ≤ y ∧ y ≤ 4 * q + r))
      (4 * q + r - 1 + (p.head.map (PMove.toMove q) ++ bulk σ (4 * q + r + 11) (4 * q + r) p.d).length)
      (p.tail.map (PMove.toMove q)) := by
    refine segB_sound σ q (r + 3) _ p.tail (p.head.map PMove.t) 0 _ _ hsegT ?_ ?_ ?_
      (fun m hm => hres m (List.mem_append_right _ hm))
    · simp only [List.length_append, hlenH, hlenB]; omega
    · intro x hx
      simp only [okTail, Bool.or_eq_true] at hx
      rcases hx with (hx | hx) | hx
      · have h := okHead_sound hq hx; exact ⟨h.1, Or.inr ⟨h.2.1, h.2.2⟩⟩
      · have h := inB_sound hx hq
        simp only [Lin.ev, H1] at h
        refine ⟨by simp only [Lin.ev]; omega, Or.inl ?_⟩
        rw [List.map_append, List.mem_append]; right
        rw [bulk_vals hσ _ _ _ _ hJ]; left
        simp only [Lin.val, Lin.ev]; omega
      · have h := inB_sound hx hq
        simp only [Lin.ev, H1] at h
        refine ⟨by simp only [Lin.ev]; omega, Or.inl ?_⟩
        rw [List.map_append, List.mem_append]; right
        rw [bulk_vals hσ _ _ _ _ hJ]; right
        simp only [Lin.val, Lin.ev]; omega
    · intro x hx
      refine ⟨hposH x hx, Or.inl ?_⟩
      rw [List.map_append, List.mem_append]; left
      rw [eH]; exact List.mem_map_of_mem hx
  have hS : Steps σ (fun y => 1 ≤ y ∧ y ≤ 4 * q + r) (4 * q + r - 1) (ext σ r p q) := by
    unfold ext; exact Steps.append (Steps.append hH hB) hT
  have hnd : ((ext σ r p q).map Move.v).Nodup := by
    unfold ext
    rw [List.map_append, List.map_append]
    have hp : ((p.head.map (PMove.toMove q)).map Move.v ++
        (bulk σ (4 * q + r + 11) (4 * q + r) p.d).map Move.v ++
        (p.tail.map (PMove.toMove q)).map Move.v).Perm
        ((bulk σ (4 * q + r + 11) (4 * q + r) p.d).map Move.v ++
          ((p.head.map (PMove.toMove q)).map Move.v ++ (p.tail.map (PMove.toMove q)).map Move.v)) := by
      rw [List.append_assoc]; exact List.perm_append_comm_assoc _ _ _
    refine hp.nodup_iff.2 ?_
    rw [List.nodup_append]
    refine ⟨bulk_nodup hσ _ _ _ hJ, hHT.nodup_iff.2 (pool_nodup r c p.d q hJ), ?_⟩
    intro a ha b hb hab
    subst hab
    rw [bulk_vals hσ _ _ _ _ hJ] at ha
    have hb' := (hHT.mem_iff).1 hb
    rw [pool_vals r c p.d q hJ] at hb'
    simp only [eE] at hb'
    omega
  have hmem : ∀ y, y ∈ (ext σ r p q).map Move.v ↔ 4 * q + r < y ∧ y ≤ 3 * (4 * q + r) + c := by
    intro y
    unfold ext
    simp only [List.map_append, List.mem_append]
    have hpv := pool_vals r c p.d q hJ y
    have hbv := bulk_vals hσ (4 * q + r + 11) (4 * q + r) p.d y hJ
    have hht := (hHT.mem_iff (a := y))
    simp only [List.mem_append, eE] at hht hpv
    constructor
    · rintro ((h | h) | h)
      · have := hpv.1 (hht.1 (Or.inl h)); omega
      · have := hbv.1 h; omega
      · have := hpv.1 (hht.1 (Or.inr h)); omega
    · intro h
      by_cases hb : (4 * q + r + 1 + p.d ≤ y ∧ y ≤ 3 * (4 * q + r) / 2 - 3) ∨
          (3 * (4 * q + r) / 2 + 2 + 2 * ((4 * q + r) % 2) + 3 * p.d ≤ y ∧ y ≤ 3 * (4 * q + r) - 8)
      · exact Or.inl (Or.inr (hbv.2 hb))
      · rcases hht.2 (hpv.2 (by omega)) with h' | h'
        · exact Or.inl (Or.inl h')
        · exact Or.inr h'
  have hlen : (ext σ r p q).length = 3 * (4 * q + r) + c - (4 * q + r) := by
    unfold ext
    simp only [List.length_append, List.length_map, hhl, htl, bulk_length, lastJ]
    omega
  exact perfect_extend hpre (by omega) hS hnd hmem hlen

end Orders

