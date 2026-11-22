import Mathlib.Data.List.GetD
import Mathlib.Data.Nat.Bitwise

/-!
# LeanGale: a kernel-checkable refutation certificate for graph colourings

A state is one natural number: bits `8u … 8u+7` hold the colour domain of vertex `u` (bit `8u + i` set:
colour `i` still possible), so the kernel updates a whole state with a few GMP operations. Adjacency is
given by spread masks: `sp.getD v 0` has bit `8w` set iff `v` and `w` are adjacent (`SpreadOK` checks this
against an edge predicate `E`). A certificate is a preorder token list of a refutation tree, run on a stack
of states (`run`):
* `br v k`: split the top state into `c v = k` (`assign`: `dom v := {k}`, colour `k` removed from the
  neighbours of `v`) and `c v ≠ k` (bit `8v + k` cleared);
* `unit v k`: `dom v = {k}`, so assign `v ↦ k` (no split);
* `cut v a b`: every domain is symmetric in the colours `a ≠ b`, so without loss of generality `c v ≠ b`
  (swapping `a` and `b` in a colouring keeps it proper and inside the domains);
* `wipe u`: `dom u = ∅`; the top state is closed.
`exists_mono`: if `run` accepts from full domains, every colouring with `K ≤ 8` colours has two adjacent
vertices of the same colour.
-/

namespace LeanGale.Cert

/-- clear bit `j` of `s` -/
def clr (s j : ℕ) : ℕ := s ^^^ (s &&& 2 ^ j)

theorem testBit_clr (s j i : ℕ) : (clr s j).testBit i = (s.testBit i && !decide (j = i)) := by
  rw [clr, Nat.testBit_xor, Nat.testBit_and, Nat.testBit_two_pow]
  cases s.testBit i <;> cases decide (j = i) <;> rfl

/-- `v ↦ k`: `dom v := {k}` and colour `k` is removed from the domains of the neighbours of `v`
(`F = 2^(8N) − 1`) -/
def assign (sp : List ℕ) (F v k s : ℕ) : ℕ :=
  ((s &&& (F ^^^ (sp.getD v 0 <<< k))) &&& (F ^^^ (255 <<< (8 * v)))) ||| 2 ^ (8 * v + k)

/-- every domain (of the vertices `0 … N-1`) is symmetric in the colours `a`, `b` -/
def symm (N a b s : ℕ) : Bool := (List.range N).all fun u => s.testBit (8 * u + a) == s.testBit (8 * u + b)

/-- tokens of a refutation certificate -/
inductive Tok where
  | br (v k : ℕ)
  | unit (v k : ℕ)
  | cut (v a b : ℕ)
  | wipe (u : ℕ)

/-- run a certificate on a stack of states; `true` iff every state is closed -/
def run (sp : List ℕ) (N K : ℕ) : List Tok → List ℕ → Bool
  | [], st => st.isEmpty
  | .wipe u :: ts, s :: st => (decide (u < N) && (s >>> (8 * u)) % 256 == 0) && run sp N K ts st
  | .br v k :: ts, s :: st =>
      (decide (v < N) && decide (k < 8)) &&
        run sp N K ts (assign sp (2 ^ (8 * N) - 1) v k s :: clr s (8 * v + k) :: st)
  | .unit v k :: ts, s :: st =>
      (decide (v < N) && decide (k < 8) && (s >>> (8 * v)) % 256 == 2 ^ k) &&
        run sp N K ts (assign sp (2 ^ (8 * N) - 1) v k s :: st)
  | .cut v a b :: ts, s :: st =>
      (decide (v < N) && a != b && decide (a < K) && decide (b < K) && symm N a b s) &&
        run sp N K ts (clr s (8 * v + b) :: st)
  | _ :: _, [] => false

/-- `c` is a proper colouring for the spread masks -/
def Proper (sp : List ℕ) (c : ℕ → ℕ) : Prop := ∀ v w, (sp.getD v 0).testBit (8 * w) = true → c v ≠ c w

/-- every vertex's colour lies in its domain -/
def Inv (N : ℕ) (c : ℕ → ℕ) (s : ℕ) : Prop := ∀ u < N, s.testBit (8 * u + c u) = true

/-- some proper colouring with `K` colours lies inside the domains of `s` -/
def Good (sp : List ℕ) (N K s : ℕ) : Prop := ∃ c : ℕ → ℕ, (∀ u, c u < K) ∧ Proper sp c ∧ Inv N c s

/-- spread masks have bits only at slot starts below `8N` -/
def SpreadBits (sp : List ℕ) (N : ℕ) : Prop := ∀ v j, (sp.getD v 0).testBit j = true → j % 8 = 0 ∧ j < 8 * N

theorem testBit_assign {sp : List ℕ} {N : ℕ} (hsp : SpreadBits sp N) {v k s u i : ℕ} (hu : u < N)
    (hi : i < 8) (hk : k < 8) :
    (assign sp (2 ^ (8 * N) - 1) v k s).testBit (8 * u + i) =
      if u = v then decide (i = k) else (s.testBit (8 * u + i) && !(decide (i = k) && (sp.getD v 0).testBit (8 * u))) := by
  have hF : (2 ^ (8 * N) - 1).testBit (8 * u + i) = true := by
    rw [Nat.testBit_two_pow_sub_one, decide_eq_true_eq]; omega
  have h255 : (255 : ℕ) = 2 ^ 8 - 1 := rfl
  simp only [assign, Nat.testBit_or, Nat.testBit_and, Nat.testBit_xor, hF, Nat.testBit_shiftLeft,
    Nat.testBit_two_pow, h255, Nat.testBit_two_pow_sub_one]
  by_cases huv : u = v
  · subst huv
    have e1 : decide (8 * u + i ≥ 8 * u) = true := by simp
    have e2 : decide (8 * u + i - 8 * u < 8) = true := by simp; omega
    rw [e1, e2]
    by_cases hik : i = k
    · subst hik; simp
    · have : decide (8 * u + k = 8 * u + i) = false := by simp; omega
      rw [this]; simp [hik]
  · rw [if_neg huv]
    have e3 : (decide (8 * u + i ≥ 8 * v) && decide (8 * u + i - 8 * v < 8)) = false := by
      by_cases h : 8 * u + i ≥ 8 * v
      · have : ¬ (8 * u + i - 8 * v < 8) := by omega
        simp [h, this]
      · simp [h]
    have e4 : decide (8 * v + k = 8 * u + i) = false := by simp; omega
    rw [e3, e4]
    by_cases hik : i = k
    · subst hik
      have e5 : decide (8 * u + i ≥ i) = true := by simp
      rw [e5, show 8 * u + i - i = 8 * u by omega]
      simp
    · have e6 : (decide (8 * u + i ≥ k) && (sp.getD v 0).testBit (8 * u + i - k)) = false := by
        by_cases h : 8 * u + i ≥ k
        · have : (sp.getD v 0).testBit (8 * u + i - k) = false := by
            by_contra hb
            rw [Bool.not_eq_false] at hb
            have := (hsp v _ hb).1
            omega
          simp only [this, Bool.and_false]
        · simp [h]
      rw [e6]
      simp [hik]

theorem slot_testBit (s v i : ℕ) (hi : i < 8) : ((s >>> (8 * v)) % 256).testBit i = s.testBit (8 * v + i) := by
  rw [show (256 : ℕ) = 2 ^ 8 from rfl, Nat.testBit_mod_two_pow, Nat.testBit_shiftRight]
  simp [hi]

theorem inv_assign {sp : List ℕ} {N : ℕ} (hsp : SpreadBits sp N) {c : ℕ → ℕ} {s v k : ℕ}
    (hc8 : ∀ u, c u < 8) (hP : Proper sp c) (hI : Inv N c s) (hv : c v = k) :
    Inv N c (assign sp (2 ^ (8 * N) - 1) v k s) := by
  intro u hu
  rw [testBit_assign hsp hu (hc8 u) (hv ▸ hc8 v)]
  by_cases huv : u = v
  · subst huv; simp [hv]
  · rw [if_neg huv, hI u hu, Bool.true_and]
    by_cases hb : (sp.getD v 0).testBit (8 * u) = true
    · have : c u ≠ k := fun h => hP v u hb (hv.trans h.symm)
      simp [this]
    · have hb' : (sp.getD v 0).testBit (8 * u) = false := by simpa using hb
      rw [hb', Bool.and_false, Bool.not_false]

theorem inv_clr {N : ℕ} {c : ℕ → ℕ} {s v b : ℕ} (hc8 : ∀ u, c u < 8) (hI : Inv N c s) (hv : c v ≠ b)
    (hb : b < 8) : Inv N c (clr s (8 * v + b)) := by
  intro u hu
  rw [testBit_clr, hI u hu, Bool.true_and, Bool.not_eq_true', decide_eq_false_iff_not]
  intro h
  have := hc8 u
  have huv : u = v := by omega
  subst huv
  exact hv (by omega)

/-- swap the colours `a` and `b` -/
def swapc (a b x : ℕ) : ℕ := if x = a then b else if x = b then a else x

theorem swapc_inj (a b : ℕ) {x y : ℕ} (h : swapc a b x = swapc a b y) : x = y := by
  unfold swapc at h
  split_ifs at h <;> omega

theorem good_cut {sp : List ℕ} {N K s v a b : ℕ} (hK : K ≤ 8) (hab : a ≠ b) (ha : a < K)
    (hb : b < K) (hs : symm N a b s = true) (hg : Good sp N K s) : Good sp N K (clr s (8 * v + b)) := by
  obtain ⟨c, hcK, hP, hI⟩ := hg
  have hc8 : ∀ u, c u < 8 := fun u => lt_of_lt_of_le (hcK u) hK
  by_cases hv : c v = b
  · have hc'K : ∀ x, swapc a b (c x) < K := fun x => by
      unfold swapc; split_ifs <;> first | exact hb | exact ha | exact hcK x
    refine ⟨fun x => swapc a b (c x), hc'K, fun x y hxy h => hP x y hxy (swapc_inj a b h), ?_⟩
    apply inv_clr (fun u => lt_of_lt_of_le (hc'K u) hK) _ _ (lt_of_lt_of_le hb hK)
    · intro u hu
      have hsym := (List.all_eq_true.mp hs) u (List.mem_range.mpr hu)
      rw [beq_iff_eq] at hsym
      have hu' := hI u hu
      show s.testBit (8 * u + swapc a b (c u)) = true
      unfold swapc
      split_ifs with h1 h2
      · rw [h1] at hu'; rw [← hsym]; exact hu'
      · rw [h2] at hu'; rw [hsym]; exact hu'
      · exact hu'
    · show swapc a b (c v) ≠ b
      unfold swapc
      rw [hv, if_neg (Ne.symm hab), if_pos rfl]
      exact hab
  · exact ⟨c, hcK, hP, inv_clr hc8 hI hv (lt_of_lt_of_le hb hK)⟩

theorem run_sound {sp : List ℕ} {N K : ℕ} (hK : K ≤ 8) (hsp : SpreadBits sp N) :
    ∀ (ts : List Tok) (st : List ℕ), run sp N K ts st = true → ∀ s ∈ st, ¬ Good sp N K s
  | [], st, h => by
      intro s hs
      have : st = [] := by simpa [run] using h
      simp [this] at hs
  | .wipe u :: ts, [], h => by simp [run] at h
  | .br v k :: ts, [], h => by simp [run] at h
  | .unit v k :: ts, [], h => by simp [run] at h
  | .cut v a b :: ts, [], h => by simp [run] at h
  | .wipe u :: ts, s :: st, h => by
      simp only [run, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at h
      obtain ⟨⟨hu, h1⟩, h2⟩ := h
      intro s' hs'
      rcases List.mem_cons.mp hs' with rfl | hs'
      · rintro ⟨c, hcK, _, hI⟩
        have hb := hI u hu
        rw [← slot_testBit _ _ _ (lt_of_lt_of_le (hcK u) hK), h1, Nat.zero_testBit] at hb
        exact Bool.false_ne_true hb
      · exact run_sound hK hsp ts st h2 s' hs'
  | .br v k :: ts, s :: st, h => by
      simp only [run, Bool.and_eq_true, decide_eq_true_eq] at h
      obtain ⟨⟨hv, hk⟩, h2⟩ := h
      have ih := run_sound hK hsp ts _ h2
      intro s' hs'
      rcases List.mem_cons.mp hs' with rfl | hs'
      · rintro ⟨c, hcK, hP, hI⟩
        have hc8 : ∀ u, c u < 8 := fun u => lt_of_lt_of_le (hcK u) hK
        by_cases hc : c v = k
        · exact ih _ List.mem_cons_self ⟨c, hcK, hP, inv_assign hsp hc8 hP hI hc⟩
        · exact ih _ (List.mem_cons_of_mem _ List.mem_cons_self) ⟨c, hcK, hP, inv_clr hc8 hI hc hk⟩
      · exact ih s' (List.mem_cons_of_mem _ (List.mem_cons_of_mem _ hs'))
  | .unit v k :: ts, s :: st, h => by
      simp only [run, Bool.and_eq_true, decide_eq_true_eq, beq_iff_eq] at h
      obtain ⟨⟨⟨hv, hk⟩, h1⟩, h2⟩ := h
      have ih := run_sound hK hsp ts _ h2
      intro s' hs'
      rcases List.mem_cons.mp hs' with rfl | hs'
      · rintro ⟨c, hcK, hP, hI⟩
        have hc8 : ∀ u, c u < 8 := fun u => lt_of_lt_of_le (hcK u) hK
        have hb := hI v hv
        rw [← slot_testBit _ _ _ (hc8 v), h1, Nat.testBit_two_pow, decide_eq_true_eq] at hb
        exact ih _ List.mem_cons_self ⟨c, hcK, hP, inv_assign hsp hc8 hP hI hb.symm⟩
      · exact ih s' (List.mem_cons_of_mem _ hs')
  | .cut v a b :: ts, s :: st, h => by
      simp only [run, Bool.and_eq_true, bne_iff_ne, ne_eq, decide_eq_true_eq] at h
      obtain ⟨⟨⟨⟨⟨_, hab⟩, ha⟩, hb⟩, hs⟩, h2⟩ := h
      have ih := run_sound hK hsp ts _ h2
      intro s' hs'
      rcases List.mem_cons.mp hs' with rfl | hs'
      · exact fun hg => ih _ List.mem_cons_self (good_cut hK hab ha hb hs hg)
      · exact ih s' (List.mem_cons_of_mem _ hs')

/-- full domains `{0, …, K-1}` for the vertices `0 … n-1` -/
def initS (K : ℕ) : ℕ → ℕ
  | 0 => 0
  | n + 1 => initS K n ||| ((2 ^ K - 1) <<< (8 * n))

theorem testBit_initS {K : ℕ} (hK : K ≤ 8) : ∀ (n u i : ℕ), i < 8 →
    (initS K n).testBit (8 * u + i) = (decide (u < n) && decide (i < K))
  | 0, u, i, _ => by simp [initS]
  | n + 1, u, i, hi => by
      rw [initS, Nat.testBit_or, testBit_initS hK n u i hi, Nat.testBit_shiftLeft, Nat.testBit_two_pow_sub_one]
      by_cases hun : u = n
      · subst hun; simp
      · by_cases h : u < n
        · have h1 : u < n + 1 := by omega
          have h2 : ¬ (8 * u + i ≥ 8 * n) := by omega
          simp [h, h1, h2]
        · have h1 : ¬ u < n + 1 := by omega
          by_cases h2 : 8 * u + i ≥ 8 * n
          · have : ¬ (8 * u + i - 8 * n < K) := by omega
            simp [h, h1, h2, this]
          · simp [h, h1, h2]

/-- **Soundness**: an accepted certificate (from full domains) excludes every proper `K`-colouring. -/
theorem no_colouring {sp : List ℕ} {N K : ℕ} {toks : List Tok} (hK : K ≤ 8) (hsp : SpreadBits sp N)
    (h : run sp N K toks [initS K N] = true) : ¬ ∃ c : ℕ → ℕ, (∀ u, c u < K) ∧ Proper sp c := by
  rintro ⟨c, hcK, hP⟩
  refine run_sound hK hsp toks _ h _ List.mem_cons_self ⟨c, hcK, hP, fun u hu => ?_⟩
  rw [testBit_initS hK N u (c u) (lt_of_lt_of_le (hcK u) hK)]
  simp [hu, hcK u]

/-- the spread masks describe a subgraph of `E` on the vertices `0 … N-1` -/
def SpreadOK (sp : List ℕ) (N : ℕ) (E : ℕ → ℕ → Bool) : Bool :=
  decide (sp.length ≤ N) && (List.range sp.length).all fun v =>
    decide (sp.getD v 0 < 2 ^ (8 * N)) && (List.range (8 * N)).all fun j =>
      !(sp.getD v 0).testBit j || (decide (j % 8 = 0) && E v (j / 8))

theorem spreadOK_sound {sp : List ℕ} {N : ℕ} {E : ℕ → ℕ → Bool} (h : SpreadOK sp N E = true)
    {v j : ℕ} (hb : (sp.getD v 0).testBit j = true) : v < N ∧ j % 8 = 0 ∧ j < 8 * N ∧ E v (j / 8) = true := by
  simp only [SpreadOK, Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true, List.mem_range,
    Bool.or_eq_true, Bool.not_eq_true'] at h
  obtain ⟨hlen, h⟩ := h
  have hv : v < sp.length := by
    by_contra hv
    rw [List.getD_eq_default _ _ (Nat.le_of_not_lt hv), Nat.zero_testBit] at hb
    exact Bool.false_ne_true hb
  obtain ⟨hlt, hj⟩ := h v hv
  have hj' : j < 8 * N := by
    by_contra hj'
    rw [Nat.testBit_lt_two_pow (lt_of_lt_of_le hlt (Nat.pow_le_pow_right (by decide : 2 > 0)
      (Nat.le_of_not_lt hj')))] at hb
    exact Bool.false_ne_true hb
  rcases hj j hj' with h1 | ⟨h1, h2⟩
  · rw [hb] at h1; exact absurd h1 (by decide)
  · exact ⟨lt_of_lt_of_le hv hlen, h1, hj', h2⟩

/-- **Colouring theorem**: every colouring of `0 … N-1` with `K ≤ 8` colours has two adjacent (in `E`)
vertices of the same colour, given spread masks inside `E` and an accepted certificate. -/
theorem exists_mono {sp : List ℕ} {K N : ℕ} {toks : List Tok} {E : ℕ → ℕ → Bool} (hK0 : 0 < K) (hK : K ≤ 8)
    (hsp : SpreadOK sp N E = true) (h : run sp N K toks [initS K N] = true)
    (col : Fin N → Fin K) : ∃ v u : Fin N, E v u = true ∧ col v = col u := by
  by_contra hne
  simp only [not_exists, not_and] at hne
  have hbits : SpreadBits sp N := fun v j hb => by
    obtain ⟨_, h1, h2, _⟩ := spreadOK_sound hsp hb
    exact ⟨h1, h2⟩
  apply no_colouring hK hbits h
  refine ⟨fun u => if hu : u < N then (col ⟨u, hu⟩ : ℕ) else 0, fun u => ?_, fun v w hb => ?_⟩
  · dsimp only
    split_ifs
    · exact (col _).isLt
    · exact hK0
  · obtain ⟨hv, _, hj, hE⟩ := spreadOK_sound hsp hb
    have hw : w < N := by omega
    rw [show 8 * w / 8 = w by omega] at hE
    dsimp only
    rw [dif_pos hv, dif_pos hw]
    intro heq
    exact hne ⟨v, hv⟩ ⟨w, hw⟩ hE (Fin.ext heq)

end LeanGale.Cert
