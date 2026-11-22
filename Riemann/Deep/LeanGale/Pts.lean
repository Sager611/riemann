import Riemann.Deep.Gale7.Lean.Crit
import Riemann.Deep.LeanGale.Cert

/-!
# LeanGale: explicit finite point sets with integer coordinates

`PtSet`: `N` points of ℝⁿ given by integer coordinate lists (each coordinate × 10⁸), the squared diameter
`D2` and a threshold `dc2` (both × 10¹⁶), a pair `(iD, jD)` attaining `D2` and a colouring `col`.
`PtSet.ok K` (a Bool, evaluated by the kernel): the lengths are right, distinct points differ, every squared
distance is ≤ `D2` with equality at `(iD, jD)`, and `col` uses `K` colours with every same-colour pair at
squared distance ≤ `dc2`. Then (Gale7's `sqd`, `diam`): `card_X`, `diam_X`, `upper`; the lower bound
`lower_of_cert` comes from a refutation certificate (Cert.lean).
-/

namespace LeanGale

open Gale7 Cert

/-- squared distance of integer coordinate lists -/
def sqZ (a b : List ℤ) : ℤ := (List.zipWith (fun x y => (x - y) ^ 2) a b).sum

/-- a coordinate list as a point of ℝⁿ (coordinates × 10⁻⁸) -/
noncomputable def toPt (n : ℕ) (a : List ℤ) : Fin n → ℝ := fun t => ((a.getD t 0 : ℤ) : ℝ) / 10 ^ 8

theorem sum_sq_getD : ∀ (a b : List ℤ) (n : ℕ), a.length = n → b.length = n →
    ∑ t : Fin n, (((a.getD t 0 : ℤ) : ℝ) - ((b.getD t 0 : ℤ) : ℝ)) ^ 2 = ((sqZ a b : ℤ) : ℝ)
  | [], [], 0, _, _ => by simp [sqZ]
  | x :: xs, y :: ys, n + 1, ha, hb => by
      rw [Fin.sum_univ_succ]
      simp only [Fin.val_zero, List.getD_cons_zero, Fin.val_succ, List.getD_cons_succ]
      rw [sum_sq_getD xs ys n (by simpa using ha) (by simpa using hb)]
      simp only [sqZ, List.zipWith_cons_cons, List.sum_cons]
      push_cast
      ring
  | [], _ :: _, _, ha, hb => by simp at ha hb; omega
  | _ :: _, [], _, ha, hb => by simp at ha hb; omega
  | [], [], _ + 1, ha, _ => by simp at ha
  | _ :: _, _ :: _, 0, ha, _ => by simp at ha

theorem sqd_toPt {n : ℕ} {a b : List ℤ} (ha : a.length = n) (hb : b.length = n) :
    sqd (toPt n a) (toPt n b) = (sqZ a b : ℝ) / 10 ^ 16 := by
  rw [sqd, ← sum_sq_getD a b n ha hb, div_eq_mul_inv, Finset.sum_mul]
  refine Finset.sum_congr rfl fun t _ => ?_
  simp only [toPt]
  ring

/-- `P i j` for all `i, j < N` -/
def allPairs (N : ℕ) (P : ℕ → ℕ → Bool) : Bool :=
  (List.range N).all fun i => (List.range N).all fun j => P i j

theorem allPairs_sound {N : ℕ} {P : ℕ → ℕ → Bool} (h : allPairs N P = true) {i j : ℕ} (hi : i < N)
    (hj : j < N) : P i j = true := by
  simp only [allPairs, List.all_eq_true, List.mem_range] at h
  exact h i hi j hj

/-- an explicit finite point set (see the module docstring) -/
structure PtSet where
  N : ℕ
  n : ℕ
  pts : List (List ℤ)
  D2 : ℤ
  dc2 : ℤ
  iD : ℕ
  jD : ℕ
  col : List ℕ

namespace PtSet

variable (S : PtSet)

/-- squared distance (× 10¹⁶) of the points `i`, `j` -/
def sq (i j : ℕ) : ℤ := sqZ (S.pts.getD i []) (S.pts.getD j [])

/-- edges of G(≥ dc²) -/
def E (i j : ℕ) : Bool := decide (S.dc2 ≤ S.sq i j)

/-- the kernel-checked facts -/
def ok (K : ℕ) : Bool :=
  (S.pts.length == S.N) && S.pts.all (fun p => p.length == S.n) && decide (0 < S.dc2) &&
  decide (S.iD < S.N) && decide (S.jD < S.N) && (S.sq S.iD S.jD == S.D2) &&
  allPairs S.N (fun i j => decide (S.sq i j ≤ S.D2) && (i == j || decide (0 < S.sq i j))) &&
  (S.col.length == S.N) && S.col.all (fun c => decide (c < K)) &&
  allPairs S.N (fun i j => S.col.getD i 0 != S.col.getD j 0 || decide (S.sq i j ≤ S.dc2))

/-- the points -/
noncomputable def pt (i : Fin S.N) : Fin S.n → ℝ := toPt S.n (S.pts.getD i [])

/-- the point set -/
noncomputable def X : Finset (Fin S.n → ℝ) := Finset.univ.image S.pt

variable {S} {K : ℕ}

theorem ok_parts (h : S.ok K = true) :
    S.pts.length = S.N ∧ (∀ p ∈ S.pts, p.length = S.n) ∧ 0 < S.dc2 ∧ S.iD < S.N ∧ S.jD < S.N ∧
      S.sq S.iD S.jD = S.D2 ∧
      allPairs S.N (fun i j => decide (S.sq i j ≤ S.D2) && (i == j || decide (0 < S.sq i j))) = true ∧
      S.col.length = S.N ∧ (∀ c ∈ S.col, c < K) ∧
      allPairs S.N (fun i j => S.col.getD i 0 != S.col.getD j 0 || decide (S.sq i j ≤ S.dc2)) = true := by
  simp only [ok, Bool.and_eq_true, beq_iff_eq, List.all_eq_true, decide_eq_true_eq] at h
  obtain ⟨⟨⟨⟨⟨⟨⟨⟨⟨h1, h2⟩, h3⟩, h4⟩, h5⟩, h6⟩, h7⟩, h8⟩, h9⟩, h10⟩ := h
  exact ⟨h1, h2, h3, h4, h5, h6, h7, h8, h9, h10⟩

theorem len_getD (h : S.ok K = true) (i : Fin S.N) : (S.pts.getD i []).length = S.n := by
  obtain ⟨h1, h2, _⟩ := ok_parts h
  apply h2
  rw [List.getD_eq_getElem _ _ (by rw [h1]; exact i.isLt)]
  exact List.getElem_mem _

theorem sqd_pt (h : S.ok K = true) (i j : Fin S.N) :
    sqd (S.pt i) (S.pt j) = (S.sq i j : ℝ) / 10 ^ 16 :=
  sqd_toPt (len_getD h i) (len_getD h j)

theorem pt_injective (h : S.ok K = true) : Function.Injective S.pt := by
  obtain ⟨_, _, _, _, _, _, h7, _⟩ := ok_parts h
  intro i j hij
  by_contra hne
  have := allPairs_sound h7 i.isLt j.isLt
  simp only [Bool.and_eq_true, decide_eq_true_eq, Bool.or_eq_true, beq_iff_eq] at this
  rcases this.2 with h1 | h1
  · exact hne (Fin.ext h1)
  · have h2 := sqd_pt h i j
    rw [hij, sqd_self] at h2
    have : (0 : ℝ) < (S.sq i j : ℝ) / 10 ^ 16 := div_pos (by exact_mod_cast h1) (by positivity)
    linarith

/-- the set has `N` points -/
theorem card_X (h : S.ok K = true) : S.X.card = S.N := by
  rw [X, Finset.card_image_of_injective _ (pt_injective h), Finset.card_univ, Fintype.card_fin]

theorem mem_X {x : Fin S.n → ℝ} : x ∈ S.X ↔ ∃ i, S.pt i = x := by
  simp [X]

/-- every squared distance is ≤ D² -/
theorem sqd_le (h : S.ok K = true) {x y : Fin S.n → ℝ} (hx : x ∈ S.X) (hy : y ∈ S.X) :
    sqd x y ≤ (S.D2 : ℝ) / 10 ^ 16 := by
  obtain ⟨i, rfl⟩ := mem_X.mp hx
  obtain ⟨j, rfl⟩ := mem_X.mp hy
  obtain ⟨_, _, _, _, _, _, h7, _⟩ := ok_parts h
  have := allPairs_sound h7 i.isLt j.isLt
  simp only [Bool.and_eq_true, decide_eq_true_eq] at this
  rw [sqd_pt h, div_le_div_iff_of_pos_right (by positivity)]
  exact_mod_cast this.1

/-- the diameter: √(D² / 10¹⁶) -/
theorem diam_X (h : S.ok K = true) : diam S.X = Real.sqrt ((S.D2 : ℝ) / 10 ^ 16) := by
  obtain ⟨_, _, _, h4, h5, h6, _⟩ := ok_parts h
  apply le_antisymm
  · exact diam_le (Real.sqrt_nonneg _) fun x hx y hy => Real.sqrt_le_sqrt (sqd_le h hx hy)
  · have hx : S.pt ⟨S.iD, h4⟩ ∈ S.X := mem_X.mpr ⟨_, rfl⟩
    have hy : S.pt ⟨S.jD, h5⟩ ∈ S.X := mem_X.mpr ⟨_, rfl⟩
    have := le_diam hx hy
    rw [eucl, sqd_pt h] at this
    simp only at this
    rw [h6] at this
    exact this

/-- the colouring `col` as a colouring of ℝⁿ (0 off the set) -/
noncomputable def colF (K : ℕ) (hK : 0 < K) (x : Fin S.n → ℝ) : Fin K := by
  classical
  exact if hx : ∃ i, S.pt i = x then ⟨S.col.getD (Classical.choose hx).val 0 % K, Nat.mod_lt _ hK⟩ else ⟨0, hK⟩

/-- **Upper bound**: a colouring of the set with `K` colours whose same-colour pairs are at squared
distance ≤ dc² -/
theorem upper (h : S.ok K = true) (hK : 0 < K) :
    ∃ c : (Fin S.n → ℝ) → Fin K, ∀ x ∈ S.X, ∀ y ∈ S.X, c x = c y → sqd x y ≤ (S.dc2 : ℝ) / 10 ^ 16 := by
  classical
  obtain ⟨h1, _, _, _, _, _, _, h8, h9, h10⟩ := ok_parts h
  refine ⟨S.colF K hK, fun x hx y hy hxy => ?_⟩
  obtain ⟨i, rfl⟩ := mem_X.mp hx
  obtain ⟨j, rfl⟩ := mem_X.mp hy
  have key : ∀ k : Fin S.N, S.colF K hK (S.pt k) = ⟨S.col.getD k 0 % K, Nat.mod_lt _ hK⟩ := by
    intro k
    have hk : ∃ i, S.pt i = S.pt k := ⟨k, rfl⟩
    unfold colF
    rw [dif_pos hk]
    have e := pt_injective h (Classical.choose_spec hk)
    simp only [e]
  rw [key, key, Fin.mk.injEq] at hxy
  have hc : ∀ k : Fin S.N, S.col.getD k 0 < K := fun k => by
    apply h9
    rw [List.getD_eq_getElem _ _ (by rw [h8]; exact k.isLt)]
    exact List.getElem_mem _
  rw [Nat.mod_eq_of_lt (hc i), Nat.mod_eq_of_lt (hc j)] at hxy
  have := allPairs_sound h10 i.isLt j.isLt
  simp only [Bool.or_eq_true, bne_iff_ne, ne_eq, decide_eq_true_eq] at this
  rcases this with h11 | h11
  · exact absurd hxy h11
  · rw [sqd_pt h, div_le_div_iff_of_pos_right (by positivity)]
    exact_mod_cast h11

/-- **Lower bound from a certificate**: every colouring of ℝⁿ with `K ≤ 8` colours has two points of the
set of the same colour at squared distance ≥ dc² -/
theorem lower_of_cert (h : S.ok K = true) (hK0 : 0 < K) (hK : K ≤ 8) {sp : List ℕ} {toks : List Tok}
    (hsp : SpreadOK sp S.N S.E = true) (hrun : run sp S.N K toks [initS K S.N] = true)
    (c : (Fin S.n → ℝ) → Fin K) :
    ∃ x ∈ S.X, ∃ y ∈ S.X, c x = c y ∧ (S.dc2 : ℝ) / 10 ^ 16 ≤ sqd x y := by
  obtain ⟨v, u, hE, hc⟩ := exists_mono hK0 hK hsp hrun (fun i => c (S.pt i))
  refine ⟨S.pt v, mem_X.mpr ⟨v, rfl⟩, S.pt u, mem_X.mpr ⟨u, rfl⟩, hc, ?_⟩
  rw [sqd_pt h, div_le_div_iff_of_pos_right (by positivity)]
  simp only [E, decide_eq_true_eq] at hE
  exact_mod_cast hE

end PtSet

/-- the vertices of H(p, q) -/
def vxList (p q : ℕ) : List (Vx p q) :=
  (List.finRange p).map Vx.P ++ (List.finRange q).map Vx.Q ++
    ((List.finRange p).flatMap fun i => (List.finRange q).map fun j => Vx.A i j) ++
    ((List.finRange p).flatMap fun i => (List.finRange q).map fun j => Vx.B i j) ++
    ((List.finRange p).flatMap fun i => (List.finRange q).map fun j => Vx.C i j)

theorem mem_vxList {p q : ℕ} (v : Vx p q) : v ∈ vxList p q := by
  cases v <;> simp [vxList]

namespace PtSet

variable {S : PtSet} {K : ℕ}

/-- **Lower bound from labels**: if the labels §idx§ send every edge of H(m₁+1, m₂+1) to a pair at
squared distance ≥ dc² and §K < m₁ + m₂ + 2§, every colouring of ℝⁿ with §K§ colours has two points of the
set of the same colour at squared distance ≥ dc² (the colouring lemma, Gale7's §Vx.not_colourable§). -/
theorem lower_of_labels (h : S.ok K = true) {m₁ m₂ : ℕ} (hK : K < m₁ + 1 + (m₂ + 1))
    (idx : Vx (m₁ + 1) (m₂ + 1) → ℕ)
    (hidx : (vxList (m₁ + 1) (m₂ + 1)).all (fun v => decide (idx v < S.N)) = true)
    (hH : (vxList (m₁ + 1) (m₂ + 1)).all (fun v => (vxList (m₁ + 1) (m₂ + 1)).all fun w =>
      !adjD (orb v) (orb w) (decide (ix1 v = ix1 w)) (decide (ix2 v = ix2 w)) ||
        decide (S.dc2 ≤ S.sq (idx v) (idx w))) = true)
    (c : (Fin S.n → ℝ) → Fin K) :
    ∃ x ∈ S.X, ∃ y ∈ S.X, c x = c y ∧ (S.dc2 : ℝ) / 10 ^ 16 ≤ sqd x y := by
  have hi : ∀ v, idx v < S.N := fun v => by
    have := (List.all_eq_true.mp hidx) v (mem_vxList v)
    simpa using this
  let f : Vx (m₁ + 1) (m₂ + 1) → Fin S.N := fun v => ⟨idx v, hi v⟩
  obtain ⟨v, w, hvw, hc⟩ := Vx.not_colourable (K := Fin K) (by rw [Fintype.card_fin]; exact hK)
    (fun v => c (S.pt (f v)))
  refine ⟨S.pt (f v), mem_X.mpr ⟨_, rfl⟩, S.pt (f w), mem_X.mpr ⟨_, rfl⟩, hc, ?_⟩
  have h1 := (List.all_eq_true.mp ((List.all_eq_true.mp hH) v (mem_vxList v))) w (mem_vxList w)
  rw [adj_adjD hvw] at h1
  simp only [Bool.not_true, Bool.false_or, decide_eq_true_eq] at h1
  rw [sqd_pt h, div_le_div_iff_of_pos_right (by positivity)]
  exact_mod_cast h1

end PtSet

end LeanGale
