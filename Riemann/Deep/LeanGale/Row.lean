import Riemann.Deep.Gale7.Lean.Crit
import Riemann.Deep.LeanGale.Splx

/-!
# LeanGale: one row of the family table, checked in integers

`Params`: the 8 parameters (r_A, r_B, a, b, c, d, e, f) of family_params.txt as integers over 10¹², and two
bounds `tZ`, `DZ`. `clsZ m₁ m₂ p` = m₁ m₂ 10²⁴ × (Gale7's rational class value `clsQ`) (`clsQ_eq`).
`rowOK m₁ m₂ p` (a Bool, evaluated by the kernel for every row): every class value ≤ `DZ`, every edge class
of H (`adjD`) ≥ `tZ`, every class of two distinct vertices (`neD`) > 0, and the exact σ test
`sigmaOK (m₁ + m₂) tZ DZ` (σₙ² · DZ < tZ). `rowOK_sound`: then the family member Vₙ (`famRow`, any unit
simplices U, W) has no partition into n + 1 parts of diameter ≤ σₙ · diam (Gale7's `GaleSplit`).
-/

namespace LeanGale

open Gale7 Vx

/-- a row of family_params.txt: the parameters × 10¹² and the bounds `tZ ≤` edge classes, classes `≤ DZ`
(both × m₁ m₂ 10²⁴) -/
structure Params where
  rA : ℕ
  rB : ℕ
  a : ℕ
  b : ℕ
  c : ℕ
  d : ℕ
  e : ℕ
  f : ℕ
  tZ : ℕ
  DZ : ℕ
  deriving Inhabited

/-- first-block coefficient of an orbit × 10¹²: (r_A, 0, a, −c, e) -/
def k1Z (p : Params) : Orb → ℤ
  | .P => p.rA
  | .Q => 0
  | .A => p.a
  | .B => -(p.c : ℤ)
  | .C => p.e

/-- second-block coefficient of an orbit × 10¹²: (0, r_B, −b, d, f) -/
def k2Z (p : Params) : Orb → ℤ
  | .P => 0
  | .Q => p.rB
  | .A => -(p.b : ℤ)
  | .B => p.d
  | .C => p.f

/-- m₁ m₂ 10²⁴ × the class value of a descriptor -/
def clsZ (m₁ m₂ : ℕ) (p : Params) (o o' : Orb) (e₁ e₂ : Bool) : ℤ :=
  (m₂ : ℤ) * ((m₁ : ℤ) * (k1Z p o ^ 2 + k1Z p o' ^ 2) - 2 * k1Z p o * k1Z p o' * (if e₁ then (m₁ : ℤ) else -1))
  + (m₁ : ℤ) * ((m₂ : ℤ) * (k2Z p o ^ 2 + k2Z p o' ^ 2) - 2 * k2Z p o * k2Z p o' * (if e₂ then (m₂ : ℤ) else -1))

/-- a parameter as a rational number -/
def prm (x : ℕ) : ℚ := (x : ℚ) / 10 ^ 12

theorem k1Q_prm (p : Params) (o : Orb) :
    k1Q (prm p.rA) (prm p.a) (prm p.c) (prm p.e) o = (k1Z p o : ℚ) / 10 ^ 12 := by
  cases o <;> simp [k1Q, k1Z, prm, neg_div]

theorem k2Q_prm (p : Params) (o : Orb) :
    k2Q (prm p.rB) (prm p.b) (prm p.d) (prm p.f) o = (k2Z p o : ℚ) / 10 ^ 12 := by
  cases o <;> simp [k2Q, k2Z, prm, neg_div]

theorem clsQ_eq {m₁ m₂ : ℕ} (hm₁ : m₁ ≠ 0) (hm₂ : m₂ ≠ 0) (p : Params) (o o' : Orb) (e₁ e₂ : Bool) :
    clsQ m₁ m₂ (prm p.rA) (prm p.rB) (prm p.a) (prm p.b) (prm p.c) (prm p.d) (prm p.e) (prm p.f) o o' e₁ e₂ =
      (clsZ m₁ m₂ p o o' e₁ e₂ : ℚ) / ((m₁ : ℚ) * m₂ * 10 ^ 24) := by
  have h1 : (m₁ : ℚ) ≠ 0 := Nat.cast_ne_zero.mpr hm₁
  have h2 : (m₂ : ℚ) ≠ 0 := Nat.cast_ne_zero.mpr hm₂
  simp only [clsQ, clsZ, k1Q_prm, k2Q_prm]
  cases e₁ <;> cases e₂ <;> simp only [if_true, if_false, Bool.false_eq_true] <;> push_cast <;>
    field_simp <;> ring

/-- the five orbits -/
def orbs : List Orb := [.P, .Q, .A, .B, .C]

/-- all 100 descriptors -/
def descs : List (Orb × Orb × Bool × Bool) :=
  orbs.flatMap fun o => orbs.flatMap fun o' =>
    [true, false].flatMap fun e₁ => [true, false].map fun e₂ => (o, o', e₁, e₂)

theorem mem_descs (o o' : Orb) (e₁ e₂ : Bool) : (o, o', e₁, e₂) ∈ descs := by
  cases o <;> cases o' <;> cases e₁ <;> cases e₂ <;> decide

/-- exact test σₙ² · DZ < tZ (σₙ² = (n+1)/(n+2) for even n, 1/2 + √((n−1)/(n+3))/2 for odd n) -/
def sigmaOK (n tZ DZ : ℕ) : Bool :=
  if n % 2 = 0 then decide ((n + 1) * DZ < (n + 2) * tZ)
  else decide (DZ < 2 * tZ) && decide ((n - 1) * DZ ^ 2 < (n + 3) * (2 * tZ - DZ) ^ 2)

/-- the three class checks for one descriptor -/
def descOK (m₁ m₂ : ℕ) (p : Params) (x : Orb × Orb × Bool × Bool) : Bool :=
  decide (clsZ m₁ m₂ p x.1 x.2.1 x.2.2.1 x.2.2.2 ≤ p.DZ) &&
  (!adjD x.1 x.2.1 x.2.2.1 x.2.2.2 || decide ((p.tZ : ℤ) ≤ clsZ m₁ m₂ p x.1 x.2.1 x.2.2.1 x.2.2.2)) &&
  (!neD x.1 x.2.1 x.2.2.1 x.2.2.2 || decide (0 < clsZ m₁ m₂ p x.1 x.2.1 x.2.2.1 x.2.2.2))

/-- the row check (kernel-evaluated for every row of the table) -/
def rowOK (m₁ m₂ : ℕ) (p : Params) : Bool :=
  decide (1 ≤ m₁) && decide (1 ≤ m₂) && sigmaOK (m₁ + m₂) p.tZ p.DZ && descs.all (descOK m₁ m₂ p)

/-- consecutive rows for n = start, start + 1, … with the split ⌊n/2⌋ + ⌈n/2⌉ -/
def chunkOK : ℕ → List Params → Bool
  | _, [] => true
  | n, p :: ps => rowOK (n / 2) (n - n / 2) p && chunkOK (n + 1) ps

theorem chunkOK_append (n : ℕ) (ps qs : List Params) :
    chunkOK n (ps ++ qs) = (chunkOK n ps && chunkOK (n + ps.length) qs) := by
  induction ps generalizing n with
  | nil => simp [chunkOK]
  | cons p ps ih =>
      simp only [List.cons_append, chunkOK, ih, List.length_cons, Bool.and_assoc]
      rw [show n + 1 + ps.length = n + (ps.length + 1) by omega]

theorem chunkOK_getD {n : ℕ} {ps : List Params} (h : chunkOK n ps = true) {i : ℕ} (hi : i < ps.length) :
    rowOK ((n + i) / 2) ((n + i) - (n + i) / 2) (ps.getD i default) = true := by
  induction ps generalizing n i with
  | nil => exact absurd hi (Nat.not_lt_zero _)
  | cons p ps ih =>
      simp only [chunkOK, Bool.and_eq_true] at h
      cases i with
      | zero => simpa using h.1
      | succ i =>
          have := ih h.2 (i := i) (by simpa using hi)
          rw [show n + 1 + i = n + (i + 1) by omega] at this
          simpa using this

/-- the family member of a row (NOTE §2), parameters `prm` (rational, cast to ℝ) -/
noncomputable def famRow {m₁ m₂ : ℕ} (U : Fin (m₁ + 1) → Fin m₁ → ℝ) (W : Fin (m₂ + 1) → Fin m₂ → ℝ)
    (p : Params) : Vx (m₁ + 1) (m₂ + 1) → Fin (m₁ + m₂) → ℝ :=
  famPt U W (prm p.rA) (prm p.rB) (prm p.a) (prm p.b) (prm p.c) (prm p.d) (prm p.e) (prm p.f)

theorem sqd_famRow {m₁ m₂ : ℕ} (hm₁ : m₁ ≠ 0) (hm₂ : m₂ ≠ 0) {U : Fin (m₁ + 1) → Fin m₁ → ℝ}
    {W : Fin (m₂ + 1) → Fin m₂ → ℝ} (hU : IsUnitSimplex U) (hW : IsUnitSimplex W) (p : Params)
    (v w : Vx (m₁ + 1) (m₂ + 1)) :
    sqd (famRow U W p v) (famRow U W p w) =
      (clsZ m₁ m₂ p (orb v) (orb w) (decide (ix1 v = ix1 w)) (decide (ix2 v = ix2 w)) : ℝ) /
        ((m₁ : ℝ) * m₂ * 10 ^ 24) := by
  rw [famRow, sqd_famPt_cls hU hW, clsQ_eq hm₁ hm₂]
  push_cast
  ring

/-- σₙ² · D < t from the exact integer test -/
theorem sigma_sq_mul_lt {n tZ DZ : ℕ} (hn : 1 ≤ n) (h : sigmaOK n tZ DZ = true) :
    sigma n ^ 2 * (DZ : ℝ) < tZ := by
  unfold sigmaOK at h
  split_ifs at h with hpar
  · rw [decide_eq_true_eq] at h
    rw [sigma_even_sq hpar]
    have h' : ((n : ℝ) + 1) * DZ < ((n : ℝ) + 2) * tZ := by exact_mod_cast h
    rw [div_mul_eq_mul_div, div_lt_iff₀ (by positivity)]
    linarith
  · rw [Bool.and_eq_true, decide_eq_true_eq, decide_eq_true_eq] at h
    obtain ⟨h1, h2⟩ := h
    rw [sigma_odd_sq (by omega) hn]
    have hle : DZ ≤ 2 * tZ := h1.le
    have h2' : ((n : ℝ) - 1) * (DZ : ℝ) ^ 2 < ((n : ℝ) + 3) * (2 * (tZ : ℝ) - DZ) ^ 2 := by
      have := h2
      have e1 : (((n - 1 : ℕ)) : ℝ) = (n : ℝ) - 1 := by rw [Nat.cast_sub hn]; simp
      have e2 : (((2 * tZ - DZ : ℕ)) : ℝ) = 2 * (tZ : ℝ) - DZ := by rw [Nat.cast_sub hle]; push_cast; ring
      have h3 : (((n - 1) * DZ ^ 2 : ℕ) : ℝ) < (((n + 3) * (2 * tZ - DZ) ^ 2 : ℕ) : ℝ) := by exact_mod_cast this
      push_cast [e1, e2] at h3
      rw [← e1, ← e2]
      exact_mod_cast h3
    have hu : 0 < 2 * (tZ : ℝ) - DZ := by
      have : (DZ : ℝ) < 2 * tZ := by exact_mod_cast h1
      linarith
    have hs : Real.sqrt (((n : ℝ) - 1) / ((n : ℝ) + 3)) * DZ < 2 * (tZ : ℝ) - DZ := by
      rcases Nat.eq_zero_or_pos DZ with hD | hD
      · have hD0 : (DZ : ℝ) = 0 := by exact_mod_cast hD
        rw [hD0, mul_zero, sub_zero]
        linarith
      · have hD' : (0 : ℝ) < DZ := by exact_mod_cast hD
        rw [← lt_div_iff₀ hD', Real.sqrt_lt' (div_pos hu hD'), div_pow,
          div_lt_div_iff₀ (by positivity) (by positivity)]
        nlinarith
    nlinarith

/-- **Row soundness**: an accepted row gives a family member with no partition into n + 1 parts of
diameter ≤ σₙ · diam (n = m₁ + m₂), for any unit simplices U, W. -/
theorem rowOK_sound {m₁ m₂ : ℕ} {p : Params} (h : rowOK m₁ m₂ p = true)
    {U : Fin (m₁ + 1) → Fin m₁ → ℝ} {W : Fin (m₂ + 1) → Fin m₂ → ℝ}
    (hU : IsUnitSimplex U) (hW : IsUnitSimplex W) :
    ¬ GaleSplit (m₁ + m₂) (Finset.univ.image (famRow U W p)) := by
  simp only [rowOK, Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true] at h
  obtain ⟨⟨⟨hm₁, hm₂⟩, hσ⟩, hall⟩ := h
  have hM : (0 : ℝ) < (m₁ : ℝ) * m₂ * 10 ^ 24 := by
    have : (0 : ℝ) < m₁ := by exact_mod_cast hm₁
    have : (0 : ℝ) < m₂ := by exact_mod_cast hm₂
    positivity
  have hd : ∀ v w : Vx (m₁ + 1) (m₂ + 1), descOK m₁ m₂ p
      (orb v, orb w, decide (ix1 v = ix1 w), decide (ix2 v = ix2 w)) = true :=
    fun v w => hall _ (mem_descs _ _ _ _)
  refine not_galeSplit_of_H (p := m₁ + 1) (q := m₂ + 1) (by omega) (famRow U W p)
    (t := (p.tZ : ℝ) / ((m₁ : ℝ) * m₂ * 10 ^ 24)) (D2 := (p.DZ : ℝ) / ((m₁ : ℝ) * m₂ * 10 ^ 24))
    (by positivity) (fun v w => ?_) (fun v w hvw => ?_) ?_
  · rw [sqd_famRow (by omega) (by omega) hU hW, div_le_div_iff_of_pos_right hM]
    have := hd v w
    simp only [descOK, Bool.and_eq_true, decide_eq_true_eq] at this
    exact_mod_cast this.1.1
  · rw [sqd_famRow (by omega) (by omega) hU hW, div_le_div_iff_of_pos_right hM]
    have := hd v w
    simp only [descOK, Bool.and_eq_true, decide_eq_true_eq, Bool.or_eq_true, Bool.not_eq_true'] at this
    rcases this.1.2 with h1 | h1
    · rw [adj_adjD hvw] at h1; exact absurd h1 (by decide)
    · exact_mod_cast h1
  · rw [← mul_div_assoc, div_lt_div_iff_of_pos_right hM]
    exact sigma_sq_mul_lt (by omega) hσ

/-- distinct vertices give distinct points -/
theorem famRow_injective {m₁ m₂ : ℕ} {p : Params} (h : rowOK m₁ m₂ p = true)
    {U : Fin (m₁ + 1) → Fin m₁ → ℝ} {W : Fin (m₂ + 1) → Fin m₂ → ℝ}
    (hU : IsUnitSimplex U) (hW : IsUnitSimplex W) : Function.Injective (famRow U W p) := by
  simp only [rowOK, Bool.and_eq_true, decide_eq_true_eq, List.all_eq_true] at h
  obtain ⟨⟨⟨hm₁, hm₂⟩, _⟩, hall⟩ := h
  intro v w hvw
  by_contra hne
  have hd := hall _ (mem_descs (orb v) (orb w) (decide (ix1 v = ix1 w)) (decide (ix2 v = ix2 w)))
  simp only [descOK, Bool.and_eq_true, decide_eq_true_eq, Bool.or_eq_true, Bool.not_eq_true'] at hd
  have hpos : (0 : ℝ) < sqd (famRow U W p v) (famRow U W p w) := by
    rw [sqd_famRow (by omega) (by omega) hU hW]
    rcases hd.2 with h1 | h1
    · rw [ne_neD hne] at h1; exact absurd h1 (by decide)
    · have hM : (0 : ℝ) < (m₁ : ℝ) * m₂ * 10 ^ 24 := by
        have : (0 : ℝ) < m₁ := by exact_mod_cast hm₁
        have : (0 : ℝ) < m₂ := by exact_mod_cast hm₂
        positivity
      exact div_pos (by exact_mod_cast h1) hM
  rw [hvw, sqd_self] at hpos
  exact lt_irrefl _ hpos

end LeanGale
