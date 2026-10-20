/-
Riemann/Deep/Yoshida/AReflect.lean — rational checkers for affine-cell certificates (`Affine.lean`).

* `chkA`: walks consecutive cells; `2πμ` and `J` are bounded once per grid point (`muLowQ`, `jLowQ`)
  and each cell checks `s − (A + B r) ≤ W` at its two ends with the second-order loss
  `δ = (mH − mL)(hi − lo)/(4(hi + lo)) + E_j (hi − lo)²`;
* moments of affine cells on a dyadic grid as integer sums;
* `weilPos_of_acert`: the certificate gives Weil positivity on every window `L < log 2`, `L ≤ 2a`.
-/
import Riemann.Deep.Yoshida.Affine
import Riemann.Deep.Yoshida.Reflect

noncomputable section

open Finset

namespace Fayed.Yoshida.AReflect

open Fayed.Yoshida Fayed.Yoshida.Reflect

/-- `log π < 1.14474` (`exp` partial sum against `π < 3.1416`). -/
theorem log_pi_lt' : Real.log Real.pi < 1.14474 := by
  rw [Real.log_lt_iff_lt_exp Real.pi_pos]
  have h := Real.sum_le_exp_of_nonneg (x := 1.14474) (by norm_num) 13
  have hs : (3.1416 : ℝ) < ∑ i ∈ range 13, (1.14474 : ℝ) ^ i / (i.factorial : ℝ) := by
    norm_num [Finset.sum_range_succ, Nat.factorial]
  linarith [Real.pi_lt_d4]

/-! ### Cells -/

def muQ (P : Params) (x : ℚ) : ℚ := muLowQ P.g P.lp P.M P.Mp x

/-- One affine cell, given the end bounds. -/
def cellA (P : Params) (Ej lo hi A B mL jL mH jH : ℚ) : Bool :=
  decide (0 ≤ lo) && decide (lo < hi) && decide (mL ≤ mH) && decide (0 ≤ A + B * lo) &&
    decide (0 ≤ A + B * hi) &&
    decide (P.s - (A + B * lo) ≤ mL + jL - ((mH - mL) * ((hi - lo) / (4 * (hi + lo))) + Ej * (hi - lo) ^ 2)) &&
    decide (P.s - (A + B * hi) ≤ mH + jH - ((mH - mL) * ((hi - lo) / (4 * (hi + lo))) + Ej * (hi - lo) ^ 2))

/-- Walk the cells from `x`; `mL, jL` bound `2πμ`, `J` at `x`. -/
def chkA (P : Params) (Ej : ℚ) : ℚ → ℚ → ℚ → List (ℚ × ℚ × ℚ × ℚ) → ℚ → Bool
  | x, _, _, [], R => decide (x = R)
  | x, mL, jL, c :: cs, R =>
    decide (c.1 = x) && cellA P Ej c.1 c.2.1 c.2.2.1 c.2.2.2 mL jL (muQ P c.2.1) (jLowQ P c.2.1)
      && chkA P Ej c.2.1 (muQ P c.2.1) (jLowQ P c.2.1) cs R

def toR4 (cs : List (ℚ × ℚ × ℚ × ℚ)) : List Cell4 :=
  cs.map fun c => ((c.1 : ℝ), (c.2.1 : ℝ), (c.2.2.1 : ℝ), (c.2.2.2 : ℝ))

theorem cellA_sound {P : Params} (hP : Valid P) {Ej : ℚ}
    (hE : 4 / 3 * (P.a : ℝ) ^ 3 * Real.cosh P.a ≤ Ej) {lo hi A B mL jL mH jH : ℚ}
    (hmL : (mL : ℝ) ≤ 2 * Real.pi * Zeta23.mu lo) (hjL : (jL : ℝ) ≤ J P.a lo)
    (hmH : (mH : ℝ) ≤ 2 * Real.pi * Zeta23.mu hi) (hjH : (jH : ℝ) ≤ J P.a hi)
    (h : cellA P Ej lo hi A B mL jL mH jH = true) :
    AffOK (P.a : ℝ) (P.s : ℝ) ((lo : ℝ), (hi : ℝ), (A : ℝ), (B : ℝ)) := by
  unfold cellA at h
  simp only [Bool.and_eq_true, decide_eq_true_eq] at h
  obtain ⟨⟨⟨⟨⟨⟨h0, h1⟩, h2⟩, h3⟩, h4⟩, h5⟩, h6⟩ := h
  have c (x y : ℚ) (hxy : x ≤ y) : (x : ℝ) ≤ (y : ℝ) := by exact_mod_cast hxy
  have h5' := c _ _ h5
  have h6' := c _ _ h6
  push_cast at h5' h6'
  exact affOK_of hP.ha.le (by exact_mod_cast h0) (by exact_mod_cast h1) hmL hmH hjL hjH hE
    (by exact_mod_cast h2) (by have := c _ _ h3; push_cast at this; exact this)
    (by have := c _ _ h4; push_cast at this; exact this) h5' h6'

theorem chkA_sound {P : Params} (hP : Valid P) {Ej : ℚ}
    (hE : 4 / 3 * (P.a : ℝ) ^ 3 * Real.cosh P.a ≤ Ej) :
    ∀ (cs : List (ℚ × ℚ × ℚ × ℚ)) (x mL jL R : ℚ), 0 ≤ x →
    (mL : ℝ) ≤ 2 * Real.pi * Zeta23.mu x → (jL : ℝ) ≤ J P.a x → chkA P Ej x mL jL cs R = true →
    Chain4 (x : ℝ) (toR4 cs) (R : ℝ) ∧ ∀ c ∈ toR4 cs, AffOK (P.a : ℝ) (P.s : ℝ) c
  | [], x, mL, jL, R, _, _, _, h => by
    simp only [chkA, decide_eq_true_eq] at h
    simp [toR4, Chain4, h]
  | c :: cs, x, mL, jL, R, hx, hmL, hjL, h => by
    simp only [chkA, Bool.and_eq_true, decide_eq_true_eq] at h
    obtain ⟨⟨hcx, hcell⟩, hrest⟩ := h
    have hcell' := hcell
    unfold cellA at hcell'
    simp only [Bool.and_eq_true, decide_eq_true_eq] at hcell'
    have hhi : 0 ≤ c.2.1 := le_trans hcell'.1.1.1.1.1.1 hcell'.1.1.1.1.1.2.le
    have hmH := muLowQ_le hP.M hP.g hP.lp c.2.1
    have hjH := jLowQ_le hP hhi
    rw [← hcx] at hmL hjL
    have ih := chkA_sound hP hE cs c.2.1 _ _ R hhi hmH hjH hrest
    refine ⟨?_, ?_⟩
    · simp only [toR4, List.map_cons, Chain4]
      exact ⟨by rw [hcx], ih.1⟩
    · intro d hd
      simp only [toR4, List.map_cons, List.mem_cons] at hd
      rcases hd with rfl | hd
      · exact cellA_sound hP hE hmL hjL hmH hjH hcell
      · exact ih.2 d hd

theorem chkA_append {P : Params} {Ej : ℚ} :
    ∀ (cs ds : List (ℚ × ℚ × ℚ × ℚ)) (x mL jL y R : ℚ), cs ≠ [] →
    chkA P Ej x mL jL cs y = true → chkA P Ej y (muQ P y) (jLowQ P y) ds R = true →
    chkA P Ej x mL jL (cs ++ ds) R = true
  | [], _, _, _, _, _, _, h, _, _ => absurd rfl h
  | [c], ds, x, mL, jL, y, R, _, h1, h2 => by
    simp only [chkA, Bool.and_eq_true, decide_eq_true_eq] at h1
    obtain ⟨⟨hx, hc⟩, hy⟩ := h1
    simp only [List.cons_append, List.nil_append, chkA, Bool.and_eq_true, decide_eq_true_eq]
    rw [← hy] at h2
    exact ⟨⟨hx, hc⟩, h2⟩
  | c :: c' :: cs, ds, x, mL, jL, y, R, _, h1, h2 => by
    simp only [chkA, Bool.and_eq_true, decide_eq_true_eq] at h1
    obtain ⟨⟨hx, hc⟩, hrest⟩ := h1
    have := chkA_append (c' :: cs) ds c.2.1 (muQ P c.2.1) (jLowQ P c.2.1) y R
      (List.cons_ne_nil _ _) (by simpa only [chkA, Bool.and_eq_true, decide_eq_true_eq] using hrest) h2
    simp only [List.cons_append, chkA, Bool.and_eq_true, decide_eq_true_eq] at this ⊢
    exact ⟨⟨hx, hc⟩, this⟩

def kQ4 (a : ℚ) (K : ℕ) (cs : List (ℚ × ℚ × ℚ × ℚ)) : Bool :=
  cs.all fun c => decide (c.2.2.1 = 0 ∧ c.2.2.2 = 0) || decide (c.2.1 * a ≤ ((K : ℚ) + 1) / 2)

theorem kQ4_sound {a : ℚ} {K : ℕ} {cs : List (ℚ × ℚ × ℚ × ℚ)} (h : kQ4 a K cs = true) :
    ∀ c ∈ toR4 cs, (c.2.2.1 = 0 ∧ c.2.2.2 = 0) ∨ c.2.1 * (a : ℝ) ≤ (K + 1) / 2 := by
  intro c hc
  unfold toR4 at hc
  obtain ⟨c', hc', rfl⟩ := List.mem_map.1 hc
  have := List.all_eq_true.1 h c' hc'
  simp only [Bool.or_eq_true, decide_eq_true_eq] at this
  rcases this with h1 | h1
  · left; simp [h1.1, h1.2]
  · right
    have : ((c'.2.1 * a : ℚ) : ℝ) ≤ ((((K : ℚ) + 1) / 2 : ℚ) : ℝ) := by exact_mod_cast h1
    push_cast at this
    exact this

/-! ### Moments -/

def gmQ (A B lo hi : ℚ) (p : ℕ) : ℚ :=
  A * ((hi ^ (p + 1) - lo ^ (p + 1)) / ((p : ℚ) + 1)) + B * ((hi ^ (p + 2) - lo ^ (p + 2)) / ((p : ℚ) + 2))

def affMomQ (c : ℚ × ℚ × ℚ × ℚ) (p : ℕ) : ℚ :=
  gmQ c.2.2.1 c.2.2.2 c.1 c.2.1 p + gmQ c.2.2.1 (-c.2.2.2) (-c.2.1) (-c.1) p

def momA : List (ℚ × ℚ × ℚ × ℚ) → ℕ → ℚ
  | [], _ => 0
  | c :: cs, p => affMomQ c p + momA cs p

theorem cast_momA : ∀ (cs : List (ℚ × ℚ × ℚ × ℚ)) (p : ℕ), ((momA cs p : ℚ) : ℝ) = affSum (toR4 cs) p
  | [], p => by simp [momA, affSum, toR4]
  | c :: cs, p => by
    rw [momA, Rat.cast_add, cast_momA cs p]
    simp only [toR4, List.map_cons, affSum]
    congr 1
    unfold affMomQ affMom gmQ gm
    push_cast
    ring

/-- Cells on a dyadic grid: `(Lo, Hi, TA, TB)` ↦ `(Lo/2^e, Hi/2^e, TA/2^f, TB/2^f)`. -/
def toQ4 (e f : ℕ) (c : ℕ × ℕ × ℤ × ℤ) : ℚ × ℚ × ℚ × ℚ :=
  ((c.1 : ℚ) / 2 ^ e, (c.2.1 : ℚ) / 2 ^ e, (c.2.2.1 : ℚ) / 2 ^ f, (c.2.2.2 : ℚ) / 2 ^ f)

def momAN : List (ℕ × ℕ × ℤ × ℤ) → ℕ → ℤ
  | [], _ => 0
  | c :: cs, p => c.2.2.1 * ((c.2.1 : ℤ) ^ (p + 1) - (c.1 : ℤ) ^ (p + 1)) + momAN cs p

def momBN : List (ℕ × ℕ × ℤ × ℤ) → ℕ → ℤ
  | [], _ => 0
  | c :: cs, p => c.2.2.2 * ((c.2.1 : ℤ) ^ (p + 2) - (c.1 : ℤ) ^ (p + 2)) + momBN cs p

def momAF (cs : List (ℕ × ℕ × ℤ × ℤ)) (e f p : ℕ) : ℚ :=
  if p % 2 = 0 then
    2 * ((momAN cs p : ℚ) / (((p : ℚ) + 1) * 2 ^ (e * (p + 1)) * 2 ^ f)
      + (momBN cs p : ℚ) / (((p : ℚ) + 2) * 2 ^ (e * (p + 2)) * 2 ^ f))
  else 0

lemma affMomQ_grid (e f : ℕ) (c : ℕ × ℕ × ℤ × ℤ) (p : ℕ) :
    affMomQ (toQ4 e f c) p = if p % 2 = 0 then
      2 * (((c.2.2.1 * ((c.2.1 : ℤ) ^ (p + 1) - (c.1 : ℤ) ^ (p + 1)) : ℤ) : ℚ)
          / (((p : ℚ) + 1) * 2 ^ (e * (p + 1)) * 2 ^ f)
        + ((c.2.2.2 * ((c.2.1 : ℤ) ^ (p + 2) - (c.1 : ℤ) ^ (p + 2)) : ℤ) : ℚ)
          / (((p : ℚ) + 2) * 2 ^ (e * (p + 2)) * 2 ^ f))
    else 0 := by
  unfold affMomQ gmQ toQ4
  simp only
  push_cast
  split_ifs with hpar
  · have hodd : Odd (p + 1) := by
      rw [Nat.odd_add_one, Nat.not_odd_iff_even]
      exact Nat.even_iff.2 hpar
    have hev : Even (p + 2) := by
      rw [show p + 2 = p + 1 + 1 by ring, Nat.even_add_one, Nat.not_even_iff_odd]
      exact hodd
    rw [hodd.neg_pow, hodd.neg_pow, hev.neg_pow, hev.neg_pow, div_pow, div_pow, div_pow, div_pow,
      pow_mul, pow_mul]
    field_simp
    ring
  · have hev : Even (p + 1) := by
      rw [Nat.even_add_one, Nat.not_even_iff_odd]
      exact Nat.odd_iff.2 (by omega)
    have hodd : Odd (p + 2) := by
      rw [show p + 2 = p + 1 + 1 by ring, Nat.odd_add_one, Nat.not_odd_iff_even]
      exact hev
    rw [hev.neg_pow, hev.neg_pow, hodd.neg_pow, hodd.neg_pow]
    ring

theorem momA_grid (e f p : ℕ) : ∀ cs : List (ℕ × ℕ × ℤ × ℤ),
    momA (cs.map (toQ4 e f)) p = momAF cs e f p
  | [] => by simp [momA, momAF, momAN, momBN]
  | c :: cs => by
    simp only [List.map_cons, momA]
    rw [momA_grid e f p cs, affMomQ_grid]
    unfold momAF
    simp only [momAN, momBN]
    split_ifs with hp
    · push_cast
      ring
    · simp

def momCheckA (cs : List (ℕ × ℕ × ℤ × ℤ)) (e f K : ℕ) (momD : List ℚ) : Bool :=
  (List.range (2 * K + 1)).all fun p => decide (momD.getD p 0 = momAF cs e f p)

/-! ### The rational certificate theorem -/

/-- **Rational affine certificate ⇒ Weil positivity on every window `L < log 2`, `L ≤ 2a`.** -/
theorem weilPos_of_acert (P : Params) (hP : Valid P) (Ej : ℚ)
    (cN : List (ℕ × ℕ × ℤ × ℤ)) (e f : ℕ) (R : ℚ) (K : ℕ) (momD : List ℚ)
    (Lm Gm LGm Mm U MU Dm Wm : List (List ℚ)) (α ε β : ℚ)
    (hcells : chkA P Ej 0 (muQ P 0) (jLowQ P 0) (cN.map (toQ4 e f)) R = true)
    (htail : tailQ P R = true)
    (hK : kQ4 P.a K (cN.map (toQ4 e f)) = true) (hmom : momCheckA cN e f K momD = true)
    (hGm : glQ K P.a Lm Gm = true) (hLGm : tmulQ K Lm Gm LGm = true)
    (hM : mQ K α Lm LGm Mm momD = true) (hMU : mulQ K Mm U MU = true)
    (hD : tmulQ K U MU Dm = true) (hUW : invQ K U Wm = true) (hdom : domQ K Dm = true)
    (hnum : decide (0 < ε ∧ 0 ≤ α ∧ 0 ≤ P.s ∧ 4 / 3 * P.a ^ 3 * P.chhi ≤ Ej
      ∧ 8 * P.a ^ (2 * K + 1) * momD.getD (2 * K) 0 / ((K.factorial : ℚ)) ^ 2 ≤ β
      ∧ (1 + ε) * α + (1 + 1 / ε) * β ≤ 2 * piLo * P.s) = true)
    {L : ℝ} (hL : L < Real.log 2) (hL1 : L ≤ 2 * (P.a : ℝ)) :
    Riemann.Deep.Weil.WeilPos L := by
  simp only [decide_eq_true_eq] at hnum
  obtain ⟨hε, hα, hs, hEj, hβ, hn⟩ := hnum
  set cs := cN.map (toQ4 e f) with hcs
  have hE : 4 / 3 * (P.a : ℝ) ^ 3 * Real.cosh P.a ≤ Ej := by
    have h1 : ((4 / 3 * P.a ^ 3 * P.chhi : ℚ) : ℝ) ≤ (Ej : ℝ) := by exact_mod_cast hEj
    push_cast at h1
    have ha := hP.ha
    have h2 : (P.a : ℝ) ^ 3 * Real.cosh P.a ≤ (P.a : ℝ) ^ 3 * P.chhi :=
      mul_le_mul_of_nonneg_left hP.chhi (by positivity)
    nlinarith
  have hmom' : ∀ p ≤ 2 * K, ((momD.getD p 0 : ℚ) : ℝ) = affSum (toR4 cs) p := by
    intro p hp
    have := List.all_eq_true.1 hmom p (List.mem_range.2 (by omega))
    simp only [decide_eq_true_eq] at this
    rw [this, ← momA_grid e f p cN, ← hcs, cast_momA]
  have hGm' := all_range₂ hGm
  have hM' := all_range₂ hM
  have hpsd := hpsd_of K (α : ℝ) (entR Lm) (entR Gm) (entR LGm) (entR Mm) (Ghat (P.a : ℝ))
    (fun p => affSum (toR4 cs) p)
    (fun k hk j hj => by
      unfold entR
      rw [hGm' k hk j hj]
      push_cast
      simp_rw [cast_ghatQ])
    (tmulQ_sound hLGm)
    (fun i hi j hj => by
      unfold entR
      rw [hM' i hi j hj]
      push_cast
      rw [hmom' (i + j) (by omega)])
    (qf_nonneg_of_cert K (entR Mm) (entR U) (entR Wm) (entR MU) (entR Dm) (mulQ_sound hMU)
      (tmulQ_sound hD) (invQ_sound hUW) (domQ_sound hdom))
  have hβ' : 8 * (P.a : ℝ) ^ (2 * K + 1) * affSum (toR4 cs) (2 * K) / (K.factorial : ℝ) ^ 2 ≤ β := by
    rw [← hmom' (2 * K) le_rfl]
    have : ((8 * P.a ^ (2 * K + 1) * momD.getD (2 * K) 0 / ((K.factorial : ℚ)) ^ 2 : ℚ) : ℝ) ≤ (β : ℝ) := by
      exact_mod_cast hβ
    push_cast at this
    exact this
  have hn' : (1 + (ε : ℝ)) * α + (1 + 1 / (ε : ℝ)) * β ≤ 2 * Real.pi * P.s := by
    have h1 : (((1 + ε) * α + (1 + 1 / ε) * β : ℚ) : ℝ) ≤ ((2 * piLo * P.s : ℚ) : ℝ) := by
      exact_mod_cast hn
    push_cast at h1
    have h2 := piLo_lt
    have h3 : (0 : ℝ) ≤ P.s := by exact_mod_cast hs
    nlinarith
  have h0 := chkA_sound hP hE cs 0 _ _ R le_rfl (by have := muLowQ_le hP.M hP.g hP.lp 0; simpa [muQ] using this)
    (by simpa using jLowQ_le hP (le_refl (0 : ℚ))) hcells
  rw [Rat.cast_zero] at h0
  exact weilPos_of_affine hP.ha (toR4 cs) h0.1 h0.2 (tailQ_sound hP htail) K (kQ4_sound hK) (entR Lm)
    hpsd (by exact_mod_cast hε) (by exact_mod_cast hα) hβ' hn' hL hL1

end Fayed.Yoshida.AReflect
