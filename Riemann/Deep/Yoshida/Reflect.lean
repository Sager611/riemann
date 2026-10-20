/-
Riemann/Deep/Yoshida/Reflect.lean — rational certificate checkers (evaluated by `decide +kernel`)
and their soundness, for the cell and moment certificates of `Moment.lean`.

* `muLowQ`: rounded-down digamma partial sum + telescoping tail (`two_pi_mu_ge_tel`);
* `cosLowQ`/`sinLowQ`: range reduction by `2kπ` (`π` to 20 digits) + alternating Taylor bounds;
* `cellQ`: one cell `W ≥ s − τ` on `lo ≤ |r| < hi`; `chainQ`, `tailQ`, `momQ`;
* `psd_of_cert`: `xᵀMx ≥ 0` from `UᵀMU` diagonally dominant and `UW = I`.
All checks are pure rational arithmetic; no `native_decide`.
-/
import Riemann.Deep.Yoshida.Moment
import Riemann.Deep.Yoshida.MuTail

noncomputable section

open Finset

namespace Fayed.Yoshida.Reflect

open Fayed.Yoshida Fayed.Yoshida.Taylor

/-- Rational literal from numerator and denominator (cheap to elaborate in long data lists). -/
def Q (n : ℤ) (d : ℕ) : ℚ := (n : ℚ) / (d : ℚ)

/-! ### Rounding -/

def D40 : ℚ := 1099511627776

/-- Round down to a multiple of `2⁻⁴⁰`. -/
def flD (x : ℚ) : ℚ := ((⌊x * D40⌋ : ℤ) : ℚ) / D40

lemma flD_le (x : ℚ) : flD x ≤ x := by
  unfold flD
  rw [div_le_iff₀ (by norm_num [D40])]
  exact Int.floor_le _

lemma cast_flD_le (x : ℚ) : ((flD x : ℚ) : ℝ) ≤ (x : ℝ) := by exact_mod_cast flD_le x

/-! ### μ -/

def fQ (t2 : ℚ) (n : ℕ) : ℚ := 1 / ((n : ℚ) + 1) - ((n : ℚ) + 1 + 1 / 4) / (((n : ℚ) + 1 + 1 / 4) ^ 2 + t2)

def h1Q (n : ℕ) : ℚ := (1 / 4) / ((n : ℚ) + 5 / 4)

def h2Q (t2 : ℚ) (n : ℕ) : ℚ := (t2 / 2) / (((n : ℚ) + 3 / 4) ^ 2 + t2)

def sumFl (t2 : ℚ) : ℕ → ℚ
  | 0 => 0
  | n + 1 => sumFl t2 n + flD (fQ t2 n)

def muLowQ (g lp : ℚ) (M Mp : ℕ) (r : ℚ) : ℚ :=
  -g - lp - (1 / 4) / ((1 / 4) ^ 2 + (r / 2) ^ 2) + sumFl ((r / 2) ^ 2) M
    + (h1Q M - h1Q Mp + (if 1 / 2 ≤ (r / 2) ^ 2 then h2Q ((r / 2) ^ 2) M - h2Q ((r / 2) ^ 2) Mp else 0))

lemma cast_fQ (r : ℚ) (n : ℕ) : ((fQ ((r / 2) ^ 2) n : ℚ) : ℝ) = fT (r : ℝ) n := by
  unfold fQ fT
  push_cast
  ring

lemma sumFl_le (r : ℚ) : ∀ M : ℕ, ((sumFl ((r / 2) ^ 2) M : ℚ) : ℝ) ≤ ∑ n ∈ range M, fT (r : ℝ) n
  | 0 => by simp [sumFl]
  | n + 1 => by
    rw [sumFl, Finset.sum_range_succ]
    push_cast
    have h1 := sumFl_le r n
    have h2 := cast_flD_le (fQ ((r / 2) ^ 2) n)
    rw [cast_fQ] at h2
    linarith

theorem muLowQ_le {g lp : ℚ} {M Mp : ℕ} (hM : M ≤ Mp)
    (hg : Real.eulerMascheroniConstant ≤ (g : ℝ)) (hl : Real.log Real.pi ≤ (lp : ℝ)) (r : ℚ) :
    ((muLowQ g lp M Mp r : ℚ) : ℝ) ≤ 2 * Real.pi * Zeta23.mu (r : ℝ) := by
  have h := two_pi_mu_ge_tel (r : ℝ) M Mp hM hg hl
  refine le_trans ?_ h
  unfold muLowQ
  have hs := sumFl_le r M
  have e1 : ∀ n : ℕ, ((h1Q n : ℚ) : ℝ) = h1R n := by
    intro n; unfold h1Q h1R; push_cast; ring
  have e2 : ∀ n : ℕ, ((h2Q ((r / 2) ^ 2) n : ℚ) : ℝ) = h2R (((r : ℝ) / 2) ^ 2) n := by
    intro n; unfold h2Q h2R; push_cast; ring
  have hc : (1 / 2 ≤ (r / 2) ^ 2) ↔ (1 / 2 ≤ ((r : ℝ) / 2) ^ 2) := by
    rw [← Rat.cast_le (K := ℝ)]
    push_cast
    rfl
  push_cast
  rw [e1, e1]
  by_cases hq : (1 : ℚ) / 2 ≤ (r / 2) ^ 2
  · rw [ite_cond_eq_true _ _ (eq_true hq), ite_cond_eq_true _ _ (eq_true (hc.1 hq))]
    push_cast
    rw [e2, e2]
    linarith
  · rw [ite_cond_eq_false _ _ (eq_false hq), ite_cond_eq_false _ _ (eq_false (fun h => hq (hc.2 h)))]
    push_cast
    linarith

/-! ### Trigonometric lower bounds -/

def cosTQ (n : ℕ) (x : ℚ) : ℚ := ∑ k ∈ range n, (-1 : ℚ) ^ k * x ^ (2 * k) / ((2 * k).factorial : ℚ)

def sinTQ (n : ℕ) (x : ℚ) : ℚ :=
  ∑ k ∈ range n, (-1 : ℚ) ^ k * x ^ (2 * k + 1) / ((2 * k + 1).factorial : ℚ)

lemma cast_cosTQ (n : ℕ) (x : ℚ) : ((cosTQ n x : ℚ) : ℝ) = cosT n (x : ℝ) := by
  unfold cosTQ cosT; push_cast; rfl

lemma cast_sinTQ (n : ℕ) (x : ℚ) : ((sinTQ n x : ℚ) : ℝ) = sinT n (x : ℝ) := by
  unfold sinTQ sinT; push_cast; rfl

def piLo : ℚ := 314159265358979323846 / 100000000000000000000

def piHi : ℚ := 314159265358979323847 / 100000000000000000000

lemma piLo_lt : ((piLo : ℚ) : ℝ) < Real.pi := by
  have h := Real.pi_gt_d20
  unfold piLo
  push_cast
  norm_num at h ⊢
  linarith

lemma pi_lt_piHi : Real.pi < ((piHi : ℚ) : ℝ) := by
  have h := Real.pi_lt_d20
  unfold piHi
  push_cast
  norm_num at h ⊢
  linarith

/-- Number of periods removed. -/
def redK (θ : ℚ) : ℕ := (⌊θ / (2 * piLo) + 1 / 2⌋).toNat

/-- Taylor level: `cosT 10`, `sinT 10` (lower) and `sinT 9` (upper); `|y| ≤ π` after reduction. -/
def TM : ℕ := 4

def cosLowQ (θ : ℚ) : ℚ :=
  flD (cosTQ (2 * TM + 2) |θ - 2 * (redK θ : ℚ) * piLo|) - 2 * (redK θ : ℚ) * (piHi - piLo)

def sinLowQ (θ : ℚ) : ℚ :=
  (if 0 ≤ θ - 2 * (redK θ : ℚ) * piLo then flD (sinTQ (2 * TM + 2) (θ - 2 * (redK θ : ℚ) * piLo))
    else flD (-sinTQ (2 * TM + 1) (-(θ - 2 * (redK θ : ℚ) * piLo)))) - 2 * (redK θ : ℚ) * (piHi - piLo)

/-- `|u − y| ≤ 2k(π̄ − π̲)` for `u = θ − 2kπ`, `y = θ − 2kπ̲`. -/
lemma red_close (θ : ℚ) :
    |((θ : ℝ) - (redK θ : ℝ) * (2 * Real.pi)) - (((θ - 2 * (redK θ : ℚ) * piLo : ℚ)) : ℝ)|
      ≤ (((2 * (redK θ : ℚ) * (piHi - piLo) : ℚ)) : ℝ) := by
  have h1 := piLo_lt
  have h2 := pi_lt_piHi
  have hk : (0 : ℝ) ≤ (redK θ : ℝ) := Nat.cast_nonneg _
  push_cast
  rw [abs_le]
  constructor <;> nlinarith

theorem cosLowQ_le (θ : ℚ) : ((cosLowQ θ : ℚ) : ℝ) ≤ Real.cos (θ : ℝ) := by
  have hr := red_close θ
  have hc := Real.abs_cos_sub_cos_le ((θ : ℝ) - (redK θ : ℝ) * (2 * Real.pi))
    (((θ - 2 * (redK θ : ℚ) * piLo : ℚ)) : ℝ)
  rw [Real.cos_sub_nat_mul_two_pi] at hc
  have h4 := abs_le.1 (hc.trans hr)
  have hT := cos_ge TM (abs_nonneg (((θ - 2 * (redK θ : ℚ) * piLo : ℚ)) : ℝ))
  rw [Real.cos_abs] at hT
  have hf := cast_flD_le (cosTQ (2 * TM + 2) |θ - 2 * (redK θ : ℚ) * piLo|)
  rw [cast_cosTQ, Rat.cast_abs] at hf
  unfold cosLowQ
  rw [Rat.cast_sub]
  linarith [h4.1]

theorem sinLowQ_le (θ : ℚ) : ((sinLowQ θ : ℚ) : ℝ) ≤ Real.sin (θ : ℝ) := by
  have hr := red_close θ
  have hc := Real.abs_sin_sub_sin_le ((θ : ℝ) - (redK θ : ℝ) * (2 * Real.pi))
    (((θ - 2 * (redK θ : ℚ) * piLo : ℚ)) : ℝ)
  rw [Real.sin_sub_nat_mul_two_pi] at hc
  have h4 := abs_le.1 (hc.trans hr)
  have hy : ∀ y : ℚ, ((if 0 ≤ y then flD (sinTQ (2 * TM + 2) y) else flD (-sinTQ (2 * TM + 1) (-y)) : ℚ) : ℝ)
      ≤ Real.sin (y : ℝ) := by
    intro y
    split_ifs with h
    · have hf := cast_flD_le (sinTQ (2 * TM + 2) y)
      rw [cast_sinTQ] at hf
      have hT := sin_ge TM (show (0 : ℝ) ≤ (y : ℝ) by exact_mod_cast h)
      linarith
    · have hf := cast_flD_le (-sinTQ (2 * TM + 1) (-y))
      rw [Rat.cast_neg, cast_sinTQ, Rat.cast_neg] at hf
      have hy0 : (0 : ℝ) ≤ -((y : ℚ) : ℝ) := by
        have : y < 0 := lt_of_not_ge h
        have : ((y : ℚ) : ℝ) < 0 := by exact_mod_cast this
        linarith
      have hT := (bounds TM).1 _ hy0
      rw [Real.sin_neg] at hT
      linarith
  have hy' := hy (θ - 2 * (redK θ : ℚ) * piLo)
  unfold sinLowQ
  rw [Rat.cast_sub]
  linarith [h4.1]

/-! ### Cells -/

/-- Certificate parameters (all rational). -/
structure Params where
  a : ℚ
  s : ℚ
  shlo : ℚ
  shhi : ℚ
  chlo : ℚ
  chhi : ℚ
  Lam : ℚ
  g : ℚ
  lp : ℚ
  M : ℕ
  Mp : ℕ

/-- The real facts behind the parameters. -/
structure Valid (P : Params) : Prop where
  ha : 0 < (P.a : ℝ)
  shlo : (P.shlo : ℝ) ≤ Real.sinh P.a
  shhi : Real.sinh P.a ≤ P.shhi
  chlo : (P.chlo : ℝ) ≤ Real.cosh P.a
  chhi : Real.cosh P.a ≤ P.chhi
  lam : 8 * (P.a : ℝ) ^ 2 * Real.cosh P.a ≤ P.Lam
  g : Real.eulerMascheroniConstant ≤ P.g
  lp : Real.log Real.pi ≤ P.lp
  M : P.M ≤ P.Mp

def jLowQ (P : Params) (p : ℚ) : ℚ :=
  (fun cl sl : ℚ => 4 * (1 / 2 * (if 0 ≤ cl then P.shlo * cl else P.shhi * cl)
    + p * (if 0 ≤ sl then P.chlo * sl else P.chhi * sl)) / (1 / 4 + p ^ 2))
    (cosLowQ (2 * P.a * p)) (sinLowQ (2 * P.a * p))

theorem jLowQ_le {P : Params} (hP : Valid P) {p : ℚ} (hp : 0 ≤ p) :
    ((jLowQ P p : ℚ) : ℝ) ≤ J (P.a : ℝ) (p : ℝ) := by
  have hc := cosLowQ_le (2 * P.a * p)
  have hs := sinLowQ_le (2 * P.a * p)
  push_cast at hc hs
  unfold jLowQ
  set cl := cosLowQ (2 * P.a * p)
  set sl := sinLowQ (2 * P.a * p)
  have hp' : (0 : ℝ) ≤ (p : ℝ) := by exact_mod_cast hp
  have hX : (((if 0 ≤ cl then P.shlo * cl else P.shhi * cl : ℚ)) : ℝ) ≤ Real.sinh P.a * cl := by
    split_ifs with h
    · push_cast
      exact mul_le_mul_of_nonneg_right hP.shlo (by exact_mod_cast h)
    · push_cast
      exact mul_le_mul_of_nonpos_right hP.shhi (by
        have : cl < 0 := lt_of_not_ge h
        exact_mod_cast this.le)
  have hY : (((if 0 ≤ sl then P.chlo * sl else P.chhi * sl : ℚ)) : ℝ) ≤ Real.cosh P.a * sl := by
    split_ifs with h
    · push_cast
      exact mul_le_mul_of_nonneg_right hP.chlo (by exact_mod_cast h)
    · push_cast
      exact mul_le_mul_of_nonpos_right hP.chhi (by
        have : sl < 0 := lt_of_not_ge h
        exact_mod_cast this.le)
  have := J_ge_of hP.ha.le hp' hc hs hX hY
  simp only
  push_cast
  push_cast at this
  exact this

/-- One cell check. -/
def cellQ (P : Params) (c : ℚ × ℚ × ℚ) : Bool :=
  decide (0 ≤ c.1) && decide (c.1 ≤ c.2.1) && decide (0 ≤ c.2.2) &&
    decide (P.s - c.2.2 ≤ muLowQ P.g P.lp P.M P.Mp c.1 + jLowQ P ((c.1 + c.2.1) / 2)
      - P.Lam * ((c.2.1 - c.1) / 2))

def toR (cs : List (ℚ × ℚ × ℚ)) : List (ℝ × ℝ × ℝ) :=
  cs.map fun c => ((c.1 : ℝ), (c.2.1 : ℝ), (c.2.2 : ℝ))

theorem cellQ_sound {P : Params} (hP : Valid P) {c : ℚ × ℚ × ℚ} (h : cellQ P c = true) :
    CellOK (P.a : ℝ) (P.s : ℝ) ((c.1 : ℝ), (c.2.1 : ℝ), (c.2.2 : ℝ)) := by
  unfold cellQ at h
  simp only [Bool.and_eq_true, decide_eq_true_eq] at h
  obtain ⟨⟨⟨h0, h1⟩, h2⟩, h3⟩ := h
  have hp : 0 ≤ (c.1 + c.2.1) / 2 := by linarith
  have hm := muLowQ_le hP.M hP.g hP.lp c.1
  have hj := jLowQ_le hP hp
  have h3' : ((P.s - c.2.2 : ℚ) : ℝ) ≤ ((muLowQ P.g P.lp P.M P.Mp c.1 + jLowQ P ((c.1 + c.2.1) / 2)
      - P.Lam * ((c.2.1 - c.1) / 2) : ℚ) : ℝ) := by exact_mod_cast h3
  push_cast at h3' hj
  refine cellOK_of (p := ((c.1 : ℝ) + c.2.1) / 2) (d := ((c.2.1 : ℝ) - c.1) / 2) hP.ha.le
    (by exact_mod_cast h0) (by exact_mod_cast h1) (by exact_mod_cast h2) (by linarith) (by linarith)
    hm hj hP.lam ?_
  linarith

theorem cells_sound {P : Params} (hP : Valid P) {cs : List (ℚ × ℚ × ℚ)} (h : cs.all (cellQ P) = true) :
    ∀ c ∈ toR cs, CellOK (P.a : ℝ) (P.s : ℝ) c := by
  intro c hc
  unfold toR at hc
  obtain ⟨c', hc', rfl⟩ := List.mem_map.1 hc
  exact cellQ_sound hP (List.all_eq_true.1 h c' hc')

def chainQ : ℚ → List (ℚ × ℚ × ℚ) → ℚ → Bool
  | x, [], R => decide (x = R)
  | x, c :: cs, R => decide (c.1 = x) && chainQ c.2.1 cs R

theorem chainQ_sound : ∀ (cs : List (ℚ × ℚ × ℚ)) (x R : ℚ), chainQ x cs R = true →
    Chain (x : ℝ) (toR cs) (R : ℝ)
  | [], x, R, h => by
    simp only [chainQ, decide_eq_true_eq] at h
    simp [toR, Chain, h]
  | c :: cs, x, R, h => by
    simp only [chainQ, Bool.and_eq_true, decide_eq_true_eq] at h
    have ih := chainQ_sound cs c.2.1 R h.2
    simp only [toR, List.map_cons, Chain]
    exact ⟨by rw [h.1], ih⟩

def tailQ (P : Params) (R : ℚ) : Bool :=
  decide (1 ≤ R) && decide (P.s ≤ muLowQ P.g P.lp P.M P.Mp R - 4 * (P.shhi / (2 * R ^ 2) + P.chhi / R))

theorem tailQ_sound {P : Params} (hP : Valid P) {R : ℚ} (h : tailQ P R = true) :
    ∀ r : ℝ, (R : ℝ) ≤ |r| → (P.s : ℝ) ≤ W (P.a : ℝ) r := by
  unfold tailQ at h
  simp only [Bool.and_eq_true, decide_eq_true_eq] at h
  have h2 : ((P.s : ℚ) : ℝ) ≤ ((muLowQ P.g P.lp P.M P.Mp R - 4 * (P.shhi / (2 * R ^ 2) + P.chhi / R) : ℚ) : ℝ) := by
    exact_mod_cast h.2
  push_cast at h2
  exact tail_of hP.ha.le (by exact_mod_cast h.1) (muLowQ_le hP.M hP.g hP.lp R) hP.shhi hP.chhi h2

def kQ (a : ℚ) (K : ℕ) (cs : List (ℚ × ℚ × ℚ)) : Bool :=
  cs.all fun c => decide (c.2.2 = 0) || decide (c.2.1 * a ≤ ((K : ℚ) + 1) / 2)

theorem kQ_sound {a : ℚ} {K : ℕ} {cs : List (ℚ × ℚ × ℚ)} (h : kQ a K cs = true) :
    ∀ c ∈ toR cs, c.2.2 = 0 ∨ c.2.1 * (a : ℝ) ≤ (K + 1) / 2 := by
  intro c hc
  unfold toR at hc
  obtain ⟨c', hc', rfl⟩ := List.mem_map.1 hc
  have := List.all_eq_true.1 h c' hc'
  simp only [Bool.or_eq_true, decide_eq_true_eq] at this
  rcases this with h1 | h1
  · left; simp [h1]
  · right
    have : ((c'.2.1 * a : ℚ) : ℝ) ≤ ((((K : ℚ) + 1) / 2 : ℚ) : ℝ) := by exact_mod_cast h1
    push_cast at this
    exact this

/-! ### Moments -/

def cellMomQ (lo hi : ℚ) (p : ℕ) : ℚ :=
  ((hi ^ (p + 1) - lo ^ (p + 1)) + ((-lo) ^ (p + 1) - (-hi) ^ (p + 1))) / ((p : ℚ) + 1)

def momQ : List (ℚ × ℚ × ℚ) → ℕ → ℚ
  | [], _ => 0
  | c :: cs, p => (if c.2.2 = 0 then 0 else c.2.2 * cellMomQ c.1 c.2.1 p) + momQ cs p

theorem cast_momQ : ∀ (cs : List (ℚ × ℚ × ℚ)) (p : ℕ), ((momQ cs p : ℚ) : ℝ) = momSum (toR cs) p
  | [], p => by simp [momQ, momSum, toR]
  | c :: cs, p => by
    simp only [momQ, toR, List.map_cons, momSum]
    have ih := cast_momQ cs p
    simp only [toR] at ih
    rw [← ih]
    split_ifs with h
    · simp [h]
    · unfold cellMomQ cellMom
      push_cast
      ring

/-! ### Matrices (lists of rows, ℕ-indexed) -/

def ent (A : List (List ℚ)) (i j : ℕ) : ℚ := (A.getD i []).getD j 0

def entR (A : List (List ℚ)) (i j : ℕ) : ℝ := ((ent A i j : ℚ) : ℝ)

/-- Check `C = A·B` on `K × K`. -/
def mulQ (K : ℕ) (A B C : List (List ℚ)) : Bool :=
  (List.range K).all fun i => (List.range K).all fun j =>
    decide (ent C i j = ∑ l ∈ range K, ent A i l * ent B l j)

/-- Check `C = Aᵀ·B` on `K × K`. -/
def tmulQ (K : ℕ) (A B C : List (List ℚ)) : Bool :=
  (List.range K).all fun i => (List.range K).all fun j =>
    decide (ent C i j = ∑ l ∈ range K, ent A l i * ent B l j)

theorem mulQ_sound {K : ℕ} {A B C : List (List ℚ)} (h : mulQ K A B C = true) :
    ∀ i < K, ∀ j < K, entR C i j = ∑ l ∈ range K, entR A i l * entR B l j := by
  intro i hi j hj
  unfold mulQ at h
  have := List.all_eq_true.1 (List.all_eq_true.1 h i (List.mem_range.2 hi)) j (List.mem_range.2 hj)
  simp only [decide_eq_true_eq] at this
  unfold entR
  rw [this]
  push_cast
  rfl

theorem tmulQ_sound {K : ℕ} {A B C : List (List ℚ)} (h : tmulQ K A B C = true) :
    ∀ i < K, ∀ j < K, entR C i j = ∑ l ∈ range K, entR A l i * entR B l j := by
  intro i hi j hj
  unfold tmulQ at h
  have := List.all_eq_true.1 (List.all_eq_true.1 h i (List.mem_range.2 hi)) j (List.mem_range.2 hj)
  simp only [decide_eq_true_eq] at this
  unfold entR
  rw [this]
  push_cast
  rfl

/-- Check `U·W = I`. -/
def invQ (K : ℕ) (U W : List (List ℚ)) : Bool :=
  (List.range K).all fun i => (List.range K).all fun j =>
    decide (∑ l ∈ range K, ent U i l * ent W l j = if i = j then 1 else 0)

theorem invQ_sound {K : ℕ} {U W : List (List ℚ)} (h : invQ K U W = true) :
    ∀ i < K, ∀ j < K, ∑ l ∈ range K, entR U i l * entR W l j = if i = j then 1 else 0 := by
  intro i hi j hj
  unfold invQ at h
  have := List.all_eq_true.1 (List.all_eq_true.1 h i (List.mem_range.2 hi)) j (List.mem_range.2 hj)
  simp only [decide_eq_true_eq] at this
  unfold entR
  have h2 : ((∑ l ∈ range K, ent U i l * ent W l j : ℚ) : ℝ) = (((if i = j then 1 else 0 : ℚ)) : ℝ) := by
    rw [this]
  push_cast at h2
  rw [h2]
  split_ifs <;> simp

/-- Check symmetric diagonal dominance. -/
def domQ (K : ℕ) (D : List (List ℚ)) : Bool :=
  (List.range K).all fun i =>
    decide (∑ l ∈ range K, (if l = i then 0 else (|ent D i l| + |ent D l i|) / 2) ≤ ent D i i)

theorem domQ_sound {K : ℕ} {D : List (List ℚ)} (h : domQ K D = true) :
    ∀ i < K, ∑ l ∈ range K, (if l = i then 0 else (|entR D i l| + |entR D l i|) / 2) ≤ entR D i i := by
  intro i hi
  unfold domQ at h
  have := List.all_eq_true.1 h i (List.mem_range.2 hi)
  simp only [decide_eq_true_eq] at this
  unfold entR
  have h2 : ((∑ l ∈ range K, (if l = i then 0 else (|ent D i l| + |ent D l i|) / 2) : ℚ) : ℝ)
      ≤ ((ent D i i : ℚ) : ℝ) := by exact_mod_cast this
  push_cast at h2
  refine le_trans (le_of_eq ?_) h2
  refine Finset.sum_congr rfl fun l _ => ?_
  split_ifs <;> push_cast <;> rfl

/-! ### Real PSD lemmas -/

theorem qf_nonneg_of_dom (K : ℕ) (D : ℕ → ℕ → ℝ)
    (hD : ∀ i < K, ∑ l ∈ range K, (if l = i then 0 else (|D i l| + |D l i|) / 2) ≤ D i i)
    (y : ℕ → ℝ) : 0 ≤ ∑ i ∈ range K, ∑ j ∈ range K, D i j * y i * y j := by
  have hpt : ∀ i j, (if j = i then D i i * y i ^ 2 else -(|D i j| * (y i ^ 2 + y j ^ 2) / 2))
      ≤ D i j * y i * y j := by
    intro i j
    split_ifs with h
    · subst h; ring_nf; rfl
    · have h1 : -(|D i j| * |y i * y j|) ≤ D i j * (y i * y j) := by
        have := neg_abs_le (D i j * (y i * y j))
        rw [abs_mul] at this
        exact this
      have h2 : |y i * y j| ≤ (y i ^ 2 + y j ^ 2) / 2 := by
        rw [abs_mul]
        nlinarith [sq_nonneg (|y i| - |y j|), sq_abs (y i), sq_abs (y j)]
      have h3 := mul_le_mul_of_nonneg_left h2 (abs_nonneg (D i j))
      nlinarith
  have hsum : ∑ i ∈ range K, ∑ j ∈ range K,
      (if j = i then D i i * y i ^ 2 else -(|D i j| * (y i ^ 2 + y j ^ 2) / 2))
      ≤ ∑ i ∈ range K, ∑ j ∈ range K, D i j * y i * y j :=
    Finset.sum_le_sum fun i _ => Finset.sum_le_sum fun j _ => hpt i j
  refine le_trans ?_ hsum
  -- rewrite the lower bound as Σ_i y_i² (D_ii − Σ_{l≠i} (|D_il| + |D_li|)/2)
  have e : ∑ i ∈ range K, ∑ j ∈ range K,
      (if j = i then D i i * y i ^ 2 else -(|D i j| * (y i ^ 2 + y j ^ 2) / 2))
      = ∑ i ∈ range K, y i ^ 2 * (D i i
          - ∑ l ∈ range K, (if l = i then 0 else (|D i l| + |D l i|) / 2)) := by
    have e1 : ∀ i j, (if j = i then D i i * y i ^ 2 else -(|D i j| * (y i ^ 2 + y j ^ 2) / 2))
        = (if j = i then D i i * y i ^ 2 else 0) - (if j = i then 0 else |D i j| * y i ^ 2 / 2)
          - (if j = i then 0 else |D i j| * y j ^ 2 / 2) := by
      intro i j; split_ifs <;> ring
    simp_rw [e1, Finset.sum_sub_distrib]
    have e2 : ∑ i ∈ range K, ∑ j ∈ range K, (if j = i then (0 : ℝ) else |D i j| * y j ^ 2 / 2)
        = ∑ i ∈ range K, ∑ j ∈ range K, (if i = j then (0 : ℝ) else |D j i| * y i ^ 2 / 2) := by
      rw [Finset.sum_comm]
    rw [e2]
    rw [← Finset.sum_sub_distrib, ← Finset.sum_sub_distrib]
    refine Finset.sum_congr rfl fun i hi => ?_
    rw [Finset.sum_ite_eq' (range K) i (fun _ => D i i * y i ^ 2), if_pos hi, mul_sub,
      Finset.mul_sum]
    rw [show D i i * y i ^ 2 - ∑ j ∈ range K, (if j = i then (0 : ℝ) else |D i j| * y i ^ 2 / 2)
        - ∑ j ∈ range K, (if i = j then (0 : ℝ) else |D j i| * y i ^ 2 / 2)
        = D i i * y i ^ 2 - ∑ j ∈ range K, ((if j = i then (0 : ℝ) else |D i j| * y i ^ 2 / 2)
          + (if i = j then (0 : ℝ) else |D j i| * y i ^ 2 / 2)) by rw [Finset.sum_add_distrib]; ring]
    rw [mul_comm (y i ^ 2) (D i i)]
    congr 1
    refine Finset.sum_congr rfl fun l _ => ?_
    by_cases h : l = i
    · subst h; simp
    · rw [if_neg h, if_neg (Ne.symm h), if_neg h]
      ring
  rw [e]
  refine Finset.sum_nonneg fun i hi => mul_nonneg (sq_nonneg _) ?_
  linarith [hD i (List.mem_range.1 (by simpa using hi))]

/-- Reorder a four-fold finite sum `Σ_i Σ_j Σ_k Σ_l` as `Σ_k Σ_l Σ_i Σ_j`. -/
lemma sum4_comm (K : ℕ) (F : ℕ → ℕ → ℕ → ℕ → ℝ) :
    ∑ i ∈ range K, ∑ j ∈ range K, ∑ k ∈ range K, ∑ l ∈ range K, F i j k l
      = ∑ k ∈ range K, ∑ l ∈ range K, ∑ i ∈ range K, ∑ j ∈ range K, F i j k l := by
  calc ∑ i ∈ range K, ∑ j ∈ range K, ∑ k ∈ range K, ∑ l ∈ range K, F i j k l
      = ∑ i ∈ range K, ∑ k ∈ range K, ∑ j ∈ range K, ∑ l ∈ range K, F i j k l :=
        Finset.sum_congr rfl fun i _ => Finset.sum_comm
    _ = ∑ k ∈ range K, ∑ i ∈ range K, ∑ j ∈ range K, ∑ l ∈ range K, F i j k l := Finset.sum_comm
    _ = ∑ k ∈ range K, ∑ i ∈ range K, ∑ l ∈ range K, ∑ j ∈ range K, F i j k l :=
        Finset.sum_congr rfl fun k _ => Finset.sum_congr rfl fun i _ => Finset.sum_comm
    _ = ∑ k ∈ range K, ∑ l ∈ range K, ∑ i ∈ range K, ∑ j ∈ range K, F i j k l :=
        Finset.sum_congr rfl fun k _ => Finset.sum_comm

/-- Congruence: `xᵀMx = yᵀ(UᵀMU)y` with `x = Uy`, `y = Wx`, and dominance. -/
theorem qf_nonneg_of_cert (K : ℕ) (M U Wm MU D : ℕ → ℕ → ℝ)
    (hMU : ∀ i < K, ∀ j < K, MU i j = ∑ l ∈ range K, M i l * U l j)
    (hD : ∀ i < K, ∀ j < K, D i j = ∑ l ∈ range K, U l i * MU l j)
    (hUW : ∀ i < K, ∀ j < K, ∑ l ∈ range K, U i l * Wm l j = if i = j then 1 else 0)
    (hdom : ∀ i < K, ∑ l ∈ range K, (if l = i then 0 else (|D i l| + |D l i|) / 2) ≤ D i i)
    (x : ℕ → ℝ) : 0 ≤ ∑ i ∈ range K, ∑ j ∈ range K, M i j * x i * x j := by
  set y : ℕ → ℝ := fun j => ∑ l ∈ range K, Wm j l * x l
  have hx : ∀ i ∈ range K, ∑ k ∈ range K, U i k * y k = x i := by
    intro i hi
    have hi' := List.mem_range.1 (by simpa using hi)
    simp only [y, Finset.mul_sum]
    rw [Finset.sum_comm]
    have e : ∀ l ∈ range K, ∑ k ∈ range K, U i k * (Wm k l * x l) = (if i = l then 1 else 0) * x l := by
      intro l hl
      rw [← hUW i hi' l (List.mem_range.1 (by simpa using hl)), Finset.sum_mul]
      refine Finset.sum_congr rfl fun k _ => ?_
      ring
    rw [Finset.sum_congr rfl e]
    simp [Finset.sum_ite_eq, hi]
  clear_value y
  have h0 := qf_nonneg_of_dom K D hdom y
  refine le_of_le_of_eq h0 ?_
  have e1 : ∑ i ∈ range K, ∑ j ∈ range K, M i j * x i * x j
      = ∑ i ∈ range K, ∑ j ∈ range K, M i j * (∑ k ∈ range K, U i k * y k)
          * (∑ l ∈ range K, U j l * y l) := by
    refine Finset.sum_congr rfl fun i hi => Finset.sum_congr rfl fun j hj => ?_
    rw [hx i hi, hx j hj]
  rw [e1]
  have e2 : ∀ i j, M i j * (∑ k ∈ range K, U i k * y k) * (∑ l ∈ range K, U j l * y l)
      = ∑ k ∈ range K, ∑ l ∈ range K, U i k * M i j * U j l * y k * y l := by
    intro i j
    rw [mul_assoc, Finset.sum_mul_sum, Finset.mul_sum]
    refine Finset.sum_congr rfl fun k _ => ?_
    rw [Finset.mul_sum]
    refine Finset.sum_congr rfl fun l _ => ?_
    ring
  simp_rw [e2]
  rw [sum4_comm K (fun i j k l => U i k * M i j * U j l * y k * y l)]
  refine Finset.sum_congr rfl fun k hk => Finset.sum_congr rfl fun l hl => ?_
  rw [hD k (List.mem_range.1 (by simpa using hk)) l (List.mem_range.1 (by simpa using hl)),
    Finset.sum_mul, Finset.sum_mul]
  refine Finset.sum_congr rfl fun i hi => ?_
  rw [hMU i (List.mem_range.1 (by simpa using hi)) l (List.mem_range.1 (by simpa using hl)),
    Finset.mul_sum, Finset.sum_mul, Finset.sum_mul]
  refine Finset.sum_congr rfl fun j _ => ?_
  ring

/-- The moment-PSD hypothesis of `weilPos_of_moment` from `xᵀMx ≥ 0`,
`M = α(Lm + Lmᵀ − LmᵀĜLm) − H`. -/
theorem hpsd_of (K : ℕ) (α : ℝ) (Lm Gm LGm M G : ℕ → ℕ → ℝ) (H : ℕ → ℝ)
    (hGm : ∀ k < K, ∀ j < K, Gm k j = ∑ l ∈ range K, G k l * Lm l j)
    (hLGm : ∀ i < K, ∀ j < K, LGm i j = ∑ k ∈ range K, Lm k i * Gm k j)
    (hM : ∀ i < K, ∀ j < K, M i j = α * (Lm i j + Lm j i - LGm i j) - H (i + j))
    (hpos : ∀ x : ℕ → ℝ, 0 ≤ ∑ i ∈ range K, ∑ j ∈ range K, M i j * x i * x j) (x : ℕ → ℝ) :
    ∑ m ∈ range K, ∑ k ∈ range K, x m * x k * H (m + k)
      ≤ α * (2 * ∑ m ∈ range K, (∑ i ∈ range K, Lm m i * x i) * x m
        - ∑ m ∈ range K, ∑ k ∈ range K,
            ((∑ i ∈ range K, Lm m i * x i) * (∑ i ∈ range K, Lm k i * x i)) * G m k) := by
  have h := hpos x
  have eM : ∑ i ∈ range K, ∑ j ∈ range K, M i j * x i * x j
      = α * (∑ i ∈ range K, ∑ j ∈ range K, Lm i j * x i * x j
          + ∑ i ∈ range K, ∑ j ∈ range K, Lm j i * x i * x j
          - ∑ i ∈ range K, ∑ j ∈ range K, LGm i j * x i * x j)
        - ∑ i ∈ range K, ∑ j ∈ range K, x i * x j * H (i + j) := by
    have e0 : ∑ i ∈ range K, ∑ j ∈ range K, M i j * x i * x j
        = ∑ i ∈ range K, ∑ j ∈ range K, (α * (Lm i j * x i * x j) + α * (Lm j i * x i * x j)
          - α * (LGm i j * x i * x j) - x i * x j * H (i + j)) := by
      refine Finset.sum_congr rfl fun i hi => Finset.sum_congr rfl fun j hj => ?_
      rw [hM i (List.mem_range.1 (by simpa using hi)) j (List.mem_range.1 (by simpa using hj))]
      ring
    rw [e0, mul_sub, mul_add]
    simp only [Finset.mul_sum, Finset.sum_add_distrib, Finset.sum_sub_distrib]
  have e1 : ∑ i ∈ range K, ∑ j ∈ range K, Lm j i * x i * x j
      = ∑ i ∈ range K, ∑ j ∈ range K, Lm i j * x i * x j := by
    rw [Finset.sum_comm]
    refine Finset.sum_congr rfl fun i _ => Finset.sum_congr rfl fun j _ => ?_
    ring
  have e2 : ∑ i ∈ range K, ∑ j ∈ range K, Lm i j * x i * x j
      = ∑ m ∈ range K, (∑ i ∈ range K, Lm m i * x i) * x m := by
    refine Finset.sum_congr rfl fun m _ => ?_
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl fun i _ => ?_
    ring
  have e3 : ∑ i ∈ range K, ∑ j ∈ range K, LGm i j * x i * x j
      = ∑ m ∈ range K, ∑ k ∈ range K,
          ((∑ i ∈ range K, Lm m i * x i) * (∑ i ∈ range K, Lm k i * x i)) * G m k := by
    have e4 : ∀ i ∈ range K, ∀ j ∈ range K, LGm i j * x i * x j
        = ∑ m ∈ range K, ∑ k ∈ range K, Lm m i * G m k * Lm k j * x i * x j := by
      intro i hi j hj
      rw [hLGm i (List.mem_range.1 (by simpa using hi)) j (List.mem_range.1 (by simpa using hj)),
        Finset.sum_mul, Finset.sum_mul]
      refine Finset.sum_congr rfl fun m hm => ?_
      rw [hGm m (List.mem_range.1 (by simpa using hm)) j (List.mem_range.1 (by simpa using hj)),
        Finset.mul_sum, Finset.sum_mul, Finset.sum_mul]
      refine Finset.sum_congr rfl fun k _ => ?_
      ring
    rw [Finset.sum_congr rfl fun i hi => Finset.sum_congr rfl fun j hj => e4 i hi j hj]
    rw [sum4_comm K (fun i j m k => Lm m i * G m k * Lm k j * x i * x j)]
    refine Finset.sum_congr rfl fun m _ => Finset.sum_congr rfl fun k _ => ?_
    rw [Finset.sum_mul_sum, Finset.sum_mul]
    refine Finset.sum_congr rfl fun i _ => ?_
    rw [Finset.sum_mul]
    refine Finset.sum_congr rfl fun j _ => ?_
    ring
  rw [eM, e1, e2, e3] at h
  linarith

/-! ### The rational certificate theorem -/

def sgn2Q (m : ℕ) : ℚ := (-1) ^ (m / 2)

def ghatQ (a : ℚ) (m k : ℕ) : ℚ :=
  sgn2Q m * sgn2Q k * ((a ^ (m + k + 1) - (-a) ^ (m + k + 1)) / ((m + k : ℕ) + 1 : ℚ))
    / ((m.factorial : ℚ) * k.factorial)

lemma cast_ghatQ (a : ℚ) (m k : ℕ) : ((ghatQ a m k : ℚ) : ℝ) = Ghat (a : ℝ) m k := by
  unfold ghatQ Ghat sgn2Q sgn2
  push_cast
  ring

def momCheck (cs : List (ℚ × ℚ × ℚ)) (K : ℕ) (momD : List ℚ) : Bool :=
  (List.range (2 * K + 1)).all fun p => decide (momD.getD p 0 = momQ cs p)

def glQ (K : ℕ) (a : ℚ) (Lm Gm : List (List ℚ)) : Bool :=
  (List.range K).all fun k => (List.range K).all fun j =>
    decide (ent Gm k j = ∑ l ∈ range K, ghatQ a k l * ent Lm l j)

def mQ (K : ℕ) (α : ℚ) (Lm LGm Mm : List (List ℚ)) (momD : List ℚ) : Bool :=
  (List.range K).all fun i => (List.range K).all fun j =>
    decide (ent Mm i j = α * (ent Lm i j + ent Lm j i - ent LGm i j) - momD.getD (i + j) 0)

/-- Cells on a dyadic grid: `(Lo, Hi, T)` ↦ `(Lo/2^e, Hi/2^e, T/2^f)`. -/
def toQc (e f : ℕ) (c : ℕ × ℕ × ℕ) : ℚ × ℚ × ℚ :=
  ((c.1 : ℚ) / 2 ^ e, (c.2.1 : ℚ) / 2 ^ e, (c.2.2 : ℚ) / 2 ^ f)

/-- Integer moment sum `Σ T (Hi^{p+1} − Lo^{p+1})`. -/
def momN : List (ℕ × ℕ × ℕ) → ℕ → ℕ
  | [], _ => 0
  | c :: cs, p => c.2.2 * (c.2.1 ^ (p + 1) - c.1 ^ (p + 1)) + momN cs p

/-- Closed form of `momQ` on a dyadic grid. -/
def momF (cs : List (ℕ × ℕ × ℕ)) (e f p : ℕ) : ℚ :=
  if p % 2 = 0 then 2 * (momN cs p : ℚ) / (((p : ℚ) + 1) * 2 ^ (e * (p + 1)) * 2 ^ f) else 0

lemma cellMomQ_grid (e : ℕ) (Lo Hi : ℕ) (p : ℕ) (h : Lo ≤ Hi) :
    cellMomQ ((Lo : ℚ) / 2 ^ e) ((Hi : ℚ) / 2 ^ e) p
      = if p % 2 = 0 then 2 * (((Hi ^ (p + 1) - Lo ^ (p + 1) : ℕ) : ℚ)) / (((p : ℚ) + 1) * 2 ^ (e * (p + 1)))
        else 0 := by
  unfold cellMomQ
  have hp : Lo ^ (p + 1) ≤ Hi ^ (p + 1) := Nat.pow_le_pow_left h _
  rw [Nat.cast_sub hp]
  push_cast
  split_ifs with hpar
  · have hodd : Odd (p + 1) := by
      rw [Nat.odd_add_one, Nat.not_odd_iff_even]
      exact Nat.even_iff.2 hpar
    rw [hodd.neg_pow, hodd.neg_pow, div_pow, div_pow, pow_mul]
    field_simp
    ring
  · have hev : Even (p + 1) := by
      rw [Nat.even_add_one, Nat.not_even_iff_odd]
      exact Nat.odd_iff.2 (by omega)
    rw [hev.neg_pow, hev.neg_pow]
    ring

theorem momQ_grid (e f p : ℕ) : ∀ cs : List (ℕ × ℕ × ℕ), (∀ c ∈ cs, c.1 ≤ c.2.1) →
    momQ (cs.map (toQc e f)) p = momF cs e f p
  | [], _ => by simp [momQ, momF, momN]
  | c :: cs, h => by
    simp only [List.mem_cons, forall_eq_or_imp] at h
    simp only [List.map_cons, momQ]
    rw [momQ_grid e f p cs h.2]
    have hc := cellMomQ_grid e c.1 c.2.1 p h.1
    have key : ∀ x y : ℚ, (if x = 0 then 0 else x * y) = x * y := by
      intro x y
      split_ifs with hx
      · rw [hx, zero_mul]
      · rfl
    rw [key]
    simp only [toQc]
    rw [hc]
    unfold momF
    simp only [momN]
    split_ifs with hp
    · push_cast
      field_simp
    · simp

def momCheckN (cs : List (ℕ × ℕ × ℕ)) (e f K : ℕ) (momD : List ℚ) : Bool :=
  (List.range (2 * K + 1)).all fun p => decide (momD.getD p 0 = momF cs e f p)

theorem momCheckN_sound {cs : List (ℕ × ℕ × ℕ)} {e f K : ℕ} {momD : List ℚ}
    (hle : ∀ c ∈ cs, c.1 ≤ c.2.1) (h : momCheckN cs e f K momD = true) :
    momCheck (cs.map (toQc e f)) K momD = true := by
  unfold momCheckN at h
  unfold momCheck
  rw [List.all_eq_true] at h ⊢
  intro p hp
  have := h p hp
  simp only [decide_eq_true_eq] at this ⊢
  rw [this, momQ_grid e f p cs hle]

def leQ (cs : List (ℕ × ℕ × ℕ)) : Bool := cs.all fun c => decide (c.1 ≤ c.2.1)

theorem leQ_sound {cs : List (ℕ × ℕ × ℕ)} (h : leQ cs = true) : ∀ c ∈ cs, c.1 ≤ c.2.1 := by
  intro c hc
  have := List.all_eq_true.1 h c hc
  simpa using this

lemma all_range₂ {K : ℕ} {q : ℕ → ℕ → Prop} [∀ i j, Decidable (q i j)]
    (h : ((List.range K).all fun i => (List.range K).all fun j => decide (q i j)) = true) :
    ∀ i < K, ∀ j < K, q i j := by
  intro i hi j hj
  have := List.all_eq_true.1 (List.all_eq_true.1 h i (List.mem_range.2 hi)) j (List.mem_range.2 hj)
  simpa using this

/-- **Rational certificate ⇒ Weil positivity on the window `2a`.** -/
theorem weilPos_of_qcert (P : Params) (hP : Valid P) (hlog : 2 * (P.a : ℝ) < Real.log 2)
    (cs : List (ℚ × ℚ × ℚ)) (R : ℚ) (K : ℕ) (momD : List ℚ)
    (Lm Gm LGm Mm U MU Dm Wm : List (List ℚ)) (α ε β : ℚ)
    (hcells : cs.all (cellQ P) = true) (hchain : chainQ 0 cs R = true) (htail : tailQ P R = true)
    (hK : kQ P.a K cs = true) (hmom : momCheck cs K momD = true)
    (hGm : glQ K P.a Lm Gm = true) (hLGm : tmulQ K Lm Gm LGm = true)
    (hM : mQ K α Lm LGm Mm momD = true) (hMU : mulQ K Mm U MU = true)
    (hD : tmulQ K U MU Dm = true) (hUW : invQ K U Wm = true) (hdom : domQ K Dm = true)
    (hnum : decide (0 < ε ∧ 0 ≤ α ∧ 0 ≤ P.s
      ∧ 8 * P.a ^ (2 * K + 1) * momD.getD (2 * K) 0 / ((K.factorial : ℚ)) ^ 2 ≤ β
      ∧ (1 + ε) * α + (1 + 1 / ε) * β ≤ 2 * piLo * P.s) = true) :
    Riemann.Deep.Weil.WeilPos (2 * (P.a : ℝ)) := by
  simp only [decide_eq_true_eq] at hnum
  obtain ⟨hε, hα, hs, hβ, hn⟩ := hnum
  have hmom' : ∀ p ≤ 2 * K, ((momD.getD p 0 : ℚ) : ℝ) = momSum (toR cs) p := by
    intro p hp
    have := List.all_eq_true.1 hmom p (List.mem_range.2 (by omega))
    simp only [decide_eq_true_eq] at this
    rw [this, cast_momQ]
  have hGm' := all_range₂ hGm
  have hM' := all_range₂ hM
  have hpsd := hpsd_of K (α : ℝ) (entR Lm) (entR Gm) (entR LGm) (entR Mm) (Ghat (P.a : ℝ))
    (fun p => momSum (toR cs) p)
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
  have hβ' : 8 * (P.a : ℝ) ^ (2 * K + 1) * momSum (toR cs) (2 * K) / (K.factorial : ℝ) ^ 2 ≤ β := by
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
  have hch := chainQ_sound cs 0 R hchain
  rw [Rat.cast_zero] at hch
  exact weilPos_of_moment hP.ha hlog (toR cs) hch (cells_sound hP hcells)
    (tailQ_sound hP htail) K (kQ_sound hK) (entR Lm) hpsd (by exact_mod_cast hε)
    (by exact_mod_cast hα) hβ' hn'

end Fayed.Yoshida.Reflect
