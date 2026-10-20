/-
Riemann/Deep/Yoshida/MuTail.lean — sharper lower bounds for `2πμ` and `γ`.

* `two_pi_mu_ge_tel`: the digamma partial sum plus a telescoping lower bound for the tail
  (`f_n ≥ h₁(n) − h₁(n+1) + [t² ≥ 1/2](h₂(n) − h₂(n+1))`, `h₁ = (1/4)/(n+5/4)`,
  `h₂ = (t²/2)/((n+3/4)² + t²)`, `t = r/2`), cut at any `Mp ≥ M`.
* `euler_lt_of_harmonic`: `γ < H_n − log n` with `H_n` bounded by a rounded sum.
-/
import Riemann.Deep.Yoshida.MuBounds

noncomputable section

open Finset

namespace Fayed.Yoshida

open Zeta23

/-- Digamma series term at `t = r/2`. -/
def fT (r : ℝ) (n : ℕ) : ℝ :=
  1 / ((n : ℝ) + 1) - ((n : ℝ) + 1 + 1 / 4) / (((n : ℝ) + 1 + 1 / 4) ^ 2 + (r / 2) ^ 2)

def h1R (n : ℕ) : ℝ := (1 / 4) / ((n : ℝ) + 5 / 4)

def h2R (t2 : ℝ) (n : ℕ) : ℝ := (t2 / 2) / (((n : ℝ) + 3 / 4) ^ 2 + t2)

lemma fT_ge (r : ℝ) (n : ℕ) :
    h1R n - h1R (n + 1) + (if 1 / 2 ≤ (r / 2) ^ 2 then h2R ((r / 2) ^ 2) n - h2R ((r / 2) ^ 2) (n + 1)
      else 0) ≤ fT r n := by
  unfold fT h1R h2R
  set t2 := (r / 2) ^ 2 with ht2
  have ht0 : 0 ≤ t2 := by positivity
  have hn : (0 : ℝ) ≤ n := Nat.cast_nonneg n
  push_cast
  set x : ℝ := (n : ℝ) + 5 / 4 with hx
  have hx0 : 5 / 4 ≤ x := by linarith
  have e1 : (n : ℝ) + 1 + 1 / 4 = x := by rw [hx]; ring
  have e2 : (n : ℝ) + 1 + 5 / 4 = x + 1 := by rw [hx]; ring
  have e3 : (n : ℝ) + 3 / 4 = x - 1 / 2 := by rw [hx]; ring
  have e4 : (n : ℝ) + 1 + 3 / 4 = x + 1 / 2 := by rw [hx]; ring
  have e5 : (n : ℝ) + 1 = x - 1 / 4 := by rw [hx]; ring
  rw [e1, e2, e3, e4, e5]
  -- first part: 1/(x - 1/4) - 1/x ≥ (1/4)/x - (1/4)/(x+1)
  have k1 : (1 / 4) / x - (1 / 4) / (x + 1) ≤ 1 / (x - 1 / 4) - 1 / x := by
    rw [div_sub_div _ _ (by positivity) (by positivity), div_sub_div _ _ (by linarith) (by positivity)]
    rw [div_le_div_iff₀ (by positivity) (by nlinarith)]
    nlinarith
  have hsplit : 1 / (x - 1 / 4) - x / (x ^ 2 + t2) = (1 / (x - 1 / 4) - 1 / x) + (1 / x - x / (x ^ 2 + t2)) := by
    ring
  rw [hsplit]
  have k0 : 0 ≤ 1 / x - x / (x ^ 2 + t2) := by
    rw [div_sub_div _ _ (by positivity) (by positivity)]
    apply div_nonneg _ (by positivity)
    nlinarith
  split_ifs with h
  · have k2 : t2 / 2 / ((x - 1 / 2) ^ 2 + t2) - t2 / 2 / ((x + 1 / 2) ^ 2 + t2)
        ≤ 1 / x - x / (x ^ 2 + t2) := by
      have hA : 0 < (x - 1 / 2) ^ 2 + t2 := by nlinarith
      have hB : 0 < (x + 1 / 2) ^ 2 + t2 := by nlinarith
      have hC : 0 < x ^ 2 + t2 := by nlinarith
      have l1 : t2 / 2 / ((x - 1 / 2) ^ 2 + t2) - t2 / 2 / ((x + 1 / 2) ^ 2 + t2)
          = t2 * x / (((x - 1 / 2) ^ 2 + t2) * ((x + 1 / 2) ^ 2 + t2)) := by
        field_simp
        ring
      have l2 : 1 / x - x / (x ^ 2 + t2) = t2 / (x * (x ^ 2 + t2)) := by
        field_simp
        ring
      rw [l1, l2, div_le_div_iff₀ (by positivity) (by positivity)]
      have key : x ^ 2 * (x ^ 2 + t2) ≤ ((x - 1 / 2) ^ 2 + t2) * ((x + 1 / 2) ^ 2 + t2) := by
        nlinarith [sq_nonneg x, mul_nonneg (sq_nonneg x) (by linarith : (0 : ℝ) ≤ t2 - 1 / 2)]
      have := mul_le_mul_of_nonneg_left key (by positivity : (0 : ℝ) ≤ t2)
      nlinarith
    linarith
  · linarith

lemma sum_Ico_fT_ge (r : ℝ) (M : ℕ) : ∀ Mp, M ≤ Mp →
    h1R M - h1R Mp + (if 1 / 2 ≤ (r / 2) ^ 2 then h2R ((r / 2) ^ 2) M - h2R ((r / 2) ^ 2) Mp else 0)
      ≤ ∑ n ∈ Ico M Mp, fT r n := by
  intro Mp hMp
  induction Mp, hMp using Nat.le_induction with
  | base => simp
  | succ Mp hMp ih =>
    rw [Finset.sum_Ico_succ_top hMp]
    have := fT_ge r Mp
    split_ifs at this ih ⊢ <;> linarith

/-- **Digamma partial sum with a telescoping tail.** -/
theorem two_pi_mu_ge_tel (r : ℝ) (M Mp : ℕ) (hM : M ≤ Mp) {g l : ℝ}
    (hg : Real.eulerMascheroniConstant ≤ g) (hl : Real.log Real.pi ≤ l) :
    -g - l - (1 / 4) / ((1 / 4) ^ 2 + (r / 2) ^ 2) + ∑ n ∈ range M, fT r n
      + (h1R M - h1R Mp + (if 1 / 2 ≤ (r / 2) ^ 2 then h2R ((r / 2) ^ 2) M - h2R ((r / 2) ^ 2) Mp
        else 0)) ≤ 2 * Real.pi * mu r := by
  have hπ := Real.pi_pos
  rw [MuFields.mu_eq, MuFields.re_digamma_vertical (a := 1 / 4) (by norm_num) (by norm_num)]
  have hnn : ∀ n : ℕ, 0 ≤ fT r n := by
    intro n
    unfold fT
    have hα : (0 : ℝ) < (n : ℝ) + 1 + 1 / 4 := by positivity
    have h1 : ((n : ℝ) + 1 + 1 / 4) / (((n : ℝ) + 1 + 1 / 4) ^ 2 + (r / 2) ^ 2)
        ≤ 1 / ((n : ℝ) + 1 + 1 / 4) := by
      rw [div_le_div_iff₀ (by positivity) hα]
      nlinarith [sq_nonneg (r / 2)]
    have h2 : 1 / ((n : ℝ) + 1 + 1 / 4) ≤ 1 / ((n : ℝ) + 1) :=
      one_div_le_one_div_of_le (by positivity) (by linarith)
    linarith
  have hsumm := MuFields.summable_re_terms (a := 1 / 4) (by norm_num) (by norm_num) (r / 2)
  have hsum := hsumm.sum_le_tsum (range Mp) (fun n _ => hnn n)
  have hsplit := Finset.sum_range_add_sum_Ico (fun n => fT r n) hM
  have htel := sum_Ico_fT_ge r M Mp hM
  have e : 2 * Real.pi * (1 / (2 * Real.pi) * (-Real.eulerMascheroniConstant
      - 1 / 4 / ((1 / 4) ^ 2 + (r / 2) ^ 2) + ∑' n : ℕ, (1 / ((n : ℝ) + 1)
        - ((n : ℝ) + 1 + 1 / 4) / (((n : ℝ) + 1 + 1 / 4) ^ 2 + (r / 2) ^ 2)))
      - Real.log Real.pi / (2 * Real.pi))
      = -Real.eulerMascheroniConstant - 1 / 4 / ((1 / 4) ^ 2 + (r / 2) ^ 2)
        + ∑' n : ℕ, (1 / ((n : ℝ) + 1) - ((n : ℝ) + 1 + 1 / 4) / (((n : ℝ) + 1 + 1 / 4) ^ 2
          + (r / 2) ^ 2)) - Real.log Real.pi := by
    field_simp
  rw [e]
  have hs' : ∑ n ∈ range Mp, fT r n = ∑ n ∈ range Mp, (1 / ((n : ℝ) + 1)
      - ((n : ℝ) + 1 + 1 / 4) / (((n : ℝ) + 1 + 1 / 4) ^ 2 + (r / 2) ^ 2)) := rfl
  have ht' : ∑' n : ℕ, fT r n = ∑' n : ℕ, (1 / ((n : ℝ) + 1)
      - ((n : ℝ) + 1 + 1 / 4) / (((n : ℝ) + 1 + 1 / 4) ^ 2 + (r / 2) ^ 2)) := rfl
  linarith

/-! ### A sharper bound for `γ` -/

/-- `Σ_{i<n} ⌈D/(i+1)⌉`. -/
def hsumN (D : ℕ) : ℕ → ℕ
  | 0 => 0
  | n + 1 => hsumN D n + (D + n) / (n + 1)

lemma harmonic_le_hsumN (D : ℕ) (hD : 0 < D) : ∀ n : ℕ, (harmonic n : ℝ) ≤ (hsumN D n : ℝ) / D
  | 0 => by simp [hsumN]
  | n + 1 => by
    rw [harmonic_succ, hsumN]
    push_cast
    have ih := harmonic_le_hsumN D hD n
    have hq : D ≤ (D + n) / (n + 1) * (n + 1) := by
      have h1 := Nat.div_add_mod (D + n) (n + 1)
      have h2 := Nat.mod_lt (D + n) (Nat.succ_pos n)
      nlinarith
    have hD' : (0 : ℝ) < D := by exact_mod_cast hD
    have hq' : (D : ℝ) ≤ (((D + n) / (n + 1) : ℕ) : ℝ) * ((n : ℝ) + 1) := by exact_mod_cast hq
    have k : ((n : ℝ) + 1)⁻¹ ≤ (((D + n) / (n + 1) : ℕ) : ℝ) / D := by
      rw [inv_eq_one_div, div_le_div_iff₀ (by positivity) hD']
      linarith
    have e : ((hsumN D n : ℝ) + (((D + n) / (n + 1) : ℕ) : ℝ)) / D
        = (hsumN D n : ℝ) / D + (((D + n) / (n + 1) : ℕ) : ℝ) / D := by ring
    rw [e]
    push_cast at k
    linarith

set_option maxRecDepth 100000 in
theorem hsumN_val : hsumN 1099511627776 16384 = 11304416284840 := by decide +kernel

/-- `γ < 0.57725` (from `γ < H_{2^14} − 14 log 2`). -/
theorem euler_lt' : Real.eulerMascheroniConstant < 0.57725 := by
  have h := Real.eulerMascheroniConstant_lt_eulerMascheroniSeq' 16384
  have e : Real.eulerMascheroniSeq' 16384 = (harmonic 16384 : ℝ) - 14 * Real.log 2 := by
    rw [Real.eulerMascheroniSeq', if_neg (by norm_num)]
    rw [show ((16384 : ℕ) : ℝ) = 2 ^ 14 by norm_num, Real.log_pow]
    push_cast
    ring
  rw [e] at h
  have hH := harmonic_le_hsumN 1099511627776 (by norm_num) 16384
  rw [hsumN_val] at hH
  generalize (harmonic 16384 : ℝ) = H at h hH
  have hl := Real.log_two_gt_d9
  norm_num at hH hl
  linarith

end Fayed.Yoshida
