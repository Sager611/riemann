import Riemann.Deep.FirstZero.NextLib

/-!
# Minimality of the first zero, part 1: certified Taylor models of ζ along segments

`Mf N s` is the order-8 Euler–Maclaurin main term (as in `em_main_eqN`); `M1f`, `M2f` are its first
two derivatives, written explicitly. On a segment s = p + u·d (u ∈ [0,h], ‖d‖ = 1):

  ‖ζ(p + u d) − (A + u d B)‖ ≤ ε + r_A + r_B u + C u²            (`seg_model`)

from a ball A for Mf(p), a ball B for M1f(p), a crude bound C for ‖M2f‖ on the segment
(`norm_M2f_le`) and the Euler–Maclaurin remainder ε (`em_err_le`). `re_pos_of_model` turns a
model plus an anchor c into Re(ζ·c̄) > 0 on the whole segment; `arg_le_*`/`le_arg_*` bound
arguments from balls.
-/

open Complex Set

namespace Fayed.FirstZero.Min
open Fayed.FirstZero

/-! ## Powers -/

lemma cpow_gen {x : ℝ} (hx : 0 < x) (σ t : ℝ) :
    (x : ℂ) ^ (-((σ : ℂ) + (t : ℂ) * I)) =
      ((x ^ (-σ) : ℝ) : ℂ) * cexp (((-(t * Real.log x)) : ℝ) * I) := by
  have hx0 : (x : ℂ) ≠ 0 := by exact_mod_cast hx.ne'
  rw [cpow_def_of_ne_zero hx0, ← ofReal_log hx.le, Real.rpow_def_of_pos hx, ofReal_exp,
    ← Complex.exp_add]
  congr 1
  push_cast
  ring

lemma rpow_ball {n : ℝ} (hn : 0 < n) (p q : ℕ) (hq : 0 < q) {c ε : ℝ} (hε : 0 ≤ ε) (hc : 0 < c - ε)
    (h1 : (c - ε) ^ q * n ^ p ≤ 1) (h2 : 1 ≤ (c + ε) ^ q * n ^ p) :
    |n ^ (-((p : ℝ) / q)) - c| ≤ ε := by
  set x := n ^ (-((p : ℝ) / q)) with hx
  have hx0 : 0 < x := Real.rpow_pos_of_pos hn _
  have hq' : (q : ℝ) ≠ 0 := by positivity
  have hnp : 0 < n ^ p := pow_pos hn p
  have hxx : x ^ q * n ^ p = 1 := by
    rw [hx, ← Real.rpow_natCast, ← Real.rpow_mul hn.le,
      show -((p : ℝ) / q) * (q : ℝ) = -(p : ℝ) by field_simp,
      Real.rpow_neg hn.le, Real.rpow_natCast, inv_mul_cancel₀ hnp.ne']
  have l1 : c - ε ≤ x := by
    by_contra hcon
    push Not at hcon
    have : x ^ q < (c - ε) ^ q := pow_lt_pow_left₀ hcon hx0.le (by omega)
    nlinarith [mul_lt_mul_of_pos_right this hnp]
  have l2 : x ≤ c + ε := by
    by_contra hcon
    push Not at hcon
    have : (c + ε) ^ q < x ^ q := pow_lt_pow_left₀ hcon (by linarith) (by omega)
    nlinarith [mul_lt_mul_of_pos_right this hnp]
  rw [abs_le]; constructor <;> linarith

/-! ## Derivatives of the Pochhammer factors `gc` -/

/-- d/ds gc s j. -/
noncomputable def gc1 (s : ℂ) : ℕ → ℂ
  | 0 => 0
  | j + 1 => gc1 s j * (-s - 1 - j) - gc s j

/-- d²/ds² gc s j. -/
noncomputable def gc2 (s : ℂ) : ℕ → ℂ
  | 0 => 0
  | j + 1 => gc2 s j * (-s - 1 - j) - 2 * gc1 s j

lemma hasDerivAt_lin (s : ℂ) (j : ℕ) : HasDerivAt (fun z : ℂ => -z - 1 - (j : ℂ)) (-1) s :=
  ((hasDerivAt_id' s).neg.sub_const 1).sub_const (j : ℂ)

lemma hasDerivAt_gc (s : ℂ) (j : ℕ) : HasDerivAt (fun z => gc z j) (gc1 s j) s := by
  induction j with
  | zero => simp only [gc, gc1]; exact hasDerivAt_const s 1
  | succ j ih =>
    simp only [gc, gc1]
    convert ih.mul (hasDerivAt_lin s j) using 1
    ring

lemma hasDerivAt_gc1 (s : ℂ) (j : ℕ) : HasDerivAt (fun z => gc1 z j) (gc2 s j) s := by
  induction j with
  | zero => simp only [gc1, gc2]; exact hasDerivAt_const s 0
  | succ j ih =>
    simp only [gc1, gc2]
    convert (ih.mul (hasDerivAt_lin s j)).sub (hasDerivAt_gc s j) using 1
    ring

/-- Crude bounds for gc, gc1, gc2 on ‖s‖ ≤ R. -/
noncomputable def G0 (R : ℝ) : ℕ → ℝ
  | 0 => 1
  | j + 1 => G0 R j * (R + 1 + j)

noncomputable def G1 (R : ℝ) : ℕ → ℝ
  | 0 => 0
  | j + 1 => G1 R j * (R + 1 + j) + G0 R j

noncomputable def G2 (R : ℝ) : ℕ → ℝ
  | 0 => 0
  | j + 1 => G2 R j * (R + 1 + j) + 2 * G1 R j

lemma gc_bounds {s : ℂ} {R : ℝ} (hR : ‖s‖ ≤ R) (j : ℕ) :
    ‖gc s j‖ ≤ G0 R j ∧ ‖gc1 s j‖ ≤ G1 R j ∧ ‖gc2 s j‖ ≤ G2 R j := by
  induction j with
  | zero => simp [gc, gc1, gc2, G0, G1, G2]
  | succ j ih =>
    obtain ⟨h0, h1, h2⟩ := ih
    have hf : ‖-s - 1 - (j : ℂ)‖ ≤ R + 1 + j := by
      rw [show -s - 1 - (j : ℂ) = -(s + ((1 + j : ℕ) : ℂ)) by push_cast; ring, norm_neg]
      refine (norm_add_le _ _).trans ?_
      rw [Complex.norm_natCast]
      push_cast
      try linarith
    have hG0 : 0 ≤ G0 R j := (norm_nonneg _).trans h0
    have hG1 : 0 ≤ G1 R j := (norm_nonneg _).trans h1
    refine ⟨?_, ?_, ?_⟩
    · simp only [gc, G0, norm_mul]
      exact mul_le_mul h0 hf (norm_nonneg _) hG0
    · simp only [gc1, G1]
      calc ‖gc1 s j * (-s - 1 - j) - gc s j‖ ≤ ‖gc1 s j‖ * ‖-s - 1 - (j : ℂ)‖ + ‖gc s j‖ := by
            refine (norm_sub_le _ _).trans ?_; rw [norm_mul]
        _ ≤ G1 R j * (R + 1 + j) + G0 R j := by
            linarith [mul_le_mul h1 hf (norm_nonneg _) hG1]
    · simp only [gc2, G2]
      calc ‖gc2 s j * (-s - 1 - j) - 2 * gc1 s j‖ ≤ ‖gc2 s j‖ * ‖-s - 1 - (j : ℂ)‖ + 2 * ‖gc1 s j‖ := by
            refine (norm_sub_le _ _).trans ?_; rw [norm_mul, norm_mul]; norm_num
        _ ≤ G2 R j * (R + 1 + j) + 2 * G1 R j := by
            linarith [mul_le_mul h2 hf (norm_nonneg _) ((norm_nonneg _).trans h2)]

/-! ## The Bernoulli part of `EMKN` and its derivatives -/

noncomputable def Bs (N : ℕ) (s : ℂ) : ℂ :=
  1 / 12 * gc s 0 / (N : ℂ) - 1 / 720 * gc s 2 / (N : ℂ) ^ 3 + 1 / 30240 * gc s 4 / (N : ℂ) ^ 5
    - 1 / 1209600 * gc s 6 / (N : ℂ) ^ 7

noncomputable def Bs1 (N : ℕ) (s : ℂ) : ℂ :=
  1 / 12 * gc1 s 0 / (N : ℂ) - 1 / 720 * gc1 s 2 / (N : ℂ) ^ 3 + 1 / 30240 * gc1 s 4 / (N : ℂ) ^ 5
    - 1 / 1209600 * gc1 s 6 / (N : ℂ) ^ 7

noncomputable def Bs2 (N : ℕ) (s : ℂ) : ℂ :=
  1 / 12 * gc2 s 0 / (N : ℂ) - 1 / 720 * gc2 s 2 / (N : ℂ) ^ 3 + 1 / 30240 * gc2 s 4 / (N : ℂ) ^ 5
    - 1 / 1209600 * gc2 s 6 / (N : ℂ) ^ 7

lemma EMKN_eq (N : ℕ) (s : ℂ) : EMKN N s = 1 - (N : ℂ) / (1 - s) - 1 / 2 + s * Bs N s := rfl

lemma hasDerivAt_Bs (N : ℕ) (s : ℂ) : HasDerivAt (Bs N) (Bs1 N s) s := by
  have h := (((((hasDerivAt_gc s 0).const_mul (1 / 12 : ℂ)).div_const (N : ℂ)).sub
    (((hasDerivAt_gc s 2).const_mul (1 / 720 : ℂ)).div_const ((N : ℂ) ^ 3))).add
    (((hasDerivAt_gc s 4).const_mul (1 / 30240 : ℂ)).div_const ((N : ℂ) ^ 5))).sub
    (((hasDerivAt_gc s 6).const_mul (1 / 1209600 : ℂ)).div_const ((N : ℂ) ^ 7))
  exact h

lemma hasDerivAt_Bs1 (N : ℕ) (s : ℂ) : HasDerivAt (Bs1 N) (Bs2 N s) s := by
  have h := (((((hasDerivAt_gc1 s 0).const_mul (1 / 12 : ℂ)).div_const (N : ℂ)).sub
    (((hasDerivAt_gc1 s 2).const_mul (1 / 720 : ℂ)).div_const ((N : ℂ) ^ 3))).add
    (((hasDerivAt_gc1 s 4).const_mul (1 / 30240 : ℂ)).div_const ((N : ℂ) ^ 5))).sub
    (((hasDerivAt_gc1 s 6).const_mul (1 / 1209600 : ℂ)).div_const ((N : ℂ) ^ 7))
  exact h

noncomputable def EMKN1 (N : ℕ) (s : ℂ) : ℂ := -(N : ℂ) / (1 - s) ^ 2 + (Bs N s + s * Bs1 N s)

noncomputable def EMKN2 (N : ℕ) (s : ℂ) : ℂ :=
  -(2 * (N : ℂ)) / (1 - s) ^ 3 + (2 * Bs1 N s + s * Bs2 N s)

lemma hasDerivAt_EMKN (N : ℕ) {s : ℂ} (hs : s ≠ 1) : HasDerivAt (EMKN N) (EMKN1 N s) s := by
  have h1s : (1 : ℂ) - s ≠ 0 := sub_ne_zero.mpr (Ne.symm hs)
  have hq : HasDerivAt (fun z : ℂ => (N : ℂ) / (1 - z)) ((N : ℂ) / (1 - s) ^ 2) s := by
    have := (hasDerivAt_const s (N : ℂ)).div ((hasDerivAt_id' s).const_sub 1) h1s
    convert this using 1
    field_simp
    ring
  have hp : HasDerivAt (fun z : ℂ => z * Bs N z) (Bs N s + s * Bs1 N s) s := by
    convert (hasDerivAt_id' s).mul (hasDerivAt_Bs N s) using 1
    ring
  have h := (((hasDerivAt_const s (1 : ℂ)).sub hq).sub_const (1 / 2 : ℂ)).add hp
  convert h using 1
  · funext z; rfl
  · simp only [EMKN1]; ring

lemma hasDerivAt_EMKN1 (N : ℕ) {s : ℂ} (hs : s ≠ 1) : HasDerivAt (EMKN1 N) (EMKN2 N s) s := by
  have h1s : (1 : ℂ) - s ≠ 0 := sub_ne_zero.mpr (Ne.symm hs)
  have hq : HasDerivAt (fun z : ℂ => -(N : ℂ) / (1 - z) ^ 2) (-(2 * (N : ℂ)) / (1 - s) ^ 3) s := by
    have h2 : HasDerivAt (fun z : ℂ => (1 - z) ^ 2) (2 * (1 - s) * (-1)) s := by
      convert ((hasDerivAt_id' s).const_sub 1).pow 2 using 1
      simp
    have := (hasDerivAt_const s (-(N : ℂ))).div h2 (pow_ne_zero 2 h1s)
    convert this using 1
    field_simp
    ring
  have hp : HasDerivAt (fun z : ℂ => z * Bs1 N z) (Bs1 N s + s * Bs2 N s) s := by
    convert (hasDerivAt_id' s).mul (hasDerivAt_Bs1 N s) using 1
    ring
  have h := hq.add ((hasDerivAt_Bs N s).add hp)
  exact h.congr_deriv (by simp only [EMKN2]; ring)

/-! ## The Euler–Maclaurin main term and its derivatives -/

noncomputable def lg (n : ℕ) : ℂ := (Real.log n : ℂ)

noncomputable def Mf (N : ℕ) (s : ℂ) : ℂ :=
  1 + (∑ k ∈ Finset.range (N - 2), nps s (k + 2)) + nps s N * EMKN N s

noncomputable def M1f (N : ℕ) (s : ℂ) : ℂ :=
  (∑ k ∈ Finset.range (N - 2), -lg (k + 2) * nps s (k + 2))
    + nps s N * (EMKN1 N s - lg N * EMKN N s)

noncomputable def M2f (N : ℕ) (s : ℂ) : ℂ :=
  (∑ k ∈ Finset.range (N - 2), -lg (k + 2) * (-lg (k + 2) * nps s (k + 2)))
    + nps s N * (EMKN2 N s - 2 * lg N * EMKN1 N s + lg N ^ 2 * EMKN N s)

lemma hasDerivAt_nps {n : ℕ} (hn : n ≠ 0) (s : ℂ) :
    HasDerivAt (fun z => nps z n) (-lg n * nps s n) s := by
  have hneg : HasDerivAt (fun z : ℂ => -z) (-1) s := hasDerivAt_neg' s
  have h : HasDerivAt (fun z : ℂ => (n : ℂ) ^ (-z)) ((n : ℂ) ^ (-s) * Complex.log n * (-1)) s :=
    hneg.const_cpow (Or.inl (Nat.cast_ne_zero.mpr hn))
  have hlog : Complex.log (n : ℂ) = lg n := by
    rw [lg, Complex.ofReal_log (Nat.cast_nonneg n), Complex.ofReal_natCast]
  exact h.congr_deriv (by rw [hlog, nps]; ring)

lemma hasDerivAt_Mf {N : ℕ} (hN : 2 ≤ N) {s : ℂ} (hs : s ≠ 1) :
    HasDerivAt (Mf N) (M1f N s) s := by
  have hsum : HasDerivAt (fun z => ∑ k ∈ Finset.range (N - 2), nps z (k + 2))
      (∑ k ∈ Finset.range (N - 2), -lg (k + 2) * nps s (k + 2)) s :=
    HasDerivAt.fun_sum fun k _ => hasDerivAt_nps (by omega) s
  have h := ((hasDerivAt_const s (1 : ℂ)).add hsum).add
    ((hasDerivAt_nps (n := N) (by omega) s).mul (hasDerivAt_EMKN N hs))
  exact h.congr_deriv (by simp only [M1f]; ring)

lemma hasDerivAt_M1f {N : ℕ} (hN : 2 ≤ N) {s : ℂ} (hs : s ≠ 1) :
    HasDerivAt (M1f N) (M2f N s) s := by
  have hsum : HasDerivAt (fun z => ∑ k ∈ Finset.range (N - 2), -lg (k + 2) * nps z (k + 2))
      (∑ k ∈ Finset.range (N - 2), -lg (k + 2) * (-lg (k + 2) * nps s (k + 2))) s :=
    HasDerivAt.fun_sum fun k _ => (hasDerivAt_nps (by omega) s).const_mul _
  have hE : HasDerivAt (fun z => EMKN1 N z - lg N * EMKN N z)
      (EMKN2 N s - lg N * EMKN1 N s) s :=
    (hasDerivAt_EMKN1 N hs).sub ((hasDerivAt_EMKN N hs).const_mul _)
  have h := hsum.add ((hasDerivAt_nps (n := N) (by omega) s).mul hE)
  exact h.congr_deriv (by simp only [M2f]; ring)

/-! ## Crude bound for M2f -/

noncomputable def B0b (N : ℕ) (R : ℝ) : ℝ :=
  1 / 12 * G0 R 0 / N + 1 / 720 * G0 R 2 / (N : ℝ) ^ 3 + 1 / 30240 * G0 R 4 / (N : ℝ) ^ 5
    + 1 / 1209600 * G0 R 6 / (N : ℝ) ^ 7

noncomputable def B1b (N : ℕ) (R : ℝ) : ℝ :=
  1 / 12 * G1 R 0 / N + 1 / 720 * G1 R 2 / (N : ℝ) ^ 3 + 1 / 30240 * G1 R 4 / (N : ℝ) ^ 5
    + 1 / 1209600 * G1 R 6 / (N : ℝ) ^ 7

noncomputable def B2b (N : ℕ) (R : ℝ) : ℝ :=
  1 / 12 * G2 R 0 / N + 1 / 720 * G2 R 2 / (N : ℝ) ^ 3 + 1 / 30240 * G2 R 4 / (N : ℝ) ^ 5
    + 1 / 1209600 * G2 R 6 / (N : ℝ) ^ 7

lemma norm_four_le {g : ℕ → ℂ} {G : ℕ → ℝ} (hg : ∀ j, ‖g j‖ ≤ G j) (N : ℕ) (hN : 0 < N) :
    ‖1 / 12 * g 0 / (N : ℂ) - 1 / 720 * g 2 / (N : ℂ) ^ 3 + 1 / 30240 * g 4 / (N : ℂ) ^ 5
      - 1 / 1209600 * g 6 / (N : ℂ) ^ 7‖ ≤
    1 / 12 * G 0 / N + 1 / 720 * G 2 / (N : ℝ) ^ 3 + 1 / 30240 * G 4 / (N : ℝ) ^ 5
      + 1 / 1209600 * G 6 / (N : ℝ) ^ 7 := by
  have hN' : (0 : ℝ) < N := by exact_mod_cast hN
  have e : ∀ (c : ℝ) (j k : ℕ), 0 ≤ c → ‖(c : ℂ) * g j / (N : ℂ) ^ k‖ ≤ c * G j / (N : ℝ) ^ k := by
    intro c j k hc
    rw [norm_div, norm_mul, norm_pow, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hc,
      Complex.norm_natCast]
    gcongr
    exact hg j
  have e1 := e (1 / 12) 0 1 (by norm_num)
  have e2 := e (1 / 720) 2 3 (by norm_num)
  have e3 := e (1 / 30240) 4 5 (by norm_num)
  have e4 := e (1 / 1209600) 6 7 (by norm_num)
  push_cast at e1 e2 e3 e4
  rw [pow_one] at e1
  calc _ ≤ ‖1 / 12 * g 0 / (N : ℂ)‖ + ‖1 / 720 * g 2 / (N : ℂ) ^ 3‖ + ‖1 / 30240 * g 4 / (N : ℂ) ^ 5‖
        + ‖1 / 1209600 * g 6 / (N : ℂ) ^ 7‖ := by
        refine (norm_sub_le _ _).trans ?_
        gcongr
        refine (norm_add_le _ _).trans ?_
        gcongr
        exact norm_sub_le _ _
    _ ≤ _ := by
        rw [pow_one] at e1
        linarith

lemma norm_Bs_le (N : ℕ) (hN : 0 < N) {s : ℂ} {R : ℝ} (hR : ‖s‖ ≤ R) :
    ‖Bs N s‖ ≤ B0b N R ∧ ‖Bs1 N s‖ ≤ B1b N R ∧ ‖Bs2 N s‖ ≤ B2b N R :=
  ⟨norm_four_le (fun j => (gc_bounds hR j).1) N hN,
   norm_four_le (fun j => (gc_bounds hR j).2.1) N hN,
   norm_four_le (fun j => (gc_bounds hR j).2.2) N hN⟩

noncomputable def K0b (N : ℕ) (R D : ℝ) : ℝ := 1 / 2 + N / D + R * B0b N R
noncomputable def K1b (N : ℕ) (R D : ℝ) : ℝ := N / D ^ 2 + (B0b N R + R * B1b N R)
noncomputable def K2b (N : ℕ) (R D : ℝ) : ℝ := 2 * N / D ^ 3 + (2 * B1b N R + R * B2b N R)

lemma norm_EMKN_le {N : ℕ} (hN : 0 < N) {s : ℂ} {R D : ℝ} (hR : ‖s‖ ≤ R) (hD0 : 0 < D)
    (hD : D ≤ ‖1 - s‖) :
    ‖EMKN N s‖ ≤ K0b N R D ∧ ‖EMKN1 N s‖ ≤ K1b N R D ∧ ‖EMKN2 N s‖ ≤ K2b N R D := by
  obtain ⟨b0, b1, b2⟩ := norm_Bs_le N hN hR
  have hR0 : 0 ≤ R := (norm_nonneg _).trans hR
  have hN0 : (0 : ℝ) ≤ N := Nat.cast_nonneg N
  have hd : ∀ k : ℕ, (c : ℝ) → 0 ≤ c → ‖(c : ℂ) / (1 - s) ^ k‖ ≤ c / D ^ k := by
    intro k c hc
    rw [norm_div, norm_pow, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hc]
    gcongr
  have q1 := hd 1 N hN0
  have q2 := hd 2 N hN0
  have q3 := hd 3 (2 * N) (by positivity)
  push_cast at q1 q2 q3
  rw [pow_one, pow_one] at q1
  have m0 : ‖s * Bs N s‖ ≤ R * B0b N R := by
    rw [norm_mul]; exact mul_le_mul hR b0 (norm_nonneg _) hR0
  have m1 : ‖s * Bs1 N s‖ ≤ R * B1b N R := by
    rw [norm_mul]; exact mul_le_mul hR b1 (norm_nonneg _) hR0
  have m2 : ‖s * Bs2 N s‖ ≤ R * B2b N R := by
    rw [norm_mul]; exact mul_le_mul hR b2 (norm_nonneg _) hR0
  refine ⟨?_, ?_, ?_⟩
  · rw [EMKN_eq, show (1 : ℂ) - N / (1 - s) - 1 / 2 + s * Bs N s
        = (1 / 2 : ℂ) - N / (1 - s) + s * Bs N s by ring]
    calc _ ≤ ‖(1 / 2 : ℂ)‖ + ‖(N : ℂ) / (1 - s)‖ + ‖s * Bs N s‖ :=
          (norm_add_le _ _).trans (by gcongr; exact norm_sub_le _ _)
      _ ≤ 1 / 2 + N / D + R * B0b N R := by
          have : ‖(1 / 2 : ℂ)‖ = 1 / 2 := by norm_num
          linarith
  · calc ‖EMKN1 N s‖ ≤ ‖-(N : ℂ) / (1 - s) ^ 2‖ + (‖Bs N s‖ + ‖s * Bs1 N s‖) :=
          (norm_add_le _ _).trans (by gcongr; exact norm_add_le _ _)
      _ ≤ K1b N R D := by
          rw [neg_div, norm_neg]
          unfold K1b; linarith
  · calc ‖EMKN2 N s‖ ≤ ‖-(2 * (N : ℂ)) / (1 - s) ^ 3‖ + (‖2 * Bs1 N s‖ + ‖s * Bs2 N s‖) :=
          (norm_add_le _ _).trans (by gcongr; exact norm_add_le _ _)
      _ ≤ K2b N R D := by
          rw [neg_div, norm_neg, norm_mul, show ‖(2 : ℂ)‖ = 2 by norm_num]
          unfold K2b; linarith

lemma norm_nps_le {n : ℕ} (hn : 0 < n) {s : ℂ} {σ0 : ℝ} (hσ : σ0 ≤ s.re) :
    ‖nps s n‖ ≤ (n : ℝ) ^ (-σ0) := by
  rw [nps, norm_natCast_cpow_of_pos hn, neg_re]
  exact Real.rpow_le_rpow_of_exponent_le (by exact_mod_cast hn) (by linarith)

lemma norm_lg (n : ℕ) : ‖lg n‖ = Real.log n := by
  rw [lg, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg (Real.log_natCast_nonneg n)]

lemma norm_M2f_le {N : ℕ} (hN : 2 ≤ N) {s : ℂ} {R D σ0 S v L : ℝ} (hR : ‖s‖ ≤ R) (hD0 : 0 < D)
    (hD : D ≤ ‖1 - s‖) (hσ : σ0 ≤ s.re)
    (hS : ∑ k ∈ Finset.range (N - 2), Real.log ((k + 2 : ℕ) : ℝ) ^ 2 * ((k + 2 : ℕ) : ℝ) ^ (-σ0) ≤ S)
    (hv : (N : ℝ) ^ (-σ0) ≤ v) (hL : Real.log N ≤ L) :
    ‖M2f N s‖ ≤ S + v * (K2b N R D + 2 * L * K1b N R D + L ^ 2 * K0b N R D) := by
  obtain ⟨k0, k1, k2⟩ := norm_EMKN_le (by omega : 0 < N) hR hD0 hD
  have hL0 : 0 ≤ Real.log N := Real.log_natCast_nonneg N
  have hsum : ‖∑ k ∈ Finset.range (N - 2), -lg (k + 2) * (-lg (k + 2) * nps s (k + 2))‖ ≤ S := by
    refine (norm_sum_le _ _).trans (le_trans (Finset.sum_le_sum fun k _ => ?_) hS)
    rw [norm_mul, norm_neg, norm_mul, norm_neg, norm_lg, ← mul_assoc, ← sq]
    exact mul_le_mul_of_nonneg_left (norm_nps_le (by omega) hσ) (sq_nonneg _)
  have hnN : ‖nps s N‖ ≤ v := (norm_nps_le (by omega) hσ).trans hv
  have hK0 : 0 ≤ K0b N R D := (norm_nonneg _).trans k0
  have hK1 : 0 ≤ K1b N R D := (norm_nonneg _).trans k1
  have hK2 : 0 ≤ K2b N R D := (norm_nonneg _).trans k2
  have hinner : ‖EMKN2 N s - 2 * lg N * EMKN1 N s + lg N ^ 2 * EMKN N s‖ ≤
      K2b N R D + 2 * L * K1b N R D + L ^ 2 * K0b N R D := by
    calc _ ≤ ‖EMKN2 N s‖ + ‖2 * lg N * EMKN1 N s‖ + ‖lg N ^ 2 * EMKN N s‖ :=
          (norm_add_le _ _).trans (by gcongr; exact norm_sub_le _ _)
      _ = ‖EMKN2 N s‖ + 2 * Real.log N * ‖EMKN1 N s‖ + Real.log N ^ 2 * ‖EMKN N s‖ := by
          rw [norm_mul, norm_mul, norm_mul, norm_pow, norm_lg, show ‖(2 : ℂ)‖ = 2 by norm_num]
      _ ≤ _ := by
          have a1 : Real.log N * ‖EMKN1 N s‖ ≤ L * K1b N R D :=
            mul_le_mul hL k1 (norm_nonneg _) (hL0.trans hL)
          have a2 : Real.log N ^ 2 * ‖EMKN N s‖ ≤ L ^ 2 * K0b N R D :=
            mul_le_mul (pow_le_pow_left₀ hL0 hL 2) k0 (norm_nonneg _) (sq_nonneg _)
          linarith
  calc ‖M2f N s‖ ≤ S + ‖nps s N‖ * ‖EMKN2 N s - 2 * lg N * EMKN1 N s + lg N ^ 2 * EMKN N s‖ := by
        rw [M2f]; refine (norm_add_le _ _).trans ?_; rw [norm_mul]; linarith
    _ ≤ _ := by
        gcongr
        exact (norm_nonneg _).trans hnN

/-! ## Euler–Maclaurin remainder -/

lemma zeta_sub_Mf {N : ℕ} (hN : 2 ≤ N) {s : ℂ} (hs : 0 < s.re) (hs1 : s ≠ 1) :
    ‖riemannZeta s - Mf N s‖ ≤
      ‖s‖ * (‖gc s 7‖ / 1000000 * ((N : ℝ) ^ (-(s.re + 7)) / (s.re + 7))) := by
  have hsz : s ≠ 0 := by intro h; rw [h] at hs; simp at hs
  have hEM := zeta_EM8 hs hs1 (N := N) (by omega)
  rw [em_main_eqN hN s hsz hs1] at hEM
  exact hEM

lemma em_err_le {N : ℕ} (hN : 2 ≤ N) {s : ℂ} {σ0 R v ε : ℝ} (c : ℕ → ℝ) (hσ0 : 0 < σ0)
    (hσ : σ0 ≤ s.re) (hs1 : s ≠ 1) (hR : ‖s‖ ≤ R) (hc : ∀ j : ℕ, j < 7 → ‖s + 1 + (j : ℂ)‖ ≤ c j)
    (hv : (N : ℝ) ^ (-σ0) ≤ v)
    (hε : R * ((∏ j ∈ Finset.range 7, c j) / 1000000 * (v / (N : ℝ) ^ 7 / (σ0 + 7))) ≤ ε) :
    ‖riemannZeta s - Mf N s‖ ≤ ε := by
  have hs : 0 < s.re := lt_of_lt_of_le hσ0 hσ
  refine (zeta_sub_Mf hN hs hs1).trans (le_trans ?_ hε)
  have hN1 : (1 : ℝ) ≤ N := by exact_mod_cast (by omega : 1 ≤ N)
  have hN0 : (0 : ℝ) < N := by linarith
  have hpow : (N : ℝ) ^ (-(s.re + 7)) ≤ v / (N : ℝ) ^ 7 := by
    calc (N : ℝ) ^ (-(s.re + 7)) ≤ (N : ℝ) ^ (-(σ0 + 7)) :=
          Real.rpow_le_rpow_of_exponent_le hN1 (by linarith)
      _ = (N : ℝ) ^ (-σ0) / (N : ℝ) ^ 7 := by
          rw [show -(σ0 + 7) = -σ0 + -((7 : ℕ) : ℝ) by push_cast; ring, Real.rpow_add hN0,
            Real.rpow_neg hN0.le ((7 : ℕ) : ℝ), Real.rpow_natCast, div_eq_mul_inv]
      _ ≤ v / (N : ℝ) ^ 7 := by gcongr
  have hg := norm_gc_le s c 7 hc
  have hv0 : 0 ≤ v := le_trans (Real.rpow_nonneg hN0.le _) hv
  have hc0 : 0 ≤ ∏ j ∈ Finset.range 7, c j := (norm_nonneg _).trans hg
  calc ‖s‖ * (‖gc s 7‖ / 1000000 * ((N : ℝ) ^ (-(s.re + 7)) / (s.re + 7)))
      ≤ R * ((∏ j ∈ Finset.range 7, c j) / 1000000 * ((v / (N : ℝ) ^ 7) / (σ0 + 7))) := by
        have hX : ‖gc s 7‖ / 1000000 * ((N : ℝ) ^ (-(s.re + 7)) / (s.re + 7)) ≤
            (∏ j ∈ Finset.range 7, c j) / 1000000 * ((v / (N : ℝ) ^ 7) / (σ0 + 7)) :=
          mul_le_mul (div_le_div_of_nonneg_right hg (by norm_num))
            (div_le_div₀ (by positivity) hpow (by linarith) (by linarith))
            (by positivity) (by positivity)
        exact mul_le_mul hR hX (by positivity) ((norm_nonneg _).trans hR)

/-- ‖s + a‖ from bounds on Re s and |Im s|. -/
lemma norm_shift_le {s : ℂ} {x1 Y a c : ℝ} (ha : 0 ≤ a) (hre0 : 0 ≤ s.re) (hre1 : s.re ≤ x1)
    (him : |s.im| ≤ Y) (hc : 0 ≤ c) (h : (x1 + a) ^ 2 + Y ^ 2 ≤ c ^ 2) : ‖s + a‖ ≤ c := by
  refine norm_le_of_sq_le hc ?_
  have h1 : (s + a).re = s.re + a := by simp
  have h2 : (s + a).im = s.im := by simp
  rw [h1, h2]
  have : (s.re + a) ^ 2 ≤ (x1 + a) ^ 2 := by
    apply pow_le_pow_left₀ (by linarith); linarith
  have : s.im ^ 2 ≤ Y ^ 2 := by
    rw [← sq_abs]; exact pow_le_pow_left₀ (abs_nonneg _) him 2
  linarith

/-! ## Second-order Taylor bound along a segment -/

lemma hasDerivAt_seg (p d : ℂ) (u : ℝ) : HasDerivAt (fun v : ℝ => p + (v : ℂ) * d) d u := by
  have := (((hasDerivAt_id u).ofReal_comp).mul_const d).const_add p
  simpa using this

lemma taylor2 {f f1 f2 : ℂ → ℂ} {p d : ℂ} (hd : ‖d‖ = 1) {h C : ℝ}
    (hf : ∀ u ∈ Icc (0 : ℝ) h, HasDerivAt f (f1 (p + u * d)) (p + u * d))
    (hf1 : ∀ u ∈ Icc (0 : ℝ) h, HasDerivAt f1 (f2 (p + u * d)) (p + u * d))
    (hC : ∀ u ∈ Icc (0 : ℝ) h, ‖f2 (p + u * d)‖ ≤ C) :
    ∀ x ∈ Icc (0 : ℝ) h, ‖f (p + x * d) - f p - x * d * f1 p‖ ≤ C * x ^ 2 := by
  intro x hx
  have hC0 : 0 ≤ C := (norm_nonneg _).trans (hC 0 ⟨le_rfl, hx.1.trans hx.2⟩)
  have sub : ∀ u ∈ Icc (0 : ℝ) x, u ∈ Icc (0 : ℝ) h := fun u hu => ⟨hu.1, hu.2.trans hx.2⟩
  have hG : ∀ u ∈ Icc (0 : ℝ) x, HasDerivWithinAt (fun v : ℝ => f1 (p + v * d))
      (f2 (p + u * d) * d) (Icc 0 x) u := fun u hu =>
    ((hf1 u (sub u hu)).comp u (hasDerivAt_seg p d u)).hasDerivWithinAt
  have hGb0 := norm_image_sub_le_of_norm_deriv_le_segment' hG (C := C) (fun u hu => by
      rw [norm_mul, hd, mul_one]; exact hC u (sub u (Ico_subset_Icc_self hu)))
  have hGb : ∀ u ∈ Icc (0 : ℝ) x, ‖f1 (p + u * d) - f1 p‖ ≤ C * u := by
    intro u hu
    have := hGb0 u hu
    simpa using this
  have hF : ∀ u ∈ Icc (0 : ℝ) x, HasDerivWithinAt (fun v : ℝ => f (p + v * d) - v * d * f1 p)
      (f1 (p + u * d) * d - d * f1 p) (Icc 0 x) u := by
    intro u hu
    have h1 := (hf u (sub u hu)).comp u (hasDerivAt_seg p d u)
    have h2 : HasDerivAt (fun v : ℝ => (v : ℂ) * d * f1 p) (d * f1 p) u := by
      have := (((hasDerivAt_id u).ofReal_comp).mul_const d).mul_const (f1 p)
      simpa using this
    exact (h1.sub h2).hasDerivWithinAt
  have hFb := norm_image_sub_le_of_norm_deriv_le_segment' hF (C := C * x) (fun u hu => by
      rw [show f1 (p + ↑u * d) * d - d * f1 p = (f1 (p + ↑u * d) - f1 p) * d by ring, norm_mul,
        hd, mul_one]
      calc _ ≤ C * u := hGb u (Ico_subset_Icc_self hu)
        _ ≤ C * x := by gcongr; exact hu.2.le) x ⟨hx.1, le_rfl⟩
  have e : f (p + x * d) - f p - x * d * f1 p =
      (f (p + x * d) - x * d * f1 p) - (f (p + ((0 : ℝ) : ℂ) * d) - ((0 : ℝ) : ℂ) * d * f1 p) := by
    simp only [ofReal_zero, zero_mul, add_zero, sub_zero]; ring
  rw [e]
  calc _ ≤ C * x * (x - 0) := hFb
    _ = C * x ^ 2 := by ring

/-- **Segment model**: ζ(p + u d) is within ε + r_A + r_B u + C u² of A + u d B. -/
lemma seg_model {N : ℕ} (hN : 2 ≤ N) {p d : ℂ} (hd : ‖d‖ = 1) {h : ℝ}
    (hs1 : ∀ u ∈ Icc (0 : ℝ) h, p + u * d ≠ 1) {A B : ℂ} {rA rB C ε : ℝ}
    (hA : ‖Mf N p - A‖ ≤ rA) (hB : ‖M1f N p - B‖ ≤ rB)
    (hC : ∀ u ∈ Icc (0 : ℝ) h, ‖M2f N (p + u * d)‖ ≤ C)
    (hε : ∀ u ∈ Icc (0 : ℝ) h, ‖riemannZeta (p + u * d) - Mf N (p + u * d)‖ ≤ ε) :
    ∀ u ∈ Icc (0 : ℝ) h, ‖riemannZeta (p + u * d) - (A + u * d * B)‖ ≤ ε + rA + rB * u + C * u ^ 2 := by
  intro u hu
  have ht := taylor2 hd (fun v hv => hasDerivAt_Mf hN (hs1 v hv))
    (fun v hv => hasDerivAt_M1f hN (hs1 v hv)) hC u hu
  have e : riemannZeta (p + u * d) - (A + u * d * B) =
      (riemannZeta (p + u * d) - Mf N (p + u * d)) + (Mf N (p + u * d) - Mf N p - u * d * M1f N p)
        + (Mf N p - A) + (u * d) * (M1f N p - B) := by ring
  rw [e]
  have hu0 : ‖((u : ℝ) : ℂ) * d‖ = u := by
    rw [norm_mul, hd, mul_one, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hu.1]
  calc _ ≤ ‖riemannZeta (p + u * d) - Mf N (p + u * d)‖ + ‖Mf N (p + u * d) - Mf N p - u * d * M1f N p‖
        + ‖Mf N p - A‖ + ‖(u * d) * (M1f N p - B)‖ := by
        refine (norm_add_le _ _).trans ?_
        gcongr
        refine (norm_add_le _ _).trans ?_
        gcongr
        exact norm_add_le _ _
    _ ≤ ε + C * u ^ 2 + rA + u * rB := by
        rw [norm_mul, hu0]
        gcongr
        · exact hε u hu
        · exact hu.1
    _ = ε + rA + rB * u + C * u ^ 2 := by ring

/-- A model with error E and an anchor c with Re((model)·c̄) > ‖c‖E at both ends puts ζ·c̄ in
the right half plane along the whole segment. -/
lemma re_pos_of_model {p d : ℂ} {h : ℝ} (hh : 0 < h) {A B c : ℂ} {E cn : ℝ}
    (hmod : ∀ u ∈ Icc (0 : ℝ) h, ‖riemannZeta (p + u * d) - (A + u * d * B)‖ ≤ E)
    (hcn : ‖c‖ ≤ cn)
    (h0 : cn * E < (A * (starRingEnd ℂ) c).re)
    (h1 : cn * E < ((A + h * d * B) * (starRingEnd ℂ) c).re) :
    ∀ u ∈ Icc (0 : ℝ) h, 0 < (riemannZeta (p + u * d) * (starRingEnd ℂ) c).re := by
  intro u hu
  have hE : 0 ≤ E := (norm_nonneg _).trans (hmod u hu)
  set t := u / h with ht
  have ht0 : 0 ≤ t := div_nonneg hu.1 hh.le
  have ht1 : t ≤ 1 := (div_le_one hh).mpr hu.2
  have hlin : ((A + u * d * B) * (starRingEnd ℂ) c).re =
      (1 - t) * (A * (starRingEnd ℂ) c).re + t * ((A + h * d * B) * (starRingEnd ℂ) c).re := by
    have : (A + u * d * B) * (starRingEnd ℂ) c =
        ((1 - t : ℝ) : ℂ) * (A * (starRingEnd ℂ) c) + ((t : ℝ) : ℂ) * ((A + h * d * B) * (starRingEnd ℂ) c) := by
      have hh' : (h : ℂ) ≠ 0 := by exact_mod_cast hh.ne'
      rw [ht]; push_cast; field_simp; ring
    rw [this, add_re, re_ofReal_mul, re_ofReal_mul]
  have hM : cn * E < ((A + u * d * B) * (starRingEnd ℂ) c).re := by
    rw [hlin]
    rcases le_total (A * (starRingEnd ℂ) c).re ((A + h * d * B) * (starRingEnd ℂ) c).re with hle | hle
    · nlinarith [mul_le_mul_of_nonneg_left hle ht0]
    · nlinarith [mul_le_mul_of_nonneg_left hle (sub_nonneg.mpr ht1)]
  have herr : ‖(riemannZeta (p + u * d) - (A + u * d * B)) * (starRingEnd ℂ) c‖ ≤ E * cn := by
    rw [norm_mul, Complex.norm_conj]
    exact mul_le_mul (hmod u hu) hcn (norm_nonneg _) hE
  have hre := Complex.re_le_norm ((riemannZeta (p + u * d) - (A + u * d * B)) * (starRingEnd ℂ) c)
  have e : riemannZeta (p + u * d) * (starRingEnd ℂ) c =
      (A + u * d * B) * (starRingEnd ℂ) c + (riemannZeta (p + u * d) - (A + u * d * B)) * (starRingEnd ℂ) c := by
    ring
  rw [e, add_re]
  have : -(E * cn) ≤ ((riemannZeta (p + u * d) - (A + u * d * B)) * (starRingEnd ℂ) c).re := by
    have := Complex.abs_re_le_norm ((riemannZeta (p + u * d) - (A + u * d * B)) * (starRingEnd ℂ) c)
    rw [abs_le] at this; linarith [this.1]
  nlinarith

/-! ## Arguments from balls -/

lemma arg_le_of_ball_pos {w : ℂ} {a b r U : ℝ} (hw : ‖w - (a + b * I)‖ ≤ r) (har : r < a)
    (hb : 0 ≤ b + r) (hU : (b + r) / (a - r) ≤ U) : arg w ≤ U := by
  have hre : a - r ≤ w.re := by
    have := Complex.abs_re_le_norm (w - (a + b * I)); rw [abs_le] at this; simp at this; linarith
  have him : w.im ≤ b + r := by
    have := Complex.abs_im_le_norm (w - (a + b * I)); rw [abs_le] at this; simp at this; linarith
  have hre0 : 0 < w.re := by linarith
  have hθ := Complex.abs_arg_lt_pi_div_two_iff.mpr (Or.inl hre0)
  rw [abs_lt] at hθ
  have hq : 0 ≤ (b + r) / (a - r) := div_nonneg hb (by linarith)
  rcases lt_or_ge (arg w) 0 with hneg | hpos
  · linarith
  · have ht := Real.le_tan hpos hθ.2
    rw [Complex.tan_arg] at ht
    refine ht.trans (le_trans ?_ hU)
    rcases lt_or_ge w.im 0 with hi | hi
    · exact (div_neg_of_neg_of_pos hi hre0).le.trans hq
    · calc w.im / w.re ≤ (b + r) / w.re := by gcongr
        _ ≤ (b + r) / (a - r) := by gcongr

lemma arg_le_of_ball_neg {w : ℂ} {a b r M U : ℝ} (hw : ‖w - (a + b * I)‖ ≤ r) (har : r < a)
    (hb : b + r < 0) (hM : 0 ≤ M - r) (hM2 : a ^ 2 + b ^ 2 ≤ (M - r) ^ 2) (hU : (b + r) / M ≤ U) :
    arg w ≤ U := by
  have hre : a - r ≤ w.re := by
    have := Complex.abs_re_le_norm (w - (a + b * I)); rw [abs_le] at this; simp at this; linarith
  have him : w.im ≤ b + r := by
    have := Complex.abs_im_le_norm (w - (a + b * I)); rw [abs_le] at this; simp at this; linarith
  have hre0 : 0 < w.re := by linarith
  have hw0 : w ≠ 0 := by intro h; rw [h] at hre0; simp at hre0
  have hn0 : 0 < ‖w‖ := norm_pos_iff.mpr hw0
  have hc : ‖((a : ℂ) + b * I)‖ ≤ M - r := norm_le_of_sq_le hM (by simpa using hM2)
  have hnw : ‖w‖ ≤ M := by
    calc ‖w‖ = ‖(w - (a + b * I)) + (a + b * I)‖ := by ring_nf
      _ ≤ r + (M - r) := (norm_add_le _ _).trans (add_le_add hw hc)
      _ = M := by ring
  have hθneg : arg w < 0 := Complex.arg_neg_iff.mpr (by linarith)
  have hs := Real.sin_le (x := -arg w) (by linarith)
  rw [Real.sin_neg, Complex.sin_arg] at hs
  have h1 : arg w ≤ w.im / ‖w‖ := by linarith
  refine h1.trans (le_trans ?_ hU)
  calc w.im / ‖w‖ ≤ (b + r) / ‖w‖ := by gcongr
    _ ≤ (b + r) / M := by
        rw [div_le_div_iff₀ hn0 (by linarith)]
        nlinarith

lemma conj_ball {w : ℂ} {a b r : ℝ} (hw : ‖w - (a + b * I)‖ ≤ r) :
    ‖(starRingEnd ℂ) w - (a + (-b : ℝ) * I)‖ ≤ r := by
  have : (starRingEnd ℂ) w - (a + (-b : ℝ) * I) = (starRingEnd ℂ) (w - (a + b * I)) := by
    simp [map_sub, map_add, Complex.conj_ofReal]
  rw [this, Complex.norm_conj]; exact hw

lemma arg_conj_of_re_pos {w : ℂ} (hw : 0 < w.re) : arg ((starRingEnd ℂ) w) = -arg w := by
  rw [Complex.arg_conj, if_neg]
  intro h
  have := Complex.abs_arg_lt_pi_div_two_iff.mpr (Or.inl hw)
  rw [h, abs_of_pos Real.pi_pos] at this
  linarith [Real.pi_pos]

lemma le_arg_of_ball_pos {w : ℂ} {a b r L : ℝ} (hw : ‖w - (a + b * I)‖ ≤ r) (har : r < a)
    (hb : 0 ≤ r - b) (hL : L ≤ (b - r) / (a - r)) : L ≤ arg w := by
  have hre : a - r ≤ w.re := by
    have := Complex.abs_re_le_norm (w - (a + b * I)); rw [abs_le] at this; simp at this; linarith
  have h := arg_le_of_ball_pos (conj_ball hw) har (by linarith) (U := -L) (by
    rw [show -b + r = -(b - r) by ring, neg_div]; linarith)
  rw [arg_conj_of_re_pos (by linarith)] at h
  linarith

lemma le_arg_of_ball_neg {w : ℂ} {a b r M L : ℝ} (hw : ‖w - (a + b * I)‖ ≤ r) (har : r < a)
    (hb : r - b < 0) (hM : 0 ≤ M - r) (hM2 : a ^ 2 + b ^ 2 ≤ (M - r) ^ 2) (hL : L ≤ (b - r) / M) :
    L ≤ arg w := by
  have hre : a - r ≤ w.re := by
    have := Complex.abs_re_le_norm (w - (a + b * I)); rw [abs_le] at this; simp at this; linarith
  have h := arg_le_of_ball_neg (conj_ball hw) har (by linarith) hM (by nlinarith) (U := -L) (by
    rw [show -b + r = -(b - r) by ring, neg_div]; linarith)
  rw [arg_conj_of_re_pos (by linarith)] at h
  linarith

/-- Ball for z · conj(c), c = cr + ci·i exact. -/
lemma ball_mul_conj {z : ℂ} {a b r cr ci cn : ℝ} (hz : ‖z - (a + b * I)‖ ≤ r)
    (hc : ‖((cr : ℂ) + ci * I)‖ ≤ cn) :
    ‖z * (starRingEnd ℂ) ((cr : ℂ) + ci * I) - (((a * cr + b * ci : ℝ) : ℂ) + ((b * cr - a * ci : ℝ) : ℂ) * I)‖
      ≤ r * cn := by
  have e : z * (starRingEnd ℂ) ((cr : ℂ) + ci * I) - (((a * cr + b * ci : ℝ) : ℂ) + ((b * cr - a * ci : ℝ) : ℂ) * I)
      = (z - (a + b * I)) * (starRingEnd ℂ) ((cr : ℂ) + ci * I) := by
    simp only [map_add, map_mul, Complex.conj_ofReal, Complex.conj_I]
    push_cast
    ring_nf
    rw [I_sq]
    ring
  rw [e, norm_mul, Complex.norm_conj]
  exact mul_le_mul hz hc (norm_nonneg _) ((norm_nonneg _).trans hz)

end Fayed.FirstZero.Min
