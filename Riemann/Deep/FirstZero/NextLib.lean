import Riemann.Deep.FirstZero.FirstZero

/-!
# General machinery for certified zeros of ζ on the critical line

Generalizes the N = 6 / height-14 tools of Crit.lean, GammaRatio.lean, FirstZero.lean:
* `em_main_eqN`, `zeta_ballN`: order-8 Euler–Maclaurin main term for any N ≥ 2, written as
  1 + Σ_{2 ≤ n < N} n^{-s} + N^{-s}·EMKN N s;
* `Gamma_ratio_gen`, `Gammaℝ_ratio_gen`: Γ phase stability at any height ≥ 14;
* `crit_zero_of_balls`: two ζ-enclosures with a sign change of Λ ⇒ a zero in between.
-/

open Complex Set

namespace Fayed.FirstZero

/-- n^{-s}. -/
noncomputable def nps (s : ℂ) (n : ℕ) : ℂ := (n : ℂ) ^ (-s)

lemma nps_mul (s : ℂ) (a b : ℕ) : nps s (a * b) = nps s a * nps s b := by
  simp only [nps, Nat.cast_mul]
  exact natCast_mul_natCast_cpow a b (-s)

/-- The coefficient of N^{-s} in the Euler–Maclaurin main term. -/
noncomputable def EMKN (N : ℕ) (s : ℂ) : ℂ := 1 - (N : ℂ) / (1 - s) - 1 / 2 +
  s * (1 / 12 * gc s 0 / (N : ℂ) - 1 / 720 * gc s 2 / (N : ℂ) ^ 3 + 1 / 30240 * gc s 4 / (N : ℂ) ^ 5
    - 1 / 1209600 * gc s 6 / (N : ℂ) ^ 7)

lemma cpow_shiftN {N : ℕ} (hN : N ≠ 0) (s : ℂ) (j : ℕ) :
    ((((N : ℕ) : ℝ) : ℂ)) ^ (-s - 1 - (j : ℂ)) = (N : ℂ) ^ (-s) * ((N : ℂ) ^ (j + 1))⁻¹ := by
  have hN' : (N : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hN
  rw [Complex.ofReal_natCast, show -s - 1 - (j : ℂ) = -s + (-((j + 1 : ℕ) : ℂ)) by push_cast; ring,
    cpow_add _ _ hN', cpow_neg (N : ℂ) ((j + 1 : ℕ) : ℂ), cpow_natCast]

lemma em_main_eqN {N : ℕ} (hN : 2 ≤ N) (s : ℂ) (hs : s ≠ 0) (hs1 : s ≠ 1) :
    (∑ n ∈ Finset.range (N + 1), 1 / (n : ℂ) ^ s) + (-((N : ℕ) : ℂ) ^ (1 - s)) / (1 - s)
      + (-((N : ℕ) : ℂ) ^ (-s)) / 2 + s * EMPhi s ((N : ℕ) : ℝ) =
    1 + (∑ k ∈ Finset.range (N - 2), nps s (k + 2)) + nps s N * EMKN N s := by
  obtain ⟨M, rfl⟩ : ∃ M, N = M + 2 := ⟨N - 2, by omega⟩
  have hN0 : M + 2 ≠ 0 := by omega
  have hNc : ((M + 2 : ℕ) : ℂ) ≠ 0 := Nat.cast_ne_zero.mpr hN0
  have hsum : ∑ n ∈ Finset.range (M + 2 + 1), 1 / (n : ℂ) ^ s =
      1 + (∑ k ∈ Finset.range M, nps s (k + 2)) + 1 / ((M + 2 : ℕ) : ℂ) ^ s := by
    rw [Finset.sum_range_succ, Finset.sum_range_succ', Finset.sum_range_succ']
    have hterm : ∀ i ∈ Finset.range M, 1 / (((i + 1 + 1 : ℕ) : ℂ)) ^ s = nps s (i + 2) :=
      fun i _ => by rw [nps, show i + 1 + 1 = i + 2 by ring, cpow_neg, one_div]
    rw [Finset.sum_congr rfl hterm]
    simp only [Nat.cast_zero, zero_cpow hs, div_zero, add_zero, zero_add, Nat.cast_one, one_cpow,
      div_one]
    ring
  have h6' : ((M + 2 : ℕ) : ℂ) ^ (1 - s) = ((M + 2 : ℕ) : ℂ) * ((M + 2 : ℕ) : ℂ) ^ (-s) := by
    rw [sub_eq_add_neg, cpow_add _ _ hNc, cpow_one]
  have hs1' : (1 : ℂ) - s ≠ 0 := sub_ne_zero.mpr (Ne.symm hs1)
  have hinv : 1 / ((M + 2 : ℕ) : ℂ) ^ s = ((M + 2 : ℕ) : ℂ) ^ (-s) := by rw [cpow_neg, one_div]
  rw [hsum, hinv, h6']
  simp only [EMPhi, gf, cpow_shiftN hN0, nps, EMKN]
  have hM : M + 2 - 2 = M := by omega
  rw [hM]
  field_simp
  ring

lemma N_rpow_le {N : ℕ} (hN : 1 ≤ N) {c : ℝ} (hc0 : 0 < c) (hc : c ^ 2 ≤ N) :
    (N : ℝ) ^ (-((1 / 2 : ℝ) + 7)) ≤ 1 / ((N : ℝ) ^ 7 * c) := by
  have hN0 : (0 : ℝ) < N := by exact_mod_cast hN
  rw [Real.rpow_neg hN0.le, Real.rpow_add hN0, show (7 : ℝ) = ((7 : ℕ) : ℝ) by norm_num,
    Real.rpow_natCast, ← Real.sqrt_eq_rpow, one_div]
  have hs : c ≤ Real.sqrt N := Real.le_sqrt_of_sq_le hc
  apply inv_anti₀ (by positivity)
  rw [mul_comm]
  exact mul_le_mul_of_nonneg_right hs (by positivity)

/-- Turn a ball for the N-term Euler–Maclaurin main term into a ball for ζ(1/2 + it). -/
lemma zeta_ballN {N : ℕ} (hN : 2 ≤ N) {cN : ℝ} (hc0 : 0 < cN) (hc : cN ^ 2 ≤ N) (t : ℝ)
    (ht : t ≠ 0) {A B R1 R2 : ℝ}
    (hmain : ‖1 + (∑ k ∈ Finset.range (N - 2), nps ((1 / 2 : ℂ) + (t : ℂ) * I) (k + 2))
      + nps ((1 / 2 : ℂ) + (t : ℂ) * I) N * EMKN N ((1 / 2 : ℂ) + (t : ℂ) * I) - (A + B * I)‖ ≤ R1)
    (hem : ‖(1 / 2 : ℂ) + (t : ℂ) * I‖ * (‖gc ((1 / 2 : ℂ) + (t : ℂ) * I) 7‖ / 1000000 *
      ((1 / ((N : ℝ) ^ 7 * cN)) / (1 / 2 + 7))) ≤ R2) :
    ‖riemannZeta ((1 / 2 : ℂ) + (t : ℂ) * I) - (A + B * I)‖ ≤ R1 + R2 := by
  set s : ℂ := (1 / 2 : ℂ) + (t : ℂ) * I with hsdef
  have hre : s.re = 1 / 2 := by simp [hsdef]
  have hs0 : 0 < s.re := by rw [hre]; norm_num
  have hs1 : s ≠ 1 := by
    intro h; have := congrArg Complex.im h; simp [hsdef] at this; exact ht this
  have hsz : s ≠ 0 := by
    intro h; have := congrArg Complex.re h; rw [hre] at this; norm_num at this
  have hEM := zeta_EM8 hs0 hs1 (N := N) (by omega)
  rw [em_main_eqN hN s hsz hs1] at hEM
  have hb : ‖s‖ * (‖gc s 7‖ / 1000000 * ((N : ℝ) ^ (-(s.re + 7)) / (s.re + 7))) ≤ R2 := by
    refine le_trans ?_ hem
    rw [hre]
    gcongr
    exact N_rpow_le (by omega) hc0 hc
  calc ‖riemannZeta s - (A + B * I)‖ ≤ ‖riemannZeta s - (1 + (∑ k ∈ Finset.range (N - 2), nps s (k + 2))
        + nps s N * EMKN N s)‖ + ‖(1 + (∑ k ∈ Finset.range (N - 2), nps s (k + 2))
        + nps s N * EMKN N s) - (A + B * I)‖ := norm_sub_le_norm_sub_add_norm_sub _ _ _
    _ ≤ R2 + R1 := add_le_add (hEM.trans hb) hmain
    _ = R1 + R2 := add_comm _ _

/-- **Γ moves little along a short vertical segment** (any height ≥ 7). -/
theorem Gamma_ratio_gen (z : ℂ) (hre : z.re = 1 / 4) (him : 7 ≤ z.im)
    {δ : ℝ} (hδ0 : 0 ≤ δ) (hδ : δ ≤ 1) {S : ℝ}
    (hS : 2 * (3 * δ ^ 2 + δ + 2 * δ * (3 / 4 + z.im + δ)) + 2 * (δ / 7) ≤ S) :
    ‖Gamma (z + δ * I) / Gamma z - 1‖ ≤ Real.exp S - 1 := by
  obtain ⟨z2, hz2⟩ : ∃ z2, z2 = z + δ * I := ⟨_, rfl⟩
  rw [← hz2]
  have hre2 : z2.re = 1 / 4 := by rw [hz2]; simp [hre]
  have him2 : 7 ≤ z2.im := by rw [hz2]; simp; linarith
  have hz1 : ‖z2 - 1‖ ≤ 3 / 4 + z.im + δ := by
    refine (norm_le_abs_re_add_abs_im _).trans ?_
    have h1 : (z2 - 1).re = -3 / 4 := by rw [sub_re, hre2]; norm_num
    have h2 : (z2 - 1).im = z.im + δ := by rw [hz2]; simp
    rw [h1, h2, abs_of_nonneg (show 0 ≤ z.im + δ by linarith)]
    norm_num
    linarith
  have hzre : 0 < z.re := by rw [hre]; norm_num
  have hG : Gamma z ≠ 0 := Gamma_ne_zero_of_re_pos hzre
  have hlim := (((GammaSeq_tendsto_Gamma z2).div (GammaSeq_tendsto_Gamma z) hG).sub_const 1).norm
  refine le_of_tendsto hlim (Filter.eventually_atTop.2 ⟨1, fun n hn => ?_⟩)
  obtain ⟨m, rfl⟩ : ∃ m, n = m + 1 := ⟨n - 1, by omega⟩
  rw [Pi.div_apply, GammaSeq_ratio z z2 hzre δ hz2 m]
  set K := 3 * δ ^ 2 + δ + 2 * δ * ‖z2 - 1‖ with hK
  set e : ℕ → ℂ := fun j => (1 - δ * I / (z2 + j)) *
    cexp (δ * I * ((Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) : ℝ) : ℂ)) - 1 with he_def
  have he : ∀ j, ‖e j‖ ≤ K / ((j : ℝ) + 1) ^ 2 :=
    fun j => pair_bound z2 hre2 him2 hδ0 hδ j
  have hK0 : 0 ≤ K := by positivity
  have hS2 : ∑ j ∈ Finset.range m, ‖e j‖ ≤ 2 * K := by
    calc ∑ j ∈ Finset.range m, ‖e j‖ ≤ ∑ j ∈ Finset.range m, K / ((j : ℝ) + 1) ^ 2 :=
          Finset.sum_le_sum fun j _ => he j
      _ = K * ∑ j ∈ Finset.range m, 1 / ((j : ℝ) + 1) ^ 2 := by
          rw [Finset.mul_sum]; refine Finset.sum_congr rfl fun j _ => ?_; ring
      _ ≤ K * 2 := mul_le_mul_of_nonneg_left (sum_range_inv_sq_le_two m) hK0
      _ = 2 * K := by ring
  have hX := norm_prod_one_add_sub_one_le e m
  have hprod_eq : ∏ j ∈ Finset.range m, ((1 - δ * I / (z2 + j)) *
      cexp (δ * I * ((Real.log ((j : ℝ) + 2) - Real.log ((j : ℝ) + 1) : ℝ) : ℂ))) =
      ∏ j ∈ Finset.range m, (1 + e j) := by
    refine Finset.prod_congr rfl fun j _ => ?_
    simp only [he_def]; ring
  rw [hprod_eq]
  have ha : ∀ k : ℕ, ‖-(δ * I / (z2 + k))‖ ≤ δ / 7 := by
    intro k
    have h7 := norm_shift_ge_im z2 him2 k
    rw [norm_neg, norm_div, norm_mul, Complex.norm_real, Real.norm_eq_abs, Complex.norm_I, mul_one,
      abs_of_nonneg hδ0]
    exact div_le_div_of_nonneg_left hδ0 (by norm_num) h7
  have step1 := norm_mul_one_add_sub_one_le (e := -(δ * I / (z2 + m))) hX
  have step2 := norm_mul_one_add_sub_one_le (e := -(δ * I / (z2 + ((m + 1 : ℕ) : ℂ)))) step1
  have hfin : (∏ j ∈ Finset.range m, (1 + e j)) * (1 - δ * I / (z2 + m)) *
      (1 - δ * I / (z2 + ((m + 1 : ℕ) : ℂ))) =
      ((∏ j ∈ Finset.range m, (1 + e j)) * (1 + -(δ * I / (z2 + m)))) *
        (1 + -(δ * I / (z2 + ((m + 1 : ℕ) : ℂ)))) := by ring
  rw [hfin]
  refine step2.trans ?_
  have h1 := ha m
  have h2 := ha (m + 1)
  have hKb : K ≤ 3 * δ ^ 2 + δ + 2 * δ * (3 / 4 + z.im + δ) := by
    rw [hK]; nlinarith [mul_le_mul_of_nonneg_left hz1 hδ0]
  have hsum : ∑ j ∈ Finset.range m, ‖e j‖ + ‖-(δ * I / (z2 + m))‖
      + ‖-(δ * I / (z2 + ((m + 1 : ℕ) : ℂ)))‖ ≤ S := by
    linarith
  have hmono := Real.exp_le_exp.mpr hsum
  linarith

/-- Γℝ(1/2 + i t₂) / Γℝ(1/2 + i t₁) is close to 1 when t₂ - t₁ = 2δ is small. -/
lemma Gammaℝ_ratio_gen {t1 t2 δ : ℝ} (ht1 : 14 ≤ t1) (ht2 : t2 = t1 + 2 * δ) (hδ0 : 0 ≤ δ)
    (hδ : δ ≤ 1) {S : ℝ} (hS : 2 * (3 * δ ^ 2 + δ + 2 * δ * (3 / 4 + t1 / 2 + δ)) + 2 * (δ / 7) ≤ S)
    (hS0 : 0 < S) (hS1 : S < 1) {c : ℝ}
    (hc : 11 / 5 * δ * (1 / (1 - S)) + (1 / (1 - S) - 1) ≤ c) :
    ‖Gammaℝ ((1 / 2 : ℂ) + (t2 : ℂ) * I) / Gammaℝ ((1 / 2 : ℂ) + (t1 : ℂ) * I) - 1‖ ≤ c := by
  set z : ℂ := ((1 / 2 : ℂ) + (t1 : ℂ) * I) / 2 with hz
  have hs2 : ((1 / 2 : ℂ) + (t2 : ℂ) * I) / 2 = z + δ * I := by
    rw [hz, ht2]; push_cast; ring
  have hzre : z.re = 1 / 4 := by rw [hz]; simp; norm_num
  have hzim : z.im = t1 / 2 := by rw [hz]; simp
  have hG0 := Gamma_ratio_gen z hzre (by rw [hzim]; linarith) hδ0 hδ (S := S) (by rw [hzim]; exact hS)
  have hexp : Real.exp S < 1 / (1 - S) := Real.exp_bound_div_one_sub_of_interval' hS0 hS1
  have hG : ‖Gamma (z + δ * I) / Gamma z - 1‖ ≤ 1 / (1 - S) - 1 := by linarith
  have hpi0 : ((Real.pi : ℝ) : ℂ) ≠ 0 := by exact_mod_cast Real.pi_ne_zero
  have hP : (Real.pi : ℂ) ^ (-((1 / 2 : ℂ) + (t2 : ℂ) * I) / 2) =
      (Real.pi : ℂ) ^ (-((1 / 2 : ℂ) + (t1 : ℂ) * I) / 2) *
        cexp (I * ((-(δ * Real.log Real.pi)) : ℝ)) := by
    have hPI : (Real.pi : ℂ) ^ (-((δ : ℂ) * I)) = cexp (I * ((-(δ * Real.log Real.pi)) : ℝ)) := by
      rw [cpow_def_of_ne_zero hpi0, ← ofReal_log Real.pi_pos.le]
      congr 1; push_cast; ring
    rw [show -((1 / 2 : ℂ) + (t2 : ℂ) * I) / 2 =
      -((1 / 2 : ℂ) + (t1 : ℂ) * I) / 2 + (-((δ : ℂ) * I)) by
        rw [ht2]; push_cast; ring,
      cpow_add _ _ hpi0, hPI]
  have hlogpi : Real.log Real.pi ≤ 22 / 10 := by
    have := Real.log_le_sub_one_of_pos Real.pi_pos
    have := Real.pi_lt_d2
    norm_num at *; linarith
  have hlogpi0 : 0 ≤ Real.log Real.pi := Real.log_nonneg (by linarith [Real.pi_gt_three])
  have hPb : ‖cexp (I * ((-(δ * Real.log Real.pi)) : ℝ)) - 1‖ ≤ 11 / 5 * δ := by
    refine (Real.norm_exp_I_mul_ofReal_sub_one_le).trans ?_
    rw [Real.norm_eq_abs, abs_neg, abs_of_nonneg (by positivity)]
    nlinarith
  have hGne : Gamma z ≠ 0 := Gamma_ne_zero_of_re_pos (by rw [hzre]; norm_num)
  have hPne : (Real.pi : ℂ) ^ (-((1 / 2 : ℂ) + (t1 : ℂ) * I) / 2) ≠ 0 := by
    rw [Ne, cpow_eq_zero_iff]; tauto
  have hratio : Gammaℝ ((1 / 2 : ℂ) + (t2 : ℂ) * I) / Gammaℝ ((1 / 2 : ℂ) + (t1 : ℂ) * I) =
      cexp (I * ((-(δ * Real.log Real.pi)) : ℝ)) * (Gamma (z + δ * I) / Gamma z) := by
    simp only [Gammaℝ]
    rw [hP, hs2, ← hz]
    field_simp
  rw [hratio]
  set P := cexp (I * ((-(δ * Real.log Real.pi)) : ℝ))
  set G := Gamma (z + δ * I) / Gamma z
  have e : P * G - 1 = (P - 1) * G + (G - 1) := by ring
  have hGn : ‖G‖ ≤ 1 / (1 - S) := by
    have := norm_le_norm_add_norm_sub' G 1
    simp at this; linarith
  rw [e]
  calc ‖(P - 1) * G + (G - 1)‖ ≤ ‖P - 1‖ * ‖G‖ + ‖G - 1‖ :=
        (norm_add_le _ _).trans (by rw [norm_mul])
    _ ≤ 11 / 5 * δ * (1 / (1 - S)) + (1 / (1 - S) - 1) := by gcongr
    _ ≤ c := hc

lemma re_mul_conj_neg_gen {z1 z2 ρ : ℂ} {a1 b1 r1 a2 b2 r2 R c : ℝ} (hc0 : 0 ≤ c)
    (h1 : ‖z1 - (a1 + b1 * I)‖ ≤ r1) (h2 : ‖z2 - (a2 + b2 * I)‖ ≤ r2) (hρ : ‖ρ - 1‖ ≤ c)
    (hR : (|a2| + |b2|) * r1 + (|a1| + |-b1|) * r2 + r2 * r1 + |a2 * a1 - b2 * -b1 - (a2 * a1 + b2 * b1)|
      + |a2 * -b1 + b2 * a1 - (b2 * a1 - a2 * b1)| ≤ R)
    (hnum : (a2 * a1 + b2 * b1) + R + c * (|a2 * a1 + b2 * b1| + |b2 * a1 - a2 * b1| + R) < 0) :
    (z2 * (starRingEnd ℂ) z1 * ρ).re < 0 := by
  have hc : ‖(starRingEnd ℂ) z1 - (a1 + ((-b1 : ℝ) : ℂ) * I)‖ ≤ r1 := by
    rw [← norm_conj]
    convert h1 using 2
    simp [map_sub, map_add, map_mul, Complex.conj_ofReal, Complex.conj_I]
  have hw := ball_mul (e := a2 * a1 + b2 * b1) (f := b2 * a1 - a2 * b1) (R := R) h2 hc hR
  set w := z2 * (starRingEnd ℂ) z1
  set w0 : ℂ := ((a2 * a1 + b2 * b1 : ℝ) : ℂ) + ((b2 * a1 - a2 * b1 : ℝ) : ℂ) * I
  have hre : w.re ≤ (a2 * a1 + b2 * b1) + R := by
    have := (re_le_norm (w - w0)).trans hw
    simp [w0] at this
    linarith
  have hnorm : ‖w‖ ≤ |a2 * a1 + b2 * b1| + |b2 * a1 - a2 * b1| + R := by
    have h0 : ‖w0‖ ≤ |a2 * a1 + b2 * b1| + |b2 * a1 - a2 * b1| :=
      (norm_le_abs_re_add_abs_im _).trans (le_of_eq (by simp [w0]))
    have := norm_le_norm_add_norm_sub' w w0
    linarith
  have e : (w * ρ).re = w.re + (w * (ρ - 1)).re := by
    rw [mul_sub, mul_one, sub_re]; ring
  have h3 : (w * (ρ - 1)).re ≤ ‖w‖ * c := by
    refine (re_le_norm _).trans ?_
    rw [norm_mul]
    exact mul_le_mul_of_nonneg_left hρ (norm_nonneg _)
  have h4 : ‖w‖ * c ≤ (|a2 * a1 + b2 * b1| + |b2 * a1 - a2 * b1| + R) * c :=
    mul_le_mul_of_nonneg_right hnorm hc0
  rw [e]
  linarith

/-- **Sign change of Λ from two ζ-enclosures ⇒ a zero of ζ on the critical line in between.** -/
theorem crit_zero_of_balls {t1 t2 : ℝ} (hlt : t1 < t2) {a1 b1 r1 a2 b2 r2 R c : ℝ} (hc0 : 0 ≤ c)
    (h1 : ‖riemannZeta ((1 / 2 : ℂ) + (t1 : ℂ) * I) - (a1 + b1 * I)‖ ≤ r1)
    (h2 : ‖riemannZeta ((1 / 2 : ℂ) + (t2 : ℂ) * I) - (a2 + b2 * I)‖ ≤ r2)
    (hρ : ‖Gammaℝ ((1 / 2 : ℂ) + (t2 : ℂ) * I) / Gammaℝ ((1 / 2 : ℂ) + (t1 : ℂ) * I) - 1‖ ≤ c)
    (hR : (|a2| + |b2|) * r1 + (|a1| + |-b1|) * r2 + r2 * r1 + |a2 * a1 - b2 * -b1 - (a2 * a1 + b2 * b1)|
      + |a2 * -b1 + b2 * a1 - (b2 * a1 - a2 * b1)| ≤ R)
    (hnum : (a2 * a1 + b2 * b1) + R + c * (|a2 * a1 + b2 * b1| + |b2 * a1 - a2 * b1| + R) < 0) :
    ∃ t ∈ Ioo t1 t2, riemannZeta ((1 / 2 : ℂ) + (t : ℂ) * I) = 0 := by
  set f : ℝ → ℝ := fun t => (completedRiemannZeta ((1 / 2 : ℂ) + (t : ℂ) * I)).re with hf
  have hre : ∀ t : ℝ, 0 < ((1 / 2 : ℂ) + (t : ℂ) * I).re := fun t => by simp
  have hsign : f t1 * f t2 < 0 := by
    have key := re_mul_conj_neg_gen hc0 h1 h2 hρ hR hnum
    set s1 : ℂ := (1 / 2 : ℂ) + (t1 : ℂ) * I
    set s2 : ℂ := (1 / 2 : ℂ) + (t2 : ℂ) * I
    have hG1 : Gammaℝ s1 ≠ 0 := Gammaℝ_ne_zero_of_re_pos (hre _)
    have hL1 : completedRiemannZeta s1 = riemannZeta s1 * Gammaℝ s1 := completed_eq_mul (hre _)
    have hL2 : completedRiemannZeta s2 = riemannZeta s2 * Gammaℝ s2 := completed_eq_mul (hre _)
    have hi1 := Lam_im t1
    have hi2 := Lam_im t2
    have hprod : completedRiemannZeta s2 * (starRingEnd ℂ) (completedRiemannZeta s1) =
        (riemannZeta s2 * (starRingEnd ℂ) (riemannZeta s1) * (Gammaℝ s2 / Gammaℝ s1)) *
          ((normSq (Gammaℝ s1) : ℝ) : ℂ) := by
      rw [hL1, hL2, ← mul_conj, map_mul]
      field_simp
    have hns : 0 < normSq (Gammaℝ s1) := normSq_pos.mpr hG1
    have hre_prod : (completedRiemannZeta s2 * (starRingEnd ℂ) (completedRiemannZeta s1)).re =
        f t2 * f t1 := by
      simp only [hf, mul_re, conj_re, conj_im]
      rw [hi1, hi2]; ring
    rw [hprod, re_mul_ofReal] at hre_prod
    have : f t2 * f t1 < 0 := by rw [← hre_prod]; exact mul_neg_of_neg_of_pos key hns
    linarith [mul_comm (f t1) (f t2)]
  have hcont : ContinuousOn f (Icc t1 t2) := by
    intro t _
    have hs0 : ((1 / 2 : ℂ) + (t : ℂ) * I) ≠ 0 := by
      intro h; have := congrArg Complex.re h; simp at this
    have hs1 : ((1 / 2 : ℂ) + (t : ℂ) * I) ≠ 1 := by
      intro h; have := congrArg Complex.re h; norm_num at this
    have hc1 : ContinuousAt (fun t : ℝ => (1 / 2 : ℂ) + (t : ℂ) * I) t := by fun_prop
    have hc2 := (differentiableAt_completedZeta hs0 hs1).continuousAt
    have hc3 := ContinuousAt.comp (g := completedRiemannZeta)
      (f := fun t : ℝ => (1 / 2 : ℂ) + (t : ℂ) * I) (x := t) hc2 hc1
    exact (continuous_re.continuousAt.comp hc3).continuousWithinAt
  obtain ⟨t, ht, hft⟩ : ∃ t ∈ Ioo t1 t2, f t = 0 := by
    rcases lt_or_gt_of_ne (show f t1 ≠ 0 by intro h; rw [h, zero_mul] at hsign; exact lt_irrefl _ hsign)
      with hneg | hpos
    · have h2' : 0 < f t2 := by
        by_contra h; push Not at h; nlinarith
      exact intermediate_value_Ioo hlt.le hcont ⟨hneg, h2'⟩
    · have h2' : f t2 < 0 := by
        by_contra h; push Not at h; nlinarith
      exact intermediate_value_Ioo' hlt.le hcont ⟨h2', hpos⟩
  refine ⟨t, ht, ?_⟩
  have hΛ : completedRiemannZeta ((1 / 2 : ℂ) + (t : ℂ) * I) = 0 :=
    Complex.ext hft (by simpa using Lam_im t)
  have hs0 : ((1 / 2 : ℂ) + (t : ℂ) * I) ≠ 0 := by
    intro h; have := congrArg Complex.re h; simp at this
  rw [riemannZeta_def_of_ne_zero hs0, hΛ, zero_div]

/-- log 7, log 11, log 13 to 3·10⁻⁹. -/
lemma log7_near : |Real.log 7 - 19459101491 / 10 ^ 10| ≤ 3 / 10 ^ 9 := by
  have h := Real.abs_log_sub_add_sum_range_le (x := (1 / 8 : ℝ)) (by norm_num [abs_of_pos]) 20
  have e : Real.log 7 = 3 * Real.log 2 + Real.log (1 - 1 / 8) := by
    rw [show (1 : ℝ) - 1 / 8 = 7 / 2 ^ 3 by norm_num, Real.log_div (by norm_num) (by norm_num),
      Real.log_pow]
    push_cast; ring
  have h2 := log2_near
  rw [abs_le] at h h2 ⊢
  norm_num [Finset.sum_range_succ] at h h2 ⊢
  rw [e]
  norm_num
  constructor <;> linarith [h.1, h.2, h2.1, h2.2]

lemma log11_near : |Real.log 11 - 23978952728 / 10 ^ 10| ≤ 3 / 10 ^ 9 := by
  have h := Real.abs_log_sub_add_sum_range_le (x := (1 / 12 : ℝ)) (by norm_num [abs_of_pos]) 20
  have e : Real.log 11 = 2 * Real.log 2 + Real.log 3 + Real.log (1 - 1 / 12) := by
    rw [show (1 : ℝ) - 1 / 12 = 11 / (2 ^ 2 * 3) by norm_num, Real.log_div (by norm_num) (by norm_num),
      Real.log_mul (by norm_num) (by norm_num), Real.log_pow]
    push_cast; ring
  have h2 := log2_near
  have h3 := log3_near
  rw [abs_le] at h h2 h3 ⊢
  norm_num [Finset.sum_range_succ] at h h2 h3 ⊢
  rw [e]
  norm_num
  constructor <;> linarith [h.1, h.2, h2.1, h2.2, h3.1, h3.2]

lemma log13_near : |Real.log 13 - 25649493575 / 10 ^ 10| ≤ 3 / 10 ^ 9 := by
  have h := Real.abs_log_sub_add_sum_range_le (x := (-1 / 12 : ℝ)) (by norm_num [abs_of_neg]) 20
  have e : Real.log 13 = 2 * Real.log 2 + Real.log 3 + Real.log (1 - (-1 / 12)) := by
    rw [show (1 : ℝ) - (-1 / 12) = 13 / (2 ^ 2 * 3) by norm_num, Real.log_div (by norm_num) (by norm_num),
      Real.log_mul (by norm_num) (by norm_num), Real.log_pow]
    push_cast; ring
  have h2 := log2_near
  have h3 := log3_near
  rw [abs_le] at h h2 h3 ⊢
  norm_num [Finset.sum_range_succ, abs_of_neg] at h h2 h3 ⊢
  rw [e]
  norm_num
  constructor <;> linarith [h.1, h.2, h2.1, h2.2, h3.1, h3.2]

end Fayed.FirstZero
