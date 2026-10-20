import Riemann.Deep.FirstZero.MinCount

/-!
# Minimality of the first zero, part 4: packaged certificate lemmas

Everything a generated certificate needs, with purely numeric side conditions:
nps_ball (n^{-s} at s = σ + it, σ = p/q), neg_lg_mul_ball, rpow bounds, zeta_node (ball for ζ at a
node), arg_up_* / arg_lo_* (arg of ζ·c̄ from a node ball), leftOK_seg / topOK_seg (one certified
piece of the contour).
-/

open Complex Set

namespace Fayed.FirstZero.Min
open Fayed.FirstZero

lemma nps_mul' {s : ℂ} (a b : ℕ) {n : ℕ} (h : a * b = n) : nps s n = nps s a * nps s b := by
  rw [← h, nps_mul]

lemma rpow_ball_q {n σ : ℝ} (hn : 0 < n) (p q : ℕ) (hq : 0 < q) (hσ : σ = (p : ℝ) / q) {c ε : ℝ}
    (hε : 0 ≤ ε) (hc : 0 < c - ε) (h1 : (c - ε) ^ q * n ^ p ≤ 1) (h2 : 1 ≤ (c + ε) ^ q * n ^ p) :
    |n ^ (-σ) - c| ≤ ε := by
  subst hσ; exact rpow_ball hn p q hq hε hc h1 h2

lemma rpow_le_q {n σ : ℝ} (hn : 0 < n) (p q : ℕ) (hq : 0 < q) (hσ : σ = (p : ℝ) / q) {v : ℝ}
    (hv : 0 < v) (h : 1 ≤ v ^ q * n ^ p) : n ^ (-σ) ≤ v := by
  subst hσ
  set x := n ^ (-((p : ℝ) / q)) with hx
  have hx0 : 0 < x := Real.rpow_pos_of_pos hn _
  have hq' : (q : ℝ) ≠ 0 := by positivity
  have hnp : 0 < n ^ p := pow_pos hn p
  have hxx : x ^ q * n ^ p = 1 := by
    rw [hx, ← Real.rpow_natCast, ← Real.rpow_mul hn.le,
      show -((p : ℝ) / q) * (q : ℝ) = -(p : ℝ) by field_simp,
      Real.rpow_neg hn.le, Real.rpow_natCast, inv_mul_cancel₀ hnp.ne']
  by_contra hcon
  push Not at hcon
  have : v ^ q < x ^ q := pow_lt_pow_left₀ hcon hv.le (by omega)
  nlinarith [mul_lt_mul_of_pos_right this hnp]

lemma nps_ball {n : ℕ} (hn : 0 < n) {σ t : ℝ} {cr ε a b r e f R : ℝ}
    (hpow : |(n : ℝ) ^ (-σ) - cr| ≤ ε)
    (hcis : ‖cexp (((-(t * Real.log n)) : ℝ) * I) - (a + b * I)‖ ≤ r)
    (h : (|a| + |b|) * ε + |cr| * r + ε * r + |cr * a - e| + |cr * b - f| ≤ R) :
    ‖nps ((σ : ℂ) + (t : ℂ) * I) n - (e + f * I)‖ ≤ R := by
  rw [nps, show ((n : ℕ) : ℂ) = ((n : ℝ) : ℂ) by simp, cpow_gen (by exact_mod_cast hn)]
  exact ball_real_mul hpow hcis h

lemma neg_lg_mul_ball {n : ℕ} {z : ℂ} {c ε a b r e f R : ℝ} (hx : |Real.log n - c| ≤ ε)
    (hz : ‖z - (a + b * I)‖ ≤ r)
    (h : (|a| + |b|) * ε + |-c| * r + ε * r + |-c * a - e| + |-c * b - f| ≤ R) :
    ‖-lg n * z - (e + f * I)‖ ≤ R := by
  have hx' : |(-Real.log n) - (-c)| ≤ ε := by
    rw [show -Real.log n - -c = -(Real.log n - c) by ring, abs_neg]
    exact hx
  have := ball_real_mul hx' hz h
  rwa [show ((-Real.log n : ℝ) : ℂ) * z = -lg n * z by simp [lg]] at this

lemma lg_ball {n : ℕ} {c ε : ℝ} (hx : |Real.log n - c| ≤ ε) :
    ‖lg n - ((c : ℝ) + (0 : ℝ) * I)‖ ≤ ε := by
  rw [show lg n - ((c : ℝ) + (0 : ℝ) * I) = ((Real.log n - c : ℝ) : ℂ) by simp [lg],
    Complex.norm_real, Real.norm_eq_abs]
  exact hx

lemma log_near_mul {m n : ℕ} {c1 e1 c2 e2 : ℝ} (hm : 0 < m) (hn : 0 < n)
    (h1 : |Real.log m - c1| ≤ e1) (h2 : |Real.log n - c2| ≤ e2) :
    |Real.log ((m * n : ℕ) : ℝ) - (c1 + c2)| ≤ e1 + e2 := by
  rw [Nat.cast_mul, Real.log_mul (by positivity) (by positivity)]
  calc |Real.log m + Real.log n - (c1 + c2)| = |(Real.log m - c1) + (Real.log n - c2)| := by ring_nf
    _ ≤ _ := (abs_add_le _ _).trans (add_le_add h1 h2)

lemma log_le_of_near {n : ℕ} {c ε L : ℝ} (h : |Real.log n - c| ≤ ε) (hL : c + ε ≤ L) :
    Real.log n ≤ L := by
  rw [abs_le] at h; linarith

lemma log2_near' : |Real.log ((2 : ℕ) : ℝ) - 69314718055 / 10 ^ 11| ≤ 3 / 10 ^ 10 := by
  simpa using log2_near
lemma log3_near' : |Real.log ((3 : ℕ) : ℝ) - 10986122887 / 10 ^ 10| ≤ 1 / 10 ^ 9 := by
  simpa using log3_near
lemma log5_near' : |Real.log ((5 : ℕ) : ℝ) - 16094379124 / 10 ^ 10| ≤ 1 / 10 ^ 9 := by
  simpa using log5_near
lemma log7_near' : |Real.log ((7 : ℕ) : ℝ) - 19459101491 / 10 ^ 10| ≤ 3 / 10 ^ 9 := by
  simpa using log7_near
lemma log4_near' : |Real.log ((4 : ℕ) : ℝ) - (69314718055 / 10 ^ 11 + 69314718055 / 10 ^ 11)| ≤
    3 / 10 ^ 10 + 3 / 10 ^ 10 :=
  log_near_mul (m := 2) (n := 2) (by norm_num) (by norm_num) log2_near' log2_near'
lemma log6_near' : |Real.log ((6 : ℕ) : ℝ) - (69314718055 / 10 ^ 11 + 10986122887 / 10 ^ 10)| ≤
    3 / 10 ^ 10 + 1 / 10 ^ 9 :=
  log_near_mul (m := 2) (n := 3) (by norm_num) (by norm_num) log2_near' log3_near'
lemma log8_near' : |Real.log ((8 : ℕ) : ℝ) - (69314718055 / 10 ^ 11 + (69314718055 / 10 ^ 11 +
    69314718055 / 10 ^ 11))| ≤ 3 / 10 ^ 10 + (3 / 10 ^ 10 + 3 / 10 ^ 10) :=
  log_near_mul (m := 2) (n := 4) (by norm_num) (by norm_num) log2_near' log4_near'

/-- Crude M2 bound with the Dirichlet part bounded termwise by log(N−1)²·2^{−σ₀}. -/
lemma norm_M2f_le' {N : ℕ} (hN : 2 ≤ N) {s : ℂ} {R D σ0 v v2 L L1 : ℝ} (hR : ‖s‖ ≤ R) (hD0 : 0 < D)
    (hD : D ≤ ‖1 - s‖) (hσ00 : 0 ≤ σ0) (hσ : σ0 ≤ s.re)
    (hv : (N : ℝ) ^ (-σ0) ≤ v) (hv2 : (2 : ℝ) ^ (-σ0) ≤ v2) (hL : Real.log N ≤ L)
    (hL1 : Real.log ((N - 1 : ℕ) : ℝ) ≤ L1) :
    ‖M2f N s‖ ≤ ((N - 2 : ℕ) : ℝ) * (L1 ^ 2 * v2) +
      v * (K2b N R D + 2 * L * K1b N R D + L ^ 2 * K0b N R D) := by
  refine norm_M2f_le hN hR hD0 hD hσ ?_ hv hL
  have hterm : ∀ k ∈ Finset.range (N - 2),
      Real.log ((k + 2 : ℕ) : ℝ) ^ 2 * ((k + 2 : ℕ) : ℝ) ^ (-σ0) ≤ L1 ^ 2 * v2 := by
    intro k hk
    rw [Finset.mem_range] at hk
    have hk2 : (2 : ℝ) ≤ ((k + 2 : ℕ) : ℝ) := by push_cast; linarith [(Nat.cast_nonneg k : (0:ℝ) ≤ k)]
    have hlog0 : 0 ≤ Real.log ((k + 2 : ℕ) : ℝ) := Real.log_nonneg (by linarith)
    have hlog : Real.log ((k + 2 : ℕ) : ℝ) ≤ L1 := by
      refine le_trans (Real.log_le_log (by linarith) ?_) hL1
      exact_mod_cast (by omega : k + 2 ≤ N - 1)
    have hp : ((k + 2 : ℕ) : ℝ) ^ (-σ0) ≤ v2 :=
      le_trans (Real.rpow_le_rpow_of_nonpos (by norm_num) hk2 (by linarith)) hv2
    have hp0 : 0 ≤ ((k + 2 : ℕ) : ℝ) ^ (-σ0) := Real.rpow_nonneg (by linarith) _
    exact mul_le_mul (pow_le_pow_left₀ hlog0 hlog 2) hp hp0 (sq_nonneg _)
  refine (Finset.sum_le_sum hterm).trans (le_of_eq ?_)
  rw [Finset.sum_const, Finset.card_range, nsmul_eq_mul]

/-- Node ball: ζ(s) from a ball for Mf N s and the Euler–Maclaurin remainder. -/
lemma zeta_node {N : ℕ} (hN : 2 ≤ N) {σ t : ℝ} (hσ0 : 0 < σ) (ht : 0 ≤ t) (hs1 : (σ : ℂ) + t * I ≠ 1)
    {a b r R v ε : ℝ} (hM : ‖Mf N ((σ : ℂ) + t * I) - (a + b * I)‖ ≤ r)
    (hR0 : 0 ≤ R) (hR : σ ^ 2 + t ^ 2 ≤ R ^ 2) (hv : (N : ℝ) ^ (-σ) ≤ v)
    (cj : ℕ → ℝ) (hc : ∀ j : ℕ, j < 7 → 0 ≤ cj j ∧ (σ + 1 + j) ^ 2 + t ^ 2 ≤ cj j ^ 2)
    (hε : R * ((∏ j ∈ Finset.range 7, cj j) / 1000000 * (v / (N : ℝ) ^ 7 / (σ + 7))) ≤ ε) :
    ‖riemannZeta ((σ : ℂ) + t * I) - (a + b * I)‖ ≤ r + ε := by
  have hre : ((σ : ℂ) + t * I).re = σ := by simp
  have him : ((σ : ℂ) + t * I).im = t := by simp
  have hR' : ‖(σ : ℂ) + t * I‖ ≤ R := norm_le_of_sq_le hR0 (by rw [hre, him]; exact hR)
  have hcj : ∀ j : ℕ, j < 7 → ‖(σ : ℂ) + t * I + 1 + (j : ℂ)‖ ≤ cj j := by
    intro j hj
    obtain ⟨h0, h1⟩ := hc j hj
    have := norm_shift_le (s := (σ : ℂ) + t * I) (x1 := σ) (Y := t) (a := 1 + j) (c := cj j)
      (by positivity) (by rw [hre]; linarith) (by rw [hre]) (by rw [him, abs_of_nonneg ht]) h0
      (by linarith)
    simpa [add_assoc] using this
  have he := em_err_le hN cj hσ0 (by rw [hre]) hs1 hR' hcj hv hε
  calc _ ≤ ‖riemannZeta ((σ : ℂ) + t * I) - Mf N ((σ : ℂ) + t * I)‖
        + ‖Mf N ((σ : ℂ) + t * I) - (a + b * I)‖ := norm_sub_le_norm_sub_add_norm_sub _ _ _
    _ ≤ ε + r := add_le_add he hM
    _ = r + ε := add_comm _ _

/-! ## Arg bounds for ζ·c̄ from a node ball -/

lemma norm_c_le {cr ci cn : ℝ} (hcn0 : 0 ≤ cn) (hcn : cr ^ 2 + ci ^ 2 ≤ cn ^ 2) :
    ‖((cr : ℂ) + ci * I)‖ ≤ cn := norm_le_of_sq_le hcn0 (by simpa using hcn)

lemma arg_up_pos {z : ℂ} {a b r cr ci cn U : ℝ} (hz : ‖z - (a + b * I)‖ ≤ r) (hcn0 : 0 ≤ cn)
    (hcn : cr ^ 2 + ci ^ 2 ≤ cn ^ 2) (h1 : r * cn < a * cr + b * ci)
    (h2 : 0 ≤ (b * cr - a * ci) + r * cn)
    (h3 : ((b * cr - a * ci) + r * cn) / ((a * cr + b * ci) - r * cn) ≤ U) :
    arg (z * (starRingEnd ℂ) ((cr : ℂ) + ci * I)) ≤ U :=
  arg_le_of_ball_pos (ball_mul_conj hz (norm_c_le hcn0 hcn)) h1 h2 h3

lemma arg_up_neg {z : ℂ} {a b r cr ci cn M U : ℝ} (hz : ‖z - (a + b * I)‖ ≤ r) (hcn0 : 0 ≤ cn)
    (hcn : cr ^ 2 + ci ^ 2 ≤ cn ^ 2) (h1 : r * cn < a * cr + b * ci)
    (h2 : (b * cr - a * ci) + r * cn < 0) (hM : 0 ≤ M - r * cn)
    (hM2 : (a * cr + b * ci) ^ 2 + (b * cr - a * ci) ^ 2 ≤ (M - r * cn) ^ 2)
    (h3 : ((b * cr - a * ci) + r * cn) / M ≤ U) :
    arg (z * (starRingEnd ℂ) ((cr : ℂ) + ci * I)) ≤ U :=
  arg_le_of_ball_neg (ball_mul_conj hz (norm_c_le hcn0 hcn)) h1 h2 hM hM2 h3

lemma arg_lo_pos {z : ℂ} {a b r cr ci cn L : ℝ} (hz : ‖z - (a + b * I)‖ ≤ r) (hcn0 : 0 ≤ cn)
    (hcn : cr ^ 2 + ci ^ 2 ≤ cn ^ 2) (h1 : r * cn < a * cr + b * ci)
    (h2 : 0 ≤ r * cn - (b * cr - a * ci))
    (h3 : L ≤ ((b * cr - a * ci) - r * cn) / ((a * cr + b * ci) - r * cn)) :
    L ≤ arg (z * (starRingEnd ℂ) ((cr : ℂ) + ci * I)) :=
  le_arg_of_ball_pos (ball_mul_conj hz (norm_c_le hcn0 hcn)) h1 h2 h3

lemma arg_lo_neg {z : ℂ} {a b r cr ci cn M L : ℝ} (hz : ‖z - (a + b * I)‖ ≤ r) (hcn0 : 0 ≤ cn)
    (hcn : cr ^ 2 + ci ^ 2 ≤ cn ^ 2) (h1 : r * cn < a * cr + b * ci)
    (h2 : r * cn - (b * cr - a * ci) < 0) (hM : 0 ≤ M - r * cn)
    (hM2 : (a * cr + b * ci) ^ 2 + (b * cr - a * ci) ^ 2 ≤ (M - r * cn) ^ 2)
    (h3 : L ≤ ((b * cr - a * ci) - r * cn) / M) :
    L ≤ arg (z * (starRingEnd ℂ) ((cr : ℂ) + ci * I)) :=
  le_arg_of_ball_neg (ball_mul_conj hz (norm_c_le hcn0 hcn)) h1 h2 hM hM2 h3

/-- arg ζ from a ball (anchor 1). -/
lemma arg_up_one {z : ℂ} {a b r U : ℝ} (hz : ‖z - (a + b * I)‖ ≤ r) (h1 : r < a) (h2 : 0 ≤ b + r)
    (h3 : (b + r) / (a - r) ≤ U) : arg z ≤ U := arg_le_of_ball_pos hz h1 h2 h3

/-! ## One certified piece of the contour -/

lemma conjc_re (Ar Ai cr ci : ℝ) :
    ((((Ar : ℂ) + Ai * I)) * (starRingEnd ℂ) ((cr : ℂ) + ci * I)).re = Ar * cr + Ai * ci := by
  simp [Complex.mul_re, Complex.conj_ofReal]

lemma leftOK_seg {σ a b : ℝ} (hσ0 : 0 < σ) (hσ1 : σ < 1) (ha : 0 ≤ a) (hab : a < b) {N : ℕ}
    (hN : 2 ≤ N) {Ar Ai Br Bi rA rB : ℝ}
    (hA : ‖Mf N ((σ : ℂ) + a * I) - (Ar + Ai * I)‖ ≤ rA)
    (hB : ‖M1f N ((σ : ℂ) + a * I) - (Br + Bi * I)‖ ≤ rB)
    {R v v2 L L1 C ε Dd : ℝ} (hR0 : 0 ≤ R) (hR : σ ^ 2 + b ^ 2 ≤ R ^ 2)
    (hD0 : 0 < Dd) (hD : Dd ^ 2 ≤ (1 - σ) ^ 2 + a ^ 2)
    (hv : (N : ℝ) ^ (-σ) ≤ v) (hv2 : (2 : ℝ) ^ (-σ) ≤ v2) (hL : Real.log N ≤ L)
    (hL1 : Real.log ((N - 1 : ℕ) : ℝ) ≤ L1)
    (hC : ((N - 2 : ℕ) : ℝ) * (L1 ^ 2 * v2) +
      v * (K2b N R Dd + 2 * L * K1b N R Dd + L ^ 2 * K0b N R Dd) ≤ C)
    (cj : ℕ → ℝ) (hc : ∀ j : ℕ, j < 7 → 0 ≤ cj j ∧ (σ + 1 + j) ^ 2 + b ^ 2 ≤ cj j ^ 2)
    (hε : R * ((∏ j ∈ Finset.range 7, cj j) / 1000000 * (v / (N : ℝ) ^ 7 / (σ + 7))) ≤ ε)
    {cr ci cn E : ℝ} (hcn0 : 0 ≤ cn) (hcn : cr ^ 2 + ci ^ 2 ≤ cn ^ 2)
    (hE : ε + rA + rB * (b - a) + C * (b - a) ^ 2 ≤ E)
    (h0 : cn * E < Ar * cr + Ai * ci)
    (h1 : cn * E < (Ar - (b - a) * Bi) * cr + (Ai + (b - a) * Br) * ci)
    {Ua Lb U : ℝ}
    (hUa : arg (riemannZeta ((σ : ℂ) + a * I) * (starRingEnd ℂ) ((cr : ℂ) + ci * I)) ≤ Ua)
    (hLb : Lb ≤ arg (riemannZeta ((σ : ℂ) + b * I) * (starRingEnd ℂ) ((cr : ℂ) + ci * I)))
    (hU : Ua - Lb ≤ U) : LeftOK σ a b U := by
  set h := b - a with hh
  have hh0 : 0 < h := by rw [hh]; linarith
  set p : ℂ := (σ : ℂ) + a * I with hp
  have hpu : ∀ u : ℝ, p + (u : ℂ) * I = (σ : ℂ) + ((a + u : ℝ) : ℂ) * I := by
    intro u; rw [hp]; push_cast; ring
  have hs1 : ∀ u ∈ Icc (0 : ℝ) h, p + (u : ℂ) * I ≠ 1 := by
    intro u _; rw [hpu]; exact ne_one_of_re hσ1.ne
  have hreu : ∀ u : ℝ, (p + (u : ℂ) * I).re = σ := by intro u; rw [hpu]; simp
  have himu : ∀ u : ℝ, (p + (u : ℂ) * I).im = a + u := by intro u; rw [hpu]; simp
  have hnorm : ∀ u ∈ Icc (0 : ℝ) h, ‖p + (u : ℂ) * I‖ ≤ R := by
    intro u hu
    refine norm_le_of_sq_le hR0 ?_
    rw [hreu, himu]
    have : (a + u) ^ 2 ≤ b ^ 2 := pow_le_pow_left₀ (by linarith [hu.1]) (by linarith [hu.2]) 2
    linarith
  have hC' : ∀ u ∈ Icc (0 : ℝ) h, ‖M2f N (p + (u : ℂ) * I)‖ ≤ C := by
    intro u hu
    refine le_trans (norm_M2f_le' hN (hnorm u hu) hD0 ?_ hσ0.le (by rw [hreu]) hv hv2 hL
      hL1) hC
    refine (pow_le_pow_iff_left₀ hD0.le (norm_nonneg _) two_ne_zero).mp ?_
    rw [Complex.sq_norm, Complex.normSq_apply, sub_re, one_re, sub_im, one_im, hreu, himu]
    have : a ^ 2 ≤ (a + u) ^ 2 := pow_le_pow_left₀ ha (by linarith [hu.1]) 2
    nlinarith
  have hε' : ∀ u ∈ Icc (0 : ℝ) h, ‖riemannZeta (p + (u : ℂ) * I) - Mf N (p + (u : ℂ) * I)‖ ≤ ε := by
    intro u hu
    refine em_err_le hN cj hσ0 (by rw [hreu]) (hs1 u hu) (hnorm u hu) ?_ hv hε
    intro j hj
    obtain ⟨hc0, hc1⟩ := hc j hj
    have := norm_shift_le (s := p + (u : ℂ) * I) (x1 := σ) (Y := b) (a := 1 + j) (c := cj j)
      (by positivity) (by rw [hreu]; linarith) (by rw [hreu])
      (by rw [himu, abs_of_nonneg (by linarith [hu.1])]; linarith [hu.2]) hc0 (by linarith)
    simpa [add_assoc] using this
  have hmod := seg_model hN (p := p) (d := I) (by simp) hs1 hA hB hC' hε'
  have hrB : 0 ≤ rB := (norm_nonneg _).trans hB
  have hC0 : 0 ≤ C := (norm_nonneg _).trans (hC' 0 ⟨le_rfl, hh0.le⟩)
  have hmod' : ∀ u ∈ Icc (0 : ℝ) h, ‖riemannZeta (p + (u : ℂ) * I) - (((Ar : ℂ) + Ai * I) +
      (u : ℂ) * I * ((Br : ℂ) + Bi * I))‖ ≤ E := by
    intro u hu
    refine (hmod u hu).trans (le_trans ?_ hE)
    have : rB * u ≤ rB * h := mul_le_mul_of_nonneg_left hu.2 hrB
    have : C * u ^ 2 ≤ C * h ^ 2 :=
      mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hu.1 hu.2 2) hC0
    linarith
  have hslit := re_pos_of_model hh0 hmod' (norm_c_le hcn0 hcn)
    (by rw [conjc_re]; exact h0)
    (by
      have e : (((Ar : ℂ) + Ai * I) + (h : ℂ) * I * ((Br : ℂ) + Bi * I)) =
          ((Ar - h * Bi : ℝ) : ℂ) + ((Ai + h * Br : ℝ) : ℂ) * I := by
        push_cast; ring_nf; rw [I_sq]; ring
      rw [e, conjc_re]; exact h1)
  refine leftOK_of_slit hσ1.ne hab.le (fun y hy => ?_) hUa hLb hU
  have := hslit (y - a) ⟨by linarith [hy.1], by rw [hh]; linarith [hy.2]⟩
  rwa [hpu, show a + (y - a) = y by ring] at this

lemma topOK_seg {T a b : ℝ} (hT : 0 < T) (ha : 0 < a) (hab : a < b) {N : ℕ}
    (hN : 2 ≤ N) {Ar Ai Br Bi rA rB : ℝ}
    (hA : ‖Mf N ((a : ℂ) + T * I) - (Ar + Ai * I)‖ ≤ rA)
    (hB : ‖M1f N ((a : ℂ) + T * I) - (Br + Bi * I)‖ ≤ rB)
    {R v v2 L L1 C ε : ℝ} (hR0 : 0 ≤ R) (hR : b ^ 2 + T ^ 2 ≤ R ^ 2)
    (hv : (N : ℝ) ^ (-a) ≤ v) (hv2 : (2 : ℝ) ^ (-a) ≤ v2) (hL : Real.log N ≤ L)
    (hL1 : Real.log ((N - 1 : ℕ) : ℝ) ≤ L1)
    (hC : ((N - 2 : ℕ) : ℝ) * (L1 ^ 2 * v2) +
      v * (K2b N R T + 2 * L * K1b N R T + L ^ 2 * K0b N R T) ≤ C)
    (cj : ℕ → ℝ) (hc : ∀ j : ℕ, j < 7 → 0 ≤ cj j ∧ (b + 1 + j) ^ 2 + T ^ 2 ≤ cj j ^ 2)
    (hε : R * ((∏ j ∈ Finset.range 7, cj j) / 1000000 * (v / (N : ℝ) ^ 7 / (a + 7))) ≤ ε)
    {cr ci cn E : ℝ} (hcn0 : 0 ≤ cn) (hcn : cr ^ 2 + ci ^ 2 ≤ cn ^ 2)
    (hE : ε + rA + rB * (b - a) + C * (b - a) ^ 2 ≤ E)
    (h0 : cn * E < Ar * cr + Ai * ci)
    (h1 : cn * E < (Ar + (b - a) * Br) * cr + (Ai + (b - a) * Bi) * ci)
    {Ua Lb U : ℝ}
    (hUa : arg (riemannZeta ((a : ℂ) + T * I) * (starRingEnd ℂ) ((cr : ℂ) + ci * I)) ≤ Ua)
    (hLb : Lb ≤ arg (riemannZeta ((b : ℂ) + T * I) * (starRingEnd ℂ) ((cr : ℂ) + ci * I)))
    (hU : Ua - Lb ≤ U) : TopOK T a b U := by
  set h := b - a with hh
  have hh0 : 0 < h := by rw [hh]; linarith
  set p : ℂ := (a : ℂ) + T * I with hp
  have hpu : ∀ u : ℝ, p + (u : ℂ) * 1 = ((a + u : ℝ) : ℂ) + (T : ℂ) * I := by
    intro u; rw [hp]; push_cast; ring
  have hs1 : ∀ u ∈ Icc (0 : ℝ) h, p + (u : ℂ) * 1 ≠ 1 := by
    intro u _; rw [hpu]; exact ne_one_of_im hT.ne'
  have hreu : ∀ u : ℝ, (p + (u : ℂ) * 1).re = a + u := by intro u; rw [hpu]; simp
  have himu : ∀ u : ℝ, (p + (u : ℂ) * 1).im = T := by intro u; rw [hpu]; simp
  have hnorm : ∀ u ∈ Icc (0 : ℝ) h, ‖p + (u : ℂ) * 1‖ ≤ R := by
    intro u hu
    refine norm_le_of_sq_le hR0 ?_
    rw [hreu, himu]
    have : (a + u) ^ 2 ≤ b ^ 2 := pow_le_pow_left₀ (by linarith [hu.1]) (by linarith [hu.2]) 2
    linarith
  have hC' : ∀ u ∈ Icc (0 : ℝ) h, ‖M2f N (p + (u : ℂ) * 1)‖ ≤ C := by
    intro u hu
    refine le_trans (norm_M2f_le' hN (hnorm u hu) hT ?_ ha.le (by rw [hreu]; linarith [hu.1])
      hv hv2 hL hL1) hC
    have := Complex.abs_im_le_norm (1 - (p + (u : ℂ) * 1))
    rw [sub_im, one_im, himu, zero_sub, abs_neg, abs_of_pos hT] at this
    exact this
  have hε' : ∀ u ∈ Icc (0 : ℝ) h, ‖riemannZeta (p + (u : ℂ) * 1) - Mf N (p + (u : ℂ) * 1)‖ ≤ ε := by
    intro u hu
    refine em_err_le hN cj ha (by rw [hreu]; linarith [hu.1]) (hs1 u hu) (hnorm u hu) ?_ hv hε
    intro j hj
    obtain ⟨hc0, hc1⟩ := hc j hj
    have := norm_shift_le (s := p + (u : ℂ) * 1) (x1 := b) (Y := T) (a := 1 + j) (c := cj j)
      (by positivity) (by rw [hreu]; linarith [hu.1]) (by rw [hreu]; linarith [hu.2])
      (by rw [himu, abs_of_pos hT]) hc0 (by linarith)
    simpa [add_assoc] using this
  have hmod := seg_model hN (p := p) (d := 1) (by simp) hs1 hA hB hC' hε'
  have hrB : 0 ≤ rB := (norm_nonneg _).trans hB
  have hC0 : 0 ≤ C := (norm_nonneg _).trans (hC' 0 ⟨le_rfl, hh0.le⟩)
  have hmod' : ∀ u ∈ Icc (0 : ℝ) h, ‖riemannZeta (p + (u : ℂ) * 1) - (((Ar : ℂ) + Ai * I) +
      (u : ℂ) * 1 * ((Br : ℂ) + Bi * I))‖ ≤ E := by
    intro u hu
    refine (hmod u hu).trans (le_trans ?_ hE)
    have : rB * u ≤ rB * h := mul_le_mul_of_nonneg_left hu.2 hrB
    have : C * u ^ 2 ≤ C * h ^ 2 :=
      mul_le_mul_of_nonneg_left (pow_le_pow_left₀ hu.1 hu.2 2) hC0
    linarith
  have hslit := re_pos_of_model hh0 hmod' (norm_c_le hcn0 hcn)
    (by rw [conjc_re]; exact h0)
    (by
      have e : (((Ar : ℂ) + Ai * I) + (h : ℂ) * 1 * ((Br : ℂ) + Bi * I)) =
          ((Ar + h * Br : ℝ) : ℂ) + ((Ai + h * Bi : ℝ) : ℂ) * I := by
        push_cast; ring
      rw [e, conjc_re]; exact h1)
  refine topOK_of_slit hT.ne' hab.le (fun x hx => ?_) hUa hLb hU
  have := hslit (x - a) ⟨by linarith [hx.1], by rw [hh]; linarith [hx.2]⟩
  rwa [hpu, show a + (x - a) = x by ring] at this

end Fayed.FirstZero.Min
