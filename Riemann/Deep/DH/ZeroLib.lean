import Riemann.Deep.DH.ZeroEM
import Riemann.Deep.FirstZero.MinLib

/-!
# Davenport–Heilbronn zero, part 4: numeric lemmas for the generated certificate

Logarithms by the series of log(1 - x), real powers n^{-σ} via exp((1-σ) log n / 2)^2 / n,
e^{iθ} by a degree-29 Taylor polynomial, balls for the Dirichlet blocks W s m and the partial sums
of Gmain. Every hypothesis named h is a rational inequality for norm_num.
-/

open Complex Set

namespace Fayed.DH.Zero

open Fayed.DH
open Fayed.FirstZero (nps gc ball_add ball_mul ball_real_mul sum_exp_I norm_cexp_I_sub_cexp_I)

/-! ### Logarithms -/

lemma log_step {p m : ℕ} {x cm δm c ε : ℝ} (n : ℕ) (hm0 : 0 < m) (hx : |x| < 1)
    (hpx : (p : ℝ) = m * (1 - x)) (hlm : |Real.log m - cm| ≤ δm)
    (h : δm + |x| ^ (n + 1) / (1 - |x|) +
      |cm - (∑ i ∈ Finset.range n, x ^ (i + 1) / (i + 1)) - c| ≤ ε) :
    |Real.log p - c| ≤ ε := by
  have h1x : 0 < 1 - x := by have := le_abs_self x; linarith
  have hlog : Real.log p = Real.log m + Real.log (1 - x) := by
    rw [hpx, Real.log_mul (by positivity) h1x.ne']
  have hS := Real.abs_log_sub_add_sum_range_le hx n
  set S := ∑ i ∈ Finset.range n, x ^ (i + 1) / (i + 1)
  rw [hlog]
  calc |Real.log m + Real.log (1 - x) - c|
      = |(Real.log m - cm) + (S + Real.log (1 - x)) + (cm - S - c)| := by ring_nf
    _ ≤ |Real.log m - cm| + |S + Real.log (1 - x)| + |cm - S - c| := abs_add_three _ _ _
    _ ≤ ε := by linarith

lemma log_mul_nat {a b m : ℕ} (hab : a * b = m) (ha : 0 < a) (hb : 0 < b) {c1 e1 c2 e2 c ε : ℝ}
    (h1 : |Real.log a - c1| ≤ e1) (h2 : |Real.log b - c2| ≤ e2) (h : e1 + e2 + |c1 + c2 - c| ≤ ε) :
    |Real.log m - c| ≤ ε := by
  subst hab
  have := Fayed.FirstZero.Min.log_near_mul ha hb h1 h2
  calc |Real.log ((a * b : ℕ) : ℝ) - c|
      = |(Real.log ((a * b : ℕ) : ℝ) - (c1 + c2)) + (c1 + c2 - c)| := by ring_nf
    _ ≤ |Real.log ((a * b : ℕ) : ℝ) - (c1 + c2)| + |c1 + c2 - c| := abs_add_le _ _
    _ ≤ ε := by linarith

/-- log 2, log 3, log 5 from log(8/9), log(15/16), log(24/25). -/
lemma log235 {c2 c3 c5 ε2 ε3 ε5 SA SB SC : ℝ}
    (hA : SA = ∑ i ∈ Finset.range 24, (1 / 9 : ℝ) ^ (i + 1) / (i + 1))
    (hB : SB = ∑ i ∈ Finset.range 20, (1 / 16 : ℝ) ^ (i + 1) / (i + 1))
    (hC : SC = ∑ i ∈ Finset.range 18, (1 / 25 : ℝ) ^ (i + 1) / (i + 1))
    (h2 : 3 * ((1 / 9 : ℝ) ^ 25 / (1 - 1 / 9)) + 4 * ((1 / 16 : ℝ) ^ 21 / (1 - 1 / 16)) +
      2 * ((1 / 25 : ℝ) ^ 19 / (1 - 1 / 25)) + |3 * SA + 4 * SB + 2 * SC - c2| ≤ ε2)
    (h3 : |(3 * c2 + SA) / 2 - c3| + 3 / 2 * ε2 + 1 / 2 * ((1 / 9 : ℝ) ^ 25 / (1 - 1 / 9)) ≤ ε3)
    (h5 : |4 * c2 - c3 - SB - c5| + 4 * ε2 + ε3 + (1 / 16 : ℝ) ^ 21 / (1 - 1 / 16) ≤ ε5) :
    |Real.log 2 - c2| ≤ ε2 ∧ |Real.log 3 - c3| ≤ ε3 ∧ |Real.log 5 - c5| ≤ ε5 := by
  have hA' := Real.abs_log_sub_add_sum_range_le (x := 1 / 9) (by norm_num [abs_of_pos]) 24
  have hB' := Real.abs_log_sub_add_sum_range_le (x := 1 / 16) (by norm_num [abs_of_pos]) 20
  have hC' := Real.abs_log_sub_add_sum_range_le (x := 1 / 25) (by norm_num [abs_of_pos]) 18
  rw [← hA] at hA'; rw [← hB] at hB'; rw [← hC] at hC'
  have eA : Real.log (1 - 1 / 9) = 3 * Real.log 2 - 2 * Real.log 3 := by
    rw [show (1 : ℝ) - 1 / 9 = 2 ^ 3 / 3 ^ 2 by norm_num, Real.log_div (by norm_num) (by norm_num),
      Real.log_pow, Real.log_pow]; push_cast; ring
  have eB : Real.log (1 - 1 / 16) = Real.log 3 + Real.log 5 - 4 * Real.log 2 := by
    rw [show (1 : ℝ) - 1 / 16 = 3 * 5 / 2 ^ 4 by norm_num, Real.log_div (by norm_num) (by norm_num),
      Real.log_mul (by norm_num) (by norm_num), Real.log_pow]; push_cast; ring
  have eC : Real.log (1 - 1 / 25) = 3 * Real.log 2 + Real.log 3 - 2 * Real.log 5 := by
    rw [show (1 : ℝ) - 1 / 25 = 2 ^ 3 * 3 / 5 ^ 2 by norm_num, Real.log_div (by norm_num)
      (by norm_num), Real.log_mul (by norm_num) (by norm_num), Real.log_pow, Real.log_pow]
    push_cast; ring
  rw [eA, abs_of_pos (by norm_num : (0:ℝ) < 1 / 9)] at hA'
  rw [eB, abs_of_pos (by norm_num : (0:ℝ) < 1 / 16)] at hB'
  rw [eC, abs_of_pos (by norm_num : (0:ℝ) < 1 / 25)] at hC'
  set TA := (1 / 9 : ℝ) ^ (24 + 1) / (1 - 1 / 9)
  set TB := (1 / 16 : ℝ) ^ (20 + 1) / (1 - 1 / 16)
  set TC := (1 / 25 : ℝ) ^ (18 + 1) / (1 - 1 / 25)
  clear_value TA TB TC
  rw [abs_le] at hA' hB' hC'
  obtain ⟨hA1, hA2⟩ := hA'
  obtain ⟨hB1, hB2⟩ := hB'
  obtain ⟨hC1, hC2⟩ := hC'
  have k2 : |Real.log 2 - c2| ≤ ε2 := by
    have e : Real.log 2 - c2 = -3 * (SA + (3 * Real.log 2 - 2 * Real.log 3))
        - 4 * (SB + (Real.log 3 + Real.log 5 - 4 * Real.log 2))
        - 2 * (SC + (3 * Real.log 2 + Real.log 3 - 2 * Real.log 5)) + (3 * SA + 4 * SB + 2 * SC - c2) := by
      ring
    rw [e]
    have := abs_le.mp (le_refl |3 * SA + 4 * SB + 2 * SC - c2|)
    rw [abs_le]; constructor <;> nlinarith [abs_nonneg (3 * SA + 4 * SB + 2 * SC - c2),
      le_abs_self (3 * SA + 4 * SB + 2 * SC - c2), neg_abs_le (3 * SA + 4 * SB + 2 * SC - c2)]
  have k3 : |Real.log 3 - c3| ≤ ε3 := by
    have e : Real.log 3 - c3 = (3 / 2) * (Real.log 2 - c2) - (1 / 2) * (SA + (3 * Real.log 2 - 2 * Real.log 3))
        + ((3 * c2 + SA) / 2 - c3) := by ring
    rw [e]; rw [abs_le] at k2 ⊢
    obtain ⟨k21, k22⟩ := k2
    constructor <;> linarith [le_abs_self ((3 * c2 + SA) / 2 - c3), neg_abs_le ((3 * c2 + SA) / 2 - c3)]
  have k5 : |Real.log 5 - c5| ≤ ε5 := by
    have e : Real.log 5 - c5 = 4 * (Real.log 2 - c2) - (Real.log 3 - c3)
        + (SB + (Real.log 3 + Real.log 5 - 4 * Real.log 2)) + (4 * c2 - c3 - SB - c5) := by ring
    rw [e]; rw [abs_le] at k2 k3 ⊢
    obtain ⟨k21, k22⟩ := k2
    obtain ⟨k31, k32⟩ := k3
    constructor <;> linarith [le_abs_self (4 * c2 - c3 - SB - c5), neg_abs_le (4 * c2 - c3 - SB - c5)]
  exact ⟨k2, k3, k5⟩

/-! ### Real powers -/

lemma exp_ball {y0 E ε : ℝ} (n : ℕ) (hn : 0 < n) (hy : |y0| ≤ 1)
    (h : |(∑ m ∈ Finset.range n, y0 ^ m / (m.factorial : ℝ)) - E| +
      |y0| ^ n * ((n.succ : ℝ) / ((n.factorial : ℝ) * n)) ≤ ε) :
    |Real.exp y0 - E| ≤ ε := by
  have hb := Real.exp_bound hy hn
  calc |Real.exp y0 - E|
      = |(Real.exp y0 - ∑ m ∈ Finset.range n, y0 ^ m / (m.factorial : ℝ)) +
          ((∑ m ∈ Finset.range n, y0 ^ m / (m.factorial : ℝ)) - E)| := by ring_nf
    _ ≤ _ := abs_add_le _ _
    _ ≤ ε := by linarith

lemma rpow_neg_ball {n : ℕ} (hn : 0 < n) {σ ℓ δ y0 E εE η e1 c ε : ℝ}
    (hl : |Real.log n - ℓ| ≤ δ) (hE : |Real.exp y0 - E| ≤ εE)
    (hη : |1 - σ| * δ / 2 + |(1 - σ) * ℓ / 2 - y0| ≤ η) (hη1 : η ≤ 1)
    (he1 : (|E| + εE) * (2 * η) + εE ≤ e1)
    (h : (e1 * (2 * |E| + e1) + |E ^ 2 - n * c|) / n ≤ ε) :
    |(n : ℝ) ^ (-σ) - c| ≤ ε := by
  have hn' : (0 : ℝ) < n := by exact_mod_cast hn
  set y := (1 - σ) * Real.log n / 2 with hy
  have hpow : (n : ℝ) ^ (-σ) = Real.exp y ^ 2 / n := by
    rw [Real.rpow_def_of_pos hn', ← Real.exp_nat_mul, hy]
    rw [eq_div_iff hn'.ne', show (n : ℝ) = Real.exp (Real.log n) from (Real.exp_log hn').symm]
    rw [← Real.exp_add, Real.exp_log hn']
    congr 1; push_cast; ring
  have hyy : |y - y0| ≤ η := by
    calc |y - y0| = |(1 - σ) * (Real.log n - ℓ) / 2 + ((1 - σ) * ℓ / 2 - y0)| := by
          rw [hy]; ring_nf
      _ ≤ |(1 - σ) * (Real.log n - ℓ) / 2| + |(1 - σ) * ℓ / 2 - y0| := abs_add_le _ _
      _ ≤ |1 - σ| * δ / 2 + |(1 - σ) * ℓ / 2 - y0| := by
          gcongr
          rw [abs_div, abs_mul, abs_two]
          gcongr
      _ ≤ η := hη
  have hey : |Real.exp y - E| ≤ e1 := by
    have h1 := Real.abs_exp_sub_one_le (x := y - y0) (hyy.trans hη1)
    have hE0 : |Real.exp y0| ≤ |E| + εE := by
      have := abs_sub_abs_le_abs_sub (Real.exp y0) E; linarith
    have e : Real.exp y - E = Real.exp y0 * (Real.exp (y - y0) - 1) + (Real.exp y0 - E) := by
      rw [mul_sub, ← Real.exp_add]; ring_nf
    rw [e]
    calc _ ≤ |Real.exp y0 * (Real.exp (y - y0) - 1)| + |Real.exp y0 - E| := abs_add_le _ _
      _ = |Real.exp y0| * |Real.exp (y - y0) - 1| + |Real.exp y0 - E| := by rw [abs_mul]
      _ ≤ (|E| + εE) * (2 * η) + εE := by
          gcongr
          · exact (abs_nonneg _).trans hE0
          · exact h1.trans (by linarith)
      _ ≤ e1 := he1
  have he10 : 0 ≤ e1 := (abs_nonneg _).trans hey
  have hsq : |Real.exp y ^ 2 - E ^ 2| ≤ e1 * (2 * |E| + e1) := by
    have e : Real.exp y ^ 2 - E ^ 2 = (Real.exp y - E) * ((Real.exp y - E) + 2 * E) := by ring
    rw [e, abs_mul]
    gcongr
    calc |Real.exp y - E + 2 * E| ≤ |Real.exp y - E| + |2 * E| := abs_add_le _ _
      _ ≤ e1 + 2 * |E| := by rw [abs_mul, abs_two]; linarith
      _ = 2 * |E| + e1 := by ring
  rw [hpow]
  have e : Real.exp y ^ 2 / n - c = ((Real.exp y ^ 2 - E ^ 2) + (E ^ 2 - n * c)) / n := by
    field_simp
    ring
  rw [e, abs_div, abs_of_pos hn']
  refine le_trans ?_ h
  gcongr
  exact (abs_add_le _ _).trans (by linarith)

/-! ### e^{iθ} -/

lemma theta_ball {n : ℕ} {t ℓ δ q η : ℝ} (m : ℤ) (hl : |Real.log n - ℓ| ≤ δ)
    (h : |t| * δ + |-(t * ℓ) - q - 2 * (314159265358979323846 / 10 ^ 20) * m| +
      2 * |(m : ℝ)| * (1 / 10 ^ 20) ≤ η) :
    |(-(t * Real.log n)) - (q + 2 * Real.pi * m)| ≤ η := by
  have hp1 := Real.pi_gt_d20
  have hp2 := Real.pi_lt_d20
  have hpi : |Real.pi - 314159265358979323846 / 10 ^ 20| ≤ 1 / 10 ^ 20 := by
    rw [abs_le]; constructor <;> norm_num at hp1 hp2 ⊢ <;> linarith
  have e : (-(t * Real.log n)) - (q + 2 * Real.pi * m) = -t * (Real.log n - ℓ) +
      (-(t * ℓ) - q - 2 * (314159265358979323846 / 10 ^ 20) * m) +
      (-2 * m) * (Real.pi - 314159265358979323846 / 10 ^ 20) := by ring
  rw [e]
  refine (abs_add_three _ _ _).trans (le_trans ?_ h)
  rw [abs_mul, abs_mul, abs_neg, abs_mul, abs_neg, abs_two]
  gcongr

lemma cis_ball30 {θ q η : ℝ} (m : ℤ) (hθ : |θ - (q + 2 * Real.pi * m)| ≤ η) (hq : |q| ≤ 16 / 5)
    {c s ε : ℝ}
    (h : |(∑ k ∈ Finset.range 15, (-1) ^ k * q ^ (2 * k) / ((2 * k).factorial : ℝ)) - c| +
      |(∑ k ∈ Finset.range 15, (-1) ^ k * q ^ (2 * k + 1) / ((2 * k + 1).factorial : ℝ)) - s| +
      (16 / 5) ^ 30 / (Nat.factorial 30) * 2 + η ≤ ε) :
    ‖cexp (θ * I) - (c + s * I)‖ ≤ ε := by
  have hper : cexp (θ * I) = cexp (((θ - 2 * Real.pi * m : ℝ) : ℂ) * I) := by
    rw [show (((θ - 2 * Real.pi * m : ℝ) : ℂ) * I) = (θ : ℂ) * I - (m : ℂ) * (2 * Real.pi * I) by
      push_cast; ring, Complex.exp_sub, Complex.exp_int_mul_two_pi_mul_I, div_one]
  have h1 : ‖cexp (θ * I) - cexp (q * I)‖ ≤ η := by
    rw [hper]
    refine (norm_cexp_I_sub_cexp_I _ _).trans ?_
    rw [show θ - 2 * Real.pi * m - q = θ - (q + 2 * Real.pi * m) by ring]
    exact hθ
  have hqn : ‖(q : ℂ) * I‖ = |q| := by simp
  have h2 : ‖cexp ((q : ℂ) * I) - ∑ m ∈ Finset.range 30, ((q : ℂ) * I) ^ m / (m.factorial : ℂ)‖ ≤
      (16 / 5) ^ 30 / (Nat.factorial 30) * 2 := by
    refine (Complex.exp_bound' (n := 30) (by rw [hqn]; push_cast; linarith)).trans ?_
    rw [hqn]
    gcongr
  rw [show (30 : ℕ) = 2 * 15 by norm_num, sum_exp_I] at h2
  set C := ∑ k ∈ Finset.range 15, (-1) ^ k * q ^ (2 * k) / ((2 * k).factorial : ℝ)
  set S := ∑ k ∈ Finset.range 15, (-1) ^ k * q ^ (2 * k + 1) / ((2 * k + 1).factorial : ℝ)
  have h3 : ‖((C : ℂ) + (S : ℂ) * I) - (c + s * I)‖ ≤ |C - c| + |S - s| := by
    refine (norm_le_abs_re_add_abs_im _).trans (le_of_eq ?_); simp
  calc ‖cexp (θ * I) - (c + s * I)‖ ≤ ‖cexp (θ * I) - cexp (q * I)‖ +
        ‖cexp ((q : ℂ) * I) - ((C : ℂ) + (S : ℂ) * I)‖ +
        ‖((C : ℂ) + (S : ℂ) * I) - (c + s * I)‖ := by
        refine (norm_sub_le_norm_sub_add_norm_sub (cexp (θ * I)) (cexp ((q : ℂ) * I))
          (c + s * I)).trans ?_
        have := norm_sub_le_norm_sub_add_norm_sub (cexp ((q : ℂ) * I)) ((C : ℂ) + (S : ℂ) * I)
          (c + s * I)
        linarith
    _ ≤ η + (16 / 5) ^ 30 / (Nat.factorial 30) * 2 + (|C - c| + |S - s|) := by gcongr
    _ ≤ ε := by linarith

/-! ### Balls for blocks and partial sums -/

lemma nps_one (s : ℂ) : ‖nps s 1 - ((1 : ℝ) + (0 : ℝ) * I)‖ ≤ 0 := by
  simp [nps]

lemma nps_mul_ball {s : ℂ} {a b n : ℕ} (hab : a * b = n) {a1 b1 r1 a2 b2 r2 e f R : ℝ}
    (h1 : ‖nps s a - (a1 + b1 * I)‖ ≤ r1) (h2 : ‖nps s b - (a2 + b2 * I)‖ ≤ r2)
    (h : (|a1| + |b1|) * r2 + (|a2| + |b2|) * r1 + r1 * r2 + |a1 * a2 - b1 * b2 - e| +
      |a1 * b2 + b1 * a2 - f| ≤ R) :
    ‖nps s n - (e + f * I)‖ ≤ R := by
  rw [Fayed.FirstZero.Min.nps_mul' a b hab]
  exact ball_mul h1 h2 h

lemma ball_sub {z w : ℂ} {a b r c d q : ℝ} (hz : ‖z - (a + b * I)‖ ≤ r)
    (hw : ‖w - (c + d * I)‖ ≤ q) : ‖(z - w) - (((a - c : ℝ) : ℂ) + ((b - d : ℝ) : ℂ) * I)‖ ≤ r + q := by
  have e : (z - w) - (((a - c : ℝ) : ℂ) + ((b - d : ℝ) : ℂ) * I) =
      (z - (a + b * I)) - (w - (c + d * I)) := by push_cast; ring
  rw [e]; exact (norm_sub_le _ _).trans (add_le_add hz hw)

lemma W_ball {s : ℂ} {m : ℕ} {k0 εk a1 b1 r1 a2 b2 r2 a3 b3 r3 a4 b4 r4 e f R : ℝ}
    (hk : |κ - k0| ≤ εk)
    (h1 : ‖nps s (5 * m + 1) - (a1 + b1 * I)‖ ≤ r1) (h2 : ‖nps s (5 * m + 2) - (a2 + b2 * I)‖ ≤ r2)
    (h3 : ‖nps s (5 * m + 3) - (a3 + b3 * I)‖ ≤ r3) (h4 : ‖nps s (5 * m + 4) - (a4 + b4 * I)‖ ≤ r4)
    (h : r1 + r4 + (|a2 - a3| + |b2 - b3|) * εk + |k0| * (r2 + r3) + εk * (r2 + r3) +
      |a1 - a4 + k0 * (a2 - a3) - e| + |b1 - b4 + k0 * (b2 - b3) - f| ≤ R) :
    ‖W s m - (e + f * I)‖ ≤ R := by
  have h14 := ball_sub h1 h4
  have h23 := ball_sub h2 h3
  have hk23 := ball_real_mul hk h23 (e := k0 * (a2 - a3)) (f := k0 * (b2 - b3))
    (R := (|a2 - a3| + |b2 - b3|) * εk + |k0| * (r2 + r3) + εk * (r2 + r3)) (by
      simp only [sub_self, abs_zero, add_zero]; exact le_refl _)
  have := ball_add h14 hk23 (e := e) (f := f) (R := R) (by linarith)
  simpa [W] using this

lemma Gmain_succ_ball {s : ℂ} {k : ℕ} {a b r c d q e f R : ℝ}
    (hS : ‖Gmain k s - (a + b * I)‖ ≤ r) (hW : ‖W s k - (c + d * I)‖ ≤ q)
    (h : |a + c - e| + |b + d - f| + r + q ≤ R) : ‖Gmain (k + 1) s - (e + f * I)‖ ≤ R := by
  rw [Gmain, Finset.sum_range_succ]
  exact ball_add hS hW h

lemma Gmain_zero_ball (s : ℂ) : ‖Gmain 0 s - ((0 : ℝ) + (0 : ℝ) * I)‖ ≤ 0 := by simp [Gmain]

/-- A ball with a rounded center from an exact value. -/
lemma ball_of_eq {z : ℂ} {A B a b ε : ℝ} (hz : z = A + B * I) (h : |A - a| + |B - b| ≤ ε) :
    ‖z - (a + b * I)‖ ≤ ε := by
  rw [hz]
  refine (norm_le_abs_re_add_abs_im _).trans (le_trans (le_of_eq ?_) h)
  simp

end Fayed.DH.Zero
