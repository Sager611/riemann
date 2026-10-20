import Riemann.Deep.DH.ZeroLib
import Riemann.Deep.DH.ZeroNewton

/-!
# Davenport–Heilbronn zero, part 5: κ, ball assembly, error and crude bounds, analyticity

κ_near (an enclosure of κ from square-root bounds), Gem_ball (Gem from its pieces), Ebd_le
(numeric bound on the Euler–Maclaurin error), norm_Gem_le (crude sup bound on Gem),
differentiableAt_Gem.
-/

open Complex Set

namespace Fayed.DH.Zero

open Fayed.DH
open Fayed.FirstZero (nps gc gf EMPhi ball_add ball_mul ball_real_mul norm_gc_le)

/-! ### κ -/

lemma sqrt_lo {x lo : ℝ} (hlo : 0 ≤ lo) (h : lo ^ 2 ≤ x) : lo ≤ √x :=
  calc lo = √(lo ^ 2) := (Real.sqrt_sq hlo).symm
    _ ≤ √x := Real.sqrt_le_sqrt h

lemma sqrt_hi {x hi : ℝ} (hhi : 0 ≤ hi) (h : x ≤ hi ^ 2) : √x ≤ hi :=
  calc √x ≤ √(hi ^ 2) := Real.sqrt_le_sqrt h
    _ = hi := Real.sqrt_sq hhi

lemma κ_near {a1 a2 b1 b2 k0 εk : ℝ} (ha1 : 1 < a1) (ha2 : 1 < a2) (hb2 : 0 ≤ b2) (h5lo : a1 ^ 2 ≤ 5) (h5hi : 5 ≤ a2 ^ 2)
    (hb1 : 2 ≤ b1) (hTlo : b1 ^ 2 ≤ 10 - 2 * a2) (hThi : 10 - 2 * a1 ≤ b2 ^ 2)
    (hlo : (k0 - εk) * (a2 - 1) ≤ b1 - 2) (hhi : b2 - 2 ≤ (k0 + εk) * (a1 - 1)) :
    |κ - k0| ≤ εk := by
  have r1 : a1 ≤ √5 := sqrt_lo (by linarith) h5lo
  have r2 : √5 ≤ a2 := sqrt_hi (by linarith) h5hi
  have t1 : b1 ≤ √(10 - 2 * √5) := sqrt_lo (by linarith) (by linarith)
  have t2 : √(10 - 2 * √5) ≤ b2 := sqrt_hi hb2 (by linarith)
  have hd : 0 < √5 - 1 := by linarith
  unfold κ
  rw [abs_le]
  constructor
  · rw [le_sub_iff_add_le, le_div_iff₀ hd]
    nlinarith
  · rw [sub_le_iff_le_add, div_le_iff₀ hd]
    nlinarith

/-! ### Assembling Gem from balls -/

theorem Gem_ball {s : ℂ} (hs1 : s ≠ 1) {M : ℕ}
    {k0 εk g1 g2 rg p11 p12 r1 p21 p22 r2 p31 p32 r3 p41 p42 r4 e f R : ℝ}
    (hk : |κ - k0| ≤ εk) (hG : ‖Gmain M s - (g1 + g2 * I)‖ ≤ rg)
    (h1 : ‖nps s (5 * M + 1) * Qf s ((M : ℝ) + 1 / 5) - (p11 + p12 * I)‖ ≤ r1)
    (h2 : ‖nps s (5 * M + 2) * Qf s ((M : ℝ) + 2 / 5) - (p21 + p22 * I)‖ ≤ r2)
    (h3 : ‖nps s (5 * M + 3) * Qf s ((M : ℝ) + 3 / 5) - (p31 + p32 * I)‖ ≤ r3)
    (h4 : ‖nps s (5 * M + 4) * Qf s ((M : ℝ) + 4 / 5) - (p41 + p42 * I)‖ ≤ r4)
    (h : rg + r1 + r4 + (|p21 - p31| + |p22 - p32|) * εk + |k0| * (r2 + r3) + εk * (r2 + r3) +
      |g1 + (p11 - p41 + k0 * (p21 - p31)) - e| + |g2 + (p12 - p42 + k0 * (p22 - p32)) - f| ≤ R) :
    ‖Gem M s - (e + f * I)‖ ≤ R := by
  rw [Gem_eval M hs1]
  have h14 := ball_sub h1 h4
  have h23 := ball_sub h2 h3
  have hk23 := ball_real_mul hk h23 (e := k0 * (p21 - p31)) (f := k0 * (p22 - p32))
    (R := (|p21 - p31| + |p22 - p32|) * εk + |k0| * (r2 + r3) + εk * (r2 + r3)) (by
      simp only [sub_self, abs_zero, add_zero]; exact le_refl _)
  have ht := ball_add h14 hk23 (e := p11 - p41 + k0 * (p21 - p31)) (f := p12 - p42 + k0 * (p22 - p32))
    (R := r1 + r4 + ((|p21 - p31| + |p22 - p32|) * εk + |k0| * (r2 + r3) + εk * (r2 + r3)))
    (by simp only [sub_self, abs_zero, zero_add]; linarith)
  exact ball_add hG ht (by linarith)

/-! ### The Euler–Maclaurin error bound, numerically -/

theorem Ebd_le {M : ℕ} (hM : 1 ≤ M) {s : ℂ} {σ1 v5 vM S G K E : ℝ} (hσ : σ1 ≤ s.re)
    (hσ1 : 0 ≤ σ1) (h5 : (5 : ℝ) ^ (-σ1) ≤ v5) (hMv : (M : ℝ) ^ (-(σ1 + 7)) ≤ vM)
    (hS : ‖s‖ ≤ S) (hG : ‖gc s 7‖ ≤ G) (hK : κ ≤ K)
    (h : v5 * S * (2 + 2 * K) * (G / 1000000 * (vM / (σ1 + 7))) ≤ E) : Ebd M s ≤ E := by
  have hMr : (1 : ℝ) ≤ M := by exact_mod_cast hM
  have hk := κ_nonneg
  have a1 : (5 : ℝ) ^ (-s.re) ≤ v5 :=
    (Real.rpow_le_rpow_of_exponent_le (by norm_num) (by linarith)).trans h5
  have a2 : (M : ℝ) ^ (-(s.re + 7)) ≤ vM :=
    (Real.rpow_le_rpow_of_exponent_le hMr (by linarith)).trans hMv
  have a3 : (M : ℝ) ^ (-(s.re + 7)) / (s.re + 7) ≤ vM / (σ1 + 7) := by
    have hp : 0 ≤ (M : ℝ) ^ (-(s.re + 7)) := Real.rpow_nonneg (by linarith) _
    calc (M : ℝ) ^ (-(s.re + 7)) / (s.re + 7) ≤ (M : ℝ) ^ (-(s.re + 7)) / (σ1 + 7) :=
          div_le_div_of_nonneg_left hp (by linarith) (by linarith)
      _ ≤ vM / (σ1 + 7) := div_le_div_of_nonneg_right a2 (by linarith)
  unfold Ebd
  refine le_trans ?_ h
  have p1 : 0 ≤ (5 : ℝ) ^ (-s.re) := Real.rpow_nonneg (by norm_num) _
  have p2 : 0 ≤ (M : ℝ) ^ (-(s.re + 7)) / (s.re + 7) :=
    div_nonneg (Real.rpow_nonneg (by linarith) _) (by linarith)
  have p3 : 0 ≤ ‖gc s 7‖ / 1000000 := by positivity
  have hS0 : 0 ≤ S := (norm_nonneg _).trans hS
  have hv5 : 0 ≤ v5 := p1.trans a1
  have hG0 : 0 ≤ G := (norm_nonneg _).trans hG
  have A : (5 : ℝ) ^ (-s.re) * ‖s‖ * (2 + 2 * κ) ≤ v5 * S * (2 + 2 * K) :=
    mul_le_mul (mul_le_mul a1 hS (norm_nonneg _) hv5) (by linarith) (by positivity) (by positivity)
  have B : ‖gc s 7‖ / 1000000 * ((M : ℝ) ^ (-(s.re + 7)) / (s.re + 7)) ≤
      G / 1000000 * (vM / (σ1 + 7)) :=
    mul_le_mul (div_le_div_of_nonneg_right hG (by norm_num)) a3 p2 (by positivity)
  exact mul_le_mul A B (by positivity) (by have : 0 ≤ K := hk.trans hK; positivity)

/-- ‖gc s n‖ ≤ ∏ (S + 1 + j) for ‖s‖ ≤ S. -/
lemma norm_gc_le_S {s : ℂ} {S : ℝ} (hS : ‖s‖ ≤ S) (n : ℕ) :
    ‖gc s n‖ ≤ ∏ j ∈ Finset.range n, (S + 1 + j) := by
  refine norm_gc_le s _ n fun j _ => ?_
  calc ‖s + 1 + j‖ ≤ ‖s‖ + ‖(1 : ℂ)‖ + ‖(j : ℂ)‖ := norm_add₃_le
    _ ≤ S + 1 + j := by simp only [norm_one, Complex.norm_natCast]; linarith

/-! ### Crude bounds -/

lemma norm_nps_le_one {s : ℂ} (hs : 0 ≤ s.re) {n : ℕ} (hn : 0 < n) : ‖nps s n‖ ≤ 1 := by
  rw [nps, norm_natCast_cpow_of_pos hn, neg_re]
  exact Real.rpow_le_one_of_one_le_of_nonpos (by exact_mod_cast hn) (by linarith)

lemma norm_W_le {s : ℂ} (hs : 0 ≤ s.re) (m : ℕ) : ‖W s m‖ ≤ 2 + 2 * κ := by
  have hk := κ_nonneg
  have h := fun j : ℕ => norm_nps_le_one hs (n := 5 * m + j + 1) (by omega)
  unfold W
  calc _ ≤ ‖nps s (5 * m + 1) - nps s (5 * m + 4)‖ +
        ‖(κ : ℂ) * (nps s (5 * m + 2) - nps s (5 * m + 3))‖ := norm_add_le _ _
    _ ≤ (‖nps s (5 * m + 1)‖ + ‖nps s (5 * m + 4)‖) +
        κ * (‖nps s (5 * m + 2)‖ + ‖nps s (5 * m + 3)‖) := by
        gcongr
        · exact norm_sub_le _ _
        · rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hk]
          gcongr; exact norm_sub_le _ _
    _ ≤ (1 + 1) + κ * (1 + 1) := by
        have := h 0; have := h 1; have := h 2; have := h 3
        gcongr
    _ = 2 + 2 * κ := by ring

lemma norm_Gmain_le {s : ℂ} (hs : 0 ≤ s.re) (M : ℕ) : ‖Gmain M s‖ ≤ M * (2 + 2 * κ) := by
  unfold Gmain
  calc _ ≤ ∑ m ∈ Finset.range M, ‖W s m‖ := norm_sum_le _ _
    _ ≤ ∑ m ∈ Finset.range M, (2 + 2 * κ) := Finset.sum_le_sum fun m _ => norm_W_le hs m
    _ = M * (2 + 2 * κ) := by simp; ring

lemma norm_Qf_le {s : ℂ} {x xl X S T1 G2 G4 G6 : ℝ} (hxl : 0 < xl) (hx : xl ≤ x) (hxX : x ≤ X)
    (hS : ‖s‖ ≤ S) (hT1 : 0 < T1) (hT : T1 ≤ s.im)
    (h2 : ‖gc s 2‖ ≤ G2) (h4 : ‖gc s 4‖ ≤ G4) (h6 : ‖gc s 6‖ ≤ G6) :
    ‖Qf s x‖ ≤ X / T1 + 1 / 2 +
      S * (1 / 12 / xl + 1 / 720 * G2 / xl ^ 3 + 1 / 30240 * G4 / xl ^ 5 +
        1 / 1209600 * G6 / xl ^ 7) := by
  have hx0 : 0 < x := lt_of_lt_of_le hxl hx
  have hn : ∀ k : ℕ, ‖(x : ℂ) ^ k‖ = x ^ k := fun k => by
    rw [norm_pow, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hx0]
  have hinv : ∀ k : ℕ, 1 / x ^ k ≤ 1 / xl ^ k := fun k =>
    one_div_le_one_div_of_le (pow_pos hxl k) (pow_le_pow_left₀ hxl.le hx k)
  have hs1 : T1 ≤ ‖s - 1‖ := by
    calc T1 ≤ |(s - 1).im| := by simp only [sub_im, one_im, sub_zero]; exact hT.trans (le_abs_self _)
      _ ≤ ‖s - 1‖ := Complex.abs_im_le_norm _
  have hS0 : 0 ≤ S := (norm_nonneg _).trans hS
  have g0 : ‖gc s 0‖ = 1 := by simp [Fayed.FirstZero.gc]
  have e1 : ‖(x : ℂ) / (s - 1)‖ ≤ X / T1 := by
    rw [norm_div, Complex.norm_real, Real.norm_eq_abs, abs_of_pos hx0]
    exact div_le_div₀ (by linarith) hxX hT1 hs1
  have e2 : ‖EMq s x‖ ≤ 1 / 12 / xl + 1 / 720 * G2 / xl ^ 3 + 1 / 30240 * G4 / xl ^ 5 +
      1 / 1209600 * G6 / xl ^ 7 := by
    unfold EMq
    have t : ∀ (c : ℝ) (k j : ℕ) (Gj : ℝ), 0 ≤ c → ‖gc s j‖ ≤ Gj →
        ‖(c : ℂ) * gc s j / (x : ℂ) ^ k‖ ≤ c * Gj / xl ^ k := by
      intro c k j Gj hc hG
      rw [norm_div, norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hc, hn]
      have hG0 : 0 ≤ Gj := (norm_nonneg _).trans hG
      calc c * ‖gc s j‖ / x ^ k = c * ‖gc s j‖ * (1 / x ^ k) := by ring
        _ ≤ c * Gj * (1 / xl ^ k) := by gcongr
        _ = c * Gj / xl ^ k := by ring
    have t0 := t (1 / 12) 1 0 1 (by norm_num) (by rw [g0])
    have t2 := t (1 / 720) 3 2 G2 (by norm_num) h2
    have t4 := t (1 / 30240) 5 4 G4 (by norm_num) h4
    have t6 := t (1 / 1209600) 7 6 G6 (by norm_num) h6
    push_cast at t0 t2 t4 t6
    rw [pow_one, pow_one] at t0
    calc _ ≤ ‖1 / 12 * gc s 0 / (x : ℂ)‖ + ‖1 / 720 * gc s 2 / (x : ℂ) ^ 3‖ +
          ‖1 / 30240 * gc s 4 / (x : ℂ) ^ 5‖ + ‖1 / 1209600 * gc s 6 / (x : ℂ) ^ 7‖ := by
          refine (norm_sub_le _ _).trans ?_
          gcongr
          refine (norm_add_le _ _).trans ?_
          gcongr
          exact norm_sub_le _ _
      _ ≤ _ := by
          have : 1 / 12 * 1 / xl = 1 / 12 / xl := by ring
          linarith
  unfold Qf
  calc _ ≤ ‖(x : ℂ) / (s - 1)‖ + ‖(1 / 2 : ℂ)‖ + ‖s * EMq s x‖ := norm_add₃_le
    _ ≤ X / T1 + 1 / 2 + S * (1 / 12 / xl + 1 / 720 * G2 / xl ^ 3 + 1 / 30240 * G4 / xl ^ 5 +
        1 / 1209600 * G6 / xl ^ 7) := by
        rw [norm_mul]
        have : ‖(1 / 2 : ℂ)‖ = 1 / 2 := by norm_num
        rw [this]
        gcongr

theorem norm_Gem_le {M : ℕ} (hM : 1 ≤ M) {s : ℂ} (hs : 0 ≤ s.re) (hs1 : s ≠ 1) {Q : ℝ}
    (hQ : ∀ j : ℕ, 1 ≤ j → j ≤ 4 → ‖Qf s ((M : ℝ) + j / 5)‖ ≤ Q) :
    ‖Gem M s‖ ≤ (M + Q) * (2 + 2 * κ) := by
  have hk := κ_nonneg
  rw [Gem_eval M hs1]
  have hp : ∀ j : ℕ, 1 ≤ j → j ≤ 4 →
      ‖nps s (5 * M + j) * Qf s ((M : ℝ) + j / 5)‖ ≤ Q := fun j h1 h4 => by
    rw [norm_mul]
    calc _ ≤ 1 * ‖Qf s ((M : ℝ) + j / 5)‖ :=
          mul_le_mul_of_nonneg_right (norm_nps_le_one hs (by omega)) (norm_nonneg _)
      _ ≤ Q := by rw [one_mul]; exact hQ j h1 h4
  have q1 := hp 1 le_rfl (by norm_num)
  have q2 := hp 2 (by norm_num) (by norm_num)
  have q3 := hp 3 (by norm_num) (by norm_num)
  have q4 := hp 4 (by norm_num) le_rfl
  push_cast at q1 q2 q3 q4
  have hG := norm_Gmain_le hs M
  calc _ ≤ ‖Gmain M s‖ + (‖nps s (5 * M + 1) * Qf s ((M : ℝ) + 1 / 5)‖ +
        ‖nps s (5 * M + 4) * Qf s ((M : ℝ) + 4 / 5)‖ +
        κ * (‖nps s (5 * M + 2) * Qf s ((M : ℝ) + 2 / 5)‖ +
          ‖nps s (5 * M + 3) * Qf s ((M : ℝ) + 3 / 5)‖)) := by
        refine (norm_add_le _ _).trans (add_le_add le_rfl ?_)
        refine (norm_add_le _ _).trans (add_le_add (norm_sub_le _ _) ?_)
        rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hk]
        gcongr; exact norm_sub_le _ _
    _ ≤ M * (2 + 2 * κ) + (Q + Q + κ * (Q + Q)) := by gcongr
    _ = (M + Q) * (2 + 2 * κ) := by ring

/-! ### Analyticity of Gem -/

lemma differentiable_gc (j : ℕ) : Differentiable ℂ (fun s => gc s j) := by
  induction j with
  | zero => simp only [Fayed.FirstZero.gc]; exact differentiable_const 1
  | succ j ih =>
    simp only [Fayed.FirstZero.gc]
    exact ih.mul ((differentiable_id.neg.sub (differentiable_const _)).sub (differentiable_const _))

lemma differentiable_gf (j : ℕ) {x : ℝ} (hx : 0 < x) : Differentiable ℂ (fun s => gf s j x) := by
  have hx' : (x : ℂ) ≠ 0 := ofReal_ne_zero.mpr hx.ne'
  unfold gf
  exact (differentiable_gc j).mul fun s =>
    (((differentiable_id.neg.sub (differentiable_const _)).sub
      (differentiable_const _)) s).const_cpow (Or.inl hx')

lemma differentiable_EMPhi {x : ℝ} (hx : 0 < x) : Differentiable ℂ (fun s => EMPhi s x) := by
  unfold EMPhi
  exact ((((differentiable_gf 0 hx).const_mul _).sub ((differentiable_gf 2 hx).const_mul _)).add
    ((differentiable_gf 4 hx).const_mul _)).sub ((differentiable_gf 6 hx).const_mul _)

lemma differentiableAt_Gem (M : ℕ) {s : ℂ} (hs1 : s ≠ 1) : DifferentiableAt ℂ (Gem M) s := by
  have hpos : ∀ j : ℝ, 0 < j → 0 < (M : ℝ) + j := fun j hj => by positivity
  have r : ∀ j : ℝ, 0 < j → DifferentiableAt ℂ (Sj M j) s := fun j hj =>
    (differentiableAt_Af (hpos j hj) hs1).add
      (differentiableAt_id.mul (differentiable_EMPhi (hpos j hj) s))
  unfold Gem
  refine (differentiable_Gmain M s).add (DifferentiableAt.mul ?_ ?_)
  · exact (differentiableAt_id.neg).const_cpow (Or.inl (by norm_num))
  · exact ((r _ (by norm_num)).sub (r _ (by norm_num))).add
      (((r _ (by norm_num)).sub (r _ (by norm_num))).const_mul _)

end Fayed.DH.Zero
