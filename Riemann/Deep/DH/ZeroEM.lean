import Riemann.Deep.DH.ZeroRepr

/-!
# Davenport–Heilbronn zero, part 2: order-8 Euler–Maclaurin for DH on Re s > 0

DH_sub_Gem_le: for M ≥ 1, Re s > 0, s ≠ 1,

  ‖DH s − Gem M s‖ ≤ 5^{-σ} ‖s‖ (2 + 2κ) ‖(s+1)⋯(s+7)‖ 10⁻⁶ M^{-σ-7}/(σ+7),

Gem M s = Σ_{m<M} W s m + Σ_j c_j (5M+j)^{-s} Qf(s, M + j/5) (Gem_eval), with Qf an explicit
rational function of s. The per-class sawtooth integrals are expanded with the order-8 machinery of
Riemann/Deep/FirstZero/EM.lean on the shifted pieces [M+a+i, M+a+i+1].
-/

open Complex Set Filter Topology MeasureTheory Interval

namespace Fayed.DH.Zero

open Fayed.DH
open Fayed.FirstZero (nps gc gf BP EMPhi interval_EM8 remainder_piece)

lemma saw_piece_BP (a : ℝ) (k : ℕ) (s : ℂ) :
    ∫ y in ((k : ℝ) + a)..((k : ℝ) + a) + 1, (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1)) =
      -∫ y in ((k : ℝ) + a)..((k : ℝ) + a) + 1, (BP 1 (y - ((k : ℝ) + a)) : ℂ) * gf s 0 y := by
  rw [saw_piece, ← intervalIntegral.integral_neg]
  congr 1; funext y
  simp only [BP, gf, Fayed.FirstZero.gc]
  rw [show -s - 1 - ((0 : ℕ) : ℂ) = -(s + 1) by push_cast; ring]
  push_cast; ring

lemma gf_tendsto' (s : ℂ) (hs : 0 < s.re) (j : ℕ) : Tendsto (gf s j) atTop (𝓝 0) := by
  have h := (cpow_tendsto_zero (c := -s - 1 - (j : ℂ))
    (by simp only [sub_re, neg_re, one_re, natCast_re]
        linarith [(Nat.cast_nonneg j : (0:ℝ) ≤ j)])).const_mul (gc s j)
  rw [mul_zero] at h
  exact h

/-- Shifted order-8 Euler–Maclaurin for the sawtooth integral J. -/
theorem J_EM8 {s : ℂ} (hs : 0 < s.re) {M : ℕ} {a : ℝ} (hx0 : 0 < (M : ℝ) + a) :
    ‖J a ((M : ℝ) + a) s - EMPhi s ((M : ℝ) + a)‖ ≤
      ‖gc s 7‖ / 1000000 * (((M : ℝ) + a) ^ (-(s.re + 7)) / (s.re + 7)) := by
  set f : ℝ → ℂ := fun y => (saw a y : ℂ) * (y : ℂ) ^ (-(s + 1)) with hf
  have hfi : IntegrableOn f (Ioi ((M : ℝ) + a)) := integrableOn_saw hx0 hs
  set R : ℕ → ℂ := fun i => ∫ x in (((M + i : ℕ) : ℝ) + a)..(((M + i : ℕ) : ℝ) + a) + 1,
    (BP 8 (x - (((M + i : ℕ) : ℝ) + a)) : ℂ) * gf s 7 x with hR
  have hle : ∀ i : ℕ, (M : ℝ) + a ≤ ((M + i : ℕ) : ℝ) + a := fun i => by
    push_cast; linarith [(Nat.cast_nonneg i : (0 : ℝ) ≤ i)]
  have hnpos : ∀ i : ℕ, 0 < ((M + i : ℕ) : ℝ) + a := fun i => lt_of_lt_of_le hx0 (hle i)
  have hii : ∀ p q : ℝ, (M : ℝ) + a ≤ p → p ≤ q → IntervalIntegrable f volume p q := by
    intro p q hp hpq
    rw [intervalIntegrable_iff_integrableOn_Ioc_of_le hpq]
    exact hfi.mono_set (Ioc_subset_Ioi_self.trans (Ioi_subset_Ioi hp))
  have key : ∀ m : ℕ, ∫ x in ((M : ℝ) + a)..(((M + m : ℕ) : ℝ) + a), f x =
      EMPhi s ((M : ℝ) + a) - EMPhi s (((M + m : ℕ) : ℝ) + a) + ∑ i ∈ Finset.range m, R i := by
    intro m
    induction m with
    | zero => simp
    | succ m ih =>
      have hc : (((M + (m + 1) : ℕ) : ℝ) + a) = (((M + m : ℕ) : ℝ) + a) + 1 := by push_cast; ring
      rw [hc, ← intervalIntegral.integral_add_adjacent_intervals (hii _ _ le_rfl (hle m))
        (hii _ ((((M + m : ℕ) : ℝ) + a) + 1) (hle m) (by linarith)), ih, Finset.sum_range_succ]
      have h1 := saw_piece_BP a (M + m) s
      have h2 := interval_EM8 s (hnpos m)
      simp only [hf] at h1 ⊢
      rw [h1, h2]
      simp only [hR]
      ring
  have hbound : ∀ m : ℕ, ‖∑ i ∈ Finset.range m, R i‖ ≤
      ‖gc s 7‖ / 1000000 * (((M : ℝ) + a) ^ (-(s.re + 7)) / (s.re + 7)) := by
    intro m
    set h : ℕ → ℝ := fun i => (((M + i : ℕ) : ℝ) + a) ^ (-(s.re + 7)) / (s.re + 7) with hh
    have hpiece : ∀ i, ‖R i‖ ≤ ‖gc s 7‖ / 1000000 * (h i - h (i + 1)) := by
      intro i
      have := remainder_piece s hs (hnpos i)
      have hc : (((M + (i + 1) : ℕ) : ℝ) + a) = (((M + i : ℕ) : ℝ) + a) + 1 := by push_cast; ring
      simp only [hh, hc]
      rw [← sub_div]
      exact this
    calc ‖∑ i ∈ Finset.range m, R i‖ ≤ ∑ i ∈ Finset.range m, ‖R i‖ := norm_sum_le _ _
      _ ≤ ∑ i ∈ Finset.range m, ‖gc s 7‖ / 1000000 * (h i - h (i + 1)) :=
          Finset.sum_le_sum fun i _ => hpiece i
      _ = ‖gc s 7‖ / 1000000 * (h 0 - h m) := by
          rw [← Finset.mul_sum, Finset.sum_range_sub']
      _ ≤ ‖gc s 7‖ / 1000000 * h 0 := by
          apply mul_le_mul_of_nonneg_left _ (by positivity)
          have : 0 ≤ h m := by
            simp only [hh]
            exact div_nonneg (Real.rpow_nonneg (hnpos m).le _) (by linarith)
          linarith
      _ = _ := by simp [hh]
  have hb : Tendsto (fun m : ℕ => ((M + m : ℕ) : ℝ) + a) atTop atTop := by
    refine tendsto_atTop_add_const_right _ a ?_
    exact tendsto_natCast_atTop_atTop.comp
      (tendsto_atTop_mono (fun K => Nat.le_add_left K M) tendsto_id)
  have hlim1 := intervalIntegral_tendsto_integral_Ioi _ hfi hb
  have hE : Tendsto (EMPhi s) atTop (𝓝 0) := by
    have := ((((gf_tendsto' s hs 0).const_mul (1 / 12 : ℂ)).sub
      ((gf_tendsto' s hs 2).const_mul (1 / 720 : ℂ))).add
      ((gf_tendsto' s hs 4).const_mul (1 / 30240 : ℂ))).sub
      ((gf_tendsto' s hs 6).const_mul (1 / 1209600 : ℂ))
    simp only [mul_zero, sub_zero, add_zero] at this
    exact this
  have hlim2 := hE.comp hb
  have hlim : Tendsto (fun m : ℕ => ∑ i ∈ Finset.range m, R i) atTop
      (𝓝 ((∫ x in Ioi ((M : ℝ) + a), f x) - EMPhi s ((M : ℝ) + a) + 0)) := by
    refine ((hlim1.sub_const (EMPhi s ((M : ℝ) + a))).add hlim2).congr' ?_
    filter_upwards with m
    simp only [Function.comp]
    rw [key m]; ring
  rw [add_zero] at hlim
  exact le_of_tendsto' hlim.norm hbound

/-! ### The Euler–Maclaurin main term and the error bound -/

/-- Sj M a s = A(M+a) + s Φ(M+a). -/
noncomputable def Sj (M : ℕ) (a : ℝ) (s : ℂ) : ℂ :=
  Af s ((M : ℝ) + a) + s * EMPhi s ((M : ℝ) + a)

/-- The order-8 Euler–Maclaurin approximation of DH. -/
noncomputable def Gem (M : ℕ) (s : ℂ) : ℂ :=
  Gmain M s + (5 : ℂ) ^ (-s) *
    (Sj M (1 / 5) s - Sj M (4 / 5) s + (κ : ℂ) * (Sj M (2 / 5) s - Sj M (3 / 5) s))

/-- The error bound. -/
noncomputable def Ebd (M : ℕ) (s : ℂ) : ℝ :=
  (5 : ℝ) ^ (-s.re) * ‖s‖ * (2 + 2 * κ) *
    (‖gc s 7‖ / 1000000 * ((M : ℝ) ^ (-(s.re + 7)) / (s.re + 7)))

lemma κ_nonneg : 0 ≤ κ := by
  unfold κ
  have h9 : √9 = (3 : ℝ) := by
    rw [show (9 : ℝ) = 3 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]
  have h4 : √4 = (2 : ℝ) := by
    rw [show (4 : ℝ) = 2 ^ 2 by norm_num, Real.sqrt_sq (by norm_num)]
  have h5 : √5 ≤ 3 := by
    calc √5 ≤ √9 := Real.sqrt_le_sqrt (by norm_num)
      _ = 3 := h9
  have h2 : 2 ≤ √(10 - 2 * √5) := by
    calc (2 : ℝ) = √4 := h4.symm
      _ ≤ √(10 - 2 * √5) := Real.sqrt_le_sqrt (by linarith)
  exact div_nonneg (by linarith) (by linarith [sqrt5_gt_one])

theorem DH_sub_Gem_le {M : ℕ} (hM : 1 ≤ M) {s : ℂ} (hs : 0 < s.re) (hs1 : s ≠ 1) :
    ‖DH s - Gem M s‖ ≤ Ebd M s := by
  rw [DH_eq_Grep M hs hs1]
  set E : ℝ := ‖gc s 7‖ / 1000000 * ((M : ℝ) ^ (-(s.re + 7)) / (s.re + 7)) with hE
  have hMpos : (0 : ℝ) < M := by exact_mod_cast hM
  have hj : ∀ a : ℝ, 0 ≤ a → ‖J a ((M : ℝ) + a) s - EMPhi s ((M : ℝ) + a)‖ ≤ E := by
    intro a ha
    refine (J_EM8 hs (by linarith)).trans ?_
    rw [hE]
    refine mul_le_mul_of_nonneg_left (div_le_div_of_nonneg_right ?_ (by linarith)) (by positivity)
    exact Real.rpow_le_rpow_of_nonpos hMpos (by linarith : (M : ℝ) ≤ M + a)
      (by linarith : -(s.re + 7) ≤ 0)
  have e : Grep M s - Gem M s = (5 : ℂ) ^ (-s) * s *
      ((J (1 / 5) ((M : ℝ) + 1 / 5) s - EMPhi s ((M : ℝ) + 1 / 5))
        - (J (4 / 5) ((M : ℝ) + 4 / 5) s - EMPhi s ((M : ℝ) + 4 / 5))
        + (κ : ℂ) * ((J (2 / 5) ((M : ℝ) + 2 / 5) s - EMPhi s ((M : ℝ) + 2 / 5))
          - (J (3 / 5) ((M : ℝ) + 3 / 5) s - EMPhi s ((M : ℝ) + 3 / 5)))) := by
    simp only [Grep, Gem, Rj, Sj]; ring
  rw [e, norm_mul, norm_mul]
  have h5 : ‖(5 : ℂ) ^ (-s)‖ = (5 : ℝ) ^ (-s.re) := by
    rw [show (5 : ℂ) = ((5 : ℕ) : ℂ) by norm_num, Complex.norm_natCast_cpow_of_pos (by norm_num)]
    simp
  rw [h5]
  have hk := κ_nonneg
  have hin : ‖(J (1 / 5) ((M : ℝ) + 1 / 5) s - EMPhi s ((M : ℝ) + 1 / 5))
        - (J (4 / 5) ((M : ℝ) + 4 / 5) s - EMPhi s ((M : ℝ) + 4 / 5))
        + (κ : ℂ) * ((J (2 / 5) ((M : ℝ) + 2 / 5) s - EMPhi s ((M : ℝ) + 2 / 5))
          - (J (3 / 5) ((M : ℝ) + 3 / 5) s - EMPhi s ((M : ℝ) + 3 / 5)))‖ ≤ (2 + 2 * κ) * E := by
    refine (norm_add_le _ _).trans ?_
    rw [norm_mul, Complex.norm_real, Real.norm_eq_abs, abs_of_nonneg hk]
    have a1 := (norm_sub_le _ _).trans
      (add_le_add (hj (1 / 5) (by norm_num)) (hj (4 / 5) (by norm_num)))
    have a2 := (norm_sub_le _ _).trans
      (add_le_add (hj (2 / 5) (by norm_num)) (hj (3 / 5) (by norm_num)))
    nlinarith
  unfold Ebd
  rw [← hE]
  have h50 : 0 ≤ (5 : ℝ) ^ (-s.re) := Real.rpow_nonneg (by norm_num) _
  calc (5 : ℝ) ^ (-s.re) * ‖s‖ * ‖_‖ ≤ (5 : ℝ) ^ (-s.re) * ‖s‖ * ((2 + 2 * κ) * E) :=
        mul_le_mul_of_nonneg_left hin (by positivity)
    _ = _ := by ring

/-! ### Evaluable form of the main term -/

/-- EMq s x = Φ(x) / x^{-s}. -/
noncomputable def EMq (s : ℂ) (x : ℝ) : ℂ :=
  1 / 12 * gc s 0 / (x : ℂ) - 1 / 720 * gc s 2 / (x : ℂ) ^ 3 + 1 / 30240 * gc s 4 / (x : ℂ) ^ 5
    - 1 / 1209600 * gc s 6 / (x : ℂ) ^ 7

/-- Qf s x = x/(s-1) + 1/2 + s·EMq s x. -/
noncomputable def Qf (s : ℂ) (x : ℝ) : ℂ := (x : ℂ) / (s - 1) + 1 / 2 + s * EMq s x

lemma cpow_sub_nat {x : ℝ} (hx : 0 < x) (s : ℂ) (k : ℕ) :
    (x : ℂ) ^ (-s - 1 - (k : ℂ)) = (x : ℂ) ^ (-s) / (x : ℂ) ^ (k + 1) := by
  rw [show -s - 1 - (k : ℂ) = -s - ((k + 1 : ℕ) : ℂ) by push_cast; ring,
    cpow_sub _ _ (ofReal_ne_zero.mpr hx.ne'), cpow_natCast]

lemma five_Sj (M n : ℕ) {a : ℝ} (hx : 0 < (M : ℝ) + a) (hn : (n : ℝ) = 5 * ((M : ℝ) + a))
    {s : ℂ} (hs1 : s ≠ 1) : (5 : ℂ) ^ (-s) * Sj M a s = nps s n * Qf s ((M : ℝ) + a) := by
  set x : ℝ := (M : ℝ) + a with hxdef
  have hx' : (x : ℂ) ≠ 0 := ofReal_ne_zero.mpr hx.ne'
  have hnps : nps s n = (5 : ℂ) ^ (-s) * (x : ℂ) ^ (-s) := by
    rw [nps, show (n : ℂ) = ((5 : ℝ) : ℂ) * (x : ℂ) by exact_mod_cast hn,
      mul_cpow_ofReal_nonneg (by norm_num) hx.le]
    norm_cast
  have hs1' : s - 1 ≠ 0 := sub_ne_zero.mpr hs1
  rw [hnps]
  simp only [Sj, Af, Qf, EMq, EMPhi, gf]
  rw [cpow_one_sub' hx, ← hxdef, cpow_sub_nat hx s 0, cpow_sub_nat hx s 2, cpow_sub_nat hx s 4,
    cpow_sub_nat hx s 6]
  field_simp
  ring

/-- Gem with 5^{-s}(M+a)^{-s} merged into (5M+j)^{-s}. -/
theorem Gem_eval (M : ℕ) {s : ℂ} (hs1 : s ≠ 1) :
    Gem M s = Gmain M s + (nps s (5 * M + 1) * Qf s ((M : ℝ) + 1 / 5)
      - nps s (5 * M + 4) * Qf s ((M : ℝ) + 4 / 5)
      + (κ : ℂ) * (nps s (5 * M + 2) * Qf s ((M : ℝ) + 2 / 5)
        - nps s (5 * M + 3) * Qf s ((M : ℝ) + 3 / 5))) := by
  have h : ∀ j : ℕ, 0 < j → (5 : ℂ) ^ (-s) * Sj M ((j : ℝ) / 5) s =
      nps s (5 * M + j) * Qf s ((M : ℝ) + (j : ℝ) / 5) := fun j hj =>
    five_Sj M (5 * M + j) (by have : (0 : ℝ) < j := by exact_mod_cast hj
                              positivity) (by push_cast; ring) hs1
  have h1 := h 1 (by norm_num)
  have h2 := h 2 (by norm_num)
  have h3 := h 3 (by norm_num)
  have h4 := h 4 (by norm_num)
  push_cast at h1 h2 h3 h4
  rw [Gem]
  linear_combination h1 - h4 + (κ : ℂ) * (h2 - h3)

end Fayed.DH.Zero
