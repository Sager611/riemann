module

public import Riemann.Deep.LeanLi.ITrend

/-!
# LeanLi: `|I(n) - t(n)| ≤ 4/n` (from `H_k = log k + γ + 1/(2k) + O(1/k²)`)
-/

@[expose] public section

open Real Filter Topology Set

namespace LeanLi

/-- `e_k = H_k - log k - 1/(2k)`. -/
noncomputable def eSeq (k : ℕ) : ℝ := (harmonic k : ℝ) - log k - 1 / (2 * k)

lemma f_bounds {x : ℝ} (hx : 0 ≤ x) :
    0 ≤ x / 2 + x / (2 * (1 + x)) - log (1 + x) ∧
      x / 2 + x / (2 * (1 + x)) - log (1 + x) ≤ x ^ 3 / 6 := by
  -- g(x) = x/2 + x/(2(1+x)) - log(1+x), g(0) = 0, g'(x) = x²/(2(1+x)²) ∈ [0, x²/2]
  set g : ℝ → ℝ := fun y => y / 2 + y / (2 * (1 + y)) - log (1 + y)
  have hd : ∀ y, 0 ≤ y → HasDerivAt g (y ^ 2 / (2 * (1 + y) ^ 2)) y := by
    intro y hy
    have h1 : (1 + y) ≠ 0 := by linarith
    have h2 : (2 * (1 + y)) ≠ 0 := by positivity
    have hA : HasDerivAt (fun y : ℝ => y / 2) (1 / 2) y := by
      simpa using (hasDerivAt_id y).div_const 2
    have hB : HasDerivAt (fun y : ℝ => y / (2 * (1 + y))) (1 / (2 * (1 + y) ^ 2)) y := by
      have hden : HasDerivAt (fun y : ℝ => 2 * (1 + y)) 2 y := by
        simpa using ((hasDerivAt_id y).const_add 1).const_mul 2
      convert (hasDerivAt_id' y).div hden h2 using 1
      field_simp; ring
    have hC : HasDerivAt (fun y : ℝ => log (1 + y)) (1 / (1 + y)) y := by
      have := ((hasDerivAt_id y).const_add 1).log h1
      simpa using this
    convert (hA.add hB).sub hC using 1
    field_simp; ring
  have hg0 : g 0 = 0 := by simp [g]
  have hcont : ContinuousOn g (Icc 0 x) := fun y hy => (hd y hy.1).continuousAt.continuousWithinAt
  constructor
  · -- monotone
    have hmono : MonotoneOn g (Icc 0 x) := by
      refine monotoneOn_of_deriv_nonneg (convex_Icc 0 x) hcont ?_ ?_
      · intro y hy
        rw [interior_Icc] at hy
        exact (hd y hy.1.le).differentiableAt.differentiableWithinAt
      · intro y hy
        rw [interior_Icc] at hy
        rw [(hd y hy.1.le).deriv]; positivity
    have := hmono ⟨le_refl 0, hx⟩ ⟨hx, le_refl x⟩ hx
    rw [hg0] at this
    exact this
  · -- h(y) = y³/6 - g(y) monotone
    set h : ℝ → ℝ := fun y => y ^ 3 / 6 - g y
    have hhd : ∀ y, 0 ≤ y → HasDerivAt h (y ^ 2 / 2 - y ^ 2 / (2 * (1 + y) ^ 2)) y := by
      intro y hy
      have := ((hasDerivAt_pow 3 y).div_const 6).sub (hd y hy)
      convert this using 1; push_cast; ring
    have hmono : MonotoneOn h (Icc 0 x) := by
      refine monotoneOn_of_deriv_nonneg (convex_Icc 0 x)
        (fun y hy => (hhd y hy.1).continuousAt.continuousWithinAt) ?_ ?_
      · intro y hy
        rw [interior_Icc] at hy
        exact (hhd y hy.1.le).differentiableAt.differentiableWithinAt
      · intro y hy
        rw [interior_Icc] at hy
        rw [(hhd y hy.1.le).deriv]
        have hy0 := hy.1.le
        have h1 : 1 ≤ (1 + y) ^ 2 := by nlinarith
        have : y ^ 2 / (2 * (1 + y) ^ 2) ≤ y ^ 2 / 2 := by
          apply div_le_div_of_nonneg_left (sq_nonneg y) (by norm_num) (by linarith)
        linarith
    have := hmono ⟨le_refl 0, hx⟩ ⟨hx, le_refl x⟩ hx
    simp only [h, hg0] at this
    simp only [g] at this ⊢
    linarith

lemma eSeq_succ_sub (k : ℕ) (hk : 1 ≤ k) :
    0 ≤ eSeq (k + 1) - eSeq k ∧ eSeq (k + 1) - eSeq k ≤ 1 / (6 * k ^ 3) := by
  have hk0 : (0 : ℝ) < k := by exact_mod_cast hk
  have e : eSeq (k + 1) - eSeq k =
      (1 / k) / 2 + (1 / k) / (2 * (1 + 1 / k)) - log (1 + 1 / k) := by
    unfold eSeq
    rw [harmonic_succ]
    push_cast
    have : log ((k : ℝ) + 1) = log k + log (1 + 1 / k) := by
      have e : ((k : ℝ) + 1) = k * (1 + 1 / k) := by field_simp
      rw [e, log_mul hk0.ne' (by positivity)]
    rw [this]
    field_simp
    ring
  obtain ⟨h1, h2⟩ := f_bounds (x := 1 / (k : ℝ)) (by positivity)
  rw [e]
  refine ⟨h1, h2.trans (le_of_eq ?_)⟩
  field_simp

lemma tendsto_eSeq : Tendsto eSeq atTop (𝓝 eulerMascheroniConstant) := by
  have h1 := tendsto_harmonic_sub_log
  have h2 : Tendsto (fun k : ℕ => 1 / (2 * (k : ℝ))) atTop (𝓝 0) := by
    have := (tendsto_inv_atTop_zero.comp tendsto_natCast_atTop_atTop).const_mul (1 / 2 : ℝ)
    simp only [mul_zero] at this
    refine this.congr fun k => ?_
    simp only [Function.comp]
    rw [one_div_mul_eq_div, one_div, mul_inv, mul_comm, div_eq_mul_inv]
  have := h1.sub h2
  simp only [sub_zero] at this
  exact this.congr fun k => by unfold eSeq; ring

lemma eSeq_le (n : ℕ) (hn : 1 ≤ n) : eSeq n ≤ eulerMascheroniConstant := by
  have hmono : ∀ M, eSeq n ≤ eSeq (n + M) := by
    intro M
    induction M with
    | zero => simp
    | succ M ih =>
      have := (eSeq_succ_sub (n + M) (by omega)).1
      rw [show n + (M + 1) = n + M + 1 by ring]; linarith
  apply ge_of_tendsto tendsto_eSeq
  filter_upwards [eventually_ge_atTop n] with M hM
  obtain ⟨d, rfl⟩ := Nat.exists_eq_add_of_le hM
  exact hmono d

lemma inv_cube_le (k : ℕ) (hk : 1 ≤ k) :
    1 / (6 * (k : ℝ) ^ 3) ≤ (1 / 12) * (1 / ((k : ℝ) - 1 / 2) ^ 2 - 1 / ((k : ℝ) + 1 / 2) ^ 2) := by
  have hk1 : (1 : ℝ) ≤ k := by exact_mod_cast hk
  have ha : 0 < (k : ℝ) - 1 / 2 := by linarith
  have hb : 0 < (k : ℝ) + 1 / 2 := by linarith
  have e : 1 / ((k : ℝ) - 1 / 2) ^ 2 - 1 / ((k : ℝ) + 1 / 2) ^ 2 =
      2 * k / (((k : ℝ) - 1 / 2) ^ 2 * ((k : ℝ) + 1 / 2) ^ 2) := by
    have h1 := ha.ne'
    have h2 := hb.ne'
    rw [div_sub_div _ _ (pow_ne_zero 2 h1) (pow_ne_zero 2 h2)]
    congr 1; ring
  rw [e]
  have hD : 0 < ((k : ℝ) - 1 / 2) ^ 2 * ((k : ℝ) + 1 / 2) ^ 2 := by positivity
  have hP : ((k : ℝ) - 1 / 2) ^ 2 * ((k : ℝ) + 1 / 2) ^ 2 ≤ (k : ℝ) ^ 4 := by
    have : ((k : ℝ) - 1 / 2) ^ 2 * ((k : ℝ) + 1 / 2) ^ 2 = ((k : ℝ) ^ 2 - 1 / 4) ^ 2 := by ring
    rw [this]
    have h0 : 0 ≤ (k : ℝ) ^ 2 - 1 / 4 := by nlinarith
    nlinarith
  rw [div_mul_div_comm, one_mul, div_le_div_iff₀ (by positivity) (by positivity)]
  nlinarith

lemma eSeq_ge (n : ℕ) (hn : 1 ≤ n) :
    eulerMascheroniConstant - 1 / (3 * (n : ℝ) ^ 2) ≤ eSeq n := by
  have hn1 : (1 : ℝ) ≤ n := by exact_mod_cast hn
  have hsum : ∀ M, eSeq (n + M) - eSeq n ≤
      (1 / 12) * (1 / ((n : ℝ) - 1 / 2) ^ 2 - 1 / (((n + M : ℕ) : ℝ) - 1 / 2) ^ 2) := by
    intro M
    induction M with
    | zero => simp
    | succ M ih =>
      have h := (eSeq_succ_sub (n + M) (by omega)).2
      have h2 := inv_cube_le (n + M) (by omega)
      rw [show n + (M + 1) = n + M + 1 by ring]
      push_cast at ih h h2 ⊢
      have e : ((n : ℝ) + M + 1 - 1 / 2) = ((n : ℝ) + M + 1 / 2) := by ring
      rw [e]
      linarith
  have hlim : eulerMascheroniConstant - eSeq n ≤ (1 / 12) * (1 / ((n : ℝ) - 1 / 2) ^ 2) := by
    have ht : Tendsto (fun M : ℕ => eSeq (n + M) - eSeq n) atTop
        (𝓝 (eulerMascheroniConstant - eSeq n)) :=
      (tendsto_eSeq.comp (tendsto_atTop_atTop_of_monotone (fun a b h => by omega)
        (fun b => ⟨b, by omega⟩))).sub_const _
    apply le_of_tendsto ht
    filter_upwards with M
    refine (hsum M).trans ?_
    have : 0 ≤ 1 / (((n + M : ℕ) : ℝ) - 1 / 2) ^ 2 := by positivity
    linarith
  have h3 : (1 / 12) * (1 / ((n : ℝ) - 1 / 2) ^ 2) ≤ 1 / (3 * (n : ℝ) ^ 2) := by
    have ha : 0 < (n : ℝ) - 1 / 2 := by linarith
    rw [show (1 : ℝ) / 12 * (1 / (↑n - 1 / 2) ^ 2) = 1 / (12 * (↑n - 1 / 2) ^ 2) by
      rw [one_div_mul_one_div]]
    apply one_div_le_one_div_of_le (by positivity)
    nlinarith
  linarith

/-- `H_k = log k + γ + 1/(2k) + ε_k`, `-1/(3k²) ≤ ε_k ≤ 0`. -/
lemma harmonic_fine (k : ℕ) (hk : 1 ≤ k) :
    -(1 / (3 * (k : ℝ) ^ 2)) ≤ (harmonic k : ℝ) - log k - eulerMascheroniConstant - 1 / (2 * k) ∧
      (harmonic k : ℝ) - log k - eulerMascheroniConstant - 1 / (2 * k) ≤ 0 := by
  have h1 := eSeq_le k hk
  have h2 := eSeq_ge k hk
  unfold eSeq at h1 h2
  constructor <;> linarith

lemma aux_odd (m : ℝ) (hm : 1 ≤ m) : 1 / (4 * m) + 3 / (4 * m) ≤ 4 / (2 * m + 1) := by
  rw [← add_div, div_le_div_iff₀ (by positivity) (by positivity)]
  linarith

/-- **`|I(n) - t(n)| ≤ 4/n`** for `n ≥ 1`. -/
theorem Iclosed_sub_trend_le_div (n : ℕ) (hn : 1 ≤ n) : |Iclosed n - trendT n| ≤ 4 / n := by
  have hπ : 0 < π := pi_pos
  have hlog4 : log (4 * π) = log 2 + log (2 * π) := by
    rw [show (4 : ℝ) * π = 2 * (2 * π) by ring, log_mul (by norm_num) (by positivity)]
  rcases Nat.even_or_odd' n with ⟨m, rfl | rfl⟩
  · -- n = 2m
    have hm1 : 1 ≤ m := by omega
    have hmpos : (0 : ℝ) < m := by exact_mod_cast hm1
    obtain ⟨hO, _⟩ := oddRecip_two_mul m
    obtain ⟨hP, _⟩ := card_odd_two_mul m
    obtain ⟨a1, a2⟩ := harmonic_fine (2 * m) (by omega)
    obtain ⟨b1, b2⟩ := harmonic_fine m hm1
    set d1 := (harmonic (2 * m) : ℝ) - log ((2 * m : ℕ) : ℝ) - eulerMascheroniConstant -
      1 / (2 * ((2 * m : ℕ) : ℝ)) with hd1
    set d2 := (harmonic m : ℝ) - log (m : ℝ) - eulerMascheroniConstant - 1 / (2 * (m : ℝ)) with hd2
    have hlog2m : log ((2 * m : ℕ) : ℝ) = log 2 + log m := by
      push_cast; rw [log_mul (by norm_num) hmpos.ne']
    have key : Iclosed (2 * m) - trendT (2 * m) = 2 * m * d1 - m * d2 := by
      rw [Iclosed_eq, hO, hP]
      unfold trendT
      rw [hd1, hd2, hlog2m, hlog4]
      push_cast
      field_simp
      ring
    rw [key]
    have e1 : |2 * (m : ℝ) * d1| ≤ 1 / (6 * m) := by
      rw [abs_le]
      have : (2 * (m : ℝ)) * (1 / (3 * ((2 * m : ℕ) : ℝ) ^ 2)) = 1 / (6 * m) := by
        push_cast; field_simp; ring
      constructor
      · have := mul_le_mul_of_nonneg_left a1 (by positivity : (0 : ℝ) ≤ 2 * m)
        nlinarith
      · nlinarith [mul_nonpos_iff.2 (Or.inl ⟨(by positivity : (0 : ℝ) ≤ 2 * m), a2⟩)]
    have e2 : |(m : ℝ) * d2| ≤ 1 / (3 * m) := by
      rw [abs_le]
      have : (m : ℝ) * (1 / (3 * (m : ℝ) ^ 2)) = 1 / (3 * m) := by field_simp
      constructor
      · have := mul_le_mul_of_nonneg_left b1 hmpos.le
        nlinarith
      · nlinarith [mul_nonpos_iff.2 (Or.inl ⟨hmpos.le, b2⟩)]
    calc |2 * (m : ℝ) * d1 - m * d2| ≤ |2 * (m : ℝ) * d1| + |(m : ℝ) * d2| := abs_sub _ _
      _ ≤ 1 / (6 * m) + 1 / (3 * m) := add_le_add e1 e2
      _ ≤ 4 / ((2 * m : ℕ) : ℝ) := by
        push_cast
        rw [div_add_div _ _ (by positivity) (by positivity), div_le_div_iff₀ (by positivity) (by positivity)]
        nlinarith
  · -- n = 2m + 1
    rcases Nat.eq_zero_or_pos m with h0 | hm1
    · subst h0
      have := Iclosed_sub_trend_le 1
      norm_num at this ⊢
      linarith
    have hmpos : (0 : ℝ) < m := by exact_mod_cast hm1
    obtain ⟨hO, hO'⟩ := oddRecip_two_mul m
    obtain ⟨_, hP'⟩ := card_odd_two_mul m
    obtain ⟨a1, a2⟩ := harmonic_fine (2 * m) (by omega)
    obtain ⟨b1, b2⟩ := harmonic_fine m hm1
    set d1 := (harmonic (2 * m) : ℝ) - log ((2 * m : ℕ) : ℝ) - eulerMascheroniConstant -
      1 / (2 * ((2 * m : ℕ) : ℝ)) with hd1
    set d2 := (harmonic m : ℝ) - log (m : ℝ) - eulerMascheroniConstant - 1 / (2 * (m : ℝ)) with hd2
    have hlog2m : log ((2 * m : ℕ) : ℝ) = log 2 + log m := by
      push_cast; rw [log_mul (by norm_num) hmpos.ne']
    set ℓ := log ((2 * m + 1 : ℕ) : ℝ) - log ((2 * m : ℕ) : ℝ) with hℓ
    have key : Iclosed (2 * m + 1) - trendT (2 * m + 1) =
        -(((2 * m + 1 : ℝ)) / 2 * ℓ) + 1 / 2 + (2 * m + 1) * (d1 - d2 / 2) := by
      rw [Iclosed_eq, hO', hO, hP']
      unfold trendT
      rw [hd1, hd2, hℓ, hlog2m, hlog4]
      push_cast
      field_simp
      ring
    rw [key]
    -- bounds on ℓ = log(1 + 1/(2m))
    have hℓe : ℓ = log (1 + 1 / (2 * m)) := by
      rw [hℓ, ← log_div (by positivity) (by positivity)]
      congr 1; push_cast; field_simp
    have hx : (0 : ℝ) < 1 / (2 * m) := by positivity
    have hℓ1 : ℓ ≤ 1 / (2 * m) := by
      rw [hℓe]; have := log_le_sub_one_of_pos (show (0 : ℝ) < 1 + 1 / (2 * m) by positivity); linarith
    have hℓ0 : 1 / (2 * m + 1) ≤ ℓ := by
      rw [hℓe]
      have := one_sub_inv_le_log_of_pos (show (0 : ℝ) < 1 + 1 / (2 * m) by positivity)
      have e : 1 - (1 + 1 / (2 * (m : ℝ)))⁻¹ = 1 / (2 * m + 1) := by field_simp; ring
      linarith
    have t1 : |-(((2 * m + 1 : ℝ)) / 2 * ℓ) + 1 / 2| ≤ 1 / (4 * m) := by
      rw [abs_le]
      have h1 : (2 * m + 1 : ℝ) / 2 * (1 / (2 * m + 1)) = 1 / 2 := by field_simp
      have h2 : (2 * m + 1 : ℝ) / 2 * (1 / (2 * m)) = 1 / 2 + 1 / (4 * m) := by field_simp; ring
      have g1 := mul_le_mul_of_nonneg_left hℓ0 (by positivity : (0 : ℝ) ≤ (2 * m + 1) / 2)
      have g2 := mul_le_mul_of_nonneg_left hℓ1 (by positivity : (0 : ℝ) ≤ (2 * m + 1) / 2)
      have : 0 < 1 / (4 * (m : ℝ)) := by positivity
      constructor <;> linarith
    have t2 : |(2 * m + 1 : ℝ) * (d1 - d2 / 2)| ≤ 3 / (4 * m) := by
      have ha : 1 / (3 * ((2 * m : ℕ) : ℝ) ^ 2) = 1 / (12 * (m : ℝ) ^ 2) := by
        push_cast; field_simp; ring
      rw [ha] at a1
      have hd : |d1 - d2 / 2| ≤ 1 / (4 * (m : ℝ) ^ 2) := by
        rw [abs_le]
        have : 1 / (4 * (m : ℝ) ^ 2) = 1 / (12 * m ^ 2) + 1 / (3 * m ^ 2) / 2 := by field_simp; ring
        constructor <;> linarith
      rw [abs_mul, abs_of_pos (by positivity : (0 : ℝ) < 2 * m + 1)]
      calc (2 * m + 1 : ℝ) * |d1 - d2 / 2| ≤ (2 * m + 1) * (1 / (4 * (m : ℝ) ^ 2)) :=
            mul_le_mul_of_nonneg_left hd (by positivity)
        _ ≤ 3 / (4 * m) := by
            rw [show (2 * m + 1 : ℝ) * (1 / (4 * (m : ℝ) ^ 2)) = (2 * m + 1) / (4 * m ^ 2) by ring,
              div_le_div_iff₀ (by positivity) (by positivity)]
            have : (1 : ℝ) ≤ m := by exact_mod_cast hm1
            nlinarith
    calc |-(((2 * m + 1 : ℝ)) / 2 * ℓ) + 1 / 2 + (2 * m + 1) * (d1 - d2 / 2)|
        ≤ |-(((2 * m + 1 : ℝ)) / 2 * ℓ) + 1 / 2| + |(2 * m + 1 : ℝ) * (d1 - d2 / 2)| :=
          abs_add_le _ _
      _ ≤ 1 / (4 * m) + 3 / (4 * m) := add_le_add t1 t2
      _ ≤ 4 / ((2 * m + 1 : ℕ) : ℝ) := by
        have h := aux_odd (m : ℝ) (by exact_mod_cast hm1)
        have e : (((2 * m + 1 : ℕ) : ℝ)) = 2 * (m : ℝ) + 1 := by push_cast; ring
        rw [e]; exact h

end LeanLi
