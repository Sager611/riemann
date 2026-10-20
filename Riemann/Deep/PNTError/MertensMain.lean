import Riemann.Deep.PNTError.Mertens

/-!
# Power-saving PNT, part 6: assembly of M(x) = O(x^{σ₀+ε})
-/

open Complex Set MeasureTheory intervalIntegral Filter Topology
open scoped Real Interval ArithmeticFunction.Moebius

namespace Fayed.PNTError

noncomputable section

lemma mertens_bound_aux {σ₀ σ₁ x T C B η : ℝ} (hZ : ZeroFree σ₀) (h01 : σ₀ < σ₁)
    (hσ₁0 : 0 < σ₁) (hσ₁ : σ₁ < 1) (hC : ∀ s : ℂ, σ₁ ≤ s.re → 4 ≤ |s.im| →
      ‖(riemannZeta s)⁻¹‖ ≤ C * |s.im| ^ η) (hC0 : 0 ≤ C) (hη : 0 ≤ η)
    (hB : ∀ y : ℝ, |y| ≤ 4 → ‖(riemannZeta (σ₁ + y * I))⁻¹‖ ≤ B) (hB0 : 0 ≤ B)
    (hx : 1 ≤ x) (hT : 4 ≤ T) (hfar : ∀ n : ℕ, 1 / 2 ≤ |x - n|) :
    |mertens x| ≤ 4 * x ^ 4 * (x + 1) / (π * T) * (π ^ 2 / 6) +
      1 / (2 * π) * (4 * (C * T ^ η * x ^ 4 / T) + 4 * (C * T ^ η * x ^ 4 / T)) +
      1 / (2 * π) * (x ^ σ₁ * (5 * B / σ₁ + 2 * (C * T ^ η)) * (2 * Real.log (1 + T))) := by
  have hx0 : 0 < x := by linarith
  have hT0 : 0 < T := by linarith
  have hP := perron_mu hx hT0 hfar
  have hS := shift_mu hZ h01 hσ₁0 hσ₁ hx0 hT0
  have hHb := horiz_mu hC hC0 hσ₁0.le (by linarith) hx hT (y := -T)
    (by rw [abs_neg, abs_of_pos hT0])
  have hHt := horiz_mu hC hC0 hσ₁0.le (by linarith) hx hT (y := T) (abs_of_pos hT0)
  have hV := vert_mu hC hC0 hη hB hB0 hσ₁0 hx hT
  have hc : ‖((1 / (2 * π) : ℝ) : ℂ)‖ = 1 / (2 * π) := by
    rw [Complex.norm_real, Real.norm_eq_abs, abs_of_pos (by positivity)]
  have e : ((mertens x : ℝ) : ℂ) =
      ((mertens x : ℂ) - ((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-T)..T, Fμ x ((4 : ℝ) + y * I)) +
      ((1 / (2 * π) : ℝ) : ℂ) * (I * (∫ t in σ₁..4, Fμ x (t + (-T : ℝ) * I)) -
        I * (∫ t in σ₁..4, Fμ x (t + (T : ℝ) * I))) +
      ((1 / (2 * π) : ℝ) : ℂ) * ∫ y in (-T)..T, Fμ x ((σ₁ : ℂ) + y * I) := by
    rw [← hS]
    ring
  have hnorm : |mertens x| = ‖((mertens x : ℝ) : ℂ)‖ := by
    rw [Complex.norm_real, Real.norm_eq_abs]
  rw [hnorm, e]
  refine (norm_add_le _ _).trans (add_le_add ((norm_add_le _ _).trans (add_le_add hP ?_)) ?_)
  · rw [norm_mul, hc]
    refine mul_le_mul_of_nonneg_left ((norm_sub_le _ _).trans ?_) (by positivity)
    simp only [norm_mul, Complex.norm_I, one_mul]
    exact add_le_add hHb hHt
  · rw [norm_mul, hc]
    exact mul_le_mul_of_nonneg_left hV (by positivity)

lemma num_m2 {C η x : ℝ} (hC : 0 ≤ C) (hη0 : 0 ≤ η) (hη : 5 * η ≤ 1) (hx : 1 ≤ x) :
    1 / (2 * π) * (4 * (C * (x ^ 5) ^ η * x ^ 4 / x ^ 5) + 4 * (C * (x ^ 5) ^ η * x ^ 4 / x ^ 5)) ≤
      4 * C / π := by
  have hx0 : 0 < x := by linarith
  have hπ : 0 < π := Real.pi_pos
  have hpow : (x ^ 5 : ℝ) ^ η ≤ x := by
    rw [← Real.rpow_natCast x 5, ← Real.rpow_mul hx0.le]
    calc x ^ (((5 : ℕ) : ℝ) * η) ≤ x ^ (1 : ℝ) :=
          Real.rpow_le_rpow_of_exponent_le hx (by push_cast; linarith)
      _ = x := Real.rpow_one x
  have hk : C * (x ^ 5) ^ η * x ^ 4 / x ^ 5 ≤ C := by
    rw [div_le_iff₀ (by positivity)]
    calc C * (x ^ 5) ^ η * x ^ 4 ≤ C * x * x ^ 4 :=
          mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left hpow hC) (by positivity)
      _ = C * x ^ 5 := by ring
  calc 1 / (2 * π) * (4 * (C * (x ^ 5) ^ η * x ^ 4 / x ^ 5) +
        4 * (C * (x ^ 5) ^ η * x ^ 4 / x ^ 5))
      = 4 / π * (C * (x ^ 5) ^ η * x ^ 4 / x ^ 5) := by ring
    _ ≤ 4 / π * C := mul_le_mul_of_nonneg_left hk (by positivity)
    _ = 4 * C / π := by ring

lemma num_m3 {σ₁ B C η x : ℝ} (hσ₁0 : 0 < σ₁) (hB : 0 ≤ B) (hC : 0 ≤ C) (hη0 : 0 ≤ η)
    (hx : 2 ≤ x) :
    1 / (2 * π) * (x ^ σ₁ * (5 * B / σ₁ + 2 * (C * (x ^ 5) ^ η)) * (2 * Real.log (1 + x ^ 5))) ≤
      6 / π * (5 * B / σ₁ + 2 * C) * x ^ (σ₁ + 5 * η) * Real.log x := by
  have hx0 : 0 < x := by linarith
  have hπ : 0 < π := Real.pi_pos
  have hxs : 0 < x ^ σ₁ := Real.rpow_pos_of_pos hx0 σ₁
  have hL2 : Real.log 2 ≤ Real.log x := Real.log_le_log (by norm_num) hx
  have hlog2 : 0 < Real.log 2 := Real.log_pos (by norm_num)
  have hT0 : 0 < x ^ 5 := by positivity
  have hT1 : (1 : ℝ) ≤ x ^ 5 := one_le_pow₀ (by linarith)
  have hlog1T : Real.log (1 + x ^ 5) ≤ 6 * Real.log x := by
    calc Real.log (1 + x ^ 5) ≤ Real.log (2 * x ^ 5) := Real.log_le_log (by linarith) (by linarith)
      _ = Real.log 2 + 5 * Real.log x := by
          rw [Real.log_mul (by norm_num) hT0.ne', Real.log_pow]
          norm_num
      _ ≤ 6 * Real.log x := by linarith
  have hl0 : 0 ≤ Real.log (1 + x ^ 5) := Real.log_nonneg (by linarith)
  have hpow : (x ^ 5 : ℝ) ^ η = x ^ (5 * η) := by
    rw [← Real.rpow_natCast x 5, ← Real.rpow_mul hx0.le]
    norm_num
  have hpow1 : 1 ≤ x ^ (5 * η) := Real.one_le_rpow (by linarith) (by linarith)
  have hB' : 0 ≤ 5 * B / σ₁ := div_nonneg (by linarith) hσ₁0.le
  have hK : 5 * B / σ₁ + 2 * (C * (x ^ 5) ^ η) ≤ (5 * B / σ₁ + 2 * C) * x ^ (5 * η) := by
    rw [hpow]
    nlinarith [mul_nonneg hB' (sub_nonneg.mpr hpow1)]
  have hK0 : 0 ≤ (5 * B / σ₁ + 2 * C) * x ^ (5 * η) :=
    mul_nonneg (by linarith) (by linarith)
  calc 1 / (2 * π) * (x ^ σ₁ * (5 * B / σ₁ + 2 * (C * (x ^ 5) ^ η)) * (2 * Real.log (1 + x ^ 5)))
      ≤ 1 / (2 * π) * (x ^ σ₁ * ((5 * B / σ₁ + 2 * C) * x ^ (5 * η)) * (12 * Real.log x)) := by
        apply mul_le_mul_of_nonneg_left _ (by positivity)
        exact mul_le_mul (mul_le_mul_of_nonneg_left hK hxs.le) (by linarith) (by linarith)
          (mul_nonneg hxs.le hK0)
    _ = 6 / π * (5 * B / σ₁ + 2 * C) * (x ^ σ₁ * x ^ (5 * η)) * Real.log x := by ring
    _ = 6 / π * (5 * B / σ₁ + 2 * C) * x ^ (σ₁ + 5 * η) * Real.log x := by
        rw [← Real.rpow_add hx0]

theorem mertens_halfint {σ₀ : ℝ} (hZ : ZeroFree σ₀) (h12 : 1 / 2 ≤ σ₀) {σ₁ η : ℝ}
    (h01 : σ₀ < σ₁) (hσ₁ : σ₁ < 1) (hη0 : 0 < η) (hη : 5 * η ≤ 1) :
    ∃ A B : ℝ, 0 ≤ A ∧ 0 ≤ B ∧ ∀ x : ℝ, 2 ≤ x → (∀ n : ℕ, 1 / 2 ≤ |x - n|) →
      |mertens x| ≤ A + B * x ^ (σ₁ + 5 * η) * Real.log x := by
  have h1 : σ₀ < 1 := lt_trans h01 hσ₁
  have hσ₁0 : 0 < σ₁ := by linarith
  have hπ : 0 < π := Real.pi_pos
  obtain ⟨C, hC0, hC⟩ := inv_zeta_bound hZ h12 h1 (sub_pos.mpr h01) hη0
  have hC' : ∀ s : ℂ, σ₁ ≤ s.re → 4 ≤ |s.im| → ‖(riemannZeta s)⁻¹‖ ≤ C * |s.im| ^ η :=
    fun s hs ht => hC s (by linarith) ht
  obtain ⟨B₀, hB₀0, hB₀⟩ := vert_compact_mu hZ h01 hσ₁
  refine ⟨2 * π + 4 * C / π, 6 / π * (5 * B₀ / σ₁ + 2 * C),
    add_nonneg (by positivity) (div_nonneg (by linarith) hπ.le),
    mul_nonneg (div_nonneg (by norm_num) hπ.le)
      (add_nonneg (div_nonneg (by linarith) hσ₁0.le) (by linarith)), ?_⟩
  intro x hx hfar
  have hx1 : 1 ≤ x := by linarith
  have hT4 : 4 ≤ x ^ 5 := by
    have := pow_le_pow_left₀ (by norm_num : (0 : ℝ) ≤ 2) hx 5
    norm_num at this
    linarith
  have hm := mertens_bound_aux hZ h01 hσ₁0 hσ₁ hC' hC0 hη0.le hB₀ hB₀0 hx1 hT4 hfar
  have b1 := num_b1 hx
  have b2 := num_m2 hC0 hη0.le hη hx1
  have b3 := num_m3 hσ₁0 hB₀0 hC0 hη0.le hx
  linarith only [hm, b1, b2, b3]

lemma abs_mertens_le {x : ℝ} (hx : 0 ≤ x) : |mertens x| ≤ x := by
  calc |mertens x| ≤ ∑ n ∈ Finset.Ioc 0 ⌊x⌋₊, |(μ n : ℝ)| := Finset.abs_sum_le_sum_abs _ _
    _ ≤ ∑ n ∈ Finset.Ioc 0 ⌊x⌋₊, (1 : ℝ) :=
        Finset.sum_le_sum (fun n _ => by exact_mod_cast ArithmeticFunction.abs_moebius_le_one)
    _ = ⌊x⌋₊ := by simp
    _ ≤ x := Nat.floor_le hx

/-- **Rung 3 (general form).** If ζ has no zero in Re s > σ₀ with 1/2 ≤ σ₀ < 1, then for every
ε > 0 there is C with |M(x)| ≤ C x^{σ₀+ε} for all x ≥ 2. -/
theorem mertens_error_of_zeroFree {σ₀ : ℝ} (hZ : ZeroFree σ₀) (h12 : 1 / 2 ≤ σ₀) (h1 : σ₀ < 1)
    {ε : ℝ} (hε : 0 < ε) :
    ∃ C : ℝ, ∀ x : ℝ, 2 ≤ x → |mertens x| ≤ C * x ^ (σ₀ + ε) := by
  by_cases hbig : 1 ≤ σ₀ + ε
  · refine ⟨1, fun x hx => ?_⟩
    have hxe : x ≤ x ^ (σ₀ + ε) := by
      calc x = x ^ (1 : ℝ) := (Real.rpow_one x).symm
        _ ≤ x ^ (σ₀ + ε) := Real.rpow_le_rpow_of_exponent_le (by linarith) hbig
    linarith [abs_mertens_le (show (0 : ℝ) ≤ x by linarith)]
  · rw [not_le] at hbig
    obtain ⟨A, B, hA, hB, hmain⟩ := mertens_halfint hZ h12 (σ₁ := σ₀ + ε / 2) (η := ε / 20)
      (by linarith) (by linarith) (by linarith) (by linarith)
    obtain ⟨D, hD⟩ : ∃ D : ℝ, D = A + B * (4 / ε) := ⟨_, rfl⟩
    have hD0 : 0 ≤ D := by rw [hD]; exact add_nonneg hA (mul_nonneg hB (by positivity))
    have hhalf : ∀ y : ℝ, 2 ≤ y → (∀ n : ℕ, 1 / 2 ≤ |y - n|) →
        |mertens y| ≤ D * y ^ (σ₀ + ε) := by
      intro y hy hfar
      have hy0 : 0 < y := by linarith
      have hy1 : 1 ≤ y := by linarith
      have hm := hmain y hy hfar
      have hq : 0 < ε / 4 := by linarith
      have hlog : Real.log y ≤ 4 / ε * y ^ (ε / 4) := by
        calc Real.log y ≤ y ^ (ε / 4) / (ε / 4) := Real.log_le_rpow_div hy0.le hq
          _ = 4 / ε * y ^ (ε / 4) := by field_simp
      have hexp : y ^ (σ₀ + ε / 2 + 5 * (ε / 20)) * y ^ (ε / 4) = y ^ (σ₀ + ε) := by
        rw [← Real.rpow_add hy0]
        ring_nf
      have hys : 0 ≤ y ^ (σ₀ + ε / 2 + 5 * (ε / 20)) := Real.rpow_nonneg hy0.le _
      have hy1' : 1 ≤ y ^ (σ₀ + ε) := Real.one_le_rpow hy1 (by linarith)
      calc |mertens y| ≤ A + B * y ^ (σ₀ + ε / 2 + 5 * (ε / 20)) * Real.log y := hm
        _ ≤ A * y ^ (σ₀ + ε) + B * y ^ (σ₀ + ε / 2 + 5 * (ε / 20)) * (4 / ε * y ^ (ε / 4)) := by
            apply add_le_add
            · have := mul_le_mul_of_nonneg_left hy1' hA
              linarith
            · exact mul_le_mul_of_nonneg_left hlog (mul_nonneg hB hys)
        _ = D * y ^ (σ₀ + ε) := by
            rw [hD, ← hexp]
            ring
    refine ⟨2 * D, fun x hx => ?_⟩
    have hx0 : 0 < x := by linarith
    obtain ⟨N, hN⟩ : ∃ N : ℕ, N = ⌊x⌋₊ := ⟨_, rfl⟩
    have hNx : (N : ℝ) ≤ x := by rw [hN]; exact Nat.floor_le hx0.le
    have hxN : x < N + 1 := by rw [hN]; exact Nat.lt_floor_add_one x
    have hN2 : (2 : ℝ) ≤ N := by
      have : 2 ≤ N := by rw [hN]; exact Nat.le_floor (by exact_mod_cast hx)
      exact_mod_cast this
    obtain ⟨y, hy⟩ : ∃ y : ℝ, y = N + 1 / 2 := ⟨_, rfl⟩
    have hy2 : 2 ≤ y := by rw [hy]; linarith
    have hfar : ∀ n : ℕ, 1 / 2 ≤ |y - n| := by
      intro n
      rcases le_or_gt n N with h | h
      · have : (n : ℝ) ≤ N := by exact_mod_cast h
        rw [abs_of_pos (by rw [hy]; linarith), hy]
        linarith
      · have : (N : ℝ) + 1 ≤ n := by exact_mod_cast h
        rw [abs_of_neg (by rw [hy]; linarith), hy]
        linarith
    have hfloor : ⌊y⌋₊ = N := by
      rw [Nat.floor_eq_iff (by linarith)]
      constructor <;> rw [hy] <;> linarith
    have hMxy : mertens x = mertens y := by
      unfold mertens
      rw [hfloor, ← hN]
    have hθ0 : 0 < σ₀ + ε := by linarith
    have hy2x : y ≤ 2 * x := by rw [hy]; linarith
    have hpow : y ^ (σ₀ + ε) ≤ 2 * x ^ (σ₀ + ε) := by
      calc y ^ (σ₀ + ε) ≤ (2 * x) ^ (σ₀ + ε) := Real.rpow_le_rpow (by linarith) hy2x hθ0.le
        _ = 2 ^ (σ₀ + ε) * x ^ (σ₀ + ε) := Real.mul_rpow (by norm_num) hx0.le
        _ ≤ 2 * x ^ (σ₀ + ε) := by
            apply mul_le_mul_of_nonneg_right _ (Real.rpow_nonneg hx0.le _)
            calc (2 : ℝ) ^ (σ₀ + ε) ≤ 2 ^ (1 : ℝ) :=
                  Real.rpow_le_rpow_of_exponent_le (by norm_num) hbig.le
              _ = 2 := Real.rpow_one 2
    rw [hMxy]
    calc |mertens y| ≤ D * y ^ (σ₀ + ε) := hhalf y hy2 hfar
      _ ≤ D * (2 * x ^ (σ₀ + ε)) := mul_le_mul_of_nonneg_left hpow hD0
      _ = 2 * D * x ^ (σ₀ + ε) := by ring

end

end Fayed.PNTError
