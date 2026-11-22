import Riemann.Deep.LeanNBRate.Final3
import Riemann.Deep.LeanNBRate.U2final
import Riemann.Deep.LeanNBRate.Assembly

/-!
# The upper bound for Nyman's constrained distance: assembly at a fixed large L, then L → ∞
-/

open MeasureTheory Filter Topology Set Complex
open scoped Real

noncomputable section

namespace Riemann.Deep.LeanNBRate

theorem QK_nonneg (H : ℂ → ℂ) {L T₂ B : ℝ} (hB0 : 0 ≤ B) (hT : 0 ≤ T₂) : 0 ≤ QK H L T₂ B := by
  unfold QK
  have h1 := Bc_nonneg hB0
  have h2 := ZU_nonneg H T₂
  have h3 := ZL_nonneg H T₂
  have h4 : 0 ≤ Real.exp (L * (3 / 2)) * (8 + Bc B * (64 * (T₂ + 3) + 24)) := by
    have : 0 ≤ Bc B * (64 * (T₂ + 3) + 24) := mul_nonneg h1 (by linarith)
    have := Real.exp_pos (L * (3 / 2))
    positivity
  have h5 : (0 : ℝ) ≤ (2 * π)⁻¹ := by positivity
  exact mul_nonneg (mul_nonneg h5 h4) (by linarith)

set_option maxHeartbeats 4000000 in
/-- **The estimate at a fixed large L.** -/
theorem bound_at (hP : Rudin1987_Thm9_13) (h215 : Titchmarsh_2_1_5) (hRH : RiemannHypothesis)
    {θ : ℝ} (hθ : 0 < θ) (hθ1 : θ ≤ 1) {ℓ₀ : ℝ}
    (hU2 : ∀ ℓ ≥ ℓ₀, ∀ (T : ℝ) (P : ExpSum), T ≤ ℓ ^ 7 → ZeroConfig P T →
      ‖P.eval 0‖ ^ 2 ≤ (zeroSum + θ) / ℓ * (oubil ℓ P.eval (P.dmod 0).eval P.eval (P.dmod 0).eval).re)
    {A C T₀ : ℝ} (h97 : ∀ T : ℝ, T₀ ≤ T → ∃ t ∈ Icc T (T + 1), ∀ σ ∈ Icc (-1 : ℝ) 2,
      riemannZeta (σ + t * I) ≠ 0 ∧ ‖riemannZeta (σ + t * I)‖⁻¹ ≤ C * t ^ A)
    {L : ℝ} (hL100 : 100 ≤ L) (hLT₀ : T₀ ≤ L ^ 6) (hLℓ₀ : ℓ₀ ≤ L - Real.sqrt L)
    (hL7 : L ^ 6 + 1 ≤ (L - Real.sqrt L) ^ 7) (hL5 : 4 * Real.exp 5 ≤ L ^ 3) (hLp : (pE A : ℝ) ≤ L) :
    L * nymanDist (Real.exp (-L)) ^ 2 ≤ (1 + θ) * (zeroSum + θ) * (L / (L - Real.sqrt L)) +
      3 * (1 + θ⁻¹) * (2 * π)⁻¹ * (8 * π / L + (4 * π + 4 * Zz) * (L * Qb C A L ^ 2) +
        8 * Zz * (L * Rb L ^ 2)) := by
  have hS := zeroSum_nonneg
  have hL0 : 0 < L := by linarith
  have hL1 : 1 ≤ L := by linarith
  have hsL : Real.sqrt L ≤ L / 2 := by
    have h1 : L ≤ (L / 2) ^ 2 := by nlinarith
    calc Real.sqrt L ≤ Real.sqrt ((L / 2) ^ 2) := Real.sqrt_le_sqrt h1
      _ = L / 2 := Real.sqrt_sq (by linarith)
  have hsL10 : 10 ≤ Real.sqrt L := by
    calc (10 : ℝ) = Real.sqrt (10 ^ 2) := by rw [Real.sqrt_sq (by norm_num)]
      _ ≤ Real.sqrt L := Real.sqrt_le_sqrt (by nlinarith)
  obtain ⟨ℓ, hℓdef⟩ : ∃ ℓ, ℓ = L - Real.sqrt L := ⟨_, rfl⟩
  rw [← hℓdef] at hLℓ₀ hL7 ⊢
  have hℓL2 : L / 2 ≤ ℓ := by rw [hℓdef]; linarith
  have hℓ1 : 1 ≤ ℓ := by linarith
  have hℓL : ℓ ≤ L := by rw [hℓdef]; linarith [Real.sqrt_nonneg L]
  obtain ⟨n, hndef⟩ : ∃ n : ℕ, n = ⌈L⌉₊ + 7 := ⟨_, rfl⟩
  have hn1 : L + 7 ≤ (n : ℝ) := by rw [hndef]; push_cast; linarith [Nat.le_ceil L]
  have hn2 : (n : ℝ) ≤ L + 8 := by rw [hndef]; push_cast; linarith [Nat.ceil_lt_add_one hL0.le]
  have hn7 : 7 ≤ n := by rw [hndef]; omega
  have hnpos : (0 : ℝ) < n := by linarith
  have hT₁pos : (0 : ℝ) < L ^ 2 := by positivity
  obtain ⟨a, hadef⟩ : ∃ a : ℝ, a = 2 / (n * L ^ 2) := ⟨_, rfl⟩
  have ha : 0 < a := by rw [hadef]; positivity
  have hna : (n : ℝ) * a = 2 / L ^ 2 := by rw [hadef]; field_simp
  have hDa : 2 / a = n * L ^ 2 := by rw [hadef]; field_simp
  have h2L : 2 / L ^ 2 ≤ 1 := by rw [div_le_one hT₁pos]; nlinarith
  have h2L0 : 0 ≤ 2 / L ^ 2 := by positivity
  obtain ⟨T₂, hT₂mem, hT₂good⟩ := h97 (L ^ 6) hLT₀
  have hT₂lo : L ^ 6 ≤ T₂ := hT₂mem.1
  have hT₂hi : T₂ ≤ L ^ 6 + 1 := hT₂mem.2
  have hL6 : 1 ≤ L ^ 6 := one_le_pow₀ hL1
  have hT₂1 : 1 ≤ T₂ := by linarith
  have hT₂pos : 0 < T₂ := by linarith
  obtain ⟨B, hBdef⟩ : ∃ B, B = C * T₂ ^ A := ⟨_, rfl⟩
  have hBu : ∀ x ∈ Icc (-1 : ℝ) 2, ‖(riemannZeta (x + T₂ * I))⁻¹‖ ≤ B := fun x hx => by
    rw [norm_inv, hBdef]; exact (hT₂good x hx).2
  have hBl : ∀ x ∈ Icc (-1 : ℝ) 2, ‖(riemannZeta (x + ((-T₂ : ℝ) : ℂ) * I))⁻¹‖ ≤ B := fun x hx => by
    have e : (x : ℂ) + ((-T₂ : ℝ) : ℂ) * I = (starRingEnd ℂ) ((x : ℂ) + T₂ * I) := by
      apply Complex.ext <;> simp
    rw [e, riemannZeta_conj, ← map_inv₀, Complex.norm_conj]; exact hBu x hx
  have hB0 : 0 ≤ B := (norm_nonneg _).trans (hBu 0 ⟨by norm_num, by norm_num⟩)
  have hgood : ∀ σ ∈ Icc (-1 : ℝ) 2, riemannZeta (σ + T₂ * I) ≠ 0 := fun σ hσ => (hT₂good σ hσ).1
  obtain ⟨c, hc⟩ := exists_representer (cfg T₂) (ℓ := ℓ) (by linarith)
  have hT₂7 : T₂ ≤ ℓ ^ 7 := hT₂hi.trans hL7
  have hNp' : Nself ℓ ((cfg T₂).withCoef c).eval (((cfg T₂).withCoef c).dmod 0).eval ≤ (zeroSum + θ) / ℓ :=
    representer_bound (by linarith) (cfg_zeroConfig T₂) (fun Q hQ => hU2 ℓ hLℓ₀ T₂ Q hT₂7 hQ) (by linarith) hc
  have hHZ : HZeros (Hc ℓ a n ((cfg T₂).withCoef c)) T₂ := hZeros_Hc hRH hT₂pos (by linarith) hgood hc
  obtain ⟨r, hrdef⟩ : ∃ r, r = (cfg T₂).withCoef c := ⟨_, rfl⟩
  rw [← hrdef] at hNp' hHZ
  obtain ⟨Np, hNpdef⟩ : ∃ Np, Np = Nself ℓ r.eval (r.dmod 0).eval := ⟨_, rfl⟩
  have hNp : Np ≤ (zeroSum + θ) / ℓ := by rw [hNpdef]; exact hNp'
  have hNp0 : 0 ≤ Np := by rw [hNpdef]; exact Nself_nonneg (by linarith) _ _
  have hNp1 : Np ≤ zeroSum + 1 := by
    have : (zeroSum + θ) / ℓ ≤ zeroSum + θ := div_le_self (by linarith) hℓ1
    linarith
  have hHB := hBound_Hc (by linarith : 0 ≤ ℓ) ha hn7 r
  rw [← hNpdef] at hHB
  obtain ⟨s₁, hs₁⟩ : ∃ s₁, s₁ = Real.sqrt (zeroSum + 1) := ⟨_, rfl⟩
  have hs₁0 : 0 ≤ s₁ := by rw [hs₁]; exact Real.sqrt_nonneg _
  have hsNp : Real.sqrt Np ≤ s₁ := by rw [hs₁]; exact Real.sqrt_le_sqrt hNp1
  have hsℓ : Real.sqrt (ℓ + 1) ≤ L := by
    calc Real.sqrt (ℓ + 1) ≤ Real.sqrt (L ^ 2) := Real.sqrt_le_sqrt (by nlinarith)
      _ = L := Real.sqrt_sq hL0.le
  have hM : 2 ^ n * (1 + Real.sqrt (ℓ + 1) * Real.sqrt Np * (2 / a + 1)) ≤ 2 ^ n * (1 + 2 * s₁) * L ^ 4 := by
    rw [hDa]
    have h1 : (n : ℝ) * L ^ 2 + 1 ≤ 2 * L ^ 3 := by nlinarith
    have h2 : Real.sqrt (ℓ + 1) * Real.sqrt Np * (n * L ^ 2 + 1) ≤ L * s₁ * (2 * L ^ 3) :=
      mul_le_mul (mul_le_mul hsℓ hsNp (Real.sqrt_nonneg _) hL0.le) h1 (by positivity) (by positivity)
    have h3 : 1 ≤ L ^ 4 := one_le_pow₀ hL1
    have h4 : (1 : ℝ) + L * s₁ * (2 * L ^ 3) ≤ (1 + 2 * s₁) * L ^ 4 := by nlinarith
    have h5 : (0 : ℝ) < 2 ^ n := by positivity
    calc 2 ^ n * (1 + Real.sqrt (ℓ + 1) * Real.sqrt Np * (n * L ^ 2 + 1))
        ≤ 2 ^ n * ((1 + 2 * s₁) * L ^ 4) := mul_le_mul_of_nonneg_left (by linarith) h5.le
      _ = _ := by ring
  have hβ : ℓ + n * a ≤ L := by rw [hna, hℓdef]; linarith
  have hq : (2 / a) / ((2 / a) + T₂) ≤ 2 / L ^ 3 := by
    rw [hDa, div_le_div_iff₀ (by positivity) (by positivity)]
    have h1 : (n : ℝ) * L ^ 5 ≤ 2 * L ^ 6 := by
      have h6 : (n : ℝ) ≤ 2 * L := by linarith
      have h5 : 0 ≤ L ^ 5 := by positivity
      calc (n : ℝ) * L ^ 5 ≤ 2 * L * L ^ 5 := mul_le_mul_of_nonneg_right h6 h5
        _ = 2 * L ^ 6 := by ring
    have h2 : 0 ≤ (n : ℝ) * L ^ 2 := by positivity
    have e : (n : ℝ) * L ^ 2 * L ^ 3 = n * L ^ 5 := by ring
    rw [e]
    linarith
  have hD3 : 3 + 2 * (2 / a) ≤ 5 * L ^ 3 := by
    rw [hDa]
    have h1 : (n : ℝ) * L ^ 2 ≤ (L + 8) * L ^ 2 := mul_le_mul_of_nonneg_right hn2 (sq_nonneg L)
    have h2 : 100 * L ^ 2 ≤ L ^ 3 := by
      have := mul_le_mul_of_nonneg_left hL100 (sq_nonneg L)
      nlinarith
    nlinarith
  have hZ := Z_le hHB hT₂pos hL1 hM hβ hq hD3 (by omega : n - 1 - 2 = n - 3) (by omega)
  have hBc := Bc_le hBdef hT₂1 (by linarith : T₂ ≤ 2 * L ^ 6) hL1
  have h64 : 64 * (T₂ + 3) + 24 ≤ 200 * L ^ 6 := by
    have : L ≤ L ^ 6 := le_self_pow₀ hL1 (by norm_num)
    linarith
  have hn3 : L + 3 ≤ ((n - 3 : ℕ) : ℝ) := by rw [Nat.cast_sub (by omega)]; push_cast; linarith
  have hQ := QK_le (H := Hc ℓ a n r) (by linarith : 2 ≤ L) hL5 hLp hn3 hB0 hBc h64 hZ hs₁ hT₂pos.le
  obtain ⟨M', hM'⟩ : ∃ M', M' = 1 + Real.sqrt (ℓ + 1) * Real.sqrt Np := ⟨_, rfl⟩
  obtain ⟨β', hβ'⟩ : ∃ β', β' = ℓ + 1 + n * a := ⟨_, rfl⟩
  have htz : ∀ k : ℕ, ‖Hc ℓ a n r (tz k)‖ ≤ M' * Real.exp (β' * (2 * k + 5 / 2)) := fun k => by
    have := norm_Hc_tz (by linarith : 0 ≤ ℓ) ha r k (n := n)
    rw [← hNpdef] at this; rw [hM', hβ']; exact this
  have hβ'L : β' ≤ L := by rw [hβ', hna, hℓdef]; linarith
  obtain ⟨Rt, hRtdef⟩ : ∃ Rt, Rt = M' * Real.exp (5 / 2 * (β' - L)) * (12 / 5 * Real.exp (2 * π)) :=
    ⟨_, rfl⟩
  have hRt : ∀ s : ℂ, 1 / 2 ≤ s.re → ‖Rtr (Hc ℓ a n r) L s‖ ≤ Rt := fun s hs => by
    rw [hRtdef]; exact norm_Rtr_le hL0.le hβ'L htz hs
  have hRtb : Rt ≤ Rb L := by
    rw [hRtdef]; unfold Rb KR; rw [← hs₁]
    have h0 := mul_le_mul hsℓ hsNp (Real.sqrt_nonneg _) hL0.le
    have h1 : M' ≤ (1 + s₁) * L := by rw [hM']; nlinarith
    have h2 : Real.exp (5 / 2 * (β' - L)) ≤ Real.exp 5 * Real.exp (-(5 / 2 * Real.sqrt L)) := by
      rw [← Real.exp_add, Real.exp_le_exp, hβ', hna, hℓdef]; linarith
    have hM'0 : 0 ≤ M' := by rw [hM']; positivity
    have hc' : 0 ≤ 12 / 5 * Real.exp (2 * π) := by positivity
    calc M' * Real.exp (5 / 2 * (β' - L)) * (12 / 5 * Real.exp (2 * π))
        ≤ ((1 + s₁) * L) * (Real.exp 5 * Real.exp (-(5 / 2 * Real.sqrt L))) * (12 / 5 * Real.exp (2 * π)) :=
          mul_le_mul_of_nonneg_right (mul_le_mul h1 h2 (Real.exp_pos _).le (by positivity)) hc'
      _ = _ := by ring
  have hRt0 : 0 ≤ Rt := (norm_nonneg _).trans (hRt (1 / 2 : ℂ) (by norm_num))
  have hlog : Real.log (2 * π) ≤ 7 :=
    (Real.log_le_sub_one_of_pos (by positivity)).trans (by linarith [Real.pi_le_four])
  have hLβ : (ℓ + n * a) + Real.log (2 * π) < L := by rw [hna, hℓdef]; linarith
  have hmain := nyman_le hP h215 r (by linarith : 0 ≤ ℓ) ha hT₁pos hna hL0 hHB hHZ hT₂1 hLβ hB0 hBu hBl hRt hθ
  rw [← hNpdef] at hmain
  obtain ⟨Q, hQdef⟩ : ∃ Q, Q = QK (Hc ℓ a n r) L T₂ B := ⟨_, rfl⟩
  obtain ⟨G₁, hG₁def⟩ : ∃ G₁, G₁ = G1b (Hc ℓ a n r) L T₂ B Rt := ⟨_, rfl⟩
  rw [← hQdef, ← hG₁def] at hmain
  rw [← hQdef] at hQ
  have hQ0 : 0 ≤ Q := by rw [hQdef]; exact QK_nonneg _ hB0 hT₂pos.le
  have hZ0 : 0 ≤ ZU (Hc ℓ a n r) T₂ + ZL (Hc ℓ a n r) T₂ := add_nonneg (ZU_nonneg _ _) (ZL_nonneg _ _)
  have hBc0 := Bc_nonneg hB0
  have hG₁le : G₁ ≤ Q + Rt := by
    rw [hG₁def, hQdef]; unfold G1b QK
    have h1 : Bc B ≤ 8 + Bc B * (64 * (T₂ + 3) + 24) := by nlinarith
    have h2 : (2 * π)⁻¹ * (Real.exp (L * (3 / 2)) * Bc B) * (ZU (Hc ℓ a n r) T₂ + ZL (Hc ℓ a n r) T₂) ≤
        (2 * π)⁻¹ * (Real.exp (L * (3 / 2)) * (8 + Bc B * (64 * (T₂ + 3) + 24))) *
          (ZU (Hc ℓ a n r) T₂ + ZL (Hc ℓ a n r) T₂) := by
      have hpi : (0 : ℝ) ≤ (2 * π)⁻¹ := by positivity
      exact mul_le_mul_of_nonneg_right (mul_le_mul_of_nonneg_left
        (mul_le_mul_of_nonneg_left h1 (Real.exp_pos _).le) hpi) hZ0
    linarith
  have hG₁0 : 0 ≤ G₁ := by
    rw [hG₁def]; unfold G1b
    have hpi : (0 : ℝ) ≤ (2 * π)⁻¹ := by positivity
    have := mul_nonneg (mul_nonneg hpi (mul_nonneg (Real.exp_pos (L * (3 / 2))).le hBc0)) hZ0
    linarith
  have hZz : 0 ≤ Zz := integral_nonneg fun t => by positivity
  have hLNp : L * Np ≤ (zeroSum + θ) * (L / ℓ) := by
    have h1 := (le_div_iff₀ (by linarith : (0 : ℝ) < ℓ)).mp hNp
    rw [mul_div_assoc', le_div_iff₀ (by linarith)]
    nlinarith
  have hQ2 : Q ^ 2 ≤ Qb C A L ^ 2 := pow_le_pow_left₀ hQ0 hQ 2
  have hR2 : (Rt + G₁) ^ 2 ≤ 2 * Qb C A L ^ 2 + 8 * Rb L ^ 2 := by
    have h1 : Rt + G₁ ≤ Qb C A L + 2 * Rb L := by linarith
    have h2 := pow_le_pow_left₀ (by linarith) h1 2
    nlinarith [sq_nonneg (Qb C A L - 2 * Rb L)]
  have hY : L * (8 * π / L ^ 2 + 4 * π * Q ^ 2 + (2 * Q ^ 2 + (Rt + G₁) ^ 2) * Zz) ≤
      8 * π / L + (4 * π + 4 * Zz) * (L * Qb C A L ^ 2) + 8 * Zz * (L * Rb L ^ 2) := by
    have e1 : L * (8 * π / L ^ 2) = 8 * π / L := by field_simp
    have h3 : (2 * Q ^ 2 + (Rt + G₁) ^ 2) * Zz ≤ (4 * Qb C A L ^ 2 + 8 * Rb L ^ 2) * Zz :=
      mul_le_mul_of_nonneg_right (by linarith) hZz
    have h4 : 4 * π * Q ^ 2 ≤ 4 * π * Qb C A L ^ 2 := mul_le_mul_of_nonneg_left hQ2 (by positivity)
    have h5 := mul_le_mul_of_nonneg_left h3 hL0.le
    have h6 := mul_le_mul_of_nonneg_left h4 hL0.le
    calc L * (8 * π / L ^ 2 + 4 * π * Q ^ 2 + (2 * Q ^ 2 + (Rt + G₁) ^ 2) * Zz)
        = L * (8 * π / L ^ 2) + L * (4 * π * Q ^ 2) + L * ((2 * Q ^ 2 + (Rt + G₁) ^ 2) * Zz) := by ring
      _ ≤ 8 * π / L + L * (4 * π * Qb C A L ^ 2) + L * ((4 * Qb C A L ^ 2 + 8 * Rb L ^ 2) * Zz) := by
          rw [e1]; linarith
      _ = _ := by ring
  have hc₁ : 0 ≤ 3 * (1 + θ⁻¹) * (2 * π)⁻¹ := by positivity
  calc L * nymanDist (Real.exp (-L)) ^ 2
      ≤ L * ((1 + θ) * Np + 3 * (1 + θ⁻¹) * (2 * π)⁻¹ *
          (8 * π / L ^ 2 + 4 * π * Q ^ 2 + (2 * Q ^ 2 + (Rt + G₁) ^ 2) * Zz)) :=
        mul_le_mul_of_nonneg_left hmain hL0.le
    _ = (1 + θ) * (L * Np) + 3 * (1 + θ⁻¹) * (2 * π)⁻¹ *
          (L * (8 * π / L ^ 2 + 4 * π * Q ^ 2 + (2 * Q ^ 2 + (Rt + G₁) ^ 2) * Zz)) := by ring
    _ ≤ (1 + θ) * ((zeroSum + θ) * (L / ℓ)) + 3 * (1 + θ⁻¹) * (2 * π)⁻¹ *
          (8 * π / L + (4 * π + 4 * Zz) * (L * Qb C A L ^ 2) + 8 * Zz * (L * Rb L ^ 2)) :=
        add_le_add (mul_le_mul_of_nonneg_left hLNp (by linarith)) (mul_le_mul_of_nonneg_left hY hc₁)
    _ = _ := by ring

end Riemann.Deep.LeanNBRate

