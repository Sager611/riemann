import Riemann.Deep.FasiHigham.Lean.Core
import Riemann.Deep.FasiHigham.Lean.KeyIneq

/-!
# Fasi–Higham, m = 3: the bound (c_n − R ρ^n)² ≤ K ρ^{2n} with K < R²
-/

namespace FasiLean

/-- The bound, from the decomposition and the conserved form. a = z₀ = N + d, R = R₀,
K = R² · 3a²(d² + N)/(v(v − u²)) (= (2|R₊|)², see K_note). -/
theorem core (N d a u v R : ℝ) (C : ℕ → ℝ) (hN : 1 ≤ N) (hd1 : 1 / 2 < d) (hd2 : d < 2 / 3)
    (hr : d ^ 3 + 3 * N * d = 2 * N) (ha : a = N + d) (hu : u = N - d / 2)
    (hv : v = N ^ 2 - d * N + d ^ 2 + 3 * N) (hR : R = v / (3 * (d ^ 2 + N)))
    (c0 : C 0 = 0) (c1 : C 1 = 0) (c2 : C 2 = 1)
    (hC : ∀ m, N * (N + 1) * (N + 2) * C (m + 3)
      = 3 * N * (N + 1) * C (m + 2) - 3 * N * C (m + 1) + C m) :
    0 < R ∧ R ^ 2 * (3 * a ^ 2 * (d ^ 2 + N)) / (v * (v - u ^ 2)) < R ^ 2 ∧
      ∀ n, (C (n + 2) - R * (1 / a) ^ n) ^ 2
        ≤ R ^ 2 * (3 * a ^ 2 * (d ^ 2 + N)) / (v * (v - u ^ 2)) * ((1 / a) ^ n) ^ 2 := by
  have hN0 : 0 < N := by linarith
  have ha0 : 0 < a := by rw [ha]; linarith
  have hv0 : 0 < v := by
    have : v = (N - d / 2) ^ 2 + 3 / 4 * d ^ 2 + 3 * N := by rw [hv]; ring
    rw [this]; positivity
  have hvu : v - u ^ 2 = 3 / 4 * d ^ 2 + 3 * N := by rw [hv, hu]; ring
  have hvu0 : 0 < v - u ^ 2 := by rw [hvu]; positivity
  have hva : a ^ 2 ≤ v := by
    have e : v - a ^ 2 = 3 * N * (1 - d) := by rw [hv, ha]; ring
    have : 0 ≤ 3 * N * (1 - d) := mul_nonneg (by linarith) (by linarith)
    linarith
  have hS : a ^ 2 - 2 * u * a + v = 3 * (d ^ 2 + N) := by rw [ha, hu, hv]; ring
  have hdN : 0 < d ^ 2 + N := by positivity
  have hR0 : 0 < R := by rw [hR]; positivity
  have hRS : R * (a ^ 2 - 2 * u * a + v) = v := by rw [hS, hR]; field_simp
  have hP : N * (N + 1) * (N + 2) = a * v := by rw [ha, hv]; linear_combination (-1 : ℝ) * hr
  have h2 : 3 * N * (N + 1) = 2 * u * a + v := by rw [hu, ha, hv]; ring
  have h3 : 3 * N = a + 2 * u := by rw [ha, hu]; ring
  have hρ : 1 / a * a = 1 := one_div_mul_cancel ha0.ne'
  have hP0 : N * (N + 1) * (N + 2) ≠ 0 := by positivity
  have dec := decomp N a u v R (1 / a) C hP h2 h3 hRS hρ hv0.ne' hP0 c0 c1 c2 hC
  have hf := fs_rec a u v R hv0.ne'
  have G0 : Gf u v (fs a u v R) 0 = R ^ 2 * a ^ 2 * (a ^ 2 - 2 * u * a + v) := by
    show v * (-(R * a)) ^ 2 - 2 * u * (-(R * a)) * (-(R * a ^ 2)) + (-(R * a ^ 2)) ^ 2 = _
    ring
  have key : a ^ 2 * (a ^ 2 - 2 * u * a + v) < v * (v - u ^ 2) := by
    rw [ha, hu, hv]; exact key_ineq N d hN hd1 hd2
  have hK3 : R ^ 2 * a ^ 2 * (a ^ 2 - 2 * u * a + v) = R ^ 2 * (3 * a ^ 2 * (d ^ 2 + N)) := by
    rw [hS]; ring
  refine ⟨hR0, ?_, ?_⟩
  · rw [div_lt_iff₀ (by positivity)]
    have : 3 * a ^ 2 * (d ^ 2 + N) < v * (v - u ^ 2) := by
      have e : 3 * a ^ 2 * (d ^ 2 + N) = a ^ 2 * (a ^ 2 - 2 * u * a + v) := by rw [hS]; ring
      rw [e]; exact key
    exact mul_lt_mul_of_pos_left this (by positivity)
  · intro n
    have hb := sq_le_Gf u v (fs a u v R) (n + 2)
    have hg := Gf_pow u v (fs a u v R) hf (n + 2)
    rw [G0] at hg
    have hV : 0 < v ^ n := pow_pos hv0 n
    have h1v : 1 ≤ (1 / a) ^ 2 * v := by
      rw [div_pow, one_pow, one_div_mul_eq_div, le_div_iff₀ (by positivity), one_mul]
      exact hva
    have hwV : 1 ≤ ((1 / a) ^ n) ^ 2 * v ^ n := by
      have : ((1 / a) ^ n) ^ 2 * v ^ n = ((1 / a) ^ 2 * v) ^ n := by ring
      rw [this]; exact one_le_pow₀ h1v
    have h1 : (v - u ^ 2) * fs a u v R (n + 2) ^ 2 * (v * v ^ n)
        ≤ R ^ 2 * (3 * a ^ 2 * (d ^ 2 + N)) := by
      have := mul_le_mul_of_nonneg_right hb (le_of_lt (mul_pos hv0 hV))
      calc (v - u ^ 2) * fs a u v R (n + 2) ^ 2 * (v * v ^ n)
          ≤ v * Gf u v (fs a u v R) (n + 2) * (v * v ^ n) := this
        _ = Gf u v (fs a u v R) (n + 2) * v ^ (n + 2) := by ring
        _ = _ := by rw [hg, hK3]
    have hK : fs a u v R (n + 2) ^ 2 * v ^ n
        ≤ R ^ 2 * (3 * a ^ 2 * (d ^ 2 + N)) / (v * (v - u ^ 2)) := by
      rw [le_div_iff₀ (by positivity)]
      calc fs a u v R (n + 2) ^ 2 * v ^ n * (v * (v - u ^ 2))
          = (v - u ^ 2) * fs a u v R (n + 2) ^ 2 * (v * v ^ n) := by ring
        _ ≤ _ := h1
    have ha2 : a ^ 2 * (1 / a) ^ 2 = 1 := by field_simp
    have hCf : C (n + 2) - R * (1 / a) ^ n = fs a u v R (n + 2) := by
      rw [dec (n + 2)]
      linear_combination (R * (1 / a) ^ n) * ha2
    rw [hCf]
    calc fs a u v R (n + 2) ^ 2 ≤ fs a u v R (n + 2) ^ 2 * (((1 / a) ^ n) ^ 2 * v ^ n) :=
          le_mul_of_one_le_right (sq_nonneg _) hwV
      _ = ((1 / a) ^ n) ^ 2 * (fs a u v R (n + 2) ^ 2 * v ^ n) := by ring
      _ ≤ ((1 / a) ^ n) ^ 2 * (R ^ 2 * (3 * a ^ 2 * (d ^ 2 + N)) / (v * (v - u ^ 2))) :=
          mul_le_mul_of_nonneg_left hK (sq_nonneg _)
      _ = _ := by ring

end FasiLean
