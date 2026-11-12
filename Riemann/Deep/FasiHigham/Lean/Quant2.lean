import Riemann.Deep.FasiHigham.Lean.Quant

/-!
# Fasi–Higham, m = 3: quantitative theorem (statement (2) of the brief, real form)
-/

namespace FasiLean

/-- From (x − Rw)² ≤ K w², K ≤ r², r ≥ 0, w > 0: (R − r) w ≤ x. -/
theorem lower_of_bound (x R K w r : ℝ) (hw : 0 < w) (hr0 : 0 ≤ r) (hKr : K ≤ r ^ 2)
    (hb : (x - R * w) ^ 2 ≤ K * w ^ 2) : (R - r) * w ≤ x := by
  rcases le_or_gt ((R - r) * w) x with h | h
  · exact h
  · exfalso
    have hrw : 0 ≤ r * w := mul_nonneg hr0 hw.le
    have h1 : r * w < R * w - x := by linarith
    have h2 : (r * w) * (r * w) < (R * w - x) * (R * w - x) := mul_self_lt_mul_self hrw h1
    have h3 : K * w ^ 2 ≤ r ^ 2 * w ^ 2 := mul_le_mul_of_nonneg_right hKr (by positivity)
    have e1 : (x - R * w) ^ 2 = (R * w - x) * (R * w - x) := by ring
    have e2 : r ^ 2 * w ^ 2 = (r * w) * (r * w) := by ring
    linarith

/-- QUANTITATIVE THEOREM. For N = k + 1 there are the root d ∈ (1/2, 2/3) of d³ + 3Nd = 2N, the zero
z₀ = N + d of q_{k,3}, the quadratic factor z² − 2uz + v of the complex pair (u² < v = |z_±|², z₀² < v),
the residue R₀ = −1/(z₀ q'(z₀)) > 0 and K = R₀² · 8N(N+d)²(1−d)²/(d³[(N+d)² + 3N(1−d)][3N+d²+d+1])
(NOTE.md §2's (2|R₊|)²) with 0 ≤ K < R₀², such that (c_n − R₀ z₀^{−n})² ≤ K z₀^{−2n} for all n;
hence c_n ≥ (R₀ − r) z₀^{−n} for every r ≥ 0 with K ≤ r² (r = √K = 2|R₊| < R₀). -/
theorem c_quant (k : ℕ) (N : ℝ) (hN : N = k + 1) :
    ∃ d z₀ u v R₀ K : ℝ,
      (1 / 2 < d ∧ d < 2 / 3 ∧ d ^ 3 + 3 * N * d = 2 * N) ∧
      z₀ = N + d ∧ qR k z₀ = 0 ∧
      (∀ z, N * (N + 1) * (N + 2) * qR k z = (z₀ - z) * (z ^ 2 - 2 * u * z + v)) ∧
      u ^ 2 < v ∧ z₀ ^ 2 < v ∧
      R₀ = -1 / (z₀ * qR' k z₀) ∧ 0 < R₀ ∧
      K = R₀ ^ 2 * (8 * N * (N + d) ^ 2 * (1 - d) ^ 2
        / (d ^ 3 * ((N + d) ^ 2 + 3 * N * (1 - d)) * (3 * N + d ^ 2 + d + 1))) ∧
      0 ≤ K ∧ K < R₀ ^ 2 ∧
      (∀ n, ((c k n : ℝ) - R₀ * (1 / z₀) ^ n) ^ 2 ≤ K * ((1 / z₀) ^ n) ^ 2) ∧
      (∀ n r, 0 ≤ r → K ≤ r ^ 2 → (R₀ - r) * (1 / z₀) ^ n ≤ (c k n : ℝ)) := by
  subst hN
  have hN1 : (1 : ℝ) ≤ (k : ℝ) + 1 := by
    have := Nat.cast_nonneg (α := ℝ) k
    linarith
  obtain ⟨d, hd1, hd2, hr⟩ := exists_root ((k : ℝ) + 1) hN1
  obtain ⟨hR, hK, hb⟩ := main_real k d hd1 hd2 hr
  have hN0 : (0 : ℝ) < (k : ℝ) + 1 := by linarith
  have hP0 : ((k : ℝ) + 1) * ((k : ℝ) + 1 + 1) * ((k : ℝ) + 1 + 2) ≠ 0 := by positivity
  have hv0 : 0 < vR ((k : ℝ) + 1) d := by
    have : vR ((k : ℝ) + 1) d = ((k : ℝ) + 1 - d / 2) ^ 2 + 3 / 4 * d ^ 2 + 3 * ((k : ℝ) + 1) := by
      simp only [vR]; ring
    rw [this]; positivity
  have hvu : vR ((k : ℝ) + 1) d - uR ((k : ℝ) + 1) d ^ 2 = 3 / 4 * d ^ 2 + 3 * ((k : ℝ) + 1) := by
    simp only [vR, uR]; ring
  have fact : ∀ z, ((k : ℝ) + 1) * ((k : ℝ) + 1 + 1) * ((k : ℝ) + 1 + 2) * qR k z
      = ((k : ℝ) + 1 + d - z) * (z ^ 2 - 2 * uR ((k : ℝ) + 1) d * z + vR ((k : ℝ) + 1) d) := by
    intro z
    rw [qR_eq k _ z rfl]
    simp only [uR, vR]
    linear_combination (-1 : ℝ) * hr
  refine ⟨d, (k : ℝ) + 1 + d, uR _ d, vR _ d, R0 _ d, KR _ d, ⟨hd1, hd2, hr⟩, rfl, ?_, fact, ?_, ?_,
    ?_, hR, ?_, ?_, hK, hb, ?_⟩
  · have h0 := fact ((k : ℝ) + 1 + d)
    rw [sub_self, zero_mul] at h0
    exact (mul_eq_zero.mp h0).resolve_left hP0
  · have : 0 < 3 / 4 * d ^ 2 + 3 * ((k : ℝ) + 1) := by positivity
    linarith
  · have e : vR ((k : ℝ) + 1) d - ((k : ℝ) + 1 + d) ^ 2 = 3 * ((k : ℝ) + 1) * (1 - d) := by
      simp only [vR]; ring
    have : 0 < 3 * ((k : ℝ) + 1) * (1 - d) := mul_pos (by linarith) (by linarith)
    linarith
  · have h' := qR'_eq k ((k : ℝ) + 1) ((k : ℝ) + 1 + d) rfl
    have hq' : qR' k ((k : ℝ) + 1 + d)
        = -3 * (d ^ 2 + ((k : ℝ) + 1)) / (((k : ℝ) + 1) * ((k : ℝ) + 1 + 1) * ((k : ℝ) + 1 + 2)) := by
      rw [eq_div_iff hP0, mul_comm, h']
      ring
    have hdN : 0 < d ^ 2 + ((k : ℝ) + 1) := by positivity
    rw [hq']
    simp only [R0, vR]
    field_simp
    linear_combination hr
  · simp only [KR]
    rw [mul_div_assoc, K_note _ _ hN1 hd1 hd2 hr]
  · simp only [KR]
    apply div_nonneg (by positivity)
    rw [hvu]
    positivity
  · intro n r hr0 hKr
    exact lower_of_bound _ _ _ _ r (pow_pos (by positivity) n) hr0 hKr (hb n)

end FasiLean
