import Riemann.Deep.SunComplete.Lean.Arith

/-!
# sunlean · L6Arith: the rational inequalities of the no-sieve case of Lemma 6

With `u = 2^r·s`, `K = (r+1)u` and `7K ≤ n = q − 1`, (CS) holds with `A' = 1/τ + E`, `D = nθ`.
-/

namespace SunLean

/-- Core: `n ≥ 7K`, `K = r₁u`, `1 ≤ E ≤ u`, `0 ≤ A ≤ 2u`, `n ≤ r₁D` give (CS). -/
theorem noSieve_core {n q D A E u r1 : ℚ} (hq : q = n + 1) (hn : 0 < n) (hr1 : 1 ≤ r1)
    (hK7 : 7 * (r1 * u) ≤ n) (hu1 : 1 ≤ u) (hE1 : 1 ≤ E) (hEu : E ≤ u)
    (hA2 : A ≤ 2 * u) (hA0 : 0 ≤ A) (hD : n ≤ r1 * D) (hD0 : 0 < D) :
    A < D ∧ 2 * A < q ∧ E ^ 2 * D < (D - A) * (q - 1) ^ 2 * (1 - 2 * A / q) := by
  have hr0 : 0 < r1 := by linarith
  have hu0 : 0 < u := by linarith
  have huK : u ≤ r1 * u := by
    have : 0 ≤ (r1 - 1) * u := mul_nonneg (by linarith) hu0.le
    linarith
  have hK1 : 1 ≤ r1 * u := le_trans hu1 huK
  have hD7 : 7 * u ≤ D := by
    have h1 : r1 * (7 * u) ≤ r1 * D := by linarith
    exact le_of_mul_le_mul_left h1 hr0
  have hAD : A < D := by linarith
  have h2A : 2 * A < q := by linarith
  refine ⟨hAD, h2A, ?_⟩
  have hq0 : 0 < q := by linarith
  have hY : (1 - 2 * A / q) * q = q - 2 * A := by field_simp
  have hnq : q - 1 = n := by linarith
  rw [hnq]
  -- multiply by q
  apply lt_of_mul_lt_mul_right _ hq0.le
  rw [mul_assoc ((D - A) * n ^ 2), hY]
  -- lower bound of the right side
  have hX : D * (n - 2 * (r1 * u)) ≤ n * (D - A) := by
    have h1 : n * A ≤ n * (2 * u) := mul_le_mul_of_nonneg_left hA2 hn.le
    have h2 : n * (2 * u) ≤ (r1 * D) * (2 * u) := mul_le_mul_of_nonneg_right hD (by linarith)
    linarith
  have hZ : n - 4 * (r1 * u) ≤ q - 2 * A := by linarith
  have hX0 : 0 ≤ D * (n - 2 * (r1 * u)) := mul_nonneg hD0.le (by linarith)
  have hZ0 : 0 ≤ n - 4 * (r1 * u) := by linarith
  have hR : D * (n - 2 * (r1 * u)) * n * (n - 4 * (r1 * u)) ≤ (D - A) * n ^ 2 * (q - 2 * A) := by
    have h1 : D * (n - 2 * (r1 * u)) * (n - 4 * (r1 * u)) ≤ n * (D - A) * (q - 2 * A) :=
      mul_le_mul hX hZ hZ0 (by nlinarith [hX, hX0])
    have h2 := mul_le_mul_of_nonneg_right h1 hn.le
    have e1 : D * (n - 2 * (r1 * u)) * n * (n - 4 * (r1 * u)) =
        D * (n - 2 * (r1 * u)) * (n - 4 * (r1 * u)) * n := by ring
    have e2 : (D - A) * n ^ 2 * (q - 2 * A) = n * (D - A) * (q - 2 * A) * n := by ring
    rw [e1, e2]; exact h2
  -- upper bound of the left side
  have hE2 : E ^ 2 ≤ (r1 * u) ^ 2 := pow_le_pow_left₀ (by linarith) (le_trans hEu huK) 2
  have hL : E ^ 2 * D * q ≤ (r1 * u) ^ 2 * D * q := by
    have := mul_le_mul_of_nonneg_right hE2 hD0.le
    exact mul_le_mul_of_nonneg_right this hq0.le
  -- the polynomial gap: K²(n + 1) < n(n − 2K)(n − 4K)
  have hP1 : 5 * (r1 * u) ≤ n - 2 * (r1 * u) := by linarith
  have hP2 : 3 * (r1 * u) ≤ n - 4 * (r1 * u) := by linarith
  have hP3 : 15 * (r1 * u) ^ 2 ≤ (n - 2 * (r1 * u)) * (n - 4 * (r1 * u)) := by
    have := mul_le_mul hP1 hP2 (by linarith) (by linarith)
    nlinarith [this]
  have hP4 : (r1 * u) ^ 2 * q < n * ((n - 2 * (r1 * u)) * (n - 4 * (r1 * u))) := by
    have h1 : n * (15 * (r1 * u) ^ 2) ≤ n * ((n - 2 * (r1 * u)) * (n - 4 * (r1 * u))) :=
      mul_le_mul_of_nonneg_left hP3 hn.le
    have hK2 : 0 < (r1 * u) ^ 2 := by positivity
    have hn1 : 1 ≤ n := by linarith
    have h2 : (r1 * u) ^ 2 * q ≤ (r1 * u) ^ 2 * (2 * n) :=
      mul_le_mul_of_nonneg_left (by linarith) hK2.le
    have h3 : (r1 * u) ^ 2 * (2 * n) < n * (15 * (r1 * u) ^ 2) := by nlinarith [hK2, hn]
    linarith
  have hP5 : (r1 * u) ^ 2 * D * q < D * (n - 2 * (r1 * u)) * n * (n - 4 * (r1 * u)) := by
    have := mul_lt_mul_of_pos_left hP4 hD0
    have e1 : (r1 * u) ^ 2 * D * q = D * ((r1 * u) ^ 2 * q) := by ring
    have e2 : D * (n - 2 * (r1 * u)) * n * (n - 4 * (r1 * u)) =
        D * (n * ((n - 2 * (r1 * u)) * (n - 4 * (r1 * u)))) := by ring
    rw [e1, e2]; exact this
  linarith

/-- The no-sieve case of Lemma 6 in rational form (`w = 2^r`, `r₁ = r + 1`, `θ ≥ 1/r₁`,
`E = (w − 1)s + 1`, `A' = 1/θ + E`, `D = nθ`, `0 ≤ ε ≤ 1`). -/
theorem noSieve_arith {n q w r1 s θ ε : ℚ} (hq : q = n + 1) (hn : 0 < n) (hw1 : 1 ≤ w)
    (hr1 : 1 ≤ r1) (hwr : r1 ≤ w) (hs : 1 ≤ s) (hθ : 1 ≤ r1 * θ)
    (hK : 7 * (r1 * (w * s)) ≤ n) (hε0 : 0 ≤ ε) (hε1 : ε ≤ 1) :
    q * (ε / θ + ((w - 1) * s + 1)) / (q - 1 + ((w - 1) * s + 1)) ≤ 1 / θ + ((w - 1) * s + 1) ∧
      1 / θ + ((w - 1) * s + 1) < n * θ ∧ 2 * (1 / θ + ((w - 1) * s + 1)) < q ∧
      ((w - 1) * s + 1) ^ 2 * (n * θ) <
        (n * θ - (1 / θ + ((w - 1) * s + 1))) * (q - 1) ^ 2 *
          (1 - 2 * (1 / θ + ((w - 1) * s + 1)) / q) := by
  have hθ0 : 0 < θ := by
    rcases lt_or_ge 0 θ with h | h
    · exact h
    · exfalso
      have := mul_le_mul_of_nonneg_left h (show (0 : ℚ) ≤ r1 by linarith)
      rw [mul_zero] at this
      linarith
  have hinv : 1 / θ ≤ r1 := by rw [div_le_iff₀ hθ0]; linarith
  have hinv0 : 0 < 1 / θ := by positivity
  have hεθ : ε / θ ≤ 1 / θ := div_le_div_of_nonneg_right hε1 hθ0.le
  have hεθ0 : 0 ≤ ε / θ := div_nonneg hε0 hθ0.le
  have hws : w ≤ w * s := le_mul_of_one_le_right (by linarith) hs
  have hu1 : 1 ≤ w * s := le_trans hw1 hws
  have hE1 : 1 ≤ (w - 1) * s + 1 := by
    have : 0 ≤ (w - 1) * s := mul_nonneg (by linarith) (by linarith)
    linarith
  have hEu : (w - 1) * s + 1 ≤ w * s := by
    have : (w - 1) * s = w * s - s := by ring
    linarith
  have hDn : n ≤ r1 * (n * θ) := by
    have := mul_le_mul_of_nonneg_left hθ hn.le
    have e : n * (r1 * θ) = r1 * (n * θ) := by ring
    linarith
  have hD0 : 0 < n * θ := mul_pos hn hθ0
  obtain ⟨h1, h2, h3⟩ := noSieve_core (A := 1 / θ + ((w - 1) * s + 1)) hq hn hr1 hK hu1 hE1 hEu
    (by linarith) (by linarith) hDn hD0
  refine ⟨?_, h1, h2, h3⟩
  have hpos : 0 < q - 1 + ((w - 1) * s + 1) := by linarith
  rw [div_le_iff₀ hpos]
  have hX : 0 ≤ ε / θ + ((w - 1) * s + 1) := by linarith
  have := mul_le_mul_of_nonneg_left (show q ≤ q - 1 + ((w - 1) * s + 1) by linarith) hX
  have h4 : (ε / θ + ((w - 1) * s + 1)) * (q - 1 + ((w - 1) * s + 1)) ≤
      (1 / θ + ((w - 1) * s + 1)) * (q - 1 + ((w - 1) * s + 1)) :=
    mul_le_mul_of_nonneg_right (by linarith) hpos.le
  linarith

end SunLean
