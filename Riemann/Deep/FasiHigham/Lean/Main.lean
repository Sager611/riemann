import Riemann.Deep.FasiHigham.Lean.Defs
import Riemann.Deep.FasiHigham.Lean.Bound
import Riemann.Deep.FasiHigham.Lean.Root

/-!
# Fasi–Higham, m = 3: main theorem 0 < c_n(k) for all k, n
-/

namespace FasiLean

/-- With N = k + 1 and the root d: v = |z_±|², u = Re z_±, R₀ = v/(3(d² + N)) (the residue at z₀ = N + d),
K = R₀² · 3z₀²(d² + N)/(v(v − u²)) (= (2|R₊|)², Quant.lean). -/
noncomputable def vR (N d : ℝ) : ℝ := N ^ 2 - d * N + d ^ 2 + 3 * N
noncomputable def uR (N d : ℝ) : ℝ := N - d / 2
noncomputable def R0 (N d : ℝ) : ℝ := vR N d / (3 * (d ^ 2 + N))
noncomputable def KR (N d : ℝ) : ℝ :=
  R0 N d ^ 2 * (3 * (N + d) ^ 2 * (d ^ 2 + N)) / (vR N d * (vR N d - uR N d ^ 2))

/-- The shifted coefficients, cast to ℝ, satisfy N(N+1)(N+2) C_{m+3} = 3N(N+1) C_{m+2} − 3N C_{m+1} + C_m
with N = k + 1. -/
theorem cs_rec_real (k m : ℕ) :
    ((k : ℝ) + 1) * ((k : ℝ) + 1 + 1) * ((k : ℝ) + 1 + 2) * (cs k (m + 3) : ℝ)
      = 3 * ((k : ℝ) + 1) * ((k : ℝ) + 1 + 1) * (cs k (m + 2) : ℝ)
        - 3 * ((k : ℝ) + 1) * (cs k (m + 1) : ℝ) + (cs k m : ℝ) := by
  have h1 : (k : ℝ) + 1 ≠ 0 := by positivity
  have h2 : (k : ℝ) + 2 ≠ 0 := by positivity
  have h3 : (k : ℝ) + 3 ≠ 0 := by positivity
  rw [cs_rec]
  push_cast
  field_simp
  ring

/-- From (x − Rw)² ≤ K w², K < R², R, w > 0: x > 0. -/
theorem pos_of_bound (x R K w : ℝ) (hR : 0 < R) (hw : 0 < w) (hK : K < R ^ 2)
    (hb : (x - R * w) ^ 2 ≤ K * w ^ 2) : 0 < x := by
  rcases le_or_gt x 0 with h | h
  · exfalso
    have hRw : 0 < R * w := mul_pos hR hw
    have h1 : 0 ≤ x * (x - 2 * (R * w)) := mul_nonneg_of_nonpos_of_nonpos h (by linarith)
    have h2 : K * w ^ 2 < R ^ 2 * w ^ 2 := mul_lt_mul_of_pos_right hK (by positivity)
    have e : (x - R * w) ^ 2 = x * (x - 2 * (R * w)) + R ^ 2 * w ^ 2 := by ring
    linarith
  · exact h

/-- NOTE.md §2 in real form for N = k + 1 and any root d ∈ (1/2, 2/3) of d³ + 3Nd = 2N:
0 < R₀, K < R₀², and (c_n − R₀ z₀^{−n})² ≤ K z₀^{−2n} with z₀ = N + d. -/
theorem main_real (k : ℕ) (d : ℝ) (hd1 : 1 / 2 < d) (hd2 : d < 2 / 3)
    (hr : d ^ 3 + 3 * ((k : ℝ) + 1) * d = 2 * ((k : ℝ) + 1)) :
    0 < R0 ((k : ℝ) + 1) d ∧ KR ((k : ℝ) + 1) d < R0 ((k : ℝ) + 1) d ^ 2 ∧
      ∀ n, ((c k n : ℝ) - R0 ((k : ℝ) + 1) d * (1 / ((k : ℝ) + 1 + d)) ^ n) ^ 2
        ≤ KR ((k : ℝ) + 1) d * ((1 / ((k : ℝ) + 1 + d)) ^ n) ^ 2 := by
  have hN : (1 : ℝ) ≤ (k : ℝ) + 1 := by
    have := Nat.cast_nonneg (α := ℝ) k
    linarith
  exact core ((k : ℝ) + 1) d ((k : ℝ) + 1 + d) (uR _ d) (vR _ d) (R0 _ d) (fun m => (cs k m : ℝ))
    hN hd1 hd2 hr rfl rfl rfl rfl (by norm_num [cs]) (by norm_num [cs]) (by norm_num [cs])
    (cs_rec_real k)

/-- MAIN THEOREM (Fasi–Higham, m = 3, every k): every Taylor coefficient c_n(k) of 1/q_{k,3} is positive. -/
theorem c_pos (k n : ℕ) : 0 < c k n := by
  have hN : (1 : ℝ) ≤ (k : ℝ) + 1 := by
    have := Nat.cast_nonneg (α := ℝ) k
    linarith
  obtain ⟨d, hd1, hd2, hr⟩ := exists_root ((k : ℝ) + 1) hN
  obtain ⟨hR, hK, hb⟩ := main_real k d hd1 hd2 hr
  have hw : 0 < (1 / ((k : ℝ) + 1 + d)) ^ n := pow_pos (by positivity) n
  have := pos_of_bound (c k n : ℝ) _ _ _ hR hw hK (hb n)
  exact_mod_cast this

end FasiLean
