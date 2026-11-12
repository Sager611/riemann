import Mathlib.Algebra.Order.Archimedean.Real.Basic
import Mathlib.Tactic.Linarith
import Mathlib.Tactic.Positivity
import Mathlib.Tactic.LinearCombination
import Mathlib.Tactic.Ring
import Mathlib.Tactic.FieldSimp

/-!
# Fasi–Higham, m = 3: the real-variable core (NOTE.md §2 without complex numbers)

With N = k + 1 and d the real root of d³ + 3Nd = 2N, put a = N + d (the positive zero z₀ of q),
u = N − d/2, v = N² − dN + d² + 3N (z² − 2uz + v vanishes at the complex pair z_±, |z_±|² = v),
R = v/(3(d² + N)) (the residue R₀) and ρ = 1/a. Then N(N+1)(N+2) q(z) = (a − z)(z² − 2uz + v) and
c_n = R ρ^n + f_{n+2}, where f satisfies v f_{m+2} = 2u f_{m+1} − f_m with f_0 = −R a², f_1 = −R a
(the extension of c_n − Rρ^n to n = −2, −1). G_m = v f_{m+1}² − 2u f_{m+1} f_m + f_m² obeys
v G_{m+1} = G_m, G_0 = R² a² · 3(d² + N), and (v − u²) f_m² ≤ v G_m. Hence
(c_n − Rρ^n)² ≤ K ρ^{2n} with K = R² · 3a²(d² + N)/(v(v − u²)) < R² (key inequality), so c_n > 0.
-/

namespace FasiLean

/-- f_0 = −R a², f_1 = −R a, f_{m+2} = (2u f_{m+1} − f_m)/v. -/
noncomputable def fs (a u v R : ℝ) : ℕ → ℝ
  | 0 => -(R * a ^ 2)
  | 1 => -(R * a)
  | m + 2 => (2 * u * fs a u v R (m + 1) - fs a u v R m) / v

/-- The conserved form G_m = v f_{m+1}² − 2u f_{m+1} f_m + f_m². -/
def Gf (u v : ℝ) (f : ℕ → ℝ) (m : ℕ) : ℝ := v * f (m + 1) ^ 2 - 2 * u * f (m + 1) * f m + f m ^ 2

theorem fs_rec (a u v R : ℝ) (hv : v ≠ 0) (m : ℕ) :
    v * fs a u v R (m + 2) = 2 * u * fs a u v R (m + 1) - fs a u v R m := by
  simp only [fs]
  field_simp

theorem Gf_succ (u v : ℝ) (f : ℕ → ℝ) (hf : ∀ m, v * f (m + 2) = 2 * u * f (m + 1) - f m)
    (m : ℕ) : v * Gf u v f (m + 1) = Gf u v f m := by
  simp only [Gf]
  linear_combination (v * f (m + 2) - f m) * hf m

theorem Gf_pow (u v : ℝ) (f : ℕ → ℝ) (hf : ∀ m, v * f (m + 2) = 2 * u * f (m + 1) - f m)
    (m : ℕ) : Gf u v f m * v ^ m = Gf u v f 0 := by
  induction m with
  | zero => simp
  | succ m ih =>
    rw [← ih, ← Gf_succ u v f hf m, pow_succ]
    ring

theorem sq_le_Gf (u v : ℝ) (f : ℕ → ℝ) (m : ℕ) : (v - u ^ 2) * f m ^ 2 ≤ v * Gf u v f m := by
  have h : v * Gf u v f m - (v - u ^ 2) * f m ^ 2 = (v * f (m + 1) - u * f m) ^ 2 := by
    simp only [Gf]
    ring
  nlinarith [sq_nonneg (v * f (m + 1) - u * f m)]

/-- c_m (shifted: C 0 = c_{−2}, C 1 = c_{−1}, C 2 = c_0) = R a² ρ^m + f_m. -/
theorem decomp (N a u v R ρ : ℝ) (C : ℕ → ℝ)
    (hP : N * (N + 1) * (N + 2) = a * v) (h2 : 3 * N * (N + 1) = 2 * u * a + v)
    (h3 : 3 * N = a + 2 * u) (hS : R * (a ^ 2 - 2 * u * a + v) = v) (hρ : ρ * a = 1)
    (hv : v ≠ 0) (hP0 : N * (N + 1) * (N + 2) ≠ 0)
    (c0 : C 0 = 0) (c1 : C 1 = 0) (c2 : C 2 = 1)
    (hC : ∀ m, N * (N + 1) * (N + 2) * C (m + 3)
      = 3 * N * (N + 1) * C (m + 2) - 3 * N * C (m + 1) + C m) :
    ∀ m, C m = R * a ^ 2 * ρ ^ m + fs a u v R m := by
  have hf := fs_rec a u v R hv
  have e0 : fs a u v R 0 = -(R * a ^ 2) := rfl
  have e1 : fs a u v R 1 = -(R * a) := rfl
  have key : ∀ m, C m = R * a ^ 2 * ρ ^ m + fs a u v R m
      ∧ C (m + 1) = R * a ^ 2 * ρ ^ (m + 1) + fs a u v R (m + 1)
      ∧ C (m + 2) = R * a ^ 2 * ρ ^ (m + 2) + fs a u v R (m + 2) := by
    intro m
    induction m with
    | zero =>
      refine ⟨?_, ?_, ?_⟩
      · rw [c0, e0]
        ring
      · rw [show (0 : ℕ) + 1 = 1 from rfl, c1, e1]
        linear_combination (-(R * a)) * hρ
      · rw [show (0 : ℕ) + 2 = 2 from rfl, c2]
        apply mul_left_cancel₀ hv
        have h20 := hf 0
        rw [e0, e1] at h20
        linear_combination (-(R * v * (ρ * a + 1))) * hρ - hS - h20
    | succ m ih =>
      obtain ⟨i0, i1, i2⟩ := ih
      refine ⟨i1, i2, ?_⟩
      show C (m + 3) = R * a ^ 2 * ρ ^ (m + 3) + fs a u v R (m + 3)
      apply mul_left_cancel₀ hP0
      have e2 : v * fs a u v R (m + 2) = 2 * u * fs a u v R (m + 1) - fs a u v R m := hf m
      have e3 : v * fs a u v R (m + 3) = 2 * u * fs a u v R (m + 2) - fs a u v R (m + 1) :=
        hf (m + 1)
      rw [hC m, i0, i1, i2]
      linear_combination (-(R * a ^ 2 * ρ ^ m) * (v * ρ ^ 2 - 2 * u * ρ + 1)) * hρ
        + (-(R * a ^ 2 * ρ ^ m) * ρ ^ 3 - fs a u v R (m + 3)) * hP
        + (R * a ^ 2 * ρ ^ m * ρ ^ 2 + fs a u v R (m + 2)) * h2
        + (-(R * a ^ 2 * ρ ^ m) * ρ - fs a u v R (m + 1)) * h3 - a * e3 + e2
  exact fun m => (key m).1

end FasiLean
