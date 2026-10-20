import Mathlib.Analysis.SpecialFunctions.Trigonometric.Bounds
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Analysis.Complex.Exponential

/-!
# Complex ball arithmetic with explicit rational centres

A ball is ‖z - (a + b i)‖ ≤ r with real numerals a b r. Each lemma takes a numeric side condition
that norm_num discharges; centres are rounded freely, rounding is charged to the radius.
-/

open Complex

namespace Fayed.FirstZero

lemma norm_le_abs_add_abs (z : ℂ) : ‖z‖ ≤ |z.re| + |z.im| := norm_le_abs_re_add_abs_im z

lemma ball_add {z w : ℂ} {a b r c d q e f R : ℝ} (hz : ‖z - (a + b * I)‖ ≤ r)
    (hw : ‖w - (c + d * I)‖ ≤ q) (h : |a + c - e| + |b + d - f| + r + q ≤ R) :
    ‖z + w - (e + f * I)‖ ≤ R := by
  have e1 : z + w - (e + f * I) = (z - (a + b * I)) + (w - (c + d * I)) +
      (((a + c - e : ℝ) : ℂ) + ((b + d - f : ℝ) : ℂ) * I) := by push_cast; ring
  rw [e1]
  have h3 : ‖(((a + c - e : ℝ) : ℂ) + ((b + d - f : ℝ) : ℂ) * I)‖ ≤ |a + c - e| + |b + d - f| := by
    refine (norm_le_abs_re_add_abs_im _).trans (le_of_eq ?_); simp
  calc _ ≤ ‖z - (a + b * I)‖ + ‖w - (c + d * I)‖ + ‖(((a + c - e : ℝ) : ℂ) + ((b + d - f : ℝ) : ℂ) * I)‖ :=
        (norm_add_le _ _).trans (by gcongr; exact norm_add_le _ _)
    _ ≤ r + q + (|a + c - e| + |b + d - f|) := by gcongr
    _ ≤ R := by linarith

lemma ball_mul {z w : ℂ} {a b r c d q e f R : ℝ} (hz : ‖z - (a + b * I)‖ ≤ r)
    (hw : ‖w - (c + d * I)‖ ≤ q)
    (h : (|a| + |b|) * q + (|c| + |d|) * r + r * q + |a * c - b * d - e| + |a * d + b * c - f| ≤ R) :
    ‖z * w - (e + f * I)‖ ≤ R := by
  set A : ℂ := a + b * I
  set B : ℂ := c + d * I
  have e1 : z * w - (e + f * I) = (z - A) * (w - B) + A * (w - B) + B * (z - A) +
      (((a * c - b * d - e : ℝ) : ℂ) + ((a * d + b * c - f : ℝ) : ℂ) * I) := by
    simp only [A, B]; push_cast; ring_nf; rw [I_sq]; ring
  rw [e1]
  have hA : ‖A‖ ≤ |a| + |b| := (norm_le_abs_re_add_abs_im _).trans (le_of_eq (by simp [A]))
  have hB : ‖B‖ ≤ |c| + |d| := (norm_le_abs_re_add_abs_im _).trans (le_of_eq (by simp [B]))
  have h3 : ‖(((a * c - b * d - e : ℝ) : ℂ) + ((a * d + b * c - f : ℝ) : ℂ) * I)‖ ≤
      |a * c - b * d - e| + |a * d + b * c - f| := by
    refine (norm_le_abs_re_add_abs_im _).trans (le_of_eq ?_); simp
  have hr : 0 ≤ r := (norm_nonneg _).trans hz
  have hq : 0 ≤ q := (norm_nonneg _).trans hw
  calc _ ≤ ‖(z - A) * (w - B)‖ + ‖A * (w - B)‖ + ‖B * (z - A)‖ +
        ‖(((a * c - b * d - e : ℝ) : ℂ) + ((a * d + b * c - f : ℝ) : ℂ) * I)‖ := by
        refine (norm_add_le _ _).trans ?_
        gcongr
        refine (norm_add_le _ _).trans ?_
        gcongr
        exact norm_add_le _ _
    _ ≤ r * q + (|a| + |b|) * q + (|c| + |d|) * r + (|a * c - b * d - e| + |a * d + b * c - f|) := by
        rw [norm_mul, norm_mul, norm_mul]
        gcongr
    _ ≤ R := by linarith

/-- Real scalar (known to ±ε) times a ball. -/
lemma ball_real_mul {x : ℝ} {z : ℂ} {c ε a b r e f R : ℝ} (hx : |x - c| ≤ ε)
    (hz : ‖z - (a + b * I)‖ ≤ r)
    (h : (|a| + |b|) * ε + |c| * r + ε * r + |c * a - e| + |c * b - f| ≤ R) :
    ‖(x : ℂ) * z - (e + f * I)‖ ≤ R := by
  have hx' : ‖(x : ℂ) - ((c : ℝ) + (0 : ℝ) * I)‖ ≤ ε := by
    rw [show (x : ℂ) - ((c : ℝ) + (0 : ℝ) * I) = ((x - c : ℝ) : ℂ) by push_cast; ring,
      Complex.norm_real, Real.norm_eq_abs]
    exact hx
  refine ball_mul hx' hz ?_
  have : |(0 : ℝ)| = 0 := abs_zero
  simp only [this, add_zero, zero_mul, sub_zero]
  linarith

/-- e^{iθ} from a reduced rational argument q ≈ θ - 2πm and its degree-19 Taylor polynomial. -/
lemma sum_exp_I (y : ℝ) (K : ℕ) :
    ∑ m ∈ Finset.range (2 * K), ((y : ℂ) * I) ^ m / (m.factorial : ℂ) =
      ((∑ k ∈ Finset.range K, (-1) ^ k * y ^ (2 * k) / ((2 * k).factorial : ℝ) : ℝ) : ℂ) +
      ((∑ k ∈ Finset.range K, (-1) ^ k * y ^ (2 * k + 1) / ((2 * k + 1).factorial : ℝ) : ℝ) : ℂ)
        * I := by
  induction K with
  | zero => simp
  | succ K ih =>
    rw [show 2 * (K + 1) = 2 * K + 1 + 1 by ring, Finset.sum_range_succ, Finset.sum_range_succ, ih,
      Finset.sum_range_succ, Finset.sum_range_succ]
    have h1 : ((y : ℂ) * I) ^ (2 * K) = ((-1) ^ K * y ^ (2 * K) : ℝ) := by
      rw [mul_pow, pow_mul I, I_sq]; push_cast; ring
    have h2 : ((y : ℂ) * I) ^ (2 * K + 1) = ((-1) ^ K * y ^ (2 * K + 1) : ℝ) * I := by
      rw [pow_succ, h1]; push_cast; ring
    rw [h1, h2]
    push_cast
    ring

lemma norm_cexp_I_sub_cexp_I (x y : ℝ) : ‖cexp (x * I) - cexp (y * I)‖ ≤ |x - y| := by
  have e : cexp (x * I) - cexp (y * I) = cexp (y * I) * (cexp (I * ((x - y : ℝ) : ℂ)) - 1) := by
    rw [mul_sub, ← Complex.exp_add, mul_one]; congr 2; push_cast; ring
  rw [e, norm_mul, Complex.norm_exp_ofReal_mul_I, one_mul]
  exact (Real.norm_exp_I_mul_ofReal_sub_one_le).trans (le_of_eq (Real.norm_eq_abs _))

lemma cis_ball {θ q η : ℝ} (m : ℤ) (hθ : |θ - (q + 2 * Real.pi * m)| ≤ η) (hq : |q| ≤ 16 / 5)
    {c s ε : ℝ}
    (h : |(∑ k ∈ Finset.range 10, (-1) ^ k * q ^ (2 * k) / ((2 * k).factorial : ℝ)) - c| +
      |(∑ k ∈ Finset.range 10, (-1) ^ k * q ^ (2 * k + 1) / ((2 * k + 1).factorial : ℝ)) - s| +
      (16 / 5) ^ 20 / (Nat.factorial 20) * 2 + η ≤ ε) :
    ‖cexp (θ * I) - (c + s * I)‖ ≤ ε := by
  have hper : cexp (θ * I) = cexp (((θ - 2 * Real.pi * m : ℝ) : ℂ) * I) := by
    rw [show (((θ - 2 * Real.pi * m : ℝ) : ℂ) * I) = (θ : ℂ) * I - (m : ℂ) * (2 * Real.pi * I) by
      push_cast; ring, Complex.exp_sub, Complex.exp_int_mul_two_pi_mul_I, div_one]
  have h1 : ‖cexp (θ * I) - cexp (q * I)‖ ≤ η := by
    rw [hper]
    refine (norm_cexp_I_sub_cexp_I _ _).trans ?_
    rw [show θ - 2 * Real.pi * m - q = θ - (q + 2 * Real.pi * m) by ring]
    exact hθ
  have hqn : ‖(q : ℂ) * I‖ = |q| := by simp
  have h2 : ‖cexp ((q : ℂ) * I) - ∑ m ∈ Finset.range 20, ((q : ℂ) * I) ^ m / (m.factorial : ℂ)‖ ≤
      (16 / 5) ^ 20 / (Nat.factorial 20) * 2 := by
    refine (Complex.exp_bound' (n := 20) (by rw [hqn]; push_cast; linarith)).trans ?_
    rw [hqn]
    gcongr
  rw [show (20 : ℕ) = 2 * 10 by norm_num, sum_exp_I] at h2
  set C := ∑ k ∈ Finset.range 10, (-1) ^ k * q ^ (2 * k) / ((2 * k).factorial : ℝ)
  set S := ∑ k ∈ Finset.range 10, (-1) ^ k * q ^ (2 * k + 1) / ((2 * k + 1).factorial : ℝ)
  have h3 : ‖((C : ℂ) + (S : ℂ) * I) - (c + s * I)‖ ≤ |C - c| + |S - s| := by
    refine (norm_le_abs_re_add_abs_im _).trans (le_of_eq ?_); simp
  calc ‖cexp (θ * I) - (c + s * I)‖ ≤ ‖cexp (θ * I) - cexp (q * I)‖ +
        ‖cexp ((q : ℂ) * I) - ((C : ℂ) + (S : ℂ) * I)‖ + ‖((C : ℂ) + (S : ℂ) * I) - (c + s * I)‖ := by
        refine (norm_sub_le_norm_sub_add_norm_sub (cexp (θ * I)) (cexp ((q : ℂ) * I))
          (c + s * I)).trans ?_
        have := norm_sub_le_norm_sub_add_norm_sub (cexp ((q : ℂ) * I)) ((C : ℂ) + (S : ℂ) * I)
          (c + s * I)
        linarith
    _ ≤ η + (16 / 5) ^ 20 / (Nat.factorial 20) * 2 + (|C - c| + |S - s|) := by gcongr
    _ ≤ ε := by linarith

lemma inv_sqrt_ball {p c ε : ℝ} (hp : 0 < p) (hε : 0 ≤ ε) (hc : 0 < c - ε) (h1 : (c - ε) ^ 2 * p ≤ 1)
    (h2 : 1 ≤ (c + ε) ^ 2 * p) : |(Real.sqrt p)⁻¹ - c| ≤ ε := by
  have hs : 0 < Real.sqrt p := Real.sqrt_pos.mpr hp
  set x := (Real.sqrt p)⁻¹ with hx
  have hx0 : 0 < x := inv_pos.mpr hs
  have hxx : x ^ 2 * p = 1 := by
    rw [hx, inv_pow, Real.sq_sqrt hp.le, inv_mul_cancel₀ hp.ne']
  have l1 : c - ε ≤ x := by
    by_contra hcon
    push Not at hcon
    have : x ^ 2 < (c - ε) ^ 2 := pow_lt_pow_left₀ hcon hx0.le (by norm_num)
    nlinarith
  have l2 : x ≤ c + ε := by
    by_contra hcon
    push Not at hcon
    have : (c + ε) ^ 2 < x ^ 2 := pow_lt_pow_left₀ hcon (by linarith) (by norm_num)
    nlinarith
  rw [abs_le]; constructor <;> linarith

lemma norm_le_of_sq_le {w : ℂ} {c : ℝ} (hc : 0 ≤ c) (h : w.re ^ 2 + w.im ^ 2 ≤ c ^ 2) :
    ‖w‖ ≤ c := by
  have : ‖w‖ ^ 2 ≤ c ^ 2 := by
    rw [Complex.sq_norm, Complex.normSq_apply]; nlinarith
  exact (pow_le_pow_iff_left₀ (norm_nonneg _) hc (by norm_num)).mp this

end Fayed.FirstZero
