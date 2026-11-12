import Riemann.Deep.ExploreDPP.Lean.WeakDuality

/-!
# Theorem A: the discrete-sine dual (explore34)

h = a = 1/(1+θ), k0 = 1, η_k = -sinc²(ka): P(α) = 1 - Σ_m Δ(α - m/a) ≥ |α|·1[|α| ≤ θ],
so no admissible f has J_θ(f) > 1 when λ ≤ 1/(1+θ), 0 < θ ≤ 1.
-/

open Real MeasureTheory FourierTransform Filter Asymptotics
open scoped Real

noncomputable section

namespace ExploreDPP

/-! ## Theorem A: the discrete-sine dual -/

/-- Triangle Δ(x) = max(0, 1 - |x|). -/
def tri (x : ℝ) : ℝ := max 0 (1 - |x|)

/-- η_k = -sinc²(k a) with the normalized sinc(x) = sin(πx)/(πx) = Real.sinc (π x). -/
def etaA (a : ℝ) (k : ℕ) : ℝ := -(Real.sinc (π * (k * a))) ^ 2

/-- Named hypothesis: the Fourier series of the periodized triangle (Poisson summation with
𝓕Δ = sinc²): a + 2a Σ_{k≥1} sinc²(ka) cos(2π α k a) = Σ_m Δ(α - m/a). -/
def TriangleSeries (a : ℝ) : Prop :=
  ∀ α : ℝ, a + 2 * a * ∑' n : ℕ, Real.sinc (π * ((n + 1) * a)) ^ 2 *
      Real.cos (2 * π * α * ((n + 1) * a)) = ∑' m : ℤ, tri (α - m / a)

lemma sinc_sq_le {y : ℝ} (hy : y ≠ 0) : Real.sinc y ^ 2 ≤ 1 / y ^ 2 := by
  rw [Real.sinc_of_ne_zero hy, div_pow]
  exact div_le_div_of_nonneg_right (Real.sin_sq_le_one y) (by positivity)

lemma etaA_ge (a : ℝ) (k : ℕ) : -1 ≤ etaA a k := by
  have := Real.abs_sinc_le_one (π * (k * a))
  have h2 : Real.sinc (π * (k * a)) ^ 2 ≤ 1 := by
    rw [← sq_abs]; nlinarith [abs_nonneg (Real.sinc (π * (k * a)))]
  unfold etaA; linarith

lemma summable_etaA {a : ℝ} (ha : 0 < a) : Summable (etaA a) := by
  rw [← summable_nat_add_iff 1]
  have hs : Summable (fun n : ℕ => 1 / (π * a) ^ 2 * (1 / ((n : ℝ) + 1) ^ 2)) := by
    have := (summable_nat_add_iff 1).mpr (Real.summable_one_div_nat_pow.mpr one_lt_two)
    exact (this.mul_left (1 / (π * a) ^ 2)).congr (fun n => by push_cast; ring)
  refine Summable.of_norm_bounded hs (fun n => ?_)
  have hy : π * (((n + 1 : ℕ) : ℝ) * a) ≠ 0 := by positivity
  rw [Real.norm_eq_abs, etaA, abs_neg, abs_of_nonneg (sq_nonneg _)]
  calc Real.sinc (π * (((n + 1 : ℕ) : ℝ) * a)) ^ 2
      ≤ 1 / (π * (((n + 1 : ℕ) : ℝ) * a)) ^ 2 := sinc_sq_le hy
    _ = 1 / (π * a) ^ 2 * (1 / ((n : ℝ) + 1) ^ 2) := by
      have := Real.pi_pos.ne'
      have : (n : ℝ) + 1 ≠ 0 := by positivity
      push_cast; field_simp

lemma P_etaA {a : ℝ} (hT : TriangleSeries a) (α : ℝ) :
    P a (etaA a) α = 1 - ∑' m : ℤ, tri (α - m / a) := by
  rw [← hT α]
  unfold P etaA
  have e : ∀ n : ℕ, -(Real.sinc (π * (((n + 1 : ℕ) : ℝ) * a))) ^ 2 *
      Real.cos (2 * π * α * ((n + 1) * a)) =
      -(Real.sinc (π * ((n + 1) * a)) ^ 2 * Real.cos (2 * π * α * ((n + 1) * a))) := by
    intro n; push_cast; ring
  simp only [e, tsum_neg]
  ring

lemma tri_sum_base {θ : ℝ} (hθ0 : 0 < θ) {α : ℝ} (h0 : 0 ≤ α) (h1 : α ≤ (1 + θ) / 2) :
    ∑' m : ℤ, tri (α - m / (1 / (1 + θ))) = tri α + tri (α - (1 + θ)) := by
  have hma : ∀ m : ℤ, (m : ℝ) / (1 / (1 + θ)) = m * (1 + θ) := by
    intro m; field_simp
  simp only [hma]
  rw [tsum_eq_sum (s := {0, 1})]
  · rw [Finset.sum_pair (by norm_num)]
    simp
  · intro m hm
    simp only [Finset.mem_insert, Finset.mem_singleton, not_or] at hm
    unfold tri
    apply max_eq_left
    rcases le_or_gt 2 m with h2 | h2
    · have h2' : (2 : ℝ) ≤ m := by exact_mod_cast h2
      have : α - m * (1 + θ) ≤ -1 := by nlinarith
      rw [abs_of_neg (by linarith)]; linarith
    · have h3 : m ≤ -1 := by omega
      have h3' : (m : ℝ) ≤ -1 := by exact_mod_cast h3
      have : 1 ≤ α - m * (1 + θ) := by nlinarith
      rw [abs_of_pos (by linarith)]; linarith

lemma SA_ge {θ : ℝ} (hθ0 : 0 < θ) (hθ1 : θ ≤ 1) {α : ℝ} (h0 : 0 ≤ α) (h1 : α ≤ (1 + θ) / 2) :
    (if α ≤ θ then α else 0) ≤ 1 - (tri α + tri (α - (1 + θ))) := by
  unfold tri
  rw [abs_of_nonneg h0, abs_of_neg (by linarith : α - (1 + θ) < 0)]
  simp only [max_def]
  split_ifs <;> linarith

/-- Theorem A (discrete-sine dual). Given the Fourier series of the periodized triangle at
a = 1/(1+θ): for 0 < θ ≤ 1 and λ ≤ 1/(1+θ) every admissible f has J_θ(f) ≤ 1 = f^(0). -/
theorem theoremA {θ lam : ℝ} (hθ0 : 0 < θ) (hθ1 : θ ≤ 1) (hT : TriangleSeries (1 / (1 + θ)))
    (hlam : lam ≤ 1 / (1 + θ)) {f : ℝ → ℝ} (hf : Admissible lam f) : J θ f ≤ 1 := by
  have ha0 : 0 < 1 / (1 + θ) := by positivity
  have hf' : Admissible (1 / (1 + θ)) f := hf.mono hlam
  have hhalf : 1 / (2 * (1 / (1 + θ))) = (1 + θ) / 2 := by field_simp
  refine weak_duality (h := 1 / (1 + θ)) (k0 := 1) (η := etaA (1 / (1 + θ))) ha0 (by simp)
    (summable_etaA ha0) (fun k hk1 hk => by omega) (fun k _ => etaA_ge _ k) hθ0.le
    (by rw [hhalf]; linarith) (fun α h0 h1 => ?_) hf'
  rw [hhalf] at h1
  rw [P_etaA hT, tri_sum_base hθ0 h0 h1]
  exact SA_ge hθ0 hθ1 h0 h1

end ExploreDPP
