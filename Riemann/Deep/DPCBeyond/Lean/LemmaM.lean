import Riemann.Deep.DPCBeyond.Lean.Rpow

/-!
# dpcbeyond Lemma M (new step (c)): the exponent bookkeeping of the moments

PROOF.md Lemma M. Proved here (abstract positive reals, rp = x^e from ExpLog):
* M2_bound — M2's Hölder step: ∫(n^{(a)}_m)^p ≤ meas^{1−1/w}(∫n_m^{pw})^{1/w} (Hölder, named as hypothesis hI),
  meas ≤ |B_a|/L, |B_a| ≤ C T L e^{−κa} (density bins), ∫n_m^{pw} ≤ C T (Fujii) ⇒ ≤ C T e^{−κa(1−1/w)};
* geom_mul, geom_le, M3_series — M3: Σ_a e^{θ(a+1) − κa(1−1/w)/p} has all partial sums ≤ e^θ/(1 − e^{θ−κ(1−1/w)/p})
  when pθ < κ(1 − 1/w) (convergence of a series of positive terms);
* M3_exists_w — for 0 ≤ pθ < κ some w > 1 has pθ < κ(1 − 1/w) and pθ < κ(1 − 1/w)²;
* M4_count, M4_moment — M4: #{ρ : n^{(a)}_m(x_ρ) ≥ 1} ≤ |B_a|^{1−1/w}(Σ Ñ^w)^{1/w} ≪ (|m|+2)TL e^{−κa(1−1/w)},
  then Σ_ρ n^{(a)}_m(x_ρ)^p ≪ (|m|+2) TL e^{−κa(1−1/w)²}: the exponent κ(1 − 1/w)²;
* key_ineq, inv_sq_sum, kernel_sum — M4's kernel sum: Σ_{n<M}(n+2)/(n+1)^A ≤ 4 for natural A ≥ 3, so
  Σ_{m∈ℤ}(1 + |m|)^{−A}(|m| + 2) ≤ 8 (each n = |m| ≥ 1 twice).
Named (hypotheses of these lemmas): Hölder's inequality, Fujii's moments, Lz DoubleSumZeros, the density bins.
-/

namespace DPCBLean

open Finset DPCLean

variable {ex lg : ℝ → ℝ}

lemma rp_ex (h : ExpLog ex lg) (x e : ℝ) : rp ex lg (ex x) e = ex (e * x) := by
  rw [rp_of_pos (h.pos x), lgex h]

/-- **M2**: Hölder (hI), meas ≤ |B_a|/L, |B_a| ≤ C T L e^{−κa}, ∫ n_m^{pw} ≤ C T give ∫(n^{(a)}_m)^p ≤ C T e^{−κa(1−1/w)}. -/
theorem M2_bound (h : ExpLog ex lg) {I μ B N L C T κ a w : ℝ} (hL : 0 < L) (hC : 0 < C)
    (hT : 0 < T) (hw : 1 < w) (hI : I ≤ rp ex lg μ (1 - 1 / w) * rp ex lg N (1 / w))
    (hμ : μ ≤ B / L) (hB : B ≤ C * T * L * ex (-(κ * a))) (hN : N ≤ C * T) :
    I ≤ C * T * ex (-(κ * a * (1 - 1 / w))) := by
  have hw0 : 0 < w := by linarith only [hw]
  have he1 : 0 ≤ 1 - 1 / w := by
    have h1 : 1 / w < 1 := (div_lt_one hw0).2 hw
    linarith only [h1]
  have hiw : 0 ≤ 1 / w := by positivity
  have hCT : 0 < C * T := mul_pos hC hT
  have hμ' : μ ≤ C * T * ex (-(κ * a)) := by
    have h1 : μ * L ≤ C * T * ex (-(κ * a)) * L := by
      rw [le_div_iff₀ hL] at hμ
      have e : C * T * L * ex (-(κ * a)) = C * T * ex (-(κ * a)) * L := by ring
      linarith only [hμ, hB, e]
    exact le_of_mul_le_mul_right h1 hL
  have r1 : rp ex lg μ (1 - 1 / w) ≤ rp ex lg (C * T) (1 - 1 / w) * ex (-(κ * a * (1 - 1 / w))) := by
    calc rp ex lg μ (1 - 1 / w) ≤ rp ex lg (C * T * ex (-(κ * a))) (1 - 1 / w) := rp_le_rp h hμ' he1
      _ = rp ex lg (C * T) (1 - 1 / w) * rp ex lg (ex (-(κ * a))) (1 - 1 / w) :=
          rp_mul h hCT (h.pos _) _
      _ = rp ex lg (C * T) (1 - 1 / w) * ex (-(κ * a * (1 - 1 / w))) := by
          rw [rp_ex h, show (1 - 1 / w) * -(κ * a) = -(κ * a * (1 - 1 / w)) by ring]
  have r2 : rp ex lg N (1 / w) ≤ rp ex lg (C * T) (1 / w) := rp_le_rp h hN hiw
  have e : rp ex lg (C * T) (1 - 1 / w) * rp ex lg (C * T) (1 / w) = C * T := by
    rw [← rp_add h hCT, show 1 - 1 / w + 1 / w = 1 by ring, rp_one h hCT]
  calc I ≤ rp ex lg μ (1 - 1 / w) * rp ex lg N (1 / w) := hI
    _ ≤ (rp ex lg (C * T) (1 - 1 / w) * ex (-(κ * a * (1 - 1 / w)))) * rp ex lg (C * T) (1 / w) :=
        mul_le_mul r1 r2 (rp_nonneg h _ _) (mul_nonneg (rp_nonneg h _ _) (h.pos _).le)
    _ = (rp ex lg (C * T) (1 - 1 / w) * rp ex lg (C * T) (1 / w)) * ex (-(κ * a * (1 - 1 / w))) := by
        ring
    _ = C * T * ex (-(κ * a * (1 - 1 / w))) := by rw [e]

/-- Geometric sums: (Σ_{a<M} ρ^a)(1 − ρ) = 1 − ρ^M. -/
lemma geom_mul (ρ : ℝ) (M : ℕ) : (∑ a ∈ range M, ρ ^ a) * (1 - ρ) = 1 - ρ ^ M := by
  induction M with
  | zero => simp
  | succ M ih =>
    rw [sum_range_succ, add_mul, ih, pow_succ]
    ring

lemma geom_le {ρ : ℝ} (h0 : 0 ≤ ρ) (h1 : ρ < 1) (M : ℕ) : ∑ a ∈ range M, ρ ^ a ≤ 1 / (1 - ρ) := by
  rw [le_div_iff₀ (by linarith only [h1]), geom_mul]
  have h2 := pow_nonneg h0 M
  linarith only [h2]

/-- **M3**: when pθ < κ(1 − 1/w), every partial sum of Σ_a e^{θ(a+1) − κa(1−1/w)/p} is ≤ e^θ/(1 − e^{θ−κ(1−1/w)/p})
(PROOF: ‖Off‖_p ≪ T^{1/p} Σ_a e^{θ(a+1) − κa(1−1/w)/p} ≪ T^{1/p}). -/
theorem M3_series (h : ExpLog ex lg) {θ κ w p : ℝ} (hp : 0 < p) (hpθ : p * θ < κ * (1 - 1 / w))
    (M : ℕ) : ∑ a ∈ range M, ex (θ * (a + 1) - κ * a * (1 - 1 / w) / p) ≤
      ex θ / (1 - ex (θ - κ * (1 - 1 / w) / p)) := by
  have hlt : θ < κ * (1 - 1 / w) / p := by
    rw [lt_div_iff₀ hp, mul_comm]
    exact hpθ
  have hρ1 : ex (θ - κ * (1 - 1 / w) / p) < 1 := h.lt_one (by linarith only [hlt])
  have hρ0 : 0 ≤ ex (θ - κ * (1 - 1 / w) / p) := (h.pos _).le
  have e : ∀ a ∈ range M, ex (θ * (a + 1) - κ * a * (1 - 1 / w) / p) =
      ex θ * ex (θ - κ * (1 - 1 / w) / p) ^ a := by
    intro a _
    rw [← h.nat_mul, ← h.add]
    congr 1
    ring
  rw [sum_congr rfl e, ← mul_sum]
  have h2 := mul_le_mul_of_nonneg_left (geom_le hρ0 hρ1 M) (h.pos θ).le
  rwa [mul_one_div] at h2

/-- M3/M4: for 0 ≤ pθ < κ some w > 1 has pθ < κ(1 − 1/w) and pθ < κ(1 − 1/w)² ("choose w large"). -/
theorem M3_exists_w {θ κ p : ℝ} (hκ : 0 < κ) (hpθ0 : 0 ≤ p * θ) (hpθ : p * θ < κ) :
    ∃ w : ℝ, 1 < w ∧ p * θ < κ * (1 - 1 / w) ∧ p * θ < κ * (1 - 1 / w) ^ 2 := by
  have hd : 0 < κ - p * θ := by linarith only [hpθ]
  obtain ⟨x, hx⟩ : ∃ x : ℝ, x = (κ - p * θ) / (4 * κ) := ⟨_, rfl⟩
  have hx0 : 0 < x := by rw [hx]; positivity
  have hκx : κ * x = (κ - p * θ) / 4 := by
    rw [hx, mul_div_assoc', div_eq_div_iff (by positivity) (by norm_num)]
    ring
  have hx1 : x < 1 := by
    have h1 : κ * x < κ * 1 := by rw [hκx]; linarith only [hpθ0, hκ]
    exact lt_of_mul_lt_mul_left h1 hκ.le
  refine ⟨1 / x, one_lt_one_div hx0 hx1, ?_, ?_⟩
  · rw [one_div_one_div]
    have e : κ * (1 - x) = κ - κ * x := by ring
    linarith only [e, hκx, hd]
  · rw [one_div_one_div]
    have e : κ * (1 - x) ^ 2 = κ - 2 * (κ * x) + κ * x ^ 2 := by ring
    have h2 : 0 ≤ κ * x ^ 2 := mul_nonneg hκ.le (sq_nonneg x)
    linarith only [e, hκx, hd, h2]

end DPCBLean
