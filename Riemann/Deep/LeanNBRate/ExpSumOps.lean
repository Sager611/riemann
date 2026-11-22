import Riemann.Deep.LeanNBRate.ExpSum

/-!
# Operations on exponential sums: modulation, derivative
-/

open MeasureTheory Filter Topology Set Complex

noncomputable section

namespace Riemann.Deep.LeanNBRate

namespace ExpSum

variable (r : ExpSum)

/-- Modulation by `e^{icu}`: nodes shifted by `c`. -/
def modulate (c : ℝ) : ExpSum where
  J := r.J
  node j := r.node j + c
  mult := r.mult
  coef := r.coef

theorem modulate_eval (c u : ℝ) : (r.modulate c).eval u = cexp (I * c * u) * r.eval u := by
  unfold eval modulate
  simp only
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun j _ => ?_
  rw [Finset.mul_sum]
  refine Finset.sum_congr rfl fun k _ => ?_
  have : cexp (I * ((r.node j + c : ℝ) : ℂ) * u) = cexp (I * c * u) * cexp (I * r.node j * u) := by
    rw [← Complex.exp_add]; push_cast; ring_nf
  rw [this]; ring

theorem modulate_M (c : ℝ) : (r.modulate c).M = r.M := rfl

/-- `e^{−icu} (d/du)(e^{icu} r(u))` as an exponential sum (same nodes and multiplicities). -/
def dmod (c : ℝ) : ExpSum where
  J := r.J
  node := r.node
  mult := r.mult
  coef j k := I * (c + r.node j) * r.coef j k +
    (if k + 1 < r.mult j then I * (k + 1 : ℕ) * r.coef j (k + 1) else 0)

theorem dmod_M (c : ℝ) : (r.dmod c).M = r.M := rfl

theorem dmod_node (c : ℝ) : (r.dmod c).node = r.node := rfl

/-- The termwise derivative identity behind `dmod`. -/
theorem dmod_eval_eq (c u : ℝ) : (r.dmod c).eval u =
    ∑ j, ∑ k ∈ Finset.range (r.mult j), r.coef j k *
      (I * (c + r.node j) * (I * u) ^ k + (k : ℂ) * I * (I * u) ^ (k - 1)) *
        cexp (I * r.node j * u) := by
  unfold eval dmod
  simp only
  refine Finset.sum_congr rfl fun j _ => ?_
  have hsplit : ∀ k ∈ Finset.range (r.mult j),
      (I * (c + r.node j) * r.coef j k +
          (if k + 1 < r.mult j then I * (k + 1 : ℕ) * r.coef j (k + 1) else 0)) * (I * u) ^ k *
        cexp (I * r.node j * u) =
      r.coef j k * (I * (c + r.node j) * (I * u) ^ k) * cexp (I * r.node j * u) +
        (if k + 1 < r.mult j then r.coef j (k + 1) * ((k + 1 : ℕ) * I * (I * u) ^ k) *
          cexp (I * r.node j * u) else 0) := by
    intro k _
    split_ifs <;> ring
  rw [Finset.sum_congr rfl hsplit, Finset.sum_add_distrib]
  have hshift : ∑ k ∈ Finset.range (r.mult j), (if k + 1 < r.mult j then r.coef j (k + 1) *
        ((k + 1 : ℕ) * I * (I * u) ^ k) * cexp (I * r.node j * u) else 0) =
      ∑ k ∈ Finset.range (r.mult j), r.coef j k * ((k : ℂ) * I * (I * u) ^ (k - 1)) *
        cexp (I * r.node j * u) := by
    rcases Nat.eq_zero_or_pos (r.mult j) with hm | hm
    · simp [hm]
    · obtain ⟨n, hn⟩ : ∃ n, r.mult j = n + 1 := ⟨r.mult j - 1, by omega⟩
      rw [hn, Finset.sum_range_succ' (fun k => r.coef j k * ((k : ℂ) * I * (I * u) ^ (k - 1)) *
        cexp (I * r.node j * u)), Finset.sum_range_succ]
      simp only [Nat.cast_zero, zero_mul, mul_zero, zero_mul, add_zero, lt_self_iff_false,
        if_false]
      refine Finset.sum_congr rfl fun k hk => ?_
      rw [Finset.mem_range] at hk
      rw [if_pos (by omega)]
      simp only [Nat.add_sub_cancel]
  rw [hshift, ← Finset.sum_add_distrib]
  refine Finset.sum_congr rfl fun k _ => ?_
  ring

/-- `(d/du) (e^{icu} r(u)) = e^{icu} · (r.dmod c)(u)`. -/
theorem hasDerivAt_mod_eval (c u : ℝ) :
    HasDerivAt (fun u : ℝ => cexp (I * c * u) * r.eval u)
      (cexp (I * c * u) * (r.dmod c).eval u) u := by
  have hterm : ∀ (j : Fin r.J) (k : ℕ), HasDerivAt
      (fun u : ℝ => cexp (I * c * u) * (r.coef j k * (I * u) ^ k * cexp (I * r.node j * u)))
      (cexp (I * c * u) * (r.coef j k *
        (I * (c + r.node j) * (I * u) ^ k + (k : ℂ) * I * (I * u) ^ (k - 1)) *
          cexp (I * r.node j * u))) u := by
    intro j k
    have h1 : HasDerivAt (fun u : ℝ => (I * (u : ℂ))) I u := by
      simpa using ((hasDerivAt_id (u : ℂ)).const_mul I).comp_ofReal
    have h2 : HasDerivAt (fun u : ℝ => (I * (u : ℂ)) ^ k) ((k : ℂ) * (I * u) ^ (k - 1) * I) u :=
      h1.pow k
    have h3 : HasDerivAt (fun u : ℝ => cexp (I * r.node j * u))
        (cexp (I * r.node j * u) * (I * r.node j)) u := by
      have : HasDerivAt (fun u : ℝ => (I * r.node j * (u : ℂ))) (I * r.node j) u := by
        simpa using ((hasDerivAt_id (u : ℂ)).const_mul (I * r.node j)).comp_ofReal
      exact this.cexp
    have h4 : HasDerivAt (fun u : ℝ => cexp (I * c * u)) (cexp (I * c * u) * (I * c)) u := by
      have : HasDerivAt (fun u : ℝ => (I * c * (u : ℂ))) (I * c) u := by
        simpa using ((hasDerivAt_id (u : ℂ)).const_mul (I * c)).comp_ofReal
      exact this.cexp
    have := h4.mul ((h2.const_mul (r.coef j k)).mul h3)
    convert this using 1 <;>
      first
        | (ext v; simp only [Pi.mul_apply]; ring)
        | (simp only [Pi.mul_apply]; ring)
  have hsum := HasDerivAt.fun_sum (u := Finset.univ) fun j _ =>
    HasDerivAt.fun_sum (u := Finset.range (r.mult j)) fun k _ => hterm j k
  convert hsum using 1
  · ext u; unfold eval; rw [Finset.mul_sum]; simp_rw [Finset.mul_sum]
  · rw [dmod_eval_eq, Finset.mul_sum]; simp_rw [Finset.mul_sum]

theorem hasDerivAt_eval (u : ℝ) : HasDerivAt r.eval ((r.dmod 0).eval u) u := by
  have := r.hasDerivAt_mod_eval 0 u
  simpa using this

theorem continuous_eval : Continuous r.eval := by
  unfold eval; fun_prop

end ExpSum

end Riemann.Deep.LeanNBRate
